#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2015 Puttichai Lertkultanon <L.Puttichai@gmail.com>
#
# This file is part of pymanip.
#
# pymanip is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# pymanip is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# pymanip. If not, see <http://www.gnu.org/licenses/>.

import openravepy as orpy
import numpy as np
import random
import time

from pymanip.utils import Utils
from pymanip.utils.Utils import Colorize

# cvxopt for object tracking planner
import cvxopt
import cvxopt.solvers
cvxopt.solvers.options['show_progress'] = False # disable cvxopt output

############################################################
#                      Object Tracker
############################################################
class ObjectTracker(object):
    """
    This class does the same thing as OpenRAVE's
    workspacetrajectorytracker does. Given an object trajectory (SE(3)
    trajectory), it basically calls IKFast to compute the robot's IK
    solution at each timestep.

    No-retiming is important when we plan closed-chain motions, i.e.,
    manipulation with multiple robots, since we also need motion
    synchronization.
    """
    def __init__(self, robotslist, manipulatornames, mobj):
        self.robots = robotslist
        self.manips = []
        for robot, manipname in zip(self.robots, manipulatornames):
            self.manips.append(robot.SetActiveManipulator(manipname))
            
        self.active_dofs = self.manips[0].GetArmIndices()
        
        for robot in self.robots:
            robot.SetActiveDOFs(self.active_dofs)
        self.object = mobj
        
        self._vmax = self.robots[0].GetDOFVelocityLimits()\
        [0:self.robots[0].GetActiveDOF()]
        
        self._amax = self.robots[0].GetDOFAccelerationLimits()\
        [0:self.robots[0].GetActiveDOF()]

        self._print = True
        self._hasquery = False
        self._checkclosestsolution = True
        
        
    def Track(self, lietraj, transtraj, qgrasps, timestep, q0list=None, 
              transformationslist=[]):
        """
        transformationslist provides object's transformations at every
        timestep. If transformationslist is not given, we compute it
        according to timestep. If it is given, we assume that it is
        the one computed according to timestep.

        Returns
        -------
        plannerstatus : bool
        waypointslist : list
            waypointslist is a list of N lists, where N is the number 
            of robots.
        timestamps : list
        """
        _funcname = '[ObjectTracker::Track] '

        graspedlinkslist = [self.object.GetLinks()[qgrasp[0]] 
                            for qgrasp in qgrasps]
        extentslist = [graspedlink.GetGeometries()[0].GetBoxExtents()
                       for graspedlink in graspedlinkslist]
        
        dur = lietraj.duration
        # Relative transformation between the object's COM and the COM
        # of the grasped link
        Trelslist = [np.dot(np.linalg.inv(self.object.GetTransform()),
                            graspedlink.GetTransform())
                     for graspedlink in graspedlinkslist]
        
        timestamps = np.arange(0, dur, timestep).tolist()
        if (abs(timestamps[-1] - dur) > 1e-8):
            timestamps.append(dur)

        # Compute object transformations at each time step
        if len(transformationslist) < 1:
            self._transformationslist = []
            
            M = np.eye(4) # dummy variable for storing an object's transformation
            for t in timestamps:
                M[0:3, 0:3] = lietraj.EvalRotation(t)
                M[0:3, 3] = transtraj.Eval(t)
                self._transformationslist.append(np.array(M))
                
        else:
            self._transformationslist = copy.deepcopy(transformationslist)
            
        # Check if the first solution really exists
        if q0list is None:
            q0list = []
            for (i, robot) in enumerate(self.robots):
                Tgripper = Utils.ComputeTGripper\
                (np.dot(self._transformationslist[0], Trelslist[i]), 
                 qgrasps[i], extentslist[i])
                
                q0 = self.manips[i].FindIKSolution\
                (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)
                
                if q0 is None:
                    if self._print:
                        message = Colorize\
                        ('No IK solution for robot {0} at t = 0.'.format(i), 'red')
                        print _funcname + message
                    return [False, None, None]
                else:
                    q0list.append(q0)
        
        waypointslist = []
        for q0 in q0list:
            waypointslist.append([q0])
            
        for i, t in enumerate(timestamps[1:]):
            noik = False

            for (j, robot) in enumerate(self.robots):
                Tgripper = Utils.ComputeTGripper\
                (np.dot(self._transformationslist[i], Trelslist[j]), 
                 qgrasps[j], extentslist[j])

                if self._checkclosestsolution:
                    solutionset = self.manips[j].FindIKSolutions\
                    (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)

                    if len(solutionset) < 1:
                        if self._print:
                            message = Colorize\
                            ('No IK solution for robot {0} at t = {1}'.format(j, t), 
                             'red')
                            print _funcname + message

                        return [False, None, None]
                    
                    sol = self._FindClosestIKSolution(solutionset, 
                                                      waypointslist[j][i - 1])
                else:
                    sol = self.manips[j].FindIKSolution\
                    (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)

                    if sol is None:
                        if self._print:
                            message = Colorize\
                            ('No IK solution for robot {0} at t = {1}'.format(j, t), 
                             'red')
                            print _funcname + message

                        return [False, None, None]

                waypointslist[j].append(sol)

        return [True, waypointslist, timestamps]


    def _FindClosestIKSolution(self, iksolutions, ikref):
        closestsol = iksolutions[0]
        mindist = np.sum(abs(ikref - iksolutions[0]))
        for sol in iksolutions[1:]:
            dist = np.sum(abs(ikref - sol))
            if dist < mindist:
                closestsol = sol
                mindist = dist

        return closestsol
                

# ############################################################
# #                      Object Tracker
# ############################################################
# class ObjectTracker(object):
#     """
#     This class does the same thing as OpenRAVE's
#     workspacetrajectorytracker does. Given an object trajectory (SE(3)
#     trajectory), it basically calls IKFast to compute the robot's IK
#     solution at each timestep.

#     No-retiming is important when we plan closed-chain motions, i.e.,
#     manipulation with multiple robots, since we also need motion
#     synchronization.
#     """
#     def __init__(self, robotslist, manipulatornames, mobj):
#         self.robots = robotslist
#         self.manips = []
#         for robot, manipname in zip(self.robots, manipulatornames):
#             self.manips.append(robot.SetActiveManipulator(manipname))
            
#         self.active_dofs = self.manips[0].GetArmIndices()
        
#         for robot in self.robots:
#             robot.SetActiveDOFs(self.active_dofs)
#         self.object = mobj
        
#         self._vmax = self.robots[0].GetDOFVelocityLimits()\
#         [0:self.robots[0].GetActiveDOF()]
        
#         self._amax = self.robots[0].GetDOFAccelerationLimits()\
#         [0:self.robots[0].GetActiveDOF()]

#         self._print = True
#         self._hasquery = False
#         self._checkclosestsolution = True


#     def Track(self, lietraj, transtraj, qgrasps, timestep, q0list=None,
#               qd0list=None, transformationslist=[]):
#         """
#         transformationslist provides object's transformations at every
#         timestep. If transformationslist is not given, we compute it
#         according to timestep. If it is given, we assume that it is
#         the one computed according to timestep.

#         Returns
#         -------
#         plannerstatus : bool
#         waypointslist : list
#             waypointslist is a list of N lists, where N is the number 
#             of robots.
#         timestamps : list
#         """
#         _funcname = '[ObjectTracker::Track] '

#         graspedlinkslist = [self.object.GetLinks()[qgrasp[0]] 
#                             for qgrasp in qgrasps]
#         extentslist = [graspedlink.GetGeometries()[0].GetBoxExtents()
#                        for graspedlink in graspedlinkslist]
        
#         dur = lietraj.duration
#         # Relative transformation between the object's COM and the COM
#         # of the grasped link
#         Trelslist = [np.dot(np.linalg.inv(self.object.GetTransform()),
#                             graspedlink.GetTransform())
#                      for graspedlink in graspedlinkslist]
        
#         timestamps = np.arange(0, dur, timestep).tolist()
#         if (abs(timestamps[-1] - dur) > 1e-8):
#             timestamps.append(dur)

#         # Compute object transformations at each time step
#         if len(transformationslist) < 1:
#             self._transformationslist = []
            
#             M = np.eye(4) # dummy variable for storing an object's transformation
#             for t in timestamps:
#                 M[0:3, 0:3] = lietraj.EvalRotation(t)
#                 M[0:3, 3] = transtraj.Eval(t)
#                 self._transformationslist.append(np.array(M))
                
#         else:
#             self._transformationslist = copy.deepcopy(transformationslist)
            
#         # Check if the first solution really exists
#         if q0list is None:
#             q0list = []
#             for (i, robot) in enumerate(self.robots):
#                 Tgripper = Utils.ComputeTGripper\
#                 (np.dot(self._transformationslist[0], Trelslist[i]), 
#                  qgrasps[i], extentslist[i])
                
#                 q0 = self.manips[i].FindIKSolution\
#                 (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)
                
#                 if q0 is None:
#                     if self._print:
#                         message = Colorize\
#                         ('No IK solution for robot {0} at t = 0.'.format(i), 'red')
#                         print _funcname + message
#                     return [False, None, None]
#                 else:
#                     q0list.append(q0)
        
#         waypointslist = []
#         for q0 in q0list:
#             waypointslist.append([q0])
            
#         for i, t in enumerate(timestamps[1:]):
#             noik = False

#             for (j, robot) in enumerate(self.robots):
#                 Tgripper = Utils.ComputeTGripper\
#                 (np.dot(self._transformationslist[i], Trelslist[j]), 
#                  qgrasps[j], extentslist[j])

                

#                 waypointslist[j].append(sol)

#         return [True, waypointslist, timestamps]

        

        
            
        
        
