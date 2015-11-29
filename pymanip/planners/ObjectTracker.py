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
    def __init__(self, robot, manipulatorname, mobj):
        self.robot = robot
        self.manip = self.robot.SetActiveManipulator(manipulatorname)
        self.active_dofs = self.manip.GetArmIndices()
        self.robot.SetActiveDOFs(self.active_dofs)
        self.object = mobj
        
        self._vmax = self.robot.GetDOFVelocityLimits()[0:self.robot.GetActiveDOF()]
        self._amax = self.robot.GetDOFAccelerationLimits()[0:self.robot.GetActiveDOF()]

        self._print = True
        self._hasquery = False
        self._checkclosestsolution = True
        
        
    def Track(self, lietraj, transtraj, qgrasp, timestep, q0=None, 
              transformationslist=[]):
        """
        transformationslist provides object's transformations at every
        timestep. If transformationslist is not given, we compute it
        according to timestep. If it is given, we assume that it is
        the one computed according to timestep.

        Returns
        -------
        plannerstatus : bool
        waypoints : list
        timestamps : list
        """
        _funcname = '[ObjectTracker::Track] '

        graspedlink = self.object.GetLinks()[qgrasp[0]]
        extents = graspedlink.GetGeometries()[0].GetBoxExtents()
        
        dur = lietraj.duration
        # Relative transformation between the object's COM and the COM
        # of the grasped link
        Trel = np.dot(np.linalg.inv(self.object.GetTransform()),
                      graspedlink.GetTransform())
        
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
        if q0 is None:
            Tgripper = Utils.ComputeTGripper(np.dot(self._transformationslist[0], Trel), 
                                             qgrasp, extents)
            q0 = self.manip.FindIKSolution(Tgripper,
                                           orpy.IkFilterOptions.CheckEnvCollisions)
            if q0 is None:
                if self._print:
                    message = Colorize('No IK solution at t = 0.', 'red')
                    print _funcname + message
                return [False, None, None]

        # tstartplan = time.time()
        waypoints = [q0] # containing ik solutions at each time step
        for i, t in enumerate(timestamps[1:]):
            noik = False
            Tgripper = Utils.ComputeTGripper(np.dot(self._transformationslist[i], Trel), 
                                             qgrasp, extents)

            if self._checkclosestsolution:
                solutionset = self.manip.FindIKSolutions\
                (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)

                if len(solutionset) < 1:
                    if self._print:
                        message = Colorize('No Ik solution at t = {0}.'.format(t), 'red')
                        print _funcname + message
                    return [False, None, None]
                
                sol = self._FindClosestIKSolution(solutionset, waypoints[i - 1])
            else:
                sol = self.manip.FindIKSolution(Tgripper, 
                                                orpy.IkFilterOptions.CheckEnvCollisions)
                
                if sol is None:   
                    if self._print:
                        message = Colorize('No Ik solution at t = {0}.'.format(t), 'red')
                        print _funcname + message
                    return [False, None, None]
                                
            waypoints.append(sol)
        # tendplan = time.time()
        # print Colorize('planning time = {0}'.format(tendplan - tstartplan), 'green')

        return [True, waypoints, timestamps]


    def _FindClosestIKSolution(self, iksolutions, ikref):
        closestsol = iksolutions[0]
        mindist = np.sum(abs(ikref - iksolutions[0]))
        for sol in iksolutions[1:]:
            dist = np.sum(abs(ikref - sol))
            if dist < mindist:
                closestsol = sol
                mindist = dist

        return closestsol
                


        

        
            
        
        
