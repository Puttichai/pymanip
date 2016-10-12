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
import time
import random

# For plotting manipulation graphs
from matplotlib import pyplot as plt
from pylab import ion
ion()

import BaseManipulationPlanner as bmp
from BaseManipulationPlanner import (FW, BW, REACHED, ADVANCED, TRAPPED,
                                     TRANSIT, TRANSFER, CG, CP, CGCP,
                                     CONTROLLERSLEEPDUR)
import GPGManipulationPlanner as gpgmp
from ..utils import Utils
from ..utils.Utils import Colorize
from ..utils import GraspPlacementGraph as GPG

_RNG = random.SystemRandom()


############################################################
#    Enhanced Grasp-Placement-Graph Manipulation Planner
############################################################
class EnhancedGPGMP(gpgmp.GPGManipulationPlanner):
    
    def __init__(self, robot, manipulatorname, mobj):
        super(EnhancedGPGMP, self).__init__(robot, manipulatorname, mobj)
        self._nosteering = True


    def PlanTransitTrajectory(self, qinit, T, qgrasp_goal=None, qgoal=None):
        """
        PlanTransitTrajectory plans a one-step transit trajectory from
        qinit to qgoal.
        
        Parameters
        ----------
        qinit : n-vector (initial robot configuration)
        T : 4x4 object transformation
        qgrasp_goal : list; optional
        qgoal : n-vector (goal robot configuration); optional

        Either qgrasp_goal or qgoal has to be provided.

        Returns
        -------
        plannerstatus
        trajectorystring : TOPP format
        """
        if (qgrasp_goal is None) and (qgoal is None):
            print Colorize('[PlanTransitTrajectory] Not enough information;', 
                           'red')
            print Colorize('    Either qgrasp_goal or qgoal needs to be given',
                           'red')
            return [orpy.PlannerStatus.Failed, '']
        
        # Release the grabbed object (if any)
        if (self.robot.IsGrabbing(self.object) is not None):
            self.taskmanip.ReleaseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.Release(self.object)
        
        # Set up the scene
        self.object.SetTransform(T)
        self.robot.SetActiveDOFValues(qinit)
        
        # Compute qgrasp_goal if not given
        if (qgrasp_goal is not None) and (qgoal is None):
            ibox = qgrasp_goal[0]
            graspedlink = self.object.GetLinks()[ibox]
            extents = self.boxinfos[ibox].extents
            Tgripper = Utils.ComputeTGripper(graspedlink.GetTransform(),
                                             qgrasp_goal, extents, unitscale=False)
            qgoal = self.manip.FindIKSolution\
            (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)
            if qgoal is None:
                print Colorize('[PlanTransitTrajectory] qgrasp_goal not feasible',
                               'red')
                return [orpy.PlannerStatus.Failed, '']

        if self._nosteering:
            if self._print:
                print '[PlanTransitTrajectory]',
                print Colorize('Successful', 'green')
            return [1, '']
            
        # Plan a path using Steer
        [plannerstatus, transittraj] = self.Steer(qinit, qgoal,
                                                  self._openravemaxiter,
                                                  self._postprocessing)
        
        if not (plannerstatus == 1):
            if self._print:
                print '[PlanTransitTrajectory]',
                print Colorize('OpenRAVE planner failed', 'red')
            return [plannerstatus, '']

        if self._print:
            print '[PlanTransitTrajectory]',
            print Colorize('Successful', 'green')
        
        # Convert the OpenRAVE trajectory to RampsListND
        rampslist = pu.ConvertOpenRAVETrajToRampsListND(transittraj,
                                                        self._vmax, self._amax, DELTA)
        
        return [plannerstatus, str(rampslist)]


    def PlanTransferTrajectory(self, qinit, qgrasp, Tinit, qgoal=None, Tgoal=None):
        """
        PlanTransferTrajectory plans a one-step transfer trajectory from
        qinit to qgoal.
        
        Parameters
        ----------
        qinit : n-vector (initial robot configuration)
        qgrasp : list
        Tinit : 4x4 object transformation
        qgoal : n-vector; optional
        Tgoal : 4x4 object transformation; optional

        Either qgoal or Tgoal has to be provided

        Returns
        -------
        plannerstatus
        trajectorystring : TOPP format
        """
        if (qgoal is None) and (Tgoal is None):
            print Colorize('[PlanTransferTrajectory] Not enough information;', 
                           'red')
            print Colorize('    Either qgoal or Tgoal needs to be given',
                           'red')
            return [orpy.PlannerStatus.Failed, '']

        # Release the grabbed object (if any)
        if (self.robot.IsGrabbing(self.object) is not None):
            self.taskmanip.ReleaseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.Release(self.object)

        # Compute qgoal if not given
        if (qgoal is None) and (Tgoal is not None):
            ibox = qgrasp[0]
            graspedlink = self.object.GetLinks()[ibox]
            extents = self.boxinfos[ibox].extents
            self.object.SetTransform(Tgoal)
            Tgripper = Utils.ComputeTGripper(graspedlink.GetTransform(),
                                             qgrasp, extents, unitscale=False)
            
            # Set robot reference configuration before solving for IK solutions
            self.object.SetTransform(Tinit)
            self.robot.SetActiveDOFValues(qinit)
            qgoal = self.manip.FindIKSolution\
            (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)
            if qgoal is None:
                print Colorize('[PlanTransitTrajectory] qgrasp_goal not feasible',
                               'red')
                return [orpy.PlannerStatus.Failed, '']
        
        # Set up the scene
        self.object.SetTransform(Tinit)
        self.robot.SetActiveDOFValues(qinit)
            
        # Grab the object
        holdobject = True
        isgrasping = self.CheckGrasp(qinit, Tinit, holdobject)
        if not isgrasping:
            print '[PlanTransferTrajectory]',
            print Colorize('G R A S P  E R R O R', 'red')
            raw_input() # pause

        if self._nosteering:
            # Relase the grabbed object
            self.taskmanip.ReleaseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.Release(self.object)
            
            if self._print:
                print '[PlanTransitTrajectory]',
                print Colorize('Successful', 'green')
            return [1, '']            

        # Plan a path using Steer
        [plannerstatus, transfertraj] = self.Steer(qinit, qgoal,
                                                   self._openravemaxiter,
                                                   self._postprocessing)
        
        # Relase the grabbed object
        self.taskmanip.ReleaseFingers()
        while not self.robot.GetController().IsDone():
            time.sleep(CONTROLLERSLEEPDUR)
        self.robot.Release(self.object)
        
        if not (plannerstatus == 1):
            if self._print:
                print '[PlanTransferTrajectory]',
                print Colorize('OpenRAVE planner failed', 'red')            
            return [plannerstatus, '']

        if self._print:
            print '[PlanTransferTrajectory]',
            print Colorize('Successful', 'green')
        
        # Convert the OpenRAVE trajectory to RampsListND
        rampslist = pu.ConvertOpenRAVETrajToRampsListND(transfertraj,
                                                        self._vmax, self._amax, DELTA)
        
        return [plannerstatus, str(rampslist)]


