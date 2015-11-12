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
import random
import numpy as np
import time
import copy
import TOPP

from ..utils import Heap
from ..utils import ObjectPreprocessing as op
from ..utils import parabint_utilsformanip as pu
from ..utils import Grasp as gr
from ..utils.Grasp import (pX, pY, pZ, mX, mY, mZ)
from ..utils import Utils
from ..utils.Utils import Colorize

# cvxopt for object tracking planner
import cvxopt
import cvxopt.solvers
cvxopt.solvers.options['show_progress'] = False # disable cvxopt output

############################################################
#                    Global Parameters
############################################################
# Planner parameters
FW = 0
BW = 1
REACHED = 0
ADVANCED = 1
TRAPPED = 2

# Unimanual
CG = 1
CP = 2
CGCP = 3
WORDS = ['', 'CG', 'CP', 'CGCP']
DICTIONARY = {'CG': 'CG', 
              'CP': 'CP',
              'CGCP': 'CGCP'}

# Manipulation trajectory types
TRANSIT = 0
TRANSFER = 1

EPSILON = 1e-10
DELTA = 0 # for patabint

CONTROLLERTIMEOUT = 5.
CONTROLLERSLEEPDUR = 0.001
OPENRAVEPLANNERNAME = 'birrt'


############################################################
#                Base Manipulation Planner
############################################################
class BaseManipulationPlanner(object):
    """
    This class implements basic utilities required for manipulation
    planning.
    Ready-to-use functions are
    - Steer
    - PlanTransitTrajectory
    - PlanTransferTrajectory
    - PlayManipulationTrajectory
    - PlayTransitTrajectory
    - PlayTransferTrajectory
    """
    #
    # Config
    #
    class Config(object):
        """
        Class Config

        Parameters
        ----------
        qrobot : n-vector
            A robot configuration (excluding the gripper).
        tobj : 4x4 transformation matrix
            A transformation matrix of the movable object at the current 
            configuration.
        qgrasp : list
            qgrasp must be in the format [ibox, approachingdir, slidingdir, delta].
            qgrasp must be None if the current configtype is CP.
            ibox : an index of the box being grasped (a link index).
            approachingdir : pX, pY, pZ, mX, mY, or mZ
            slidingdir : pX, pY, or, pZ
            delta : a position of the gripper along the sliding direction.
        configtype : CP, CG, or CGCP
            A type of this composite configuration.
        isurface : int
            An index of the contact surface.
            If isurface is unknown, we need to create a Config via 
            MRRTPlanner.CreateConfig
            in order to specify isurface automatically.
        """

        def __init__(self, qrobot, tobj, qgrasp, configtype, isurface):

            self.qrobot = qrobot
            self.tobj = tobj
            self.qgrasp = qgrasp
            self.type = configtype
            self.isurface = isurface

            if (configtype == CP):
                self.approachingdir = None
            else:
                self.approachingdir = qgrasp[1]


        def __str__(self):
            string = "configtype = {0}\n".format(DICTIONARY[WORDS[self.type]])
            string += "qrobot = {0}\n".format(self.qrobot)
            string += "tobj = {0}\n".format(self.tobj)
            string += "qgrasp = {0}".format(self.qgrasp)
            return string


        def Clone(self):
            qrobot = copy.copy(self.qrobot)
            tobj = copy.copy(self.tobj)
            qgrasp = copy.copy(self.qgrasp)
            configtype = copy.copy(self.type)
            isurface = copy.copy(self.isurface)
            return  BaseManipulationPlanner.Config\
            (qrobot, tobj, qgrasp, configtype, isurface)


    #
    # Vertex
    #
    class Vertex(object):
        """
        Class Vertex

        Parameters
        ----------
        config : Config
        parentindex : int
            Index of its parent
        vertextype : FW or BW
            Indicates in which tree the vertex is
        level : int
            The number of levels from the root.
            The root is at level 0.
        trajtype : TRANSIT or TRANSFER
        index : int
            The index of this vertiex in tree.verticeslist.
            This number will be assigned when added to a tree.
        id : tuple
            A vertex id is (pathlevel, isurface, approachingdir).
            id has to contain pathlevel since now the graph will possibly 
            have cycles.
        """    

        def __init__(self, config, vertextype=FW, trajtype=None, level=0):
            self.config = config
            self.parentindex = None
            self.vertextype = vertextype
            self.trajectory = None
            self.level = level
            self.trajtype = trajtype
            self.index = 0 ## to be assigned when added to the tree
            self.id = (self.level, self.config.isurface, self.config.approachingdir)


        def __str__(self):
            string = "vindex = {0}\n".format(self.index)
            string += str(self.config)
            return string


        def Clone(self):
            vnew = BaseManipulationPlanner.Vertex(self.config.Clone())

            vnew.id = copy.deepcopy(self.id)

            ## float and string are safe
            vnew.parentindex = self.parentindex
            vnew.vertextype = self.vertextype
            vnew.trajectory = self.trajectory
            vnew.level = self.level
            vnew.trajtype = self.trajtype
            vnew.index = self.index

            return vnew


    #
    # Tree
    #
    class Tree(object):

        def __init__(self, vroot, treetype=FW):
            self.verticeslist = []
            self.verticeslist.append(vroot)
            self.treetype = treetype
            self.length = 1


        def __len__(self):
            return len(self.verticeslist)


        def __getitem__(self, index):
            return self.verticeslist[index]


        def Clone(self):
            tnew = BaseManipulationPlanner.Tree(None, self.treetype)
            tnew.verticeslist = [v.Clone() for v in self.verticeslist]
            tnew.length = self.length
            return tnew


        def AddVertex(self, parentindex, vnew, trajectory, trajtype):
            parent = self.verticeslist[parentindex]
            vnew.parentindex = parentindex
            # vnew.level = parent.level + 1
            vnew.trajectory = trajectory
            vnew.trajtype = trajtype
            vnew.index = self.length
            self.verticeslist.append(vnew)
            self.length += 1


        def GenerateManipulationTrajectory(self, vindex=-1):
            trajslist = []
            trajtypeslist = []
            vertex = self.verticeslist[vindex]
            while (vertex.parentindex != None):
                parent = self.verticeslist[vertex.parentindex]
                trajslist.append(vertex.trajectory)
                trajtypeslist.append(vertex.trajtype)
                vertex = parent

            if (self.treetype == FW):
                trajslist = trajslist[::-1]
                trajtypeslist = trajtypeslist[::-1]

            return [trajslist, trajtypeslist]


        def GenerateObjectTransformationsList(self, vindex=-1):
            tobjlist = []
            vertex = self.verticeslist[vindex]
            while (vertex.parentindex != None):
                parent = self.verticeslist[vertex.parentindex]
                tobjlist.append(vertex.config.tobj)
                vertex = parent
            tobjlist.append(self.verticeslist[0].config.tobj)

            if (self.treetype == FW):
                tobjlist = tobjlist[::-1]

            return tobjlist
   
    
    #
    # BaseManipulationPlanner initialization
    #
    def __init__(self, robot, manipulatorname, mobj):
        self.robot = robot
        self.manip = self.robot.SetActiveManipulator(manipulatorname)
        self.active_dofs = self.manip.GetArmIndices()
        self.object = mobj
        self.taskmanip = orpy.interfaces.TaskManipulation(self.robot)
        self.env = self.robot.GetEnv()
        
        # Activate only the first active_dofs joints
        self.robot.SetActiveDOFs(self.active_dofs)
        self._vmax = self.robot.GetDOFVelocityLimits()[self.active_dofs]
        self._amax = self.robot.GetDOFAccelerationLimits()[self.active_dofs]
        self._RNG = random.SystemRandom()

        

        self._print = True
        self._openravemaxiter = 1000 # max iterations for OpenRAVE planners
        self._shortcutiter = 300
        self._postprocessing = False # OpenRAVE default postprocessing

        self._setsampleparameters = False
        self._hasquery = False
        self._hassolution = False
        self._runningtime = 0.0
        self._iterations = 0

        self.PreprocessObjectInfo()
        
        
    @property
    def setsampleparameters(self):
        return self._setsampleparameters
        
    
    @property
    def hasquery(self):
        return self._hasquery
    
    
    @property
    def hassolution(self):
        return self._hassolution
        
    
    @property
    def runningtime(self):
        return self._runningtime
        
    
    @property
    def iterations(self):
        return self._iterations

        
    def PreprocessObjectInfo(self):
        """
        Stores information about object's placement & grasp classes
        """
        # S contains stable contact surfaces's transformation frame
        # with respect to its COM (except S[0] which is the relative
        # transformation between COM and the first link's frame). see
        # ObjectPreprocess.py for more detail about S
        self.S = op.PlacementPreprocess(self.object)

        TCOM = np.array(self.S[0])
        # The z-axis of each surface frame is pointing out of the
        # object. However, what we need is a frame pointing into the
        # object (so that the z-axis is aligned with that of the
        # world).
        Toffset = Utils.ComputeTRot(pX, np.pi)
        # transformationset contains transformation T such that when
        # assigning self.object.SetTransform(T), self.object is
        # resting at a stable placement.
        self.transformationset = [np.dot(Toffset, np.linalg.inv(np.dot(TCOM, T)))
                                  for T in self.S[1:]]

        if self._print:
            print Colorize('Processing the object geometries', 'yellow')
        self.boxinfos = []
        for ibox in xrange(len(self.object.GetLinks())):
            ts = time.time()
            boxinfo = op.BoxInfo(self.object, ibox)
            boxinfo.GetPossibleSlidingDirections()
            boxinfo.Preprocess(self.transformationset)
            te = time.time()
            if self._print:
                print Colorize('  box {0} took {1} sec.'.format(ibox, te - ts), 
                               'yellow')
            self.boxinfos.append(boxinfo)
            
            
    def CheckGrasp(self, qrobot, Tobj, holdobject=False):
        """
        Check whether self.robot at qrobot can correctly grasp the
        object at Tobj.  If holdobject is True and the grasp is valid,
        the robot will still continue to hold the object.
        """
        isgrasping = False
        self.robot.SetActiveDOFValues(qrobot)
        self.object.SetTransform(Tobj)
        self.taskmanip.CloseFingers()
        while not self.robot.GetController().IsDone():
            time.sleep(CONTROLLERSLEEPDUR)
        if self.env.CheckCollision(self.robot, self.object):
            isgrasping = True

        # This grabbing is required regardless of holdobject.
        self.robot.Grab(self.object)
        
        if isgrasping and holdobject:
            return isgrasping
            
        self.taskmanip.ReleaseFingers()
        while not self.robot.GetController().IsDone():
            time.sleep(CONTROLLERSLEEPDUR)
        self.robot.Release(self.object)
        return isgrasping


    def ComputeGraspConfiguration(self, qgrasp, Tobj, qrobot_ref=None):
        """
        Returns a robot configuration (nearest to qrobot_ref) that can
        grasp the object, placed at Tobj, with a grasp identified by
        qgrasp.
        """
        with self.env:
            self.object.SetTransform(Tobj)
            if qrobot_ref is None:
                qrobot_ref = np.zeros(self.active_dofs)
            self.robot.SetActiveDOFValues(qrobot_ref)
            ibox = qgrasp[0]
            Tgripper = Utils.ComputeTGripper(self.object.GetLinks()[ibox].GetTransform(),
                                             qgrasp, self.boxinfos[ibox].extents)
            sol = self.manip.FindIKSolution(Tgripper, 
                                            orpy.IkFilterOptions.CheckEnvCollisions)
        if sol is None:
            return None
        if not self.CheckGrasp(sol, Tobj):
            # Invalud grasp
            return None
        return sol

    
    #
    # Local planners
    #
    def Steer(self, qinit, qgoal, nmaxiterations, postprocessing=False):
        """
        Steer is a local planner that plans a path between qinit and
        qgoal.
        Currecnt implementation uses OpenRAVE RRT.

        Parameters
        ----------
        qinit : n-vector
        qgoal : n-vector
        nmaxiterations : int
        postprocessing : bool, optional
            If True, OpenRAVE planner will run parabolic smoothing on 
            the path obtained. Otherwise, OpenRAVE planner will just 
            returns waypoints.

        Returns
        -------
        plannerstatus
        openravetraj
        """
        # Prepare planning parameters
        params = orpy.Planner.PlannerParameters()
        params.SetRobotActiveJoints(self.robot)
        params.SetInitialConfig(qinit)
        params.SetGoalConfig(qgoal)
        extraparams = ''
        extraparams += '<_nmaxiterations>' + str(nmaxiterations) + '</_nmaxiterations>'
        if not postprocessing:
            # Force no post-processing
            extraparams += '<_postprocessing></_postprocessing>'
        params.SetExtraParameters(extraparams)
        
        with self.env:
            # Intialize openrave planner
            planner = orpy.RaveCreatePlanner(self.env, OPENRAVEPLANNERNAME)
            res = planner.InitPlan(self.robot, params)
            if not res:
                # This is probably because some numerical errors??
                return [False, None]

            # Create an empty trajectory
            openravetraj = orpy.RaveCreateTrajectory(self.env, '')

            # Start planning        
            plannerstatus = planner.PlanPath(openravetraj)
        
        return [plannerstatus, openravetraj]

    
    def _PlanTransitTrajectory(self, vinit, vgoal):
        """
        PlanTransitTrajectory plans a one-step transit trajectory from
        vstart to vgoal.
        
        Parameters
        ----------
        vstart : Vertex
        vgoal : Vertex

        Returns
        -------
        plannerstatus
        trajectorystring : TOPP format
        """
        qinit = vinit.config.qrobot
        qgoal = vgoal.config.qrobot
        Tobj = vinit.config.tobj
        return self.PlanTransitTrajectory(qinit, Tobj, qgoal=qgoal)


    def _PlanTransferTrajectory(self, vinit, vgoal):
        """
        PlanTransferTrajectory plans a one-step transfer trajectory
        from vstart to vgoal.
        
        Parameters
        ----------
        vstart : Vertex
        vgoal : Vertex

        Returns
        -------
        plannerstatus
        trajectorystring : TOPP format
        """
        qinit = vinit.config.qrobot
        qgrasp = vinit.config.qgrasp
        qgoal = vgoal.config.qrobot
        Tobj = vinit.config.tobj
        return self.PlanTransferTrajectory(qinit, qgrasp, Tobj, qgoal=qgoal)


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

    
    #
    # Extension & connection utilities
    #
    def Distance(self, c0, c1, metrictype=1):
        if (metrictype == 1):
            delta_qrobot = c1.qrobot - c0.qrobot
            dist_qrobot = np.dot(delta_qrobot, delta_qrobot)
            
            ## distance in SE(3)
            p0 = c0.tobj[0:3, 3]
            p1 = c1.tobj[0:3, 3]
            R0 = c0.tobj[0:3, 0:3]
            R1 = c1.tobj[0:3, 0:3]
            
            delta_x = p1 - p0
            rvect = logvect(np.dot(R0.T, R1))

            """
            Object's rotation costs more since the robot needs to move
            a lot more when the object rotates than when it translates.
            """

            wt = 0.2 ## weight for translational distance
            wr = 0.8 ## weight for rotational distance
            
            dist_tobj = wt*np.dot(delta_x, delta_x) + wr*np.dot(rvect, rvect)
            
            return dist_qrobot + dist_tobj

        else:
            raise InputError('metrictype', metrictype)

        
    def NearestNeighborIndices(self, c_rand, treetype, newnn = None):     
        if (treetype == FW):
            vlist = self.treestart.verticeslist
        else:
            vlist = self.treeend.verticeslist
        nv = len(vlist)

        distancelist = [self.Distance(c_rand, v.config, self.metrictype) for v in vlist]
        distanceheap = Heap.Heap(distancelist)

        if (newnn == None):
            if (self.nn < 0):
                nn = nv
            else:
                nn = min(self.nn, nv)
        else:
            if (newnn < 0):
                nn = newnn
            else:
                nn = min(newnn, nv)

        nnindices = [distanceheap.ExtractMin()[0] for i in range(nn)]
        return nnindices

    
    #
    # Trajectory extraction
    #
    def GenerateFinalManipulationTrajectory(self):
        if (not self.hassolution):
            print "The planner has not found any solution so far."
            return ''
        
        [trajslist0, trajtypeslist0] = self.treestart.GenerateManipulationTrajectory()
        [trajslist1, trajtypeslist1] = self.treeend.GenerateManipulationTrajectory()
        if ((not (self.connectingtrajectorystring == '')) and 
            (self.connectingtrajectorytype is not None)):
            trajslist0.append(self.connectingtrajectorystring)
            trajtypeslist0.append(self.connectingtrajectorytype)
        trajslist = trajslist0 + trajslist1
        trajtypeslist = trajtypeslist0 + trajtypeslist1
        return [trajslist, trajtypeslist]

    
    def GenerateFinalObjectTransformationsList(self):
        if (not self.hassolution):
            print "The planner has not found any solution so far."
            return ''
        
        qobjlist0 = self.treestart.GenerateObjectTransformationsList()
        qobjlist1 = self.treeend.GenerateObjectTransformationsList()
        return qobjlist0 + qobjlist1

    
    #
    # Shortcutting
    #
    def ShortcutManipulationTrajectory(self, trajslist, trajtypeslist, Tobjslist,
                                       print_=False, plot_=False):        
        # Reset the scene
        temp = TOPP.Trajectory.PiecewisePolynomialTrajectory.FromString(trajslist[0])
        self.robot.SetActiveDOFValues(temp.Eval(0))
        if self.robot.IsGrabbing(self.object) is not None:
            self.taskmanip.ReleaseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.ReleaseAllGrabbed()
            
        ntrajs = len(trajslist)
        newtrajslist = []
        newtrajtypeslist = []
        
        # Start main shortcutting loop
        index = 0
        curtrajtype = trajtypeslist[index]
        curtobj = Tobjslist[index]
        rampslistnd = pu.RampsListND.FromString(trajslist[index])
        while (index < ntrajs - 1):
            index += 1
            if (trajtypeslist[index] == curtrajtype):
                rampslistnd.Append(pu.RampsListND.FromString(trajslist[index]))
            else:
                # Shortcut previously appended trajectories
                self.object.SetTransform(curtobj)
                self.robot.SetActiveDOFValues(rampslistnd.x0vect)
                if (curtrajtype == TRANSIT):
                    self.taskmanip.ReleaseFingers()
                    while not self.robot.GetController().IsDone():
                        time.sleep(CONTROLLERSLEEPDUR)
                    self.robot.ReleaseAllGrabbed()
                else:
                    self.taskmanip.CloseFingers()
                    while not self.robot.GetController().IsDone():
                        time.sleep(CONTROLLERSLEEPDUR)
                    self.robot.Grab(self.object)
                    
                shortcuttraj = pu.Shortcut\
                (rampslistnd, self._vmax, self._amax, DELTA, self._shortcutiter, 
                 PRINT=print_, robot=self.robot, PLOT=plot_)
                    
                newtrajslist.append(str(shortcuttraj)) # store in TOPP format
                newtrajtypeslist.append(curtrajtype)
            
                curtrajtype = trajtypeslist[index]
                curtobj = Tobjslist[index]
                rampslistnd = pu.RampsListND.FromString(trajslist[index])
            
        # Shortcut previously appended trajectories
        self.object.SetTransform(curtobj)
        self.robot.SetActiveDOFValues(rampslistnd.x0vect)
        if (curtrajtype == TRANSIT):
            self.taskmanip.ReleaseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.ReleaseAllGrabbed()
        else:
            self.taskmanip.CloseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.Grab(self.object)
            
        shortcuttraj = pu.Shortcut\
        (rampslistnd, self._vmax, self._amax, DELTA, self._shortcutiter, 
         PRINT=print_, robot=self.robot, PLOT=plot_)

        newtrajslist.append(str(shortcuttraj)) # store in TOPP format
        newtrajtypeslist.append(curtrajtype)

        # Release fingers after a TRANSFER trajectory
        if (curtrajtype == TRANSFER):
            self.taskmanip.ReleaseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.ReleaseAllGrabbed()
            
        return [newtrajslist, newtrajtypeslist]

    
    #
    # Visualization
    #
    def PlayManipulationTrajectory(self, trajslist, trajtypeslist, Tinit, Tobjslist=[],
                                   dt=0.01, timemult=1.0):
        hasref = (len(Tobjslist) > 0)
        if not hasref:
            self.object.SetTransform(Tinit)
        for i in range(len(trajslist)):
            trajstring = trajslist[i]
            trajtype = trajtypeslist[i]
            if hasref:
                T = Tobjslist[i]
            else:
                T = self.object.GetTransform()
            
            if trajtype == TRANSIT:
                self.PlayTransitTrajectory(trajstring, T, dt=dt, timemult=timemult)
            else:
                self.PlayTransferTrajectory(trajstring, T, dt=dt, timemult=timemult)
                
    
    def PlayTransitTrajectory(self, transittraj, Tobj, dt=0.01, timemult=1.0):
        topptransit = TOPP.Trajectory.PiecewisePolynomialTrajectory.FromString\
        (transittraj)

        sleeptime = dt*timemult

        # Set initial configuration
        if self.robot.IsGrabbing(self.object) is not None:
            self.taskmanip.ReleaseFingers()
            while not self.robot.GetController().IsDone():
                time.sleep(CONTROLLERSLEEPDUR)
            self.robot.ReleaseAllGrabbed()
        
        self.object.SetTransform(Tobj)
        T = np.append(np.arange(0, topptransit.duration, dt), topptransit.duration)
        for t in T:
            q = topptransit.Eval(t)
            self.robot.SetActiveDOFValues(q)
            time.sleep(sleeptime)
        
    
    def PlayTransferTrajectory(self, transfertraj, Tobj, dt=0.01, timemult=1.0):
        topptransfer = TOPP.Trajectory.PiecewisePolynomialTrajectory.FromString\
        (transfertraj)
        
        sleeptime = dt*timemult

        # Set pre-grasp configuration
        self.taskmanip.ReleaseFingers()
        while not self.robot.GetController().IsDone():
            time.sleep(CONTROLLERSLEEPDUR)
        self.robot.ReleaseAllGrabbed()

        self.object.SetTransform(Tobj)
        self.robot.SetActiveDOFValues(topptransfer.Eval(0))
        self.taskmanip.CloseFingers()
        while not self.robot.GetController().IsDone():
            time.sleep(CONTROLLERSLEEPDUR)
        self.robot.Grab(self.object)
        time.sleep(sleeptime)
        
        T = np.append(np.arange(dt, topptransfer.duration, dt), topptransfer.duration)
        for t in T:
            q = topptransfer.Eval(t)
            self.robot.SetActiveDOFValues(q)
            time.sleep(sleeptime)
            
