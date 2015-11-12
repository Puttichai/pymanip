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
                                     TRANSIT, TRANSFER, CG, CP, CGCP)
from ..utils import Utils
from ..utils.Utils import Colorize
from ..utils import GraspPlacementGraph as GPG

_RNG = random.SystemRandom()


############################################################
#       Grasp-Placement-Graph Manipulation Planner
############################################################
class GPGManipulationPlanner(bmp.BaseManipulationPlanner):
    
    def __init__(self, robot, manipulatorname, mobj):
        print Colorize('Initializing GPG Manipulation Planner . . .', 'yellow')
        super(GPGManipulationPlanner, self).__init__(robot, manipulatorname, mobj)

        print Colorize('Gathering information of vertices of GPG', 'yellow')
        self.graphvertices = []
        for isurface in xrange(len(self.S) - 1):
            for ibox in xrange(len(self.object.GetLinks())):
                for appdir in self.boxinfos[ibox].possibleapproachingdir[isurface]:
                    vnew = (isurface, appdir + 6*ibox)
                    self.graphvertices.append(vnew)

        print Colorize('Constructing Grasp-Placement Graph', 'yellow')
        self.GPG = GPG.GraspPlacementGraph(self.graphvertices)

        # Threshold value for removal of infeasible edges
        self._threshold = 15
        # Number of trials to try SampleCP & SampleCG
        self._ntrials = 20
        
        print Colorize('Initialization completed', 'yellow')


    def SetSampleParameters(self, xobjlim, yobjlim, thetaobjlim, zoffset):
        self.xobjlim = xobjlim
        self.yobjlim = yobjlim
        self.thetaobjlim = thetaobjlim
        self.zoffset = zoffset # height of the table
        self._setsampleparameters = True


    def InitQuery(self, cstart, cgoal):
        self.vstart = self.Vertex(cstart.Clone(), FW)
        self.vgoal = self.Vertex(cgoal.Clone(), BW)
        # If all paths of length k are infeasible, then try paths of
        # length k + self._increment
        if (self.vstart.id.count(None) == 0) or (self.vgoal.id.count(None) == 0):
            # Either start or goal (or both) is in GP
            self._increment = 1
        else:
            self._increment = 2
            
        self._hasquery = True
        self._hassolution = False


    def Plan(self, timeout):
        if not self._hasquery:
            print Colorize('No query for the planner yet', 'red')
            return False
        if not self._setsampleparameters:
            print Colorize('Sample parameters have not been set yet', 'red')
            return False

        # Shortest path length (from the GP Graph)
        p = self.GPG.FindShortestPathLength(self.vstart.id[1::], self.vgoal.id[1::])
        
        n = 0 # dummy variable keeping track of current path length
        self.timerecord = [] # keeps records of running time for each path length
        while (self._runningtime < timeout) and (not self._hassolution):
            ts_loop = time.time()
            
            self._k = p + self._increment*n # current shortest possible path length
            n += 1
            
            # Create new trees
            vstart = self.vstart.Clone()
            vstart.level = 0
            vgoal = self.vgoal.Clone()
            vgoal.level = self._k
            vgoal.id = (vgoal.level, ) + vgoal.id[1::]
            self.treestart = self.Tree(vstart, FW)
            self.treeend = self.Tree(vgoal, BW)

            ts_graph = time.time() # start graph processing
            if self._print:
                print Colorize('Extracting paths of length {0}. . .'.format(self._k), 
                               'yellow')
            # pathslist is a list of possible paths of length self._k
            self._pathslist = self.GPG.FindPathsOfLengthK(self.vstart.id[1::],
                                                          self.vgoal.id[1::],
                                                          self._k)
            
            if self._print:
                print Colorize('Creating direction graphs. . .', 'yellow')
            # We use Q for a directed graph
            self.QFW = GPG.CreateFWDirectedGraphFromPathsList(self._pathslist)
            self.QBW = GPG.CreateBWDirectedGraphFromPathsList(self._pathslist)

            # We need to find out what the first trajtype is (then we
            # know the rest). For simplicity, we assume that either
            # the start configuration or the goal configuration is in
            # P. This implies there is no ambiguity in the first
            # trajtype.
            self.firsttrajtype = None
            for vid in self.QFW[vstart.id]:
                if ((vstart.id[1] == vid[1]) and (not (vstart.id[2] == vid[2]))):
                    self.firsttrajtype = TRANSIT
                    break
                elif ((vstart.id[2] == vid[2]) and (not (vstart.id[1] == vid[1]))):
                    self.firsttrajtype = TRANSFER
                    break
                else:
                    # Not enough information
                    continue
            assert(self.firsttrajtype is not None)

            if self._print:
                print Colorize('Creating an edge database for storing statistics. . .',
                               'yellow')
            # Each entry stores the number of IK failure. High number
            # of IK failure can probably reflect the actual kinematic
            # infeasibility.
            self.ED = dict() # Edge database
            for key in self.QFW.graphdict:
                for val in self.QFW.graphdict[key]:
                    edge = key + val
                    if edge not in self.ED:
                        self.ED[edge] = 0

            te_graph = time.time() # end graph processing
            if self._print:
                print Colorize('  Graph processing time : {0} sec.'.\
                                   format(te_graph - ts_graph), 'yellow')
            
            # Compute weights for start and goal vertices. A vertex
            # with more weight is more likely to be sampled. This
            # weight is then roughly a bias thatt we put in the
            # exploration.
            self._weighteddist = [[], []]
            self._weighteddist[FW].append((self.treestart[0].index, 
                                           ComputeWeight(self.treestart[0].level)))
            self._weighteddist[BW].append((self.treeend[0].index, 
                                           ComputeWeight(self.treeend[0].level)))

            # Run the planner until all possible paths (in
            # self._pathslist) are declared infeasible.
            self._hassolution = self._Run(timeout - self.runningtime)
            te_loop = time.time()
            self.timerecord.append(te_loop - ts_loop)

        return self._hassolution
            

    def _Run(self, timeout):
        if self.hassolution:
            print 'The planner has already found a solution.'
            return True

        t = 0.0
        it = 0

        while (t < timeout):
            it += 1
            if self._print:
                print Colorize('({0}) iteration : {1}'.format(self._k, it), 'blue')
            
            ts = time.time()
                
            # Eliminate edges that are declared infeasible (according
            # to self._threshold)
            modified = False # True if some edges are removed
            infeasiblekeys = []
            infeasiblepaths = []
            for key in self.ED:
                if self.ED[key] > self._threshold:
                    l1 = key[0]
                    v1 = key[1:3]
                    l2 = key[3]
                    v2 = key[4:6]
                    for ipath in xrange(len(self._pathslist)):
                        path = self._pathslist[ipath]
                        if (path[l1] == v1) and (path[l2] == v2):
                            infeasiblepaths.append(ipath)
                            if key not in infeasiblekeys:
                                if self._print:
                                    print Colorize('  Eliminating edge {0}'.format(key),
                                                   'magenta')
                                infeasiblekeys.append(key)
                            modified = True

            # If some edges are removed, we also need to remove paths
            # which contain those edges.
            if modified:
                oldpathslist = self._pathslist
                self._pathslist = []
                for ipath in xrange(len(oldpathslist)):
                    if ipath not in infeasiblepaths:
                        self._pathslist.append(oldpathslist[ipath])
                if len(self._pathslist) == 0:
                    # All paths have been eliminated
                    break
                else:
                    for key in infeasiblekeys:
                        self.ED.pop(key)
                    self.QFW = GPG.CreateFWDirectedGraphFromPathsList(self._pathslist)
                    self.QBW = GPG.CreateBWDirectedGraphFromPathsList(self._pathslist)

            treedirection = np.mod(it - 1, 2)
            # Sample a vertex on the tree 'treedirection' to be
            # extended from
            index = self.SampleTree(treedirection)

            if treedirection == FW:
                # Forward extension
                status = self.ExtendFWFrom(index)
                if not (status == TRAPPED):
                    if (status == REACHED) or (self.ConnectFW() == REACHED):
                        te = time.time()
                        t += te - ts
                        self._runningtime += t
                        self._iterations += it
                        self._hassolution = True
                        if self._print:
                            print Colorize('Path found', 'green')
                            print Colorize('    Total number of iterations : {0}'.\
                                               format(self.iterations), 'green')
                            print Colorize('    Total running time : {0} sec.'.\
                                               format(self.runningtime), 'green')
                            
                        return True
            else:
                # Backward Extension
                status = self.ExtendBWFrom(index)
                if not (status == TRAPPED):
                    if (status == REACHED) or (self.ConnectBW() == REACHED):
                        te = time.time()
                        t += te - ts
                        self._runningtime += t
                        self._iterations += it
                        self._hassolution = True
                        if self._print:
                            print Colorize('Path found', 'green')
                            print Colorize('    Total number of iterations : {0}'.\
                                               format(self.iterations),
                                           'green')
                            print Colorize('    Total running time : {0} sec.'.\
                                               format(self.runningtime),
                                           'green')
                            
                        return True
                    
            te = time.time()
            t += te - ts
            
        # End while loop
        if len(self._pathslist) == 0:
            print Colorize('All possible paths of length {0} are declared infeasible'.\
                               format(self._k), 'red')
        else:
            print Colorize('Allotted time ({0} sec.) is exhausted after {1} iterations'.\
                               format(timeout, it), 'red')
        self._runningtime += t
        self._iterations += it
        return False
    
            
    def SampleTree(self, treedirection):
        return WeightedChoice(self._weighteddist[treedirection])

    
    def ExtendFWFrom(self, vindex):
        if self._print:
            print '  [ExtendFWFrom] vindex = {0}'.format(vindex)
        status = TRAPPED

        vnear = self.treestart[vindex]
        prevtype = vnear.trajtype

        # Look into the dictionary for possible extensions
        try:
            possibleextensions = self.QFW.graphdict[vnear.id]
        except KeyError:
            # Possible paths have been removed (due to kinematic infeasibility)
            oldweighteddist = self._weighteddist[FW]
            if self._print:
                print '    Sampled vertex is not on any feasible path'
                print '    Removeing vindex = {0} from treestart'.format(vindex)
            self._weighteddist[FW] = [d for d in oldweighteddist if d[0] != vindex]
            return status

        # Sample an id of a vertex to extend to
        vnext_id = _RNG.choice(possibleextensions)
        # Unpack the id
        [nextlevel, isurface, approachingdir] = vnext_id

        # The cases when vnear is either 1 step or 2 steps away from
        # the goal have to be handled separately.

        # Check whether vnear is one step away from vgoal
        onesteptogoal = (vnext_id == self.treeend[0].id)
        # Check whether vnear is two steps away from vgoal
        twostepstogoal = False
        nexttwosteps = self.QFW.graphdict[vnext_id] # set of two-step-away vertices
        if (len(nexttwosteps) == 1) and (nexttwosteps[0] == self.treeend[0].id):
            twostepstogoal = True

        if self._print:
            print '  [ExtendFWFrom] vnear = {0}; vnext = {1}'.\
            format(vnear.id, vnext_id)

        # Check whether this upcoming edge is TRANSIT or TRANSFER
        # If any of these condition is satisfied, the next trajtyps is TRANSIT
        cond1 = (prevtype == TRANSFER)
        cond2 = ((vnear.id[1] == vnext_id[1]) and (not (vnear.id[2] == vnext_id[2])))
        cond3 = ((self.firsttrajtype == TRANSIT) and (np.mod(vnear.level, 2) == 0))
        cond4 = ((self.firsttrajtype == TRANSFER) and (np.mod(vnear.level, 2) == 1))

        if cond1 or cond2 or cond3 or cond4:
            if self._print:
                print '  TRANSIT extension'
            # TRANSIT
            if onesteptogoal:
                [plannerstatus, trajectorystring] = self._PlanTransitTrajectory\
                (vnear, self.treeend[0])
                if not (plannerstatus == 1):
                    return status
                # Now we have reached the goal. 
                status = REACHED
                
                # Stack the vertex vindex on top of treestart
                self.treestart.verticeslist.append(self.treestart[vindex])
                # Stack the vertex goal on top of treeend
                self.treeend.verticeslist.append(self.treeend[0])

                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSIT
                return status
                
            if twostepstogoal:
                # We need to transit to the same grasp as in vertex goal
                qgrasp = self.treeend[0].config.qgrasp
                Tobj = vnear.config.tobj
                qref = self.treeend[0].config.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    if self._print:
                        print '    No IK solution'
                    self.ED[vnear.id + vnext_id] += 1
                    return status
            else:
                # Usual TRANSIT extension
                # Sample a new grasp
                [passed, sol, qgrasp] = self.SampleCG(vnear.config.tobj, isurface,
                                                      approachingdir)
                if not passed:
                    # Infeasible grasp is probably because it is
                    # kinematically unreachable
                    self.ED[vnear.id + vnext_id] += 1
                    return status
                
            cnew = self.Config(sol, vnear.config.tobj, qgrasp, CGCP, isurface)
            vnew = self.Vertex(cnew, FW, level=nextlevel)
            [plannerstatus, trajectorystring] = self._PlanTransitTrajectory\
            (vnear, vnew)

            if not (plannerstatus == 1):
                return status

            # Now successfully extended
            self.treestart.AddVertex(vindex, vnew, trajectorystring, TRANSIT)
            self._weighteddist[FW].append((self.treestart[-1].index, 
                                           ComputeWeight(self.treestart[-1].level)))
            status = ADVANCED
            return status

        else:
            if self._print:
                print '  TRANSFER extension'
            # TRANSFER
            if onesteptogoal:
                [plannerstatus, trajectorystring] = self._PlanTransferTrajectory\
                (vnear, self.treeend[0])
                if not (plannerstatus == 1):
                    return status
                # Now we have reached the goal. 
                status = REACHED
                
                # Stack the vertex vindex on top of treestart
                self.treestart.verticeslist.append(self.treestart[vindex])
                # Stack the vertex goal on top of treeend
                self.treeend.verticeslist.append(self.treeend[0])

                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSFER
                return status
            
            if twostepstogoal:# and (not (self.treeend[0].config.type == CP)):
                # We need to transfer the object to the same Tobj as in vertex goal
                Tobj = self.treeend[0].config.tobj
                qgrasp = vnear.config.qgrasp
                qref = vnear.config.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    if self._print:
                        print '    No IK solution'
                    self.ED[vnear.id + vnext_id] += 1
                    return status
            else:
                # Usual TRANSFER extension
                # Sample new placement
                [passed, sol, tobj] = self.SampleCP(isurface,
                                                    vnear.config.qgrasp,
                                                    vnear.config.qrobot)
                if not passed:
                    # Infeasible placement is probably because it is
                    # kinematically unreachable
                    self.ED[vnear.id + vnext_id] += 1
                    return status
                
            cnew = self.Config(sol, tobj, vnear.config.qgrasp, CGCP, isurface)
            vnew = self.Vertex(cnew, FW, level=nextlevel)
            [plannerstatus, trajectorystring] = self._PlanTransferTrajectory\
            (vnear, vnew)

            if not (plannerstatus == 1):
                return status

            # Now successfully extended
            self.treestart.AddVertex(vindex, vnew, trajectorystring, TRANSFER)
            self._weighteddist[FW].append((self.treestart[-1].index, 
                                           ComputeWeight(self.treestart[-1].level)))
            status = ADVANCED
            return status
        

    def ExtendBWFrom(self, vindex):
        if self._print:
            print '  [ExtendBWFrom] vindex = {0}'.format(vindex)
        status = TRAPPED

        vnear = self.treeend[vindex]
        prevtype = vnear.trajtype

        # Look into the dictionary for possible extensions
        try:
            possibleextensions = self.QBW.graphdict[vnear.id]
        except KeyError:
            # Possible paths have been removed (due to kinematic infeasibility)
            oldweighteddist = self._weighteddist[BW]
            if self._print:
                print '    Sampled vertex is not on any feasible path'
                print '    Removeing vindex = {0} from treeend'.format(vindex)
            self._weighteddist[BW] = [d for d in oldweighteddist if d[0] != vindex]
            return status

        # Sample an id of a vertex to extend to
        vprev_id = _RNG.choice(possibleextensions)
        # Unpack the id
        [prevlevel, isurface, approachingdir] = vprev_id

        # The cases when vnear is either 1 step or 2 steps away from
        # the start have to be handled separately.

        # Check whether vnear is one step away from vinit
        onestepfromstart = (vprev_id == self.treestart[0].id)
        # Check whether vnear is two steps away from vinit
        twostepsfromstart = False
        prevtwosteps = self.QBW.graphdict[vprev_id] # set of two-step-away vertices
        if (len(prevtwosteps) == 1) and (prevtwosteps[0] == self.treestart[0].id):
            twostepsfromstart = True

        if self._print:
            print '  [ExtendBWFrom] vprev = {0}; vnear = {1}'.\
            format(vprev_id, vnear.id)
            
        # Check whether this upcoming edge is TRANSIT or TRANSFER
        # If any of these condition is satisfied, the next trajtyps is TRANSIT
        cond1 = (prevtype == TRANSFER)
        cond2 = ((vprev_id[1] == vnear.id[1]) and (not (vprev_id[2] == vnear.id[2])))
        cond3 = ((self.firsttrajtype == TRANSIT) and (np.mod(vnear.level, 2) == 1))
        cond4 = ((self.firsttrajtype == TRANSFER) and (np.mod(vnear.level, 2) == 0))

        if cond1 or cond2 or cond3 or cond4:
            if self._print:
                print '  TRANSIT extension'
            # TRANSIT
            if onestepfromstart:
                [plannerstatus, trajectorystring] = self._PlanTransitTrajectory\
                (self.treestart[0], vnear)
                if not (plannerstatus == 1):
                    return status
                # Now we have reached the start. 
                status = REACHED
                
                # Stack the vertex start on top of treestart
                self.treestart.verticeslist.append(self.treestart[0])
                # Stack the vertex vindex on top of treeend
                self.treeend.verticeslist.append(self.treeend[vindex])

                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSIT
                return status
                
            if twostepsfromstart:# and (not (self.treestart[0].config.type == CP)):
                # We need to transit to the same grasp as in vertex start
                qgrasp = self.treestart[0].config.qgrasp
                Tobj = vnear.config.tobj
                qref = self.treestart[0].config.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    if self._print:
                        print '    No IK solution'
                    self.ED[vprev_id + vnear.id] += 1
                    return status
            else:
                # Usual TRANSIT extension
                # Sample a new grasp
                [passed, sol, qgrasp] = self.SampleCG(vnear.config.tobj, isurface,
                                                      approachingdir)
                if not passed:
                    # Infeasible grasp is probably because it is
                    # kinematically unreachable
                    self.ED[vprev_id + vnear.id] += 1
                    return status
                
            cnew = self.Config(sol, vnear.config.tobj, qgrasp, CGCP, isurface)
            vnew = self.Vertex(cnew, BW, level=prevlevel)
            [plannerstatus, trajectorystring] = self._PlanTransitTrajectory\
            (vnew, vnear)

            if not (plannerstatus == 1):
                return status

            # Now successfully extended
            self.treeend.AddVertex(vindex, vnew, trajectorystring, TRANSIT)
            self._weighteddist[BW].append((self.treeend[-1].index, 
                                           ComputeWeight(self.treeend[-1].level)))
            status = ADVANCED
            return status

        else:
            if self._print:
                print '  TRANSFER extension'
            # TRANSFER
            if onestepfromstart:
                [plannerstatus, trajectorystring] = self._PlanTransferTrajectory\
                (self.treestart[0], vnear)
                if not (plannerstatus == 1):
                    return status
                # Now we have reached the goal. 
                status = REACHED
                
                # Stack the vertex start on top of treestart
                self.treestart.verticeslist.append(self.treestart[0])
                # Stack the vertex vindex on top of treeend
                self.treeend.verticeslist.append(self.treeend[vindex])

                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSFER
                return status
            
            if twostepsfromstart:
                # We need to transfer the object to the same Tobj as in vertex goal
                Tobj = self.treestart[0].config.tobj
                qgrasp = vnear.config.qgrasp
                qref = self.treestart[0].config.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    if self._print:
                        print '    No IK solution'
                    self.ED[vprev_id + vnear.id] += 1
                    return status
            else:
                # Usual TRANSFER extension
                # Sample new placement
                [passed, sol, tobj] = self.SampleCP(isurface,
                                                    vnear.config.qgrasp,
                                                    vnear.config.qrobot)
                if not passed:
                    # Infeasible placement is probably because it is
                    # kinematically unreachable
                    self.ED[vprev_id + vnear.id] += 1
                    return status
                
            cnew = self.Config(sol, tobj, vnear.config.qgrasp, CGCP, isurface)
            vnew = self.Vertex(cnew, BW, level=prevlevel)
            [plannerstatus, trajectorystring] = self._PlanTransferTrajectory\
            (vnew, vnear)

            if not (plannerstatus == 1):
                return status

            # Now successfully extended
            self.treeend.AddVertex(vindex, vnew, trajectorystring, TRANSFER)
            self._weighteddist[BW].append((self.treeend[-1].index, 
                                           ComputeWeight(self.treeend[-1].level)))
            status = ADVANCED
            return status


    def ConnectFW(self):
        if self._print:
            print '  [ConnectFW]'
        status = TRAPPED

        vfw = self.treestart[-1] # newly added vertex
        cfw = vfw.config
        curtype = vfw.trajtype
        onestepfromvfwlist = self.QFW.graphdict[vfw.id]
            
        # Check wheter the newly added vertex on treestart is one step
        # away from goal
        cond1 = (len(onestepfromvfwlist) == 1)
        cond2 = (onestepfromvfwlist[0] == self.treeend[0].id)
        if cond1 and cond2:
            # The next vertex is vertex goal
            if curtype == TRANSIT:
                # TRANSFER
                [plannerstatus, trajectorystring] = self._PlanTransferTrajectory\
                (vfw, self.treeend[0])
                if not (plannerstatus == 1):
                    return status

                # Now we have reached the goal. 
                status = REACHED
                
                # Stack the vertex vfw on top of treestart
                self.treestart.verticeslist.append(self.treestart[vfw.index])
                # Stack the vertex goal on top of treeend
                self.treeend.verticeslist.append(self.treeend[0])

                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSFER
                return status
            else:
                # TRANSIT
                [plannerstatus, trajectorystring] = self._PlanTransitTrajectory\
                (vfw, self.treeend[0])
                if not (plannerstatus == 1):
                    return status

                # Now we have reached the goal. 
                status = REACHED
                
                # Stack the vertex vfw on top of treestart
                self.treestart.verticeslist.append(self.treestart[vfw.index])
                # Stack the vertex goal on top of treeend
                self.treeend.verticeslist.append(self.treeend[0])

                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSIT
                return status
            
        # Create a list of vertices which are two steps away from vfw
        temp = []
        for v_id in onestepfromvfwlist:
            temp += self.QFW.graphdict[v_id]
        twostepsfromvfwlist = []
        # Append only non-redundant vertices to twostepsfromvfwlist
        for v_id in temp:
            if v_id not in twostepsfromvfwlist:
                twostepsfromvfwlist.append(v_id)

        # Now examine all existing vertices on treeend which are two
        # steps away from vfw
        idslist = [v.id for v in self.treeend]
        for v_id in twostepsfromvfwlist:
            if v_id not in idslist:
                continue

            vbw_index = idslist.index(v_id)
            vbw = self.treeend[vbw_index]
            cbw = vbw.config
            
            # Now execute a two-step extension
            if curtype == TRANSIT:
                # TRANSFER --> TRANSIT
                # TRANSFER
                qgrasp = cfw.qgrasp
                Tobj = cbw.tobj
                qref = cfw.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    print '    No IK solution'
                    vnext_id = (vfw.level + 1, cbw.isurface, cfw.qgrasp[1])
                    self.ED[vfw.id + vnext_id] += 1
                    continue

                cnew = self.Config(sol, Tobj, qgrasp, CGCP, cbw.isurface)
                vnew = self.Vertex(cnew, FW, level=vfw.level + 1)
                [plannerstatus1, trajectorystring1] = self._PlanTransferTrajectory\
                (vfw, vnew)
                if (not (plannerstatus1 == 1)):
                    continue

                # Now successfully extended
                self.treestart.AddVertex(vfw.index, vnew, trajectorystring1, TRANSFER)
                self._weighteddist[FW].append((self.treestart[-1].index, 
                                               ComputeWeight(self.treestart[-1].level)))
                
                # TRANSIT
                [plannerstatus2, trajectorystring2] = self._PlanTransitTrajectory\
                (vnew, vbw)
                if (not (plannerstatus2 == 1)):
                    continue
                
                # Now successfully connected
                self.treeend.verticeslist.append(vbw)
                self.connectingtrajectorystring = trajectorystring2
                self.connectingtrajectorytype = TRANSIT
                status = REACHED
                return status
            
            else:
                # TRANSIT --> TRANSFER
                # TRANSIT
                qgrasp = cbw.qgrasp
                Tobj = cfw.tobj
                qref = cbw.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    print '    No IK solution'
                    vnext_id = (vfw.level + 1, cfw.isurface, cbw.qgrasp[1])
                    self.ED[vfw.id + vnext_id] += 1
                    continue
                
                cnew = self.Config(sol, Tobj, qgrasp, CGCP, cfw.isurface)
                vnew = self.Vertex(cnew, FW, level=vfw.level + 1)
                [plannerstatus1, trajectorystring1] = self._PlanTransitTrajectory\
                (vfw, vnew)
                if not (plannerstatus1 == 1):
                    continue
                
                # Now successfully extended
                self.treestart.AddVertex(vfw.index, vnew, trajectorystring1, TRANSIT)
                self._weighteddist[FW].append((self.treestart[-1].index, 
                                               ComputeWeight(self.treestart[-1].level)))

                # TRANSFER
                [plannerstatus2, trajectorystring2] = self._PlanTransferTrajectory\
                (vnew, vbw)
                if not (plannerstatus2 == 1):
                    continue
                
                # Now successfully extended
                self.treeend.verticeslist.append(vbw)
                self.connectingtrajectorystring = trajectorystring2
                self.connectingtrajectorytype = TRANSFER
                status = REACHED
                return status                
            
        # TRAPPED
        return status
    

    def ConnectBW(self):
        if self._print:
            print '  [ConnectBW]'
        status = TRAPPED
        
        vbw = self.treeend[-1] ## newly added vertex
        cbw = vbw.config
        curtype = vbw.trajtype
        onesteptovbwlist = self.QBW.graphdict[vbw.id]
        
        # Check whether the newly added vertex on treeend is one step
        # away from start
        cond1 = (len(onesteptovbwlist) == 1)
        cond2 = (onesteptovbwlist[0] == self.treestart[0].id)
        if cond1 and cond2:
            # The prev vertex is vertex start
            if curtype == TRANSIT:
                # TRANSFER
                [plannerstatus, trajectorystring] = self._PlanTransferTrajectory\
                (self.treestart[0], vbw)
                if not (plannerstatus == 1):
                    return status
                
                # Now we have reached the goal
                status = REACHED
                
                # Stack the vertex start on top of treestart
                self.treestart.verticeslist.append(self.treestart[0])
                # Stack the vertex vbw on top of treeend
                self.treeend.verticeslist.append(self.treeend[vbw.index])
                
                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSFER
                return status
            else:
                # TRANSIT
                [plannerstatus, trajectorystring] = self._PlanTransitTrajectory\
                (self.treestart[0], vbw)
                if (not (plannerstatus == 1)):
                    return status
                
                # Now we have reached the goal
                status = REACHED
                
                # Stack the vertex start on the top of treestart
                self.treestart.verticeslist.append(self.treestart[0])
                # Stack the vertex vbw on the top of treeend
                self.treeend.verticeslist.append(self.treeend[vbw.index])
                
                self.connectingtrajectorystring = trajectorystring
                self.connectingtrajectorytype = TRANSIT
                return status

        # Create a list of vertices which are two steps before vbw
        temp = []
        for v_id in onesteptovbwlist:
            temp += self.QBW.graphdict[v_id]
        twostepstovbwlist = []
        # Append only non-redundant vertices to twostepsfromvbwlist
        for v_id in temp:
            if v_id not in twostepstovbwlist:
                twostepstovbwlist.append(v_id)

        # Now examine all existing vertices on treestart which are two
        # steps away from vbw
        idslist = [v.id for v in self.treestart]
        for v_id in twostepstovbwlist:
            if v_id not in idslist:
                continue
            
            vfw_index = idslist.index(v_id)
            vfw = self.treestart[vfw_index]
            cfw = vfw.config
            
            # Now execute a two-step extension
            if curtype == TRANSIT:
                # TRANSIT --> TRANSFER
                # TRANFER
                qgrasp = cbw.qgrasp
                Tobj = cfw.tobj
                qref = cbw.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    print "    No IK solution"
                    vprev_id = (vbw.level - 1, cfw.isurface, cbw.qgrasp[1])
                    self.ED[vprev_id + vbw.id] += 1
                    continue

                cnew = self.Config(sol, Tobj, qgrasp, CGCP, cfw.isurface)
                vnew = self.Vertex(cnew, BW, level=vbw.level - 1)
                [plannerstatus1, trajectorystring1] = self._PlanTransferTrajectory\
                (vnew, vbw)
                if not (plannerstatus1 == 1):
                    continue
                
                # Now successfully extended
                self.treeend.AddVertex(vbw.index, vnew, trajectorystring1, TRANSFER)
                self._weighteddist[BW].append((self.treeend[-1].index, 
                                               ComputeWeight(self.treeend[-1].level)))

                # TRANSIT
                [plannerstatus2, trajectorystring2] = self._PlanTransitTrajectory\
                (vfw, vnew)
                if not (plannerstatus2 == 1):
                    continue
                
                # Now successfully extended
                self.treestart.verticeslist.append(vfw)
                self.connectingtrajectorystring = trajectorystring2
                self.connectingtrajectorytype = TRANSIT
                status = REACHED
                return status
            
            else:
                # TRANSFER --> TRANSIT
                # TRANSIT
                qgrasp = cfw.qgrasp
                Tobj = cbw.tobj
                qref = cfw.qrobot
                sol = self.ComputeGraspConfiguration(qgrasp, Tobj, qrobot_ref=qref)
                if sol is None:
                    print "    No IK solution"
                    vprev_id = (vbw.level - 1, cbw.isurface, cfw.qgrasp[1])
                    self.ED[vprev_id + vbw.id] += 1
                    continue
                
                cnew = self.Config(sol, Tobj, qgrasp, CGCP, cbw.isurface)
                vnew = self.Vertex(cnew, BW, level=vbw.level - 1)
                [plannerstatus1, trajectorystring1] = self._PlanTransitTrajectory\
                (vnew, vbw)
                if not (plannerstatus1 == 1):
                    continue
                
                # Now successfully extended
                self.treeend.AddVertex(vbw.index, vnew, trajectorystring1, TRANSIT)
                self._weighteddist[BW].append((self.treeend[-1].index, 
                                               ComputeWeight(self.treeend[-1].level)))
                
                # TRANSFER
                [plannerstatus2, trajectorystring2] = self._PlanTransferTrajectory\
                (vfw, vnew)
                if not (plannerstatus2 == 1):
                    continue
                
                # Now successfully connected
                self.treestart.verticeslist.append(vfw)
                self.connectingtrajectorystring = trajectorystring2
                self.connectingtrajectorytype = TRANSFER
                status = REACHED
                return status

        ## TRAPPED
        return status


    def SampleCP(self, isurface, qgrasp, qrobot_ref):
        """
        SampleCP samples a configuration in CP.
        
        Returns
        -------
        passed : bool
        sol : n-vector
            A configuration of a robot grasping the object with grasp
        Tobj : 4x4 transformation matrix
        """
        passed = False
        ibox = qgrasp[0]
        graspedlink = self.object.GetLinks()[ibox]
        extents = self.boxinfos[ibox].extents
        
        for _ in xrange(self._ntrials):
            xobj = _RNG.uniform(self.xobjlim[0], self.xobjlim[1])
            yobj = _RNG.uniform(self.yobjlim[0], self.yobjlim[1])
            thetaobj = _RNG.uniform(self.thetaobjlim[0], self.thetaobjlim[1])

            qobj = [xobj, yobj, thetaobj, isurface]
            Tobj = Utils.ComputeTObject(qobj, self.S, self.zoffset)
            
            if False:
                # TODO: add external implementation of placement checking
                if not self.CheckPlacement(self.object, Tobj, isurface, self.S):
                    continue
             
            self.object.SetTransform(Tobj)
            self.robot.SetDOFValues(np.zeros(7)) # HOME (with gripper open)
            objincollision = self.env.CheckCollision(self.object)
            if objincollision:
                continue
            
            # Compute IK solution for the robot
            Tgripper = Utils.ComputeTGripper(graspedlink.GetTransform(),
                                             qgrasp, extents, unitscale=False)
            with self.robot:
                self.robot.SetActiveDOFValues(qrobot_ref)
                sol = self.manip.FindIKSolution\
                (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)
            if sol is None:
                continue

            if not self.CheckGrasp(sol, Tobj):
                continue

            passed = True
            break
        
        if passed:
            return [passed, sol, Tobj]
        else:
            return [False, None, None]


    def SampleCG(self, Tobj, isurface, approachingdir):
        """
        SampleCG samples a configuration in CG.
        
        Returns
        -------
        passed : bool
        sol : n-vector
            A configuration of a robot grasping the object with grasp
        qgrasp : list
        """
        passed = False
        ibox = int(approachingdir)/6
        realapproachingdir = np.mod(approachingdir, 6)
        boxinfo = self.boxinfos[ibox]
        graspedlink = self.object.GetLinks()[ibox]
        
        self.object.SetTransform(Tobj)
        for _ in xrange(self._ntrials):
            # Sample a sliding direction
            slidingdir = _RNG.choice(boxinfo.possibleslidingdir[realapproachingdir])
            # Resample a sliding direction if necessary
            while (isurface, realapproachingdir, slidingdir) not in boxinfo.intervals:
                slidingdir = _RNG.choice(boxinfo.possibleslidingdir[realapproachingdir])

            # Sample a value for delta (where to grasp along the sliding direction)
            delta = WeightedChoice2(boxinfo.intervals[isurface, realapproachingdir, 
                                                      slidingdir])

            # Assign qgrasp
            qgrasp = [ibox, approachingdir, slidingdir, delta]
            
            # Compute Tgripper
            Tgripper = Utils.ComputeTGripper(graspedlink.GetTransform(),
                                             qgrasp, boxinfo.extents, unitscale=False)
            
            # Compute IK solution for the robot
            sol = self.manip.FindIKSolution\
            (Tgripper, orpy.IkFilterOptions.CheckEnvCollisions)
            if sol is None:
                continue

            if not self.CheckGrasp(sol, Tobj):
                continue

            passed = True
            break
        
        if passed:
            self.robot.SetActiveDOFValues(sol)
            self.object.SetTransform(Tobj)
            return [passed, sol, qgrasp]
        else:
            return [False, None, None]
        
        
    def PlotManipulationGraph(self, fignum=1, grid=False, verticesoffset=1):
        offset = verticesoffset # grasp and placement classes start from offset
        V = [v for v in self.GPG.graphdict if v[1] is not None]
        # Find the numbers of grasp and placement classes
        vgraspmax = max(V, key=lambda p: p[1])
        vgraspmin = min(V, key=lambda p: p[1])
        ngraspclasses = (vgraspmax[1] - vgraspmin[1]) + 1
        vplacementmax = max(V, key=lambda p: p[0])
        vplacementmin = min(V, key=lambda p: p[0])
        nplacementclasses = (vplacementmax[0] - vplacementmin[0]) + 1

        if grid:
            for i in xrange(nplacementclasses):
                plt.plot([i + offset, i + offset], 
                         [offset, offset + ngraspclasses - 1], 'k:')
            for i in xrange(ngraspclasses):
                plt.plot([offset, offset + nplacementclasses - 1],
                         [i + offset, i + offset], 'k:')

        # Now find the extreme vertices of each grasp and placement class
        extremevertices_placementclass = \
        np.hstack([np.ones((nplacementclasses, 1))*ngraspclasses,
                   np.zeros((nplacementclasses, 1))])
        extremevertices_graspclass = \
        np.hstack([np.ones((ngraspclasses, 1))*nplacementclasses,
                   np.zeros((ngraspclasses, 1))])

        for v in V:
            # Iterate over all vertices
            if v[1] < extremevertices_placementclass[v[0]][0]:
                extremevertices_placementclass[v[0], 0] = v[1]
            if v[1] > extremevertices_placementclass[v[0]][1]:
                extremevertices_placementclass[v[0], 1] = v[1]
                
            if v[0] < extremevertices_graspclass[v[1]][0]:
                extremevertices_graspclass[v[1], 0] = v[0]
            if v[0] > extremevertices_graspclass[v[1]][1]:
                extremevertices_graspclass[v[1], 1] = v[0]

        plt.figure(fignum)
        plt.hold(True)
        # Plot transit trajectories
        for i in xrange(nplacementclasses):
            mingraspclass = extremevertices_placementclass[i][0]
            maxgraspclass = extremevertices_placementclass[i][1]
            plt.plot([i + offset, i + offset], 
                     [mingraspclass + offset, maxgraspclass + offset], 'g-')

        # Plot transfer trajectories
        for i in xrange(ngraspclasses):
            minplacementclass = extremevertices_graspclass[i][0]
            maxplacementclass = extremevertices_graspclass[i][1]
            plt.plot([minplacementclass + offset, maxplacementclass + offset], 
                     [i + offset, i + offset], 'g-')
        
        # Plot vertices
        X = np.array([v[0] for v in V]) + offset
        Y = np.array([v[1] for v in V]) + offset
        plt.plot(X, Y, 'bo')
        plt.axis([-1 + offset, max(X) + 1, -1 + offset, max(Y) + 1])
        plt.tight_layout()


############################################################
#                       Utilities
############################################################
def ComputeWeight(l):
    """
    ComputeWeight computes a weight for an l^th-level vertex 
    """
    w = 2**(0.5*l)
    return w


def WeightedChoice(choices):
    total = sum(w for c, w in choices)
    r = random.uniform(0, total)
    upto = 0
    for c, w in choices:
        upto += w
        if upto >= r:
            return c
    assert False


def WeightedChoice2(choices):
    total = sum(w for c, w in choices)
    s = random.uniform(0, total)
    upto = 0
    for c, w in choices:
        upto += w
        if upto >= s:
            return c + (upto - s)
    assert False
