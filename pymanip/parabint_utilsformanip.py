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

from parabint import *
import random

MAX_REPEAT_SAMPLING = 100
SHORTCUT_THRESHOLD = EPSILON

## OpenRAVE dependent. not available in usual parabint library
def ConvertOpenRAVETrajToRampsListND(openravetraj, vmvect, amvect, delta = 0):
    ndof = openravetraj.GetConfigurationSpecification().\
    GetGroupFromName('joint_values').dof
    nwaypoints = openravetraj.GetNumWaypoints()
    wp0 = openravetraj.GetWaypoints(0, nwaypoints) ## 1D array of waypoints
    wp1 = []
    for i in range(nwaypoints):
        wp1.append(wp0[ndof*i:ndof*i + ndof])
    rampslist = ConvertWaypointsToRampsListND(wp1, vmvect, amvect, delta)
    return rampslist


##
def ConvertWaypointsToRampsListND(waypointslist, vmvect, amvect, delta = 0):
    ## merge collinear waypoints
    W = MergeWaypoints(waypointslist)
    nwaypoints = len(W)
    newrampslistnd = RampsListND()
    
    for i in range(nwaypoints - 1):
        q0 = W[i]
        q1 = W[i + 1]
        rampslistnd = InterpolateZeroVelND(q0, q1, vmvect, amvect, delta)
        newrampslistnd.Append(rampslistnd)
    return newrampslistnd


##
def MergeWaypoints(waypointslist):
    nwaypoints = len(waypointslist)
    newwaypointslist = []
    newwaypointslist.append(waypointslist[0])
    
    for i in range(1, nwaypoints):
        if len(newwaypointslist) >= 2:
            qprev1 = newwaypointslist[-1]
            qprev2 = newwaypointslist[-2]
            qcur = waypointslist[i]
            
            dq1 = qcur - qprev1
            dq2 = qcur - qprev2
            
            len_dq1sq = np.dot(dq1, dq1)
            len_dq2sq = np.dot(dq2, dq2)
            
            dotproduct = np.dot(dq1, dq2)
            
            if (abs(dotproduct**2 - len_dq1sq*len_dq2sq) < EPSILON):
                ## new waypoint is collinear with the previous one
                newwaypointslist.pop()
            
            newwaypointslist.append(qcur)                
        else:
            if (np.linalg.norm(waypointslist[i] - newwaypointslist[0]) > EPSILON):
                newwaypointslist.append(waypointslist[i])            
    return newwaypointslist


## ReplaceRampsListNDSegment replaces a segment from t = t0 to t = t1
## of originalrampslistnd with newsegment
def ReplaceRampsListNDSegment(originalrampslistnd, newsegment, t0, t1):
    assert(originalrampslistnd.ndof == newsegment.ndof)
    assert(t1 > t0)
    
    ndof = originalrampslistnd.ndof
    rampslistnd = []
    for j in range(ndof):
        ## replace each dof one by one
        newrampslist = RampsList()
        rampslist = originalrampslistnd[j]  
        i0, rem0 = rampslist.FindRampIndex(t0)
        i1, rem1 = rampslist.FindRampIndex(t1)
        
        ## check if t0 falls in the first ramp. 
        ## if not, insert ramp 0 to ramp i0 - 1 into newrampslist
        if i0 > 0:
            newrampslist.Append(RampsList(rampslist[0: i0]))
                
        ## remainder ramp 0
        if (abs(rem0) >= EPSILON):
            remramp0 = Ramp(rampslist[i0].v, rampslist[i0].a, rem0)
            ## initial condition has to be set because sometimes remramp0 is
            ## the beginning of the shortcut traj
            remramp0.x0 = rampslist[i0].x0
            newrampslist.Append(RampsList([remramp0]))            
        
        ## insert newsegment
        newrampslist.Append(newsegment[j])
        
        ## remainder ramp 1
        if (abs(rampslist[i1].T - rem1) >= EPSILON):
            remramp1 = Ramp(rampslist[i1].Evald(rem1), rampslist[i1].a, rampslist[i1].T - rem1)
            newrampslist.Append(RampsList([remramp1]))

        ## insert remaining ramps
        if i1 < len(rampslist) - 1:
            newrampslist.Append(RampsList(rampslist[i1 + 1: len(rampslist)]))
            
        rampslistnd.append(newrampslist)

    return RampsListND(rampslistnd)


##
def Shortcut(rampslistnd, vmvect, amvect, delta, shortcutiter, PRINT = True, robot = None, PLOT = False):
    random_number_generator = random.SystemRandom()
    ndof = rampslistnd.ndof
    nsuccessfulshortcut = 0

    if PLOT:
        from pylab import ion
        ion()
        fig = plt.figure()
        ax = fig.gca()
        ax.axis([0, shortcutiter, 0, rampslistnd.duration + 1])
        
    for it in range(shortcutiter):
        if PRINT:
            print "\n**********\niteration {0} \n**********".format(it + 1)
        dur = rampslistnd.duration
        
        ## sample two random time instants
        T = 0
        minallowedduration = max(0.04, 5*delta)
        for i in range(MAX_REPEAT_SAMPLING):
            t0 = random_number_generator.uniform(0, dur - minallowedduration)
            ## make sure that t1 - t0 > minallowedduration
            t1 = random_number_generator.uniform(t0 + minallowedduration, dur)
            T = t1 - t0
            if (T >= minallowedduration):
                break
        ## in case the above sampling fails
        if (T < minallowedduration):
            if (t0 + minallowedduration < dur):
                t1 = t0 + minallowedduration
            elif (t1 - minallowedduration > 0.0):
                t0 = t1 - minallowedduration
            else:
                t0 = 0.0
                t1 = dur
                
        ## constrain t0 and t1 not to violate the minimum switching time constraint
        i0 = bisect.bisect_left(rampslistnd.switchpointslist, t0) - 1
        i1 = bisect.bisect_left(rampslistnd.switchpointslist, t1) - 1
        ## snap t0 and t1 to the nearest switching points (if neccessary)
        ## snapping t0
        if (t0 - rampslistnd.switchpointslist[i0] < max(0.008, delta)):
            t0 = rampslistnd.switchpointslist[i0]
        ## snapping t1
        if (rampslistnd.switchpointslist[i1 + 1] - t1 < max(0.008, delta)):
            t1 = rampslistnd.switchpointslist[i1 + 1]

        T = t1 - t0

        q0 = rampslistnd.Eval(t0)
        qd0 = rampslistnd.Evald(t0)
        q1 = rampslistnd.Eval(t1)
        qd1 = rampslistnd.Evald(t1)        
        newrampslistnd = InterpolateArbitraryVelND(q0, q1, qd0, qd1, vmvect, amvect, delta)
        if not (newrampslistnd.isvoid):
            ## check the new duration
            Tnew = newrampslistnd.duration
            if (T - Tnew > SHORTCUT_THRESHOLD):
                ## check joint limits
                injointlimits = CheckJointLimits(robot, newrampslistnd)
                if (injointlimits):
                    ## check collision
                    incollision = CheckCollision(robot, newrampslistnd)
                    if (not incollision):
                        nsuccessfulshortcut += 1
                        if PRINT:
                            print "\tSuccessful Shortcut"
                        rampslistnd = ReplaceRampsListNDSegment(rampslistnd, newrampslistnd, t0, t1)
                        if PLOT:
                            ax.plot([it], [dur - (T - Tnew)], 'ro')
                            plt.draw()
                    else:
                        if PRINT:
                            print "\tIn Collision"
                        if PLOT:
                            ax.plot([it], [dur], 'ro')
                            plt.draw()
                else:
                    if PRINT:
                        print "\tNot In Joint Limits"
                    if PLOT:
                        ax.plot([it], [dur], 'ro')
                        plt.draw()
            else:
                if PRINT:
                    print "\tNot Shorter: t0 = {0}, t1 = {1}".format(t0, t1)
                    # print "\tt1 - t0 = {0}, Tnew = {1}".format(t1 - t0, Tnew)
                if PLOT:
                    ax.plot([it], [dur], 'ro')
                    plt.draw()
        else:
            if PRINT:
                print "\tInterpolation Failed"
            if PLOT:
                ax.plot([it], [dur], 'ro')
                plt.draw()

    if PRINT:
        print "Successful Shortcuts: {0}".format(nsuccessfulshortcut)

    return rampslistnd


## CheckJointLimits
def CheckJointLimits(robot, rampslistnd):
    ## user defined function
    injointlimits = False
    jointlimits = robot.GetDOFLimits()[1] ## get joint upper limits
    ndof = rampslistnd.ndof
    for i in range(ndof):
        rampslist = rampslistnd[i]
        nramps = len(rampslist)
        for j in range(nramps):
            ramp = rampslist[j]
            if (abs(ramp.a) < EPSILON):
                x_extremal = ramp.Eval(ramp.T)
            else:
                t_extremal = -ramp.v/ramp.a
                if (t_extremal < 0) or (t_extremal > ramp.T):
                    x_extremal = ramp.Eval(ramp.T)
                else:
                    x_extremal = ramp.Eval(t_extremal)
            if (abs(x_extremal) > jointlimits[i]):
                return injointlimits
    injointlimits = True
    return injointlimits


## CheckCollision checks collision for the robot along rampslistnd
def CheckCollision(robot, rampslistnd):
    ## user defined function
    env = robot.GetEnv()
    t = 0
    dt = 0.01
    incollision = False
    qgripper = robot.GetDOFValues()[6]
    while t < rampslistnd.duration:
        robot.SetDOFValues(np.append(rampslistnd.Eval(t), qgripper))
        incollision = (env.CheckCollision(robot) or robot.CheckSelfCollision())
        if (incollision):
            return incollision
        t += dt
    return incollision


def ExtractTrajectory(traj, t0, t1):
    from TOPP import Trajectory
    """ExtractTrajectory extracts a trajectory segment from t = t0 to t = t1 in traj
    """
    assert(t1 > t0)
    
    newchunkslist = []
    i0, rem0 = traj.FindChunkIndex(t0)
    i1, rem1 = traj.FindChunkIndex(t1)

    inthesamechunk = (i0 == i1)

    ## append the traj segment from t0 to the end of chunk[i0]
    newpoly_list = []
    for p in traj.chunkslist[i0].polynomialsvector:
        ## perform variable changing of p(x) = a_n(x)^n + a_(n-1)(x)^(n-1) + ...
        ## by x = y + rem0
        
        a = p.q ## coefficient vector with python convention (highest degree first)
        ## a is a poly1d object
        r = a.r ## polynomial roots
        for i in range(len(r)):
            r[i] = r[i] - rem0
        b = np.poly1d(r, True) ## reconstruct a new polynomial from roots
        ## b is a poly1d object
        b = b*a.coeffs[0] ## multiply back by a_n *** this multiplication does not commute
        
        newpoly = Trajectory.Polynomial(b.coeffs.tolist()[::-1]) ## TOPP convention is weak-term-first
        newpoly_list.append(newpoly)
        
    if (inthesamechunk):
        remchunk0 = Trajectory.Chunk(rem1 - rem0, newpoly_list)
        newchunkslist.append(remchunk0)
    else:
        remchunk0 = Trajectory.Chunk(traj.chunkslist[i0].duration - rem0, newpoly_list)
        newchunkslist.append(remchunk0)

        ## append all chunks in between chunk[i0] and chunk[i1]
        for c in traj.chunkslist[i0 + 1: i1]:
            newchunkslist.append(c)

        ## append the traj segment from the beginning of chunk[i1] to t1
        remchunk1 = Trajectory.Chunk(rem1, traj.chunkslist[i1].polynomialsvector)
        newchunkslist.append(remchunk1)

    return Trajectory.PiecewisePolynomialTrajectory(newchunkslist)
