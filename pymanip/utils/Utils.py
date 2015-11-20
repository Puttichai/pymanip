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

import numpy as np
import TOPP

# DIRECTIONS (p means plus; m means minus)
pX = 0
pY = 1
pZ = 2
mX = 3
mY = 4
mZ = 5

X = np.array([1., 0., 0.])
Y = np.array([0., 1., 0.])
Z = np.array([0., 0., 1.])

def ComputeTRot(axis, theta):
    c = np.cos(theta)
    s = np.sin(theta)
    if (axis == pX):
        T = np.array([[ 1,  0,  0,  0], 
                      [ 0,  c, -s,  0], 
                      [ 0,  s,  c,  0] ,
                      [ 0,  0,  0,  1]], dtype=float)
    elif (axis == pY):
        T = np.array([[ c,  0,  s,  0], 
                      [ 0,  1,  0,  0], 
                      [-s,  0,  c,  0] ,
                      [ 0,  0,  0,  1]], dtype=float)
    elif (axis == pZ):
        T = np.array([[ c, -s,  0,  0], 
                      [ s,  c,  0,  0], 
                      [ 0,  0,  1,  0] ,
                      [ 0,  0,  0,  1]], dtype=float)
    else:
        raise ValueError('wrong axis of rotation given')
    return T


def ComputeTTrans(axis, d):
    T = np.eye(4)
    T[axis][3] = d
    return T


def ComputeTObject(qobj, S, zoffset=0.725):
    """
    ComputeTObject computes a transformation matrix for the object
    such that it is in a stable placement.
    """
    xobj = qobj[0]
    yobj = qobj[1]
    thetaobj = qobj[2]
    isurface = qobj[3]
    
    TCOM = np.array(S[0])
    Tsurf = np.array(S[isurface + 1])

    Tp = np.eye(4)
    Tp[0][3] = xobj
    Tp[1][3] = yobj
    Tp[2][3] = zoffset

    ## rotate about x-axis first so that the z-axis is pointing downward
    Tx = ComputeTRot(pX, np.pi)
    Tz = ComputeTRot(pZ, thetaobj)
    
    T = np.dot(Tp, np.dot(Tx, Tz))
    
    Tobj = np.dot(T, np.linalg.inv(np.dot(TCOM, Tsurf)))
    
    return Tobj


def ComputeTGripper(T_obj, q_grasp, extents, unitscale=False):
    """
    ComputeTGripper computes a transformation matrix for the gripper
    given the approaching direction and the sliding direction for the
    gripper.

    convention: 
        - approaching direction and sliding direction are specified 
          in the object's local frame 
        - approaching direction k implies np.dot(zg, k) == 1 
        - sliding direction is always a positive direction (aligned with +xg)
        - the gripper can slide along the sliding axis within the range [-d, d]

    """
    # Select the last three elements of qgrasp
    # Note that qgrasp = [ibox, approachingdir, slidingdir, delta]
    approachingdir, slidingdir, delta = q_grasp[-3:]
    approachingdir = np.mod(approachingdir, 6)
    [dx, dy, dz] = extents

    if (approachingdir == pX):
        T0 = ComputeTRot(pY, np.pi/2)
        if (slidingdir == pZ):
            T1 = ComputeTRot(pZ, np.pi)
        elif (slidingdir == pY):
            T1 = ComputeTRot(pZ, np.pi/2)
    elif (approachingdir == mX):
        T0 = ComputeTRot(pY, -np.pi/2)
        if (slidingdir == pZ):
            T1 = np.eye(4)
        elif (slidingdir == pY):
            T1 = ComputeTRot(pZ, -np.pi/2)
    elif (approachingdir == pY):
        T0 = ComputeTRot(pX, -np.pi/2)
        if (slidingdir == pX):
            T1 = np.eye(4)
        elif (slidingdir == pZ):
            T1 = ComputeTRot(pZ, -np.pi/2)
    elif (approachingdir == mY):
        T0 = ComputeTRot(pX, np.pi/2)
        if (slidingdir == pX):
            T1 = np.eye(4)
        elif (slidingdir == pZ):
            T1 = ComputeTRot(pZ, np.pi/2)
    elif (approachingdir == pZ):
        T0 = np.eye(4)
        if (slidingdir == pX):
            T1 = np.eye(4)
        elif (slidingdir == pY):
            T1 = ComputeTRot(pZ, np.pi/2)            
    else:
        T0 = ComputeTRot(pX, np.pi)
        if (slidingdir == pX):
            T1 = np.eye(4)
        elif (slidingdir == pY):
            T1 = ComputeTRot(pZ, -np.pi/2)

    d = extents[slidingdir]

    if approachingdir > 2:
        T2 = ComputeTTrans(pZ, -extents[approachingdir - 3] + 0.004)
    else:
        T2 = ComputeTTrans(pZ, -extents[approachingdir] + 0.004)
            
    if unitscale:
        T3 = ComputeTTrans(pX, delta*d)
    else:
        T3 = ComputeTTrans(pX, delta)
    
    return reduce(np.dot, [T_obj, T0, T1, T2, T3])


def EnableGripper(robot):
    robot.SetActiveDOFs(np.array([0, 1, 2, 3, 4, 5, 6]))


def DisableGripper(robot):
    robot.SetActiveDOFs(np.array([0, 1, 2, 3, 4, 5]))


def ReverseTrajectory(topptrajstring):
    topptraj = TOPP.Trajectory.PiecewisePolynomialTrajectory.FromString(topptrajstring)
    newchunkslist = []
    
    for chunk in topptraj.chunkslist:
        T = chunk.duration
        newpoly_list = []
        for p in chunk.polynomialsvector:
            # Perform variable changing of p(x) = a_n(x)^n + a_(n-1)(x)^(n-1) + ...
            # by x = T - y        
            a = p.q # coefficient vector with python convention (highest degree first)
            # a is a poly1d object
            r = a.r
            newr = [T - k for k in r]            
            b = np.poly1d(newr, True) # reconstruct a new polynomial from roots
            b = b*a.coeffs[0] # multiply back by a_n
            # *** this multiplication does not commute            
            if (b(0)*a(T) < 0):
                # correct the coeffs if somehow the polynomial is flipped
                b = b*-1.0
            # TOPP convention is weak-term-first
            newpoly = TOPP.Trajectory.Polynomial(b.coeffs.tolist()[::-1])
            newpoly_list.append(newpoly)
        newchunk = TOPP.Trajectory.Chunk(chunk.duration, newpoly_list)
        newchunkslist.insert(0, newchunk)
    newtopptraj = TOPP.Trajectory.PiecewisePolynomialTrajectory(newchunkslist)
    return str(newtopptraj)
    

##################################################
colors = dict()
colors['black'] = 0
colors['red'] = 1
colors['green'] = 2
colors['yellow'] = 3
colors['blue'] = 4
colors['magenta'] = 5
colors['cyan'] = 6
colors['white'] = 7
def Colorize(string, color = 'white', bold = True):
    newstring = '\033['
    newstring += (str(int(bold)) + ';')
    newstring += ('3' + str(colors[color]))
    newstring += 'm'
    newstring += string
    newstring += '\033[0m' # reset the subsequent text back to normal
    return newstring
