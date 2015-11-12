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

from openravepy import *
from Grasp import *
import numpy as np
import copy
from scipy.spatial import ConvexHull
import Utils

X = np.array([1., 0., 0.])
Y = np.array([0., 1., 0.])
Z = np.array([0., 0., 1.])

## assume that all parts are boxes

def PlacementPreprocess(manip_object, COM = None):
    nparts = len(manip_object.GetLinks())

    with manip_object:
        manip_object.SetTransform(np.eye(4))
        if COM is None:
            ## COM and total mass        
            COM = 0
            M = 0
            for l in manip_object.GetLinks():
                m = l.GetMass()
                M += m
                COM += m*l.GetGlobalCOM()
            COM /= M

        ## Relative transformation between the object's frame and the first
        ## link's frame (the object's frame)
        ## The origin of the object's frame is set to be at the COM.
        ## The rotation is set to be the same as the first link's frame.
        TCOM = np.eye(4)
        TCOM[0][3] = COM[0]
        TCOM[1][3] = COM[1]
        TCOM[2][3] = COM[2]
        TCOMinv = np.linalg.inv(TCOM)
        
        pointsets = []
        extents = []

        manip_object.SetTransform(TCOMinv)
        for i in range(nparts):
            l = manip_object.GetLinks()[i]
            extent = l.GetGeometries()[0].GetBoxExtents()
            extent = np.around(extent, decimals = 6)
            extents.append(extent)

            com = l.GetGlobalCOM()
            [dx, dy, dz] = extent
            
            X = l.GetTransform()[0:3, 0]
            Y = l.GetTransform()[0:3, 1]
            Z = l.GetTransform()[0:3, 2]
            
            ## extreme points of the part
            pointset = [com + dx*X + dy*Y + dz*Z,
                        com + dx*X + dy*Y - dz*Z,
                        com + dx*X - dy*Y + dz*Z,
                        com + dx*X - dy*Y - dz*Z,
                        com - dx*X + dy*Y + dz*Z,
                        com - dx*X + dy*Y - dz*Z,
                        com - dx*X - dy*Y + dz*Z,
                        com - dx*X - dy*Y - dz*Z]
            pointsets += pointset

        hull = ConvexHull(pointsets, qhull_options = 'E0.001')

        ## temp contains equations describing surfaces.
        ## each equation is in the form av + b = 0
        ## however, there is redundancy in temp. 
        ## we need to remove duplicate equations
        temp = hull.equations.tolist()

        E = []
        for i in xrange(len(temp)):
            similar = False
            for e in E:
                if np.allclose(temp[i], e):
                    similar = True
                    break
            if not similar:
                E.append(temp[i])
        E = np.asarray(E)
        E = np.around(E, decimals = 6)                
        nsurfaces = len(E)

    ## Av + b <= 0 for v in the polyhedron
    A = np.array(E[0:nsurfaces, 0:3])
    b = np.array(E[0:nsurfaces, 3])

    S = []
    S.append(TCOM)
    ## S[i] is the relative transformation between the surface's
    ## coordinate and the object's cooredinate
    ## except S[0] which is TCOM

    ## NOTE: S contains only the surfaces that can be contact surfaces,
    ##       i.e., the surfaces that p lies inside.
    for i in range(nsurfaces):
        z = np.array(A[i]) ## surface normal
        d = -np.array(b[i]) ## surface offset
        p = d*z ## surface's frame origin
        if (not np.all(np.around(np.dot(A, p) + b, decimals = 6) <= 0)):
            ## p lies outside the polygon
            # print "surface {0} is an invalid contact surface".format(i)
            continue

        ## x-axis
        x = PerpendicularTo(z)
        ## y-axis
        y = np.cross(z, x)    
        ## each column of R is the axis described in the object's frame
        R = np.vstack((x, y, z)).T
        # p = np.dot(TCOM, np.append(p, 1))[0:3]
        p = np.reshape(p, (3, 1))
        T = np.vstack((np.hstack((R, p)), np.array([0., 0., 0., 1.])))
        S.append(T)

    return S


def PerpendicularTo(v):
    """ Finds an arbitrary perpendicular vector to *v*."""
    # for two vectors (x, y, z) and (a, b, c) to be perpendicular,
    # the following equation has to be fulfilled
    #     0 = ax + by + cz
    if (not (len(v) == 3)):
        raise ValueError('dimension not compatible')    
    
    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0:
        raise ValueError('zero-vector')

    # If one dimension is zero, this can be solved by setting that to
    # non-zero and the others to zero. Example: (4, 2, 0) lies in the
    # x-y-Plane, so (0, 0, 1) is orthogonal to the plane.
    if v[0] == 0:
        return np.array([1., 0., 0.])
    if v[1] == 0:
        return np.array([0., 1., 0.])
    if v[2] == 0:
        return np.array([0., 0., 1.])

    # arbitrarily set a = b = 1
    # then the equation simplifies to
    #     c = -(x + y)/z
    c = -(v[0] + v[1])/float(v[2])
    d = 1./np.sqrt(2 + abs(c)**2)
    return np.array([d, d, d*c])


class BoxInfo(object):
    DMAX = 0.0425 ## max width that the gripper can grip
    GRIPPEROFFSET = 0.08
    L = 0.35 #0.320592 (gripper length)

    def __init__(self, manip_object, linkindex):
        self.objectname = manip_object.GetName()
        self.linkindex = linkindex

        self.env = Environment()
        Clone_Bodies = 1
        self.env.Clone(manip_object.GetEnv(), Clone_Bodies)

        self.collisionchecker = RaveCreateCollisionChecker(self.env, 'ode')
        self.env.SetCollisionChecker(self.collisionchecker)
        
        ## remove all other body
        for kinbody in self.env.GetBodies():
            if not (kinbody.GetName() == self.objectname):
                self.env.Remove(kinbody)
        assert(len(self.env.GetBodies()) == 1)

        self.object = self.env.GetKinBody(self.objectname)
        self.link = self.object.GetLinks()[self.linkindex]

        ## load the gripper
        self.env.Load('../xml/gripper_aabb.xml')
        self.gripper = self.env.GetKinBody('gripper_aabb')

        ## create a floor for testing
        floor = RaveCreateKinBody(self.env, '')
        floor.InitFromBoxes(np.array([[0.0, 0.0, 0.0, 10.0, 10.0, 0.05]]))
        for geom in floor.GetLinks()[0].GetGeometries():
            geom.SetDiffuseColor(np.array([0.6, 0.6, 0.6]))
        floor.SetName('floor')
        self.env.Add(floor)
        Tfloor = np.eye(4)
        Tfloor[2][3] -= 0.050001
        floor.SetTransform(Tfloor)
        
        self.extents = self.link.GetGeometries()[0].GetBoxExtents()
        self.possibleapproachingdir = dict() # each entry depends on a contact surface
        self.possibleslidingdir = dict()     # each entry depends on an approaching dir
        self.intervals = dict()
    

    def GetPossibleSlidingDirections(self):
        """
        GetPossibleSlidingDirections returns a set containing possible
        sliding direction of the gripper for each case of approaching
        directions.
        
        return value: lib 
        
        lib[approachingdir] is a set of possible sliding direction given
        the 'approachingdir'.
        """
        objx = self.extents[0]
        objy = self.extents[1]
        objz = self.extents[2]

        ## approaching direction is +X (and -X)
        temp = []
        if (objz < self.DMAX):
            temp.append(pY)
            if (objy < self.DMAX):
                temp.append(pZ)
        elif (objy < self.DMAX):
            temp.append(pZ)
        self.possibleslidingdir[pX] = temp
        self.possibleslidingdir[mX] = temp

        ## approaching direction is +Y (and -Y)
        temp = []
        if (objz < self.DMAX):
            temp.append(pX)
            if (objx < self.DMAX):
                temp.append(pZ)
        elif (objx < self.DMAX):
            temp.append(pZ)
        self.possibleslidingdir[pY] = temp
        self.possibleslidingdir[mY] = temp

        ## approaching direction is +Z (and -Z)
        temp = []
        if (objy < self.DMAX):
            temp.append(pX)
            if (objx < self.DMAX):
                temp.append(pY)
        elif (objx < self.DMAX):
            temp.append(pY)
        self.possibleslidingdir[pZ] = temp
        self.possibleslidingdir[mZ] = temp


    def Preprocess(self, transformationset):
        """
        Preprocess examines valid approaching directions for each
        object's transformation. Also for each pair of approaching
        direction and sliding direction Preprocess examines the
        sliding range.

        transformationset contains every object's transformation T
        that results in stable configurations.
        """
        
        nsurfaces = len(transformationset)
        for isurface in xrange(nsurfaces):
            self.object.SetTransform(transformationset[isurface])
            
            plink = self.link.GetGlobalCOM()
            Tlink = self.link.GetTransform()

            xvect = np.reshape(copy.copy(Tlink[0:3, pX]), (3, ))
            yvect = np.reshape(copy.copy(Tlink[0:3, pY]), (3, ))
            zvect = np.reshape(copy.copy(Tlink[0:3, pZ]), (3, ))
            
            ## for each object's contact surface check all six
            ## surfaces of the box
            self.possibleapproachingdir[isurface] = []
            for appdir in [pX, pY, pZ, mX, mY, mZ]:
                ## for each approachingdirection                
                
                appvector = np.reshape(copy.copy(Tlink[0:3, np.mod(appdir, 3)]), (3, )) 
                # approaching vector
                if (appdir > 2):
                    appvector *= -1

                aZ = np.dot(appvector, Z)
                
                ########## CHECK I
                ## check if the approached surface is in contact
                if np.allclose(aZ, 1):
                    ## the approaching direction is aligned with Z
                    if (plink[2] - self.extents[np.mod(appdir, 3)] < self.L):
                        ## this approacing direction is invalid.
                        ## continue to the next direction.
                        continue
                
                ########## CHECK II
                ## check if the approached surface is perpendicular to the floor
                elif np.allclose(aZ, 0):
                    ax = np.dot(appvector, xvect)
                    if (np.allclose(ax, 1) or np.allclose(-ax, 1)):
                        ## approaching direction is x or -x
                        yZ = np.dot(yvect, Z)
                        zZ = np.dot(zvect, Z)
                        if (np.allclose(yZ, 1) or np.allclose(-yZ, 1)):
                            ## the local y is parallel to Z
                            if (plink[2] + self.extents[1] < self.GRIPPEROFFSET):
                                ## this approacing direction is invalid.
                                ## continue to the next direction.
                                continue
                        elif (np.allclose(zZ, 1) or np.allclose(-zZ, 1)):
                            ## the local z is parallel to Z
                            if (plink[2] + self.extents[2] < self.GRIPPEROFFSET):
                                ## this approacing direction is invalid.
                                ## continue to the next direction.
                                continue
                    ay = np.dot(appvector, yvect)
                    if (np.allclose(ay, 1) or np.allclose(-ay, 1)):
                        ## approaching direction is y or -y
                        xZ = np.dot(xvect, Z)
                        zZ = np.dot(zvect, Z)
                        if (np.allclose(xZ, 1) or np.allclose(-xZ, 1)):
                            ## the local y is parallel to Z
                            if (plink[2] + self.extents[0] < self.GRIPPEROFFSET):
                                ## this approacing direction is invalid.
                                ## continue to the next direction.
                                continue
                        elif (np.allclose(zZ, 1) or np.allclose(-zZ, 1)):
                            ## the local z is parallel to Z
                            if (plink[2] + self.extents[2] < self.GRIPPEROFFSET):
                                ## this approacing direction is invalid.
                                ## continue to the next direction.
                                continue
                    az = np.dot(appvector, zvect)
                    if (np.allclose(az, 1) or np.allclose(-az, 1)):
                        ## approaching direction is x or -x
                        xZ = np.dot(xvect, Z)
                        yZ = np.dot(yvect, Z)
                        if (np.allclose(xZ, 1) or np.allclose(-xZ, 1)):
                            ## the local x is parallel to Z
                            if (plink[2] + self.extents[0] < self.GRIPPEROFFSET):
                                ## this approacing direction is invalid.
                                ## continue to the next direction.
                                continue
                        elif (np.allclose(yZ, 1) or np.allclose(-yZ, 1)):
                            ## the local y is parallel to Z
                            if (plink[2] + self.extents[1] < self.GRIPPEROFFSET):
                                ## this approacing direction is invalid.
                                ## continue to the next direction.
                                continue

                ########## CHECK III
                elif (aZ > 0):
                    normalvector = -1.0*appvector # normal vector to the surface
                    k1 = PerpendicularTo(normalvector)
                    k2 = np.cross(normalvector, k1)
                    
                    ## normal vector is pointing to the floor
                    theta = np.pi/2 - np.arccos(aZ)
                    tantheta = np.tan(theta)
                    
                    if (np.allclose(np.dot(k1, Z), 0)):
                        ## k1 is parallel to the floor
                        k2x = np.dot(k2, xvect)
                        if (np.allclose(k2x, 1) or np.allclose(-k2x, 1)):
                            ## k2 is parallel to the local x
                            D = 2*self.extents[0]/tantheta
                            if D < self.L:
                                continue
                        k2y = np.dot(k2, yvect)
                        if (np.allclose(k2y, 1) or np.allclose(-k2y, 1)):
                            ## k2 is parallel to the local y
                            D = 2*self.extents[1]/tantheta
                            if D < self.L:
                                continue
                        k2z = np.dot(k2, zvect)
                        if (np.allclose(k2z, 1) or np.allclose(-k2z, 1)):
                            ## k2 is parallel to the local z
                            D = 2*self.extents[2]/tantheta
                            if D < self.L:
                                continue
                    elif (np.allclose(np.dot(k2, Z), 0)):
                        ## k1 is parallel to the floor
                        k1x = np.dot(k1, xvect)
                        if (np.allclose(k1x, 1) or np.allclose(-k1x, 1)):
                            ## k1 is parallel to the local x
                            D = 2*self.extents[0]/tantheta
                            if D < self.L:
                                continue
                        k1y = np.dot(k1, yvect)
                        if (np.allclose(k1y, 1) or np.allclose(-k1y, 1)):
                            ## k1 is parallel to the local y
                            D = 2*self.extents[1]/tantheta
                            if D < self.L:
                                continue
                        k1z = np.dot(k1, zvect)
                        if (np.allclose(k1z, 1) or np.allclose(-k1z, 1)):
                            ## k1 is parallel to the local z
                            D = 2*self.extents[2]/tantheta
                            if D < self.L:
                                continue
                            
                ## this approaching direction passes the three tests
                ## it is now a candidate for obtaining sliding ranges
                """
                interval = [interval_1, interval_2, ..., interval_n]
                interval_i = (start, length)
                """
                
                possibleslidingdir = self.possibleslidingdir[appdir]
                validapproachingdir = False

                if len(possibleslidingdir) == 2:
                    slidingdir1 = possibleslidingdir[0]
                    interval1 = self.ObtainSlidingRanges(Tlink, appdir, 
                                                         slidingdir1, case = 2)
                    if len(interval1) > 0:
                        self.intervals[isurface, appdir, slidingdir1] = interval1
                        validapproachingdir = True
                    
                    slidingdir2 = possibleslidingdir[1]
                    interval2 = self.ObtainSlidingRanges(Tlink, appdir, 
                                                         slidingdir2, case = 2)
                    if len(interval2) > 0:
                        self.intervals[isurface, appdir, slidingdir2] = interval2
                        validapproachingdir = True
                        
                elif len(possibleslidingdir) == 1:
                    slidingdir = possibleslidingdir[0]
                    interval = self.ObtainSlidingRanges(Tlink, appdir, slidingdir)
                    if len(interval) > 0:
                        self.intervals[isurface, appdir, slidingdir] = interval
                        validapproachingdir = True

                if not validapproachingdir:
                    continue
                
                self.possibleapproachingdir[isurface].append(appdir)
                

    
    def ObtainSlidingRanges(self, Tobj, approachingdir, slidingdir, case = 1):
        ## assume that the object is already in place        
        d = self.extents[np.mod(slidingdir, 3)]
        step = 0.04 ## each step along the sliding direction is 4 cm.
        Tstep = np.eye(4)

        qgrasp = [approachingdir, slidingdir, 0.0]
        ## Tgripper places the gripper at the middle of the object
        Tgripper = Utils.ComputeTGripper(Tobj, qgrasp, self.extents)

        interval = []

        if case == 2:
            self.gripper.SetTransform(Tgripper)
            if not (self.env.CheckCollision(self.gripper)):
                interval = [(0, 0)]

        elif case == 1:            
            domain = np.linspace(-d, d, int(2*d/step) + 1)
            Tstep[0][3] = domain[0]
            self.gripper.SetTransform(np.dot(Tgripper, Tstep))
            if (self.env.CheckCollision(self.gripper)):
                prev_in_collision = True
            else:
                prev_in_collision = False
            prevstart = domain[0]

            for i in xrange(1, len(domain)):
                Tstep[0][3] = domain[i]
                self.gripper.SetTransform(np.dot(Tgripper, Tstep))
                if (self.env.CheckCollision(self.gripper)):
                    if prev_in_collision:
                        prevstart = domain[i]
                    else:
                        if np.allclose(domain[i - 1], prevstart):
                            ## this interavl contains only one number
                            ## but we relax it by assinging length =
                            ## 0.02 (half a step) so that when we
                            ## sample this interval later, the prob of
                            ## obtaining a number here is not zero
                            interval.append((prevstart, 0.02))
                        else:
                            interval.append((prevstart, domain[i - 1] - prevstart))
                        prevstart = domain[i]
                        prev_in_collision = True
                else:
                    if (prev_in_collision):
                        prev_in_collision = False
                        prevstart = domain[i]
                        

            if not np.allclose(domain[-1], prevstart):
                interval.append((prevstart, domain[-1] - prevstart))
        
        return np.around(interval, decimals = 6).tolist()

