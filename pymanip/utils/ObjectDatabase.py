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
import time
import ObjectPreprocessing as op
import Utils
from Utils import Colorize
from Grasp import (pX, pY, pZ, mX, mY, mZ)

from os.path import expanduser, join
_home = expanduser('~')
_pymanipdir = join(_home, '.pymanip')

############################################################
#                     Object Database
############################################################
class ObjectDatabase(object):

    def __init__(self, manip_object):
        self.hash = manip_object.GetKinematicsGeometryHash()

        # S contains stable contact surfaces's transformation frame
        # with respect to its COM (except S[0] which is the relative
        # transformation between COM and the first link's frame). see
        # ObjectPreprocess.py for more detail about S
        self.S = op.PlacementPreprocess(manip_object)
        self.TCOM = np.array(self.S[0])
        
        # The z-axis of each surface frame is pointing out of the
        # object. However, what we need is a frame pointing into the
        # object (so that the z-axis is aligned with that of the
        # world).
        Toffset = Utils.ComputeTRot(pX, np.pi)

        # transformationset contains transformation T such that when
        # assigning self.object.SetTransform(T), self.object is
        # resting at a stable placement.
        self.transformationset = [np.dot(Toffset, np.linalg.inv(np.dot(self.TCOM, T)))
                                  for T in self.S[1:]]

        print Colorize('Processing the object geometries', 'yellow')

        self.boxinfos = []
        for ibox in xrange(len(manip_object.GetLinks())):
            ts = time.time()
            boxinfo = op.BoxInfo(manip_object, ibox)
            boxinfo.GetPossibleSlidingDirections()
            boxinfo.Preprocess(self.transformationset)
            te = time.time()
            print Colorize('  box {0} took {1} sec.'.format(ibox, te - ts), 'yellow')
            boxinfo.DestroyOpenRAVEObjects()
            self.boxinfos.append(boxinfo)


def SaveObjectDatabase(objectdatabase, path=_pymanipdir):
    import pickle
    filename = join(path, objectdatabase.hash + '.objectdb.pkl')
    with open(filename, 'wb') as f:
        pickle.dump(objectdatabase, f, pickle.HIGHEST_PROTOCOL)
    print 'The object database has been successfully saved to {0}'.format(filename)


def LoadObjectDatabase(objecthash, path=_pymanipdir):
    import pickle
    filename = join(path, objecthash + '.objectdb.pkl')
    try:
        with open(filename, 'rb') as f:
            objectdb = pickle.load(f)
        print 'The object database has been successfully loaded from {0}'.format(filename)
    except:
        print 'The given objectdb is not found.'
        print 'Probably it has not been created yet.'
        objectdb = None
    return objectdb
    
