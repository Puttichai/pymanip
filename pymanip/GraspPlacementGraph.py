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
import copy

# Global parameters
TRANSIT = 0
TRANSFER = 1

############################################################
#                  Grasp-Placement Graph
############################################################
class GraspPlacementGraph(object):
    
    def __init__(self, verticeslist=[], autogenerate=True):
        """
        Parameters
        ----------
        verticeslist : list, optional
            verticeslist is a list containing vertices.
            Each vertex is a pair (contactsurface, approachingdir).
            If a vertex is in CP, approachingdir is None.
        autogenerate : bool, optional
            autogenerate indicates whether or not to autogenerate 
            a graph dictionary.
        """
        self.verticeslist = copy.copy(verticeslist)
        self.graphdict = dict()
        self.nvertices = len(self.verticeslist)
        
        if autogenerate:
            self.GenerateGraphDictionary()


    def __getitem__(self, key):
        return self.graphdict[key]


    def GenerateGraphDictionary(self):
        if (self.nvertices > 0):
            # Add additional vertices (vertices in CP)
            for index in xrange(self.nvertices):
                vnew = (self.verticeslist[index][0], None)
                if vnew not in self.verticeslist:
                    self.verticeslist.append(vnew)
            self.verticeslist.sort()
            self.nvertices = len(self.verticeslist)

            # Initialize graph dictionary's entries
            for v in self.verticeslist:
                self.graphdict[v] = []

            # Add edges of the graph
            for i in xrange(self.nvertices):
                vi = self.verticeslist[i]
                vi_in_cgcp = (vi[1] is not None)
                for j in xrange(i + 1, self.nvertices):
                    vj = self.verticeslist[j]
                    vj_in_cgcp = (vj[1] is not None)
                    if ((vi[0] == vj[0]) or 
                        ((vi_in_cgcp) and (vj_in_cgcp) and (vi[1] == vj[1]))):
                        if (vj not in self.graphdict[vi]):
                            self.graphdict[vi].append(vj)
                        if (vi not in self.graphdict[vj]):
                            self.graphdict[vj].append(vi)
            
            # Add self-cycle
            for key in self.graphdict.keys():
                if key[1] is not None:
                    self.graphdict[key].append(key)
                    self.graphdict[key].sort()


    #
    # Dijkstra's algorithm
    #
    def RunDijkstra(self, vstart):
        dist = dict()
        prev = dict()
        Q = dict() ## a queue of unvisited nodes

        ## initialization
        for key in self.graphdict:
            dist[key] = np.infty
            prev[key] = None
            Q[key] = None
        dist[vstart] = 0

        while len(Q) > 0:
            [mindist, keymin] = self.ExtractMin(Q, dist)
            Q.pop(keymin)

            for v in self.graphdict[keymin]:
                alt = dist[keymin] + 1
                if alt < dist[v]:
                    dist[v] = alt
                    prev[v] = keymin

        return dist, prev


    def ExtractMin(self, Q, dist):
        ## a naive implementation of ExtractMin
        mindist = np.infty
        keymin = Q.keys()[0]
        for key in Q.keys():
            if (dist[key] < mindist):
                keymin = key
                mindist = dist[key]
        return [mindist, keymin]

    
    def FindShortestPath(self, vstart, vgoal):
        dist, prev = self.RunDijkstra(vstart)
        path = []
        curnode = copy.copy(vgoal)
        path.append(curnode)
        curdist = dist[vgoal]
        while curdist > 0:
            curnode = prev[curnode]
            curdist = dist[curnode]
            path.append(curnode)
        return path[::-1]


    def FindShortestPathLength(self, vstart, vgoal):
        dist, prev = self.RunDijkstra(vstart)
        return dist[vgoal]


    #
    # Dynamic programming for finding all paths of given length
    #
    def FindPathsOfLengthK(self, vstart, vgoal, k, sols=[], removeinfpaths=True):
        """
        FindPathsOfLengthK returns all paths connecting vstart and
        vgoal which have length k.
        """
        if len(sols) == 0:
            sols.append([vstart])

        if (k == 1):
            if vgoal in self.graphdict[vstart]:
                for path in sols:
                    path.append(vgoal)
                return sols
            else:
                return []

        newsols = []
        ## search in its neighbors
        for v in self.graphdict[vstart]:
            temp = copy.deepcopy(sols)
            for path in temp:
                path.append(v)
            paths = self.FindPathsOfLengthK(v, vgoal, k - 1, temp)
            if len(paths) == 0:
                continue

            for path in paths:
                newsols.append(path)
                
        if removeinfpaths:
            if len(newsols) == 0:
                return newsols
            else:
                return self.RemoveInfeasiblePaths(newsols)
        else:
            return newsols


    def RemoveInfeasiblePaths(self, paths0):
        """
        Parameters
        ----------
        paths0 : a list
            A list containing manipulation paths of length k.

        Returns
        -------
        sols : a list
            A list containing non-redundant manipulation paths of 
            length k.
            Note that sols \subset paths0.
        """
        paths = copy.deepcopy(paths0)
        l = len(paths[0])
        sols = []

        for path in paths:
            ## examine each path
            OK = True

            for i in xrange(1, l):
                ## examine each segment
                curnode = path[i]
                if (i == 1):
                    ## the first segment can be anything
                    prevnode = path[i - 1]
                    if (not (len(prevnode) == len(curnode))):
                        prevtype = TRANSIT
                    else:
                        if (prevnode[0] == curnode[0]):
                            prevtype = TRANSIT
                        else:
                            prevtype = TRANSFER
                    prevnode = curnode
                else:
                    ## otherwise, in order to be valid (not redundant),
                    ## the currenttype has to be different
                    if (not (len(prevnode) == len(curnode))):
                        curtype = TRANSIT
                    else:
                        if ((prevnode[0] == curnode[0]) and 
                            (prevnode[1] == curnode[1])):
                            curtype = np.mod(prevtype + 1, 2)
                        elif (prevnode[0] == curnode[0]):
                            curtype = TRANSIT
                        else:
                            curtype = TRANSFER
                    if (curtype == prevtype):
                        ## got a problem here
                        OK = False
                        break

                    ## if no problem, continue
                    prevtype = curtype
                    prevnode = curnode
            if OK:
                sols.append(path)

        return sols


############################################################
#                    Graph Utilities
############################################################
def CreateFWDirectedGraphFromPathsList(pathslist):
    """
    Gnew will contain vertices of the form (pathlevel, isurface,
    appdir) instead of (isurface, appdir)
    """
    Gnew = GraspPlacementGraph(autogenerate=False)
    k = len(pathslist[0])
    vstart = (0, ) + pathslist[0][0]
    Gnew.graphdict[vstart] = []
    Gnew.verticeslist.append(vstart)

    for path in pathslist:
        for pathlevel in xrange(1, k):
            prevkey = (pathlevel - 1, ) + path[pathlevel - 1]
            curkey = (pathlevel, ) + path[pathlevel]
            
            if (curkey not in Gnew.graphdict[prevkey]):
                Gnew.graphdict[prevkey].append(curkey)
                
            if (curkey not in Gnew.graphdict):
                Gnew.graphdict[curkey] = []
                Gnew.verticeslist.append(curkey)
                
    return Gnew


def CreateBWDirectedGraphFromPathsList(pathslist):
    """
    Gnew will contain vertices of the form (pathlevel, isurface,
    appdir) instead of (isurface, appdir)
    """
    Gnew = GraspPlacementGraph(autogenerate=False)
    k = len(pathslist[0])
    vgoal = (k - 1, ) + pathslist[0][k - 1]
    Gnew.graphdict[vgoal] = []
    Gnew.verticeslist.append(vgoal)

    for path in pathslist:
        for pathlevel in xrange(k - 1, 0, -1):
            curkey = (pathlevel, ) + path[pathlevel]
            prevkey = (pathlevel - 1, ) + path[pathlevel - 1]
            
            if (prevkey not in Gnew.graphdict[curkey]):
                Gnew.graphdict[curkey].append(prevkey)
                
            if (prevkey not in Gnew.graphdict):
                Gnew.graphdict[prevkey] = []
                Gnew.verticeslist.append(prevkey)
                
    return Gnew
