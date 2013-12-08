#!/usr/bin/env python

'''

    
    Author:          Antonia Mey
        
    Email:           antonia.mey@fu-berlin.de     
     
    Date:            03.12.2013
    
    Version:         1.0.0 Beta
    
    #=============================================================================================\n
    # COPYRIGHT NOTICE\n
    #\n
    # Written by Antonia Mey <antonia.mey@fu-berlin.de>\n
    #
    # This program is free software; you can redistribute it and/or modify it under the terms of\n
    # the GNU General Public License as published by the Free Software Foundation; either version 2\n
    # of the License, or (at your option) any later version.\n
    #\n
    # This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;\n
    # without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n
    # See the GNU General Public License for more details.\n
    # \n
    # You should have received a copy of the GNU General Public License along with this program;\n
    # if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,\n
    # Boston, MA  02110-1301, USA.\n
    #=============================================================================================\n

'''

#================================
#Imports
#================================
from MDAnalysis import *
from MDAnalysis.analysis.distances import *
import numpy as np
import argparse
from MDAnalysis import collection
import os


#=============================================================================================
# distances class definition
#=============================================================================================

class distances:
    '''
    This class extracts distances of all kinds from a trajectory and topology file\n
    prerequesit modules are: MDAnalysis

    '''

    def __init__(self):
        '''
        Constructor, does not take any arguments
        '''
        pass
    
        
    def getCAlphaDistanceTrajectory(self,_trajectoryFile, _topologyFile, skip = 2):
        '''This function returns next to nearest neighbour c-alpha distances
    
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc
           _topologyFile (str): topologyFilename.gro
    
        Returns:
           arrayList  distances array for the whole trajectory [traj_length][num_c_alpha_distances]
    
    
        Useage
    
        >>> distances = returnCAlphaDistances(_trajectoryFile='test.xtc', _topologyFile='test.gro')
        
        '''  
        u = Universe(_topologyFile,_trajectoryFile)
        print "Loading trajectory file %s"%_trajectoryFile
        protein = u.selectAtoms("protein") #Selecting only protein atoms
        nframes = u.trajectory.numframes - 1
        calphas = u.selectAtoms("name CA")
        atomnums = calphas.atoms.indices()
    
        distance_selections = []
        for k in range(len(atomnums)):
            idist = u.selectAtoms("bynum %i"%atomnums[k])
            distance_selections.append(idist)
        ndist = len(distance_selections)

        distanceListTrajectory = []
        # Loop over all frames:
        print "Loop over all frames..."
        for i in range(nframes):
            all_distances = []
            for k in range(ndist-skip):
                for l in range(k+skip,ndist):
                        result = dist(distance_selections[k],distance_selections[l])
                        all_distances.append( result[2][0])
                        
            distanceListTrajectory.append(all_distances)
            u.trajectory.next()
        print 'distance Trajectory has been obtained.'
        return distanceListTrajectory
    
    
    def getCAlphaDistanceTrajectoryParallel(self):
        print 'this method has not been implemented yet'
    
    def getHeavyAtromDistanceTrajectory(self):
        print 'this method is not implemented yet'
        
    def getCustomDistancesTrajectory(self, _topologyFile,_trjectoryFile, _atomsArray,_stop=-1):
        u = Universe(_topologyFile,_trjectoryFile)
        if self._extensionDcD(_trjectoryFile):
            collection.clear()
            selections = []
            for l in range(len(_atomsArray)):
                selectionString =self._generateSelectionString(_atomsArray[l], "bynum ")
                s = u.selectAtoms(selectionString)
                selections.append(s)
                print s 
            for s in selections:
                print s
                collection.addTimeseries(Timeseries.Distance('r',s))
            if _stop == -1:
                _stop = u.trajectory.numframes-1
            collection.compute(u.trajectory, stop=_stop)
            return collection
            
        
        
        else:
         
            print 'blub'
    
    
    def _generateSelectionString(self, selections, selectionType):
  
        finalSelection = ""
        print selections
        if len(selections)>2:
            for s in range(len(selections)-1):
                finalSelection = finalSelection+selectionType+str(selections[s])+" or "
            finalSelection =finalSelection+selectionType+str(selections[-1])
        else:
            finalSelection = finalSelection+selectionType+str(selections[0])+" or "
            finalSelection =finalSelection+selectionType+str(selections[1])
        print finalSelection
        return finalSelection
    
    
    def _extensionDcD(self, trajectoryFile):
        '''This function asserts whether the input trajectory file is in dcd format or not.
            
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
          
    
        Returns:
           boolean  true, if trajectory file is in dcd format
    
    
        Useage:
    
        >>> if == self._extensionDcD('test.dcd) 
        
        '''  
        fileName, fileExtension = os.path.splitext(trajectoryFile)
        if fileExtension == '.dcd':
            return True
        else:
            return False
        
#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
    
    #tsting the methods
    parser = argparse.ArgumentParser(description="Feature extractions from xtc files")
    parser.add_argument('-f', help='xtc inputfile');
    parser.add_argument('-c', help='coordinate file in Gromacs .gro format');
    args = parser.parse_args()

    topologyFile = args.c
    trajectoryFile = args.f
    d = distances()

    atomsArray = np.array([[423,1169]])
    #distances = d.getCAlphaDistanceTrajectory(trajectoryFile, topologyFile, 2)
    distances = d.getCustomDistancesTrajectory(topologyFile, trajectoryFile, atomsArray)
    print distances[0]
