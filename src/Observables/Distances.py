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
from Observable import observable


#=============================================================================================
# Distances class definition
#=============================================================================================

class Distances(observable):
    '''
    This class extracts Distances of all kinds from a trajectory and topology file\n
    prerequesit modules are: MDAnalysis

    '''

    def __init__(self, topologyFile, trajectoryFile):
        '''The constructor creates the MDAnalysis universe from a structure file and trajectoryFile
        
        Args:
            -topologyFile    file in pdb or gro or psf format containing the structure of molecule
            -trajectoryFile  file containing the MD trajectory from which one wishes to extract info
        
        Returns:
            none
            
        Usage:
        >>> d = distance('topologyFile.gro', 'trajectoryFile.xtc/dcd')
        '''
        observable.__init__(self, topologyFile, trajectoryFile)
        self._topologyFile = topologyFile
        self._trajectoryFile = trajectoryFile
    
        
    def get_calpha_distance_trajectory(self, skip = 2):
        '''This function returns next to nearest neighbour c-alpha Distances
    
        Args:
           skip
    
        Returns:
           arrayList  Distances array for the whole trajectory [traj_length][num_c_alpha_distances]
    
    
        Useage
    
        >>> Distances = returnCAlphaDistances(_trajectoryFile='test.xtc', _topologyFile='test.gro')
        
        '''  

        protein = self._u.selectAtoms("protein") #Selecting only protein atoms
        nframes = self._u.trajectory.numframes - 1
        calphas = self._u.selectAtoms("name CA")
        atomnums = calphas.atoms.indices()
    
        distance_selections = []
        for k in range(len(atomnums)):
            idist = self._u.selectAtoms("bynum %i"%atomnums[k])
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
            self._u.trajectory.next()
        print 'distance Trajectory has been obtained.'
        return distanceListTrajectory
    
    
    def get_calpha_distance_trajectory_parallel(self):
        print 'this method has not been implemented yet'
    
    def get_heavy_atrom_distance_trajectory(self):
        print 'this method is not implemented yet'
        
    def get_custom_distances_trajectory(self,  _atomsArray,_stop=-1):
        print self._u
        #u = Universe(_topologyFile,_trjectoryFile)
        if self.extensionDcD(self._trajectoryFile):
            collection.clear()
            selections = []
            for l in range(len(_atomsArray)):
                selectionString =self.generate_selection_string(_atomsArray[l], "bynum ")
                s = self._u.selectAtoms(selectionString)
                selections.append(s)
                print s 
            for s in selections:
                print s
                collection.addTimeseries(Timeseries.Distance('r',s))
            if _stop == -1:
                _stop = self._u.trajectory.numframes-1
            collection.compute(self._u.trajectory, stop=_stop)
            return collection
            
        
        
        else:
         
            print 'blub'
    
    
#     def _generate_selection_string(self, selections, selectionType):
#   
#         finalSelection = ""
#         print selections
#         if len(selections)>2:
#             for s in range(len(selections)-1):
#                 finalSelection = finalSelection+selectionType+str(selections[s])+" or "
#             finalSelection =finalSelection+selectionType+str(selections[-1])
#         else:
#             finalSelection = finalSelection+selectionType+str(selections[0])+" or "
#             finalSelection =finalSelection+selectionType+str(selections[1])
#         print finalSelection
#         return finalSelection
    
    
#     def _extensionDcD(self, trajectoryFile):
#         '''This function asserts whether the input trajectory file is in dcd format or not.
#             
#         Args:
#            _trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
#           
#     
#         Returns:
#            boolean  true, if trajectory file is in dcd format
#     
#     
#         Useage:
#     
#         >>> if == self._extensionDcD('test.dcd) 
#         
#         '''  
#         fileName, fileExtension = os.path.splitext(trajectoryFile)
#         if fileExtension == '.dcd':
#             return True
#         else:
#             return False
        
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
    d = Distances(topologyFile,trajectoryFile)

    atomsArray = np.array([[423,1169]])
    #Distances = d.get_calpha_distance_trajectory(trajectoryFile, topologyFile, 2)
    Distances = d.get_custom_distances_trajectory( atomsArray)
    print Distances[0]
