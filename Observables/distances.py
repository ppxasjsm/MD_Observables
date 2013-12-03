#!/usr/bin/env python

'''
    File:            /home/mi/ppxasjsm/Eclipse/MD_observables/Observables/distances.py
    Description:
    
    Author:          Antonia Mey    
    Email:           antonia.mey@fu-berlin.de      
    Date:            03.12.2013
    Version:         0.0.1 Beta
    
    #=============================================================================================
    # COPYRIGHT NOTICE
    #
    # Written by Antonia Mey <antonia.mey@fu-berlin.de>
    #
    # This program is free software; you can redistribute it and/or modify it under the terms of
    # the GNU General Public License as published by the Free Software Foundation; either version 2
    # of the License, or (at your option) any later version.
    #
    # This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    # without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    # See the GNU General Public License for more details.
    # 
    # You should have received a copy of the GNU General Public License along with this program;
    # if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
    # Boston, MA  02110-1301, USA.
    #=============================================================================================

'''

#================================
#Imports
#================================
from MDAnalysis import *
from MDAnalysis.analysis.distances import *
import numpy as np
import argparse


#=============================================================================================
# distances class definition
#=============================================================================================

class distances:
    """
    This class extracts distances of all kinds from a trajectory and topology file\n
    prerequesit modules are: MDAnalysis
    Please cite:
    [1] N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.\n 
    MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.\n 
    J. Comput. Chem. 32 (2011), 2319–2327, doi:10.1002/jcc.21787
    
    """
    def __init__(self):
        '''
        Constructor, does not take any arguments
        '''
        pass
    
        
    def returnCAlphaDistanceTrajectory(self,_trajectoryFile, _topologyFile, skip = 2):
        """This function returns next to nearest neighbour c-alpha distances
    
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc
           _topologyFile (str): topologyFilename.gro
    
        Returns:
           arrayList  distances array for the whole trajectory [traj_length][num_c_alpha_distances]
    
    
        Useage
    
        >>> distances = returnCAlphaDistances(_trajectoryFile='test.xtc', _topologyFile='test.gro')
        
        """  
        universe = Universe(_topologyFile,_trajectoryFile)
        print "Loading trajectory file %s"%_trajectoryFile
        protein = universe.selectAtoms("protein") #Selecting only protein atoms
        nframes = universe.trajectory.numframes - 1
        calphas = universe.selectAtoms("name CA")
        atomnums = calphas.atoms.indices()
    
        distance_selections = []
        for k in range(len(atomnums)):
            idist = universe.selectAtoms("bynum %i"%atomnums[k])
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
            universe.trajectory.next()
        print 'distance Trajectory has been obtained.'
        return distanceListTrajectory
    
    def ReturnHeavyAtromDistanceTrajectory(self):
        print 'this method is not implemented yet'
        
    def ReturnCustomDistances(self):
        print 'this method is not implemented yet'
    
        
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
    d.thisIsAPublicMethod()
    distances = d.returnCAlphaDistanceTrajectory(trajectoryFile, topologyFile, 2)
    print distances