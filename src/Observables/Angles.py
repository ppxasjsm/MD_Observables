#!/usr/bin/env python

'''
    File:            /home/mi/ppxasjsm/Eclipse/MD_observables/Observables/Angles.py
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

class Angles:
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
    
    def getAngle(self, _structure, _trajectory, _atomSelectionArray):
        '''This function returns the angle between three specified atoms
    
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
           _topologyFile (str): topologyFilename.gro
           _atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
    
        Returns:
           angle  angle between the three atoms specified in the _atomSelectionArray
    
    
        Useage:
    
        >>> angle = getAngle(_trajectoryFile='test.xtc/dcd', _topologyFile='test.gro', atoms)
        
        '''  
        
        u = Universe(_structure,_trajectory)
        selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
        selection = u.selectAtoms(selectionString)
        angle = selection.angle()
        return angle
        
      
        
    def getDihedralAngle(self,_structure, _trajectory, _atomSelectionArray):
        '''This function returns the dihedralangle between four specified atoms
    
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
           _topologyFile (str): topologyFilename.gro
           _atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
    
        Returns:
           angle  angle between the three atoms specified in the _atomSelectionArray
    
    
        Useage:
    
        >>> dihedral = getDihedralAngle(_trajectoryFile='test.xtc/dcd', _topologyFile='test.gro',atoms)
        
        '''  
        u = Universe(_structure,_trajectory)
        selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
        selection = u.selectAtoms(selectionString)
        dihedral = selection.dihedral()
        return dihedral
       
        
    def getDihedralTrajectory(self):
        print 'this method has not been implemented yet'
        
    def getAngleTrajectory(self):
        print 'this method has not been implemented yet'
        
    def _generateSelectionString(self, selections, selectionType):
        print selections
        finalSelection = ""
        for s in range(len(selections)-1):
            finalSelection = finalSelection+selectionType+str(selections[s])+" or "
        finalSelection =finalSelection+selectionType+str(selections[-1])
        print finalSelection
        return finalSelection
       
                
#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
    print 'Hello angle word'
    
#tsting the methods
    parser = argparse.ArgumentParser(description="Feature extractions from xtc files")
    parser.add_argument('-f', help='xtc inputfile');
    parser.add_argument('-c', help='coordinate file in Gromacs .gro format');
    args = parser.parse_args()
 
    topologyFile = args.c
    trajectoryFile = args.f
    a = Angles()
    angleArray = np.array([1169, 423, 279, 117])
    theta = a.getDihedralAngle(topologyFile,trajectoryFile,angleArray)
    print theta

        