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
import os
from MDAnalysis import collection

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
       
        
    def getDihedralTrajectory(self,_structure, _trajectory, _atomSelectionArray, _stop=-1):
        '''This function returns the dihedral angle between four specified atoms over the whole trajectory\n
        The preferred input is in dcd format. Otherwise the first trajectory frame is missing\n
        ToDo: fix issues with xtc input trajectory that dihedrals for all frames are computed.
            
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
           _topologyFile (str): topologyFilename.gro
           _atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
           _stop (int): Default =-1 (i.e. all frames)
    
        Returns:
           dihedralTraj  trajectory containing dihedrals between the four atoms specified in the _atomSelectionArray
    
    
        Useage:
    
        >>> dihedralTraj = getDihedralTrajectory(_trajectoryFile='test.xtc/dcd', _topologyFile='test.gro',atoms)
        
        '''  
        
        u = Universe(_structure,_trajectory)
        
        if self._extensionDcD(_trajectory):
            collection.clear()
            selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
            selection = u.selectAtoms(selectionString)
            collection.addTimeseries(Timeseries.Dihedral(selection))
            if _stop == -1:
                _stop = u.trajectory.numframes-1
            collection.compute(u.trajectory, stop=_stop)
            return collection[0]
           
        
               #else
        else:
            print 'Warning.....first frame is missing............................'
            print 'You have not provided a dcd file, iterating through the trajectory will be slow. '
            print 'For improved preformance please convert your input trajectory to dcd format.'
            print 'Good luck with the compuatation'
            frames = u.trajectory.numframes
            dihedrals = []
            u.trajectory.rewind()
            for t in range(frames-1):
               
                selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
                selection = u.selectAtoms(selectionString)
                dihedral = selection.dihedral()
                dihedrals.append(dihedral)
                u.trajectory.next()
            selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
            selection = u.selectAtoms(selectionString)
            dihedral = selection.dihedral()
            dihedrals.append(dihedral)
            dihedrals = np.array(dihedrals)
            dihedrals = np.radians(dihedrals)
            return dihedrals
        
    def getAngleTrajectory(self,_structure, _trajectory, _atomSelectionArray, _stop =-1):
        '''This function returns the alangle between three specified atoms over the whole trajectory\n
        The preferred input is in dcd format. Otherwise the first trajectory frame is missing\n
        ToDo: fix issues with xtc input trajectory that angles for all frames are computed.
            
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
           _topologyFile (str): topologyFilename.gro
           _atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
           _stop (int): Default =-1 (i.e. all frames)
    
        Returns:
           anglesTraj  trajectory containing angles between the three atoms specified in the _atomSelectionArray
    
    
        Useage:
    
        >>> dihedralTraj = getDihedralTrajectory(_trajectoryFile='test.xtc/dcd', _topologyFile='test.gro',atoms)
        
        '''  
        
        #if filename is a dcd file use fast method
        u = Universe(_structure,_trajectory)
        
        if self._extensionDcD(_trajectory):
            collection.clear()
            selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
            selection = u.selectAtoms(selectionString)
            collection.addTimeseries(Timeseries.Angle(selection))
            if _stop == -1:
                _stop = u.trajectory.numframes-1
            collection.compute(u.trajectory, stop=_stop)
            return collection[0][0]
            
        
        #else
        else:
            print 'Warning.......warning first frame is missing....................................'
            print 'You have not provided a dcd file, iterating through the trajectory will be slow. '
            print 'For improved preformance please convert your input trajectory to dcd format.'
            print 'Good luck with the compuatation'
            frames = u.trajectory.numframes-1
            angles = []
           
            for t in range(frames):
                
                selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
                selection = u.selectAtoms(selectionString)
                angle = selection.angle()
                angles.append(angle)
                u.trajectory.next()
            selectionString =self._generateSelectionString(_atomSelectionArray, "bynum ")
            selection = u.selectAtoms(selectionString)
            angle = selection.angle()
            angles.append(angle)
            angles = np.array(angles)
            angles = np.radians(angles)
            return angles
                
        
    def _generateSelectionString(self, selections, selectionType):
  
        finalSelection = ""
        for s in range(len(selections)-1):
            finalSelection = finalSelection+selectionType+str(selections[s])+" or "
        finalSelection =finalSelection+selectionType+str(selections[-1])
       
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
    a = Angles()
    #angleArray = np.array([1169, 423, 279, 117])
    angleArray = np.array([279,423,1169])
    #theta = a.getDihedralAngle(topologyFile,trajectoryFile,angleArray)
    traj = a.getAngleTrajectory(topologyFile, trajectoryFile, angleArray)
    print traj


        