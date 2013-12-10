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
from MDAnalysis import collection
from Observable import observable

class angles(observable):
    '''
    classdocs
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
    
    def get_angle(self, atomSelectionArray):
        '''This function returns the angle between three specified atoms
    
        Args:
           _trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
           _topologyFile (str): topologyFilename.gro
           _atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
    
        Returns:
           angle  angle between the three atoms specified in the _atomSelectionArray
    
    
        Useage:
    
        >>> angle = get_angle(_trajectoryFile='test.xtc/dcd', _topologyFile='test.gro', atoms)
        
        '''  
        
       
        selectionString =self.generate_selection_string(atomSelectionArray, "bynum ")
        selection = self._u.selectAtoms(selectionString)
        angle = selection.angle()
        return angle
        
      
        
    def get_dihedral_angle(self,atomSelectionArray):
        '''This function returns the dihedralangle between four specified atoms
    
        Args:
           atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
    
        Returns:
           angle  angle between the three atoms specified in the atomSelectionArray
    
    
        Useage:
    
        >>> dihedral = get_dihedral_angle(atoms)
        
        '''  
    
        selectionString =self.generate_selection_string(atomSelectionArray, "bynum ")
        selection = self._u.selectAtoms(selectionString)
        dihedral = selection.dihedral()
        return dihedral
       
        
    def get_dihedral_trajectory(self, atomSelectionArray, _stop=-1):
        '''This function returns the dihedral angle between four specified atoms over the whole trajectory\n
        The preferred input is in dcd format. Otherwise the first trajectory frame is missing\n
        ToDo: fix issues with xtc input trajectory that dihedrals for all frames are computed.
            
        Args:
           atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
           _stop (int): Default =-1 (i.e. all frames)
    
        Returns:
           dihedralTraj  trajectory containing dihedrals between the four atoms specified in the atomSelectionArray
    
    
        Useage:
    
        >>> dihedralTraj = get_dihedral_trajectory(_trajectoryFile='test.xtc/dcd', _topologyFile='test.gro',atoms)
        
        '''  
        
      
        
        if self.extension_DcD(self._trajectoryFile):
            collection.clear()
            selectionString =self.generate_selection_string(atomSelectionArray, "bynum ")
            selection = self._u.selectAtoms(selectionString)
            collection.addTimeseries(Timeseries.Dihedral(selection))
            if _stop == -1:
                _stop = self._u.trajectory.numframes-1
            collection.compute(self._u.trajectory, stop=_stop)
            return collection[0]
           
        
               #else
        else:
            print "Using frame iteration and not collections module"
            frames = self._u.trajectory.numframes
            dihedrals = []
            selectionString =self.generate_selection_string(atomSelectionArray, "bynum ")
            selection = self._u.selectAtoms(selectionString)
            count = 0;
            for ts in self._u.trajectory:    
                if count ==_stop:
                    break           
               
                dihedral = selection.dihedral()
                dihedrals.append(dihedral)
                count = count+1
            dihedrals.append(dihedral)
            dihedrals = np.array(dihedrals)
            dihedrals = np.radians(dihedrals)
            return dihedrals
        
    def get_angle_trajectory(self, atomSelectionArray, _stop =-1):
        '''This function returns the alangle between three specified atoms over the whole trajectory\n
        The preferred input is in dcd format. Otherwise the first trajectory frame is missing\n
        ToDo: fix issues with xtc input trajectory that angles for all frames are computed.
            
        Args:
           atomSelectionArray(numpy array): list of 3 atom indices as in the structure file
           _stop (int): Default =-1 (i.e. all frames)
    
        Returns:
           anglesTraj  trajectory containing angles between the three atoms specified in the atomSelectionArray
    
    
        Useage:
    
        >>> dihedralTraj = get_dihedral_trajectory(_trajectoryFile='test.xtc/dcd', _topologyFile='test.gro',atoms)
        
        '''  
        
        #if filename is a dcd file use fast method
      
        
        if self.extension_DcD(self._trajectoryFile):
            collection.clear()
            selectionString =self.generate_selection_string(atomSelectionArray, "bynum ")
            selection = self._u.selectAtoms(selectionString)
            collection.addTimeseries(Timeseries.Angle(selection))
            if _stop == -1:
                _stop = self._u.trajectory.numframes-1
            collection.compute(self._u.trajectory, stop=_stop)
            return collection[0][0]
            
        
        #else
        else:
            print "Using frame iteration and not collections module"
            frames = self._u.trajectory.numframes-1
            angles = []
            selectionString =self.generate_selection_string(atomSelectionArray, "bynum ")
            selection = self._u.selectAtoms(selectionString)
            count =0
            for ts in self._u.trajectory:
                if count == _stop:
                    break

                angle = selection.angle()
                angles.append(angle)
                count = count+1
            angles = np.array(angles)
            angles = np.radians(angles)
            return angles
                

                
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
    a = angles(topologyFile, trajectoryFile)
    
    angleArray = np.array([279,423,1169])
    print a.get_angle(angleArray)
    dihedralArray = np.array([1169, 423, 279, 117])
    print a.get_dihedral_angle(dihedralArray)
    print a.get_angle_trajectory(angleArray , _stop=5)
    print a.get_dihedral_trajectory(dihedralArray, _stop=5)
    
    #theta = a.get_dihedral_angle(topologyFile,trajectoryFile,angleArray)
    #traj = a.get_angle_trajectory(topologyFile, trajectoryFile, angleArray)
    #print traj


        