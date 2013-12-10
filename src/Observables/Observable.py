#!/usr/bin/env python

'''
    File:            /home/mi/ppxasjsm/Eclipse/MD_Observables/src/Observables/Observable.py
    Description:
    
    Author:          Antonia Mey    
    Email:           antonia.mey@fu-berlin.de      
    Date:            09.12.2013
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
#from MDAnalysis import *
import MDAnalysis
import os

class observable(object):
    '''
    Metaclass that allows the extraction of specific MD observables
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
        self._u =  MDAnalysis.Universe(topologyFile,trajectoryFile)
        
    def extension_DcD(self, trajectoryFile):
        '''This function asserts whether the input trajectory file is in dcd format or not.
            
        Args:
           trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
          
    
        Returns:
           boolean  true, if trajectory file is in dcd format
    
    
        Useage:
    
        >>> if == self.extension_DcD('test.dcd') 
        
        '''  
        fileName, fileExtension = os.path.splitext(trajectoryFile)
        if fileExtension == '.dcd':
            return True
        else:
            return False
        
    def generate_selection_string(self, selections, selectionType):
        ''' This method generates a selection string based on keywords for selecting atoms and residues\n
        for MDAnalysis
        
        Args: 
            -selections (int array)     numpy array containing either residue id's or atom numbers
            -selectionType (string)     either 'resid' or 'bynum' 
            
        Returns:
            slectionString (string)    string that can be used with u.selectAtoms()
            
        
        Useage:
        
        >>> selection = self.generate_selection_string(selection, 'bynum')
        
        '''
        finalSelection = ""
        
        if len(selections)>2:
            for s in range(len(selections)-1):
                finalSelection = finalSelection+selectionType+str(selections[s])+" or "
            finalSelection =finalSelection+selectionType+str(selections[-1])
        else:
            finalSelection = finalSelection+selectionType+str(selections[0])+" or "
            finalSelection =finalSelection+selectionType+str(selections[1])
   
        return finalSelection
        
        