'''
Created on 09.12.2013

@author: ppxasjsm
'''
from MDAnalysis import *
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
        self._u =  Universe(topologyFile,trajectoryFile)
        
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
        
        