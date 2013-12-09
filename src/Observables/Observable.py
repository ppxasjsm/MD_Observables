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
        
    def extensionDcD(self, trajectoryFile):
        '''This function asserts whether the input trajectory file is in dcd format or not.
            
        Args:
           trajectoryFile (str): trajectoryFilename.xtc/trajectoryFilename.dcd
          
    
        Returns:
           boolean  true, if trajectory file is in dcd format
    
    
        Useage:
    
        >>> if == self.extensionDcD('test.dcd) 
        
        '''  
        fileName, fileExtension = os.path.splitext(trajectoryFile)
        if fileExtension == '.dcd':
            return True
        else:
            return False
        
    def generate_selection_string(self, selections, selectionType):
  
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
        
        