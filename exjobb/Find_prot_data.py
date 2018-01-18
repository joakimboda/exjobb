
# coding: utf-8

# In[25]:


import Bio
from Bio.PDB import *
import sys
import re
import matplotlib.pyplot as plt
from matplotlib import cm
import math
 
         
def find_structure_params(structure):
     
    model = structure[0]
    chain = model['A'] 
    input_=[]
    count=0
    #for atom in structure.get_atoms():  If I want all atoms in a structure, depends if I want the residue
    for residue in chain:
        if residue.get_full_id()[3][0]==' ':
            for atom in residue:
                print atom.get_id()
                count=count +1
            #print residue.get_full_id()
    print count
            
           #rr_all_distance.append([residue1, residue2, distance])
        #   structure_params.append([index1, index2, distance]) #tabort residue1 o reside2
 
    return




# In[26]:


def main():
    filename_pdb = '/home/joakim/Downloads/2HIY_A.pdb'
    
    PDBobj = PDBParser()
    structure = PDBobj.get_structure(filename_pdb, filename_pdb)
    find_structure_params(structure)
    
    
    
    
if __name__ == '__main__':
  main()

