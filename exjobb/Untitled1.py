
# coding: utf-8

# In[63]:


import Bio
from Bio.PDB import *
import sys
import re
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import numpy as np

from mpl_toolkits.mplot3d import Axes3D


# In[ ]:


def calc_residue_dist(residue1, residue2,layers) :
    """Returns the C-alpha distance between two residues"""
    for atom1 in residue1:
        for atom2 in residue2:
            
                
                
                
            distance=atom1-atom2    

            if distance>17:#If the residues is to far from eachother there is no need to calculate the distance for the rest of the atoms
                return (layers)
            if distance >8:#Go to the next atom
                continue
            if distance<=0.5:
                layers.setdefault('VDW', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])

                
            if 'H' not in atom1.get_id() and 'H' not in atom2.get_id(): # DUBBELKOLLA ATT DET INTE FINNS H I ANDRA!!!!
            #if atom1.get_id() != 'H' and atom2.get_id() != 'H':    
                if atom1.get_name()=='CA' and atom2.get_name()=='CA' and distance<=8:
                    layers.setdefault('CA-CA', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                if distance<=5: #Heavy atoms distance
                    layers.setdefault('H-A', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
    return (layers)



def calc_dist_matrix(chain_one, chain_two,layers) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            layers = calc_residue_dist(residue_one, residue_two,layers)
    return (layers)


# In[ ]:


def main():
    filename_pdb = '/home/joakim/Downloads/5bk1.pdb'#'/home/joakim/Downloads/2HIY_A.pdb' #'/home/joakim/Downloads/D1A2K-a0a-merged.pdb'
    try: 
        PDBobj = PDBParser()
        structure = PDBobj.get_structure(filename_pdb, filename_pdb)
        model = structure[0]
    except IOError: 
        print 'IO Error', filename_pdb       
        while 'true':
            input1=raw_input("Error parsing PDB file! Continue(y/n)...")
            if input1.lower()=='y':
                break
            elif input1.lower()=='n':
                sys.exit()
    
    chain_used=[]
    layers={}
    layers.setdefault('VDW', []).append(['Check'])
    for chain1 in model:
        for chain2 in model:
            if chain1!=chain2 and chain2 not in chain_used:
                chain_used.append(chain1)
                layers = calc_dist_matrix(chain1, chain2,layers)
                #contact_map = dist_matrix < 12.0
    
    t=0
    print 'CA-CA'
    for i in layers['CA-CA']:
        print i
        t=t+1
    
    g=0
    print 'H-A'
    for i in layers['H-A']:
        print i
        g=g+1
    
    d=0

    print 'VDW'
    for i in layers['VDW']:
        print i
        d=d+1

    print d
    print g
    print t
#    for i in contact_map:
#        for g in i:
#            if g!= False:
#                print g

    
if __name__ == '__main__':
  main()

