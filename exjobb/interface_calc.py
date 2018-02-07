
# coding: utf-8

# In[38]:


import Bio
from Bio.PDB import *
import sys
import re
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import numpy as np
import os

import pandas
import multiprocessing


#This is for removing Bio.python warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


# In[39]:


def calc_residue_dist(residue1, residue2,layers) :
    """Returns the C-alpha distance between two residues"""
    layers_VDW=[]
    layers_CA=[]
    layers_HA=[]
    for atom1 in residue1:
        for atom2 in residue2:
            
            distance=atom1-atom2    

            if distance>20:#If the residues is to far from eachother there is no need to calculate the distance for the rest of the atoms
                return (layers)
            if distance >8:#Go to the next atom
                continue
            if distance<=0.5:
                layers_VDW.append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                #layers.setdefault('VDW', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])

                
            if 'H' not in atom1.get_id() and 'H' not in atom2.get_id(): # DUBBELKOLLA ATT DET INTE FINNS H I ANDRA!!!!
            #if atom1.get_id() != 'H' and atom2.get_id() != 'H':    
                if atom1.get_name()=='CA' and atom2.get_name()=='CA' and distance<=8: #CA-CA distance
                    layers_CA.append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                    #layers.setdefault('CA-CA', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                if distance<=5: #Heavy atoms distance
                    layers_HA.append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                    #layers.setdefault('H-A', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
    #return (layers)
    return ([layers_VDW,layers_CA,layers_HA])


def calc_dist_matrix(chain_one, chain_two,layers) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            layers.append(calc_residue_dist(residue_one, residue_two,layers))
            
    return (layers)


# In[40]:


def find_protein_value(pdb_name,value_file):
    import re
    #Re.search(r”([\w.]+)@([\w.]+)”, value_file) 


# In[41]:


def main():
    
    #Dir path
    #args = sys.argv[1:]
    args='D:\Downloads\CnM-dataset\models' #'/home/joakim/Downloads/models' #str(args[0])
    
    pdb_list=[]
    for file in os.listdir(args):
        if file.endswith(".pdb"):
            #pdb_list.append(args+'/'+file)
            pdb_list.append([args,file])
    #print pdb_list
    #filename_pdb = '/home/joakim/Downloads/*.pdb'#'/home/joakim/Downloads/2HIY_A.pdb' #'/home/joakim/Downloads/D1A2K-a0a-merged.pdb'
    
    valuefile='D:\Downloads\CnM.featuresNPqDNZ'#'/home/joakim/Downloads/CnM.featuresNPqDNZ'
    value_df=pandas.read_csv(valuefile,delim_whitespace=1)
    
    
    
    #value_file=open('/home/joakim/Downloads/CnM.featuresNPqDNZ', 'r').read().splitlines()
    #print(value_file)
    
    score=[]
    dist_matrix=[]
    for counter, filepath_pdb in enumerate(pdb_list):
        print(counter)
        filename_pdb=os.path.join(pdb_list[counter][0],pdb_list[counter][1])

        try: 
            PDBobj = PDBParser()
            structure = PDBobj.get_structure(filename_pdb, filename_pdb)
            model = structure[0]
        except IOError: 
            print('IO Error', filename_pdb)      
            while 'true':
                input1=raw_input("Error parsing PDB file! Continue(y/n)...")
                if input1.lower()=='y':
                    continue
                elif input1.lower()=='n':
                    sys.exit()

        
        score.append(value_df.loc[((value_df['#'] == pdb_list[counter][1][:-11]),'CPscore')].values[0])
                    
        chain_used=[]
        #layers={}
        #layers.setdefault('VDW', []).append(['Check'])
        #layers.setdefault('CA-CA', []).append(['Check'])
        #layers.setdefault('H-A', []).append(['Check'])
        for chain1 in model:
            for chain2 in model:
                if chain1!=chain2 and chain2 not in chain_used:
                    chain_used.append(chain1)
                    layers=[]
                    #layers = calc_dist_matrix(chain1, chain2,layers)
                    p = multiprocessing.Process(target=dist_matrix.append)
                    p.start()
                    #dist_matrix.append(calc_dist_matrix(chain1, chain2,layers))
                    #contact_map = dist_matrix < 12.0

    print('hej')
    #for counter,i in
    score=np.array(score)
    
    
    print(score[1])
    print(score.shape)
    #dist_matrix=np.array(dist_matrix)
    #print(dist_matrix.shape)
    
     #   t=0
     #   for i in layers['CA-CA']:
       #     t=t+1

     #   g=0
     #   for i in layers['H-A']:
      #      g=g+1

        #d=0

        #for i in layers['VDW']:
         #   d=d+1

        #print d-1
        #print g-1
        #print t-1
    #    for i in contact_map:
    #        for g in i:
    #            if g!= False:
    #                print g

    
if __name__ == '__main__':
    main()

