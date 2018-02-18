
# coding: utf-8

# In[41]:


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


# In[42]:


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
                #layers_VDW.append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                layers.setdefault('VDW', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])

                
            if 'H' not in atom1.get_id() and 'H' not in atom2.get_id(): # DUBBELKOLLA ATT DET INTE FINNS H I ANDRA!!!!
            #if atom1.get_id() != 'H' and atom2.get_id() != 'H':    
                if atom1.get_name()=='CA' and atom2.get_name()=='CA' and distance<=8: #CA-CA distance
                    #layers_CA.append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                    layers.setdefault('CA-CA', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                if distance<=5: #Heavy atoms distance
                    #layers_HA.append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
                    layers.setdefault('H-A', []).append([str(residue1.get_full_id()[3][1]) +'.'+ residue1.get_resname() +'-'+ str(residue2.get_full_id()[3][1]) +'.'+ residue2.get_resname(),atom1.get_name() +'-'+ atom2.get_name(), distance])
    return (layers)
    #if layers_HA==[] and layers_VDW==[] and layers_CA==[]:
     #   return
    #else:
     #   return ([layers_VDW,layers_CA,layers_HA])

    print(layers_VDW,layers_CA,layers_HA)
    return (layers_VDW,layers_CA,layers_HA)


def calc_dist_matrix(chain_one, chain_two,layers) :
    """Returns a matrix of C-alpha distances between two chains"""

    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            #contact=(calc_residue_dist(residue_one, residue_two,layers))
            layers=calc_residue_dist(residue_one, residue_two,layers)
                           

            #print('contact',contact)
            #if contact!=None:
                #layers.append(contact)
                #print(contact)
                #print('---------------------')
            #layers.append(contact)
            
    return (layers)


# In[43]:


def find_protein_value(pdb_name,value_file):
    import re
    #Re.search(r”([\w.]+)@([\w.]+)”, value_file) 


# In[44]:


def main():
    
    #Dir path
    #args = sys.argv[1:]
    pdb_dir= '/home/joakim/Downloads/models' #str(args[0]) #'D:\Downloads\CnM-dataset\models'

    while 'true':
        input1=raw_input("Save or load ")
        if input1.lower()=='-r':
            args='-r'
            break
        elif input1.lower()=='-s':
            args='-s'
            break
        else:
            break
    
    
    
    
    if not '-r' in args:
    
        pdb_list=[]
        for file in os.listdir(pdb_dir):
            if file.endswith(".pdb"):
                #pdb_list.append(args+'/'+file)
                pdb_list.append([pdb_dir,file])
        #print pdb_list
        #filename_pdb = '/home/joakim/Downloads/*.pdb'#'/home/joakim/Downloads/2HIY_A.pdb' #'/home/joakim/Downloads/D1A2K-a0a-merged.pdb'

        valuefile='/home/joakim/Downloads/CnM.featuresNPqDNZ'#'D:\Downloads\CnM.featuresNPqDNZ'
        value_df=pandas.read_csv(valuefile,delim_whitespace=1)



        #value_file=open('/home/joakim/Downloads/CnM.featuresNPqDNZ', 'r').read().splitlines()
        #print(value_file)

        score=[]
        dist_matrix=[]
        for counter, filepath_pdb in enumerate(pdb_list):
            print(counter)
            print(pdb_list[counter][1])
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

            #manager = multiprocessing.Manager()
            #return_dict = manager.list
            #jobs = []
            #job_count=0

            chain_used=[]
            layers={}
            layers.setdefault('VDW', []).append(['Check'])
            layers.setdefault('CA-CA', []).append(['Check'])
            layers.setdefault('H-A', []).append(['Check'])
            for chain1 in model:
                for chain2 in model:
                    if chain1!=chain2 and chain2 not in chain_used:
                        chain_used.append(chain1) 
                        #layers=[]
                        layers = calc_dist_matrix(chain1, chain2,layers)
                        #p = multiprocessing.Process(target=dist_matrix.append, args=(job_count,return_dict))
                        #jobs.append(p)
                        #p.start()
                        #job_count=job_count+1
                        #print(return_dict)
                        #dist_matrix.append(calc_dist_matrix(chain1, chain2,layers))

                        #contact_map = dist_matrix < 12.0
            if (len(layers['VDW']))>1:
                print('VDW',len(layers['VDW']))
                        
        #print(layers.values())
        
        #for x in dist_matrix[0]:
            #print(x)
        score=np.array(score)


    if '-s' in args and '-r' in args:
        print('Do not use save(-s) and read(-r) at the same time')
        sys.exit()
    elif '-s' in args:
        print('-Saving files...')
        np.savez_compressed('score', score)
        #np.savez_compressed('dist_matrix', dist_matrix)
        print('-Files saved as "dist_matrix.npz","score.npz"')

    elif '-r' in args:
        print('-Loading files...')
        score = np.load('score.npz')
        #dist_matrix = np.load('dist_matrix.npz')
        
        for key,array in score.items():
            score=score[key]
        #for key,array in dist_matrix.items():
            #dist_matrix=dist_matrix[key]
            
        print('-Files loaded')
    
    
    
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

