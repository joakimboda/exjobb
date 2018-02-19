
# coding: utf-8

# In[1]:


import Bio
from Bio.PDB import *
import sys
import re
#import matplotlib.pyplot as plt
#from matplotlib import cm
import math
import numpy as np
import os

import pandas
import multiprocessing

import shutil


#This is for removing Bio.python warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


# Calaculates an alpha complex from coordinates for a choosen atomtype and puts them into an dict. Dict key = name of the protein

# In[2]:


def calculate_alpha_complex(BarCollection, name, alphaprot, cwd):
    
    for key,coords in alphaprot.items():
        Bars = []
        tmpoutfile =open(cwd+'/pts_dir/tmp.csv', 'w')
        tmpoutfile.write('x1,x2,x3\n')
        
        for coord in coords:
            tmpoutfile.write(coord+'\n')
        tmpoutfile.close()
        os.system('Rscript PH_Alpha.R '+ cwd+'/pts_dir/tmp.csv '+ cwd+'/pts_dir'+'/tmp.out')
        
        tmpinfile = open(cwd+'/pts_dir/tmp.out')
        lines = tmpinfile.read().splitlines()
        
        for line in lines[1:]:
            a,b,c,d = line.split()
            if d!='Inf':
                if float(d)-float(c) >= 0.01: #If the bar is to small we don't need to use it (more data for nothing)
                    Bars.append([int(b), float(c), float(d)])
        BarCollection[name] = Bars    


# In[3]:


def calc_residue_dist(residue1, residue2, outfile, atom1_used, atom2_used, alphaprot) :

    for atom1 in residue1:
        for atom2 in residue2:
            
            distance=atom1-atom2    

            if distance>20:#If the residues is to far from eachother there is no need to calculate the distance for the rest of the atoms
                return(alphaprot)
            if distance >12:#Go to the next atom
                continue
                
            if 'H' not in atom1.get_id() and 'H' not in atom2.get_id(): # DUBBELKOLLA ATT DET INTE FINNS H I ANDRA!!!!
            #if atom1.get_id() != 'H' and atom2.get_id() != 'H':    
                
                a1_serial=str(atom1.serial_number)
                a2_serial=str(atom1.serial_number)
                if a1_serial not in atom1_used:
                    atom1_used.append(a1_serial)
                    x=str(atom1.get_coord()[0])
                    y=str(atom1.get_coord()[1])
                    z=str(atom1.get_coord()[2])
                    if 'C' in atom1.get_name():# and
                        outfile.write('0'+' '+'1'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('C', []).append(x+','+y+','+z)
                    elif 'O' in atom1.get_name():
                        outfile.write('0'+' '+'2'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('O', []).append(x+','+y+','+z)
                    elif 'N' in atom1.get_name():
                        outfile.write('0'+' '+'3'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('N', []).append(x+','+y+','+z)
                    elif 'S' in atom1.get_name():
                        outfile.write('0'+' '+'4'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('S', []).append(x+','+y+','+z)
                if a2_serial not in atom2_used:
                    atom2_used.append(a2_serial)
                    x=str(atom2.get_coord()[0])
                    y=str(atom2.get_coord()[1])
                    z=str(atom2.get_coord()[2])
                    if 'C' in atom1.get_name():
                        outfile.write('1'+' '+'1'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('C', []).append(x+','+y+','+z)
                    elif 'O' in atom1.get_name():
                        outfile.write('1'+' '+'2'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('O', []).append(x+','+y+','+z)
                    elif 'N' in atom1.get_name():
                        outfile.write('1'+' '+'3'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('N', []).append(x+','+y+','+z)
                    elif 'S' in atom1.get_name():
                        outfile.write('1'+' '+'4'+' '+x+' '+y+' '+z+'\n')
                        alphaprot.setdefault('S', []).append(x+','+y+','+z)
                
    return(alphaprot)




def calc_dist_matrix(chain_one, chain_two,outfile,alphaprot) :
    
    atom1_used=[]
    atom2_used=[]
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            
            alphaprot=calc_residue_dist(residue_one, residue_two, outfile, atom1_used, atom2_used, alphaprot)
                                    
    return(alphaprot)


# In[4]:


def main():
    
    cwd = os.getcwd()
    #Dir path
    #args = sys.argv[1:]
    
    pdb_dir= cwd+'/models' #'D:\Downloads\CnM-dataset\models'#'/home/joakim/Downloads/models' #str(args[0]) 
    
    
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
                pdb_list.append([pdb_dir,file])
      
        valuefile= cwd+'/CnM.featuresNPqDNZ'#'/home/joakim/Downloads/CnM.featuresNPqDNZ'#'D:\Downloads\CnM.featuresNPqDNZ'
        value_df=pandas.read_csv(valuefile,delim_whitespace=1)

        #Just a check to se if the data_set dir exist, in this case it removes it and all the files in it and make a new one
        if os.path.exists(cwd + "/pts_dir"):
            shutil.rmtree(cwd + "/pts_dir")
            os.makedirs(cwd + "/pts_dir")
        else:
            os.makedirs(cwd + "/pts_dir")
        


        
        score=[]
        dist_matrix=[]
        BarCollection = {}
        for counter, filepath_pdb in enumerate(pdb_list):
            name=pdb_list[counter][1][:-11]
            
            working_file=open(cwd + '/pts_dir/working_file.txt', 'a')
            working_file.write(name+'\n') #Writes what pdb have been worked on in working_file.txt
            working_file.close()
            
            print(counter+1)            
            print(name)
            
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

           
            score.append(value_df.loc[((value_df['#'] == name),'CPscore')].values[0])



            chain_used=[]
            outfile=open(cwd + "/pts_dir/"+name + '.pts', 'w')
            for chain1 in model:
                for chain2 in model:
                    alphaprot={} #Used to make temp file for alpha complex calculations, its a dict of coord for every atomtype used
                    if chain1!=chain2 and chain2 not in chain_used:
                        chain_used.append(chain1)
                        alphaprot=calc_dist_matrix(chain1, chain2, outfile,alphaprot) #Makes coord files for matlab and returns a dict for alpha complex
                        calculate_alpha_complex(BarCollection, name, alphaprot, cwd) #Calculates alpha complex
                        
            outfile.close()
    
    

    if '-s' in args and '-r' in args:
        print('Do not use save(-s) and read(-r) at the same time')
        sys.exit()
    elif '-s' in args:
        print('-Saving files...')
        #np.savez_compressed('score', score)
        np.savez_compressed('Alpha_Complex', BarCollection)
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
    
    os.remove(cwd+'/pts_dir/tmp.out')
    os.remove(cwd+'/pts_dir/tmp.csv')
    
   # print(score[1])
   # print(score.shape)
    line = 'matlab -nodisplay -nodesktop -nosplash -r '+cwd+'/bar.m'
    os.system(line)
    
if __name__ == '__main__':
    main()

