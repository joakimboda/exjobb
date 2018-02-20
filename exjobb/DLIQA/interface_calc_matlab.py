
# coding: utf-8

# In[53]:


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

ele='D'
pro=['C','N','O','S','H']
lig=['C','N','O','S','P','F','Cl','Br','I','H']

try:
    list_index_p=pro.index(atom1.get_name())
    if list_index_p==4:

    else:
        ele=str(list_index_p+1)
        outfile.write('0'+' '+ele+' '+x+' '+y+' '+z+'\n')
except:
    'odd'
    
try:
    list_index_l=pro.index(atom2.get_name())
    if list_index_p==9:
        
    else:
        ele=str(list_index_l+1)
        outfile.write('0'+' '+ele+' '+x+' '+y+' '+z+'\n')
except:
    'odd'

#print filter(lambda x: ele in x, g)


# In[38]:


def GenerateFeature_alpha_2D(protein_name, working_dir):
    
    
#if 'true':
    thr = 12.0
    rs = 0.1
    lth = int(thr/rs)
    small = 0.01

    OrderedName = np.load('temp/PerformanceOrderAlphaHand.npy')
    
    d=0
    for x in OrderedName:
        d=d+1
        print x
    print d
    X = np.zeros([16, lth, 128], float)

    InFile = open('temp/1a8i_protein_alpha.pkl')
    BarCollection = pickle.load(InFile)
    for j in range(len(OrderedName)):
        plname = OrderedName[j]; pname, lname = plname.split('_');
        f_p_i_j = np.zeros([8,lth], float)
        f_pl_i_j = np.zeros([8,lth], float)
        if 'pro_'+pname in BarCollection.keys():
            Bars = BarCollection['pro_'+pname]
            for Bar in Bars:
                if Bar[1] >= thr: continue
                if Bar[2] >= thr: Bar[2] = thr-0.00001
                bid = min(int(np.floor(Bar[1]/rs)), lth-1)
                did = min(int(np.floor(Bar[2]/rs)), lth-1)
                if Bar[0] == 0:
                    f_p_i_j[0,did] += 1.0
                    f_p_i_j[1,bid:did+1] += 1.0
                elif Bar[0] == 1:
                    f_p_i_j[2,bid] += 1.0
                    f_p_i_j[3,did] += 1.0
                    f_p_i_j[4,bid:did+1] += 1.0
                elif Bar[0] == 2:
                    f_p_i_j[5,bid] += 1.0
                    f_p_i_j[6,did] += 1.0
                    f_p_i_j[7,bid:did+1] += 1.0
        if 'com_'+plname in BarCollection.keys():
            Bars = BarCollection['com_'+plname]
            for Bar in Bars:
                if Bar[1] >= thr: continue
                if Bar[2] >= thr: Bar[2] = thr-0.00001
                bid = min(int(np.floor(Bar[1]/rs)), lth-1)
                did = min(int(np.floor(Bar[2]/rs)), lth-1)
                if Bar[0] == 0:
                    f_pl_i_j[0,did] += 1.0
                    f_pl_i_j[1,bid:did+1] += 1.0
                elif Bar[0] == 1:
                    f_pl_i_j[2,bid] += 1.0
                    f_pl_i_j[3,did] += 1.0
                    f_pl_i_j[4,bid:did+1] += 1.0
                elif Bar[0] == 2:
                    f_pl_i_j[5,bid] += 1.0
                    f_pl_i_j[6,did] += 1.0
                    f_pl_i_j[7,bid:did+1] += 1.0
        f_df_i_j = f_pl_i_j[:,:] - f_p_i_j[:,:]
        X[0:8,:,j] = f_pl_i_j[:,:]; X[8:16,:,j] = f_df_i_j[:,:]

    OutFile = open(working_dir+'/'+protein_name+'_feature_complex_alpha_2DCNN.npy', 'w')
    np.save(OutFile, X)
    OutFile.close()


# In[39]:


def element_classification(outfile,comp) :
    
    ProEleCollection = [['C'],['N'],['O'],['S'],['C','N'],['C','O'],['N','O'],['C','N','O'],['C','N','O','S']]
    LigEleCollection = [['C'],['N'],['O'],['S'],['C','N'],['C','O'],['C','S'],['N','O'],['N','S'],['O','S'],['C','N','O'],['C','N','S'],['C','O','S'],['N','O','S'],['C','N','O','S'],['C','N','O','S','P','F','Cl','Br','I']]
    
    
    PRO = comp['PROT']; LIG = comp['LIG'];   

    for epcount,ep in enumerate(ProEleCollection):
        for elcount,el in enumerate(LigEleCollection):
            
            propts = []
            for a in range(len(PRO)):
                if PRO[a][0].replace(" ","") in ep:
                    propts.append([PRO[a][0], PRO[a][1], PRO[a][2], PRO[a][3] ])
            ligpts = []
            for a in range(len(LIG)):
                if LIG[a][0].replace(" ","") in el:
                    ligpts.append([LIG[a][0], LIG[a][1], LIG[a][2], LIG[a][3] ])
            if len(propts) + len(ligpts) > 3:
                pname = ''
                for eep in ep:
                    pname = pname + eep
                lname = ''
                for eel in el:
                    lname = lname + eel
                name = 'com_'+pname+'_'+lname
                pt = propts + ligpts
                
                if name=='com_O_OS':
                    print pt

            
    return       
            
            
            
            
            
           
        
def element_clasdasdasdaassification(outfile,comp) :            
            
    ecombnr=str(epcount*100+elcount)
    ecombdictkey=str(ep)+'-'+str(el)
                    
        
    if atom1.get_name() in ep and a1_e_serial not in atom1_used:
        atom1_used.append(a1_e_serial)
        x=str(atom1.get_coord()[0])
        y=str(atom1.get_coord()[1])
        z=str(atom1.get_coord()[2])
        outfile.write('0 '+ecombnr+' '+x+' '+y+' '+z+'\n')
        alphaprot.setdefault(ecombdictkey, []).append(x+','+y+','+z)
    if atom2.get_name() in el and a2_e_serial not in atom2_used:
        atom2_used.append(a2_e_serial)
        x=str(atom2.get_coord()[0])
        y=str(atom2.get_coord()[1])
        z=str(atom2.get_coord()[2])
        outfile.write('1 '+ecombnr+' '+x+' '+y+' '+z+'\n')
        alphaprot.setdefault(ecombdictkey, []).append(x+','+y+','+z)
                
    return(alphaprot)  


# Calaculates an alpha complex from coordinates for a choosen atomtype and puts them into an dict. Dict key = name of the protein

# In[40]:


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


# In[41]:


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
                
                a1_serial=str(atom1.serial_number)+atom1.get_name()+atom2.get_name()
                a2_serial=str(atom1.serial_number)+atom1.get_name()+atom2.get_name()
                
                atom_r1=atom1.get_name()
                if 'C' in atom_r1:
                    ele=1
                elif 'H' in atom_r1: #Only used for charged rips calculation
                    ele=5
                    continue
                elif 'O' in atom_r1:
                    ele=3
                elif 'N' in atom_r1:
                    ele=2
                elif 'S' in atom_r1:
                    ele=4
                else:
                    continue
                
                atom_r2=atom2.get_name() 
                if 'C' in atom_r2:
                    ele=1
                elif 'H' in atom_r2: #Only used for charged rips calculation
                    ele=10
                    continue
                elif 'O' in atom_r2:
                    ele=3
                elif 'N' in atom_r2:
                    ele=2
                elif 'S' in atom_r2:
                    ele=4
                elif 'P' in atom_r2:
                    ele=5
                elif 'F' in atom_r2:
                    ele=6
                elif 'Cl' in atom_r2:
                    ele=7
                elif 'Br' in atom_r2:
                    ele=8
                elif 'I' in atom_r2:
                    ele=9
                else:
                    continue
                
                
                
                
                if a1_serial not in atom1_used:
                    atom1_used.append(a1_serial)
                    x=str(atom1.get_coord()[0])
                    y=str(atom1.get_coord()[1])
                    z=str(atom1.get_coord()[2])
                    #outfile.write('0'+' '+'1'+' '+x+' '+y+' '+z+'\n')
                    #alphaprot.setdefault('PROT', []).append(atom1.get_name()+','+atom2.get_name()+','+x+','+y+','+z)
                    alphaprot.setdefault('PROT', []).append([atom1.get_name(),x,y,z])
                    

   
                if a2_serial not in atom2_used:
                    atom2_used.append(a2_serial)
                    x=str(atom2.get_coord()[0])
                    y=str(atom2.get_coord()[1])
                    z=str(atom2.get_coord()[2])
                    #outfile.write('1'+' '+'1'+' '+x+' '+y+' '+z+'\n')
                    #alphaprot.setdefault('LIG', []).append(atom1.get_name()+','+atom2.get_name()+','+x+','+y+','+z)
                    alphaprot.setdefault('LIG', []).append([atom2.get_name(),x,y,z])
                    
    return(alphaprot)




def calc_dist_matrix(chain_one, chain_two,outfile,alphaprot) :
    
    atom1_used=[]
    atom2_used=[]
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            
            alphaprot=calc_residue_dist(residue_one, residue_two, outfile, atom1_used, atom2_used, alphaprot)
                                    
    return(alphaprot)


# In[42]:


def main():
    
    
    
    #python DLIQA.py -location_of_dataset -crossvalidation_file -location_of_value_file -PQR(1)
    PQR=0
    sys_args = sys.argv[1:]
    if sys_args[0]:
        pdb_dir=sys_args[0]
        sys_args = sys.argv[1:]
    if sys_args[0]:
        valuefile=sys_args[0] #################################
        sys_args = sys.argv[1:]
    if sys_args[0]:
        valuefile=sys_args[0]
        sys_args = sys.argv[1:]
    
    
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
    
    
    cwd = os.getcwd()
    
    
    pdb_dir= cwd+'/models' #'D:\Downloads\CnM-dataset\models'#'/home/joakim/Downloads/models' #str(args[0]) 
    
    

    
    
    
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

            #If PQR=1 it will try to run pdb2pqr
            if PQR==1:
                outfile=filename_pdb[:-3]+'pqr'
                line = '/software/apps/python/2.7.13/anaconda-5.0.0.1/bin/python /proj/wallner/apps/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py --ff=amber '+ filename_pdb +' '+ outfile
                os.system(line)
                filename_pdb=outfile
                
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
                
        
            #Finds the score of the file !!!!!!change to an dict
            score.append(value_df.loc[((value_df['#'] == name),'CPscore')].values[0])



            chain_used=[]
            outfile=open(cwd + "/pts_dir/"+name + '.pts', 'w')
            for chain1 in model:
                for chain2 in model:
                    alphaprot={} #Used to make temp file for alpha complex calculations, its a dict of coord for every atomtype used
                    comp={}
                    if chain1!=chain2 and chain2 not in chain_used:
                        chain_used.append(chain1)
                        comp=calc_dist_matrix(chain1, chain2, outfile,comp) #Makes coord files for matlab and returns a dict for alpha complex
                        #element_classification(outfile,comp)
                        #calculate_alpha_complex(BarCollection, name, alphaprot, cwd) #Calculates alpha complex
                        
            outfile.close()
    
    

    if '-s' in args and '-r' in args:
        print('Do not use save(-s) and read(-r) at the same time')
        sys.exit()
    elif '-s' in args:
        print('-Saving files...')
        #np.savez_compressed('score', score)
        #np.savez_compressed('Alpha_Complex', BarCollection)
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


# def calc_residue_dist(residue1, residue2, outfile, atom1_used, atom2_used, alphaprot) :
# 
#     for atom1 in residue1:
#         for atom2 in residue2:
#             
#             distance=atom1-atom2    
# 
#             if distance>20:#If the residues is to far from eachother there is no need to calculate the distance for the rest of the atoms
#                 return(alphaprot)
#             if distance >12:#Go to the next atom
#                 continue
#                 
#             if 'H' not in atom1.get_id() and 'H' not in atom2.get_id(): # DUBBELKOLLA ATT DET INTE FINNS H I ANDRA!!!!
#             #if atom1.get_id() != 'H' and atom2.get_id() != 'H':    
#                 
#                 a1_serial=str(atom1.serial_number)
#                 a2_serial=str(atom1.serial_number)
#                 if a1_serial not in atom1_used:
#                     atom1_used.append(a1_serial)
#                     x=str(atom1.get_coord()[0])
#                     y=str(atom1.get_coord()[1])
#                     z=str(atom1.get_coord()[2])
#                     if 'C' in atom1.get_name():# and
#                         outfile.write('0'+' '+'1'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('C', []).append(x+','+y+','+z)
#                     elif 'O' in atom1.get_name():
#                         outfile.write('0'+' '+'2'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('O', []).append(x+','+y+','+z)
#                     elif 'N' in atom1.get_name():
#                         outfile.write('0'+' '+'3'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('N', []).append(x+','+y+','+z)
#                     elif 'S' in atom1.get_name():
#                         outfile.write('0'+' '+'4'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('S', []).append(x+','+y+','+z)
#                 if a2_serial not in atom2_used:
#                     atom2_used.append(a2_serial)
#                     x=str(atom2.get_coord()[0])
#                     y=str(atom2.get_coord()[1])
#                     z=str(atom2.get_coord()[2])
#                     if 'C' in atom1.get_name():
#                         outfile.write('1'+' '+'1'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('C', []).append(x+','+y+','+z)
#                     elif 'O' in atom1.get_name():
#                         outfile.write('1'+' '+'2'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('O', []).append(x+','+y+','+z)
#                     elif 'N' in atom1.get_name():
#                         outfile.write('1'+' '+'3'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('N', []).append(x+','+y+','+z)
#                     elif 'S' in atom1.get_name():
#                         outfile.write('1'+' '+'4'+' '+x+' '+y+' '+z+'\n')
#                         alphaprot.setdefault('S', []).append(x+','+y+','+z)
#                 
#     return(alphaprot)
# 
# 
# 
# 
# def calc_dist_matrix(chain_one, chain_two,outfile,alphaprot) :
#     
#     atom1_used=[]
#     atom2_used=[]
#     for row, residue_one in enumerate(chain_one) :
#         for col, residue_two in enumerate(chain_two) :
#             
#             alphaprot=calc_residue_dist(residue_one, residue_two, outfile, atom1_used, atom2_used, alphaprot)
#                                     
#     return(alphaprot)
