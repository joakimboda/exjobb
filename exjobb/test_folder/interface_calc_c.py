
# coding: utf-8

# In[9]:


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


# In[10]:


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
                    layers.setdefault('dist',[]).append([atom1.get_serial_number(),atom2.get_serial_number(), distance])
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


# In[11]:


def deep_learning(protein_train,protein_target):
    from keras.models import Sequential
    from keras.layers.convolutional import Conv3D,Conv1D
    from keras.layers import Conv3D, MaxPooling3D,Activation,Reshape,Dense,AveragePooling1D,Dropout,Flatten
    from keras.layers.normalization import BatchNormalization
    from keras.optimizers import Adam
    
    import keras.backend as K #For compile
    
    print('Start training')
    seq = Sequential()

    #seq.add(Conv3D(11, 3, 3, 3, activation='relu', 
                            #border_mode='valid', name='conv1',
                            #subsample=(1, 1, 1),
                            #dim_ordering='th', 
                            #input_shape=(11,120, 120, 120)))


    seq.add(Conv1D(filters=128, kernel_size=(3), strides=(1), activation='relu',padding='valid', input_shape=(600,3)))
    
    seq.add(Conv1D(filters=128, kernel_size=(3), strides=(1), activation='relu',padding='valid', input_shape=(600,3)))

    seq.add(AveragePooling1D(pool_size=2, strides=None, padding='valid'))
    
    seq.add(Dropout(0.25))
    
    seq.add(Conv1D(filters=256, kernel_size=(3), strides=(1), activation='relu',padding='valid', input_shape=(600,3)))
    
    seq.add(Conv1D(filters=256, kernel_size=(3), strides=(1), activation='relu',padding='valid', input_shape=(600,3)))

    seq.add(AveragePooling1D(pool_size=2, strides=None, padding='valid'))
    
    seq.add(Dropout(0.25))
    
    seq.add(Flatten())
    
    seq.add(Dense(4096,activation='relu'))
    seq.add(Dense(4096,activation='relu'))
    
    seq.add(Dense(4096,activation='tanh'))
    seq.add(Dense(4096,activation='tanh'))
    
    seq.add(Dropout(0.5))

    #seq.add(Activation('linear'))
   # seq.add(Dense(256,activation='linear'))

    #seq.add(Activation('relu'))

    #seq.add(Activation('linear'))
    #seq.add(Dense(128,activation='linear'))

    #seq.add(Activation('relu'))

    #seq.add(Activation('linear'))
   # seq.add(Dense(1,activation='linear'))
    
    

    seq.summary()
    print('ready')
    def mean_pred(y_true, y_pred):
        return K.mean(y_pred)
    
    #protein_target=np.random.rand(16)
    #protein_target=np.random.rand(16,269748)
    #protein_train=np.random.rand(16,11,120,120,120)

    adam = Adam(lr=0.0003, decay=0.01)
    seq.compile(loss='mean_squared_error',
            optimizer=adam,
              metrics=['accuracy', mean_pred])

    
    seq.fit(protein_train,protein_target,
          epochs=20,
          batch_size=9)
    
    print('Training done')


# In[12]:


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
        df=pandas.read_excel('/home/joakim/Downloads/cross_val_sets.xls')
        targets={}
        for i,row in df.iterrows():
             targets.setdefault(i, []).append(row['Targets'].split())

        dir_home='/home/joakim/Downloads/CnM-dataset/MOAL_Benchmark_test/'
        datasets={}


        for dir in os.listdir(dir_home):
            for key,target in targets.items():
                if dir in str(target):
                    datasets.setdefault(key, []).append(dir)





        score=[]
        dist_matrix=[]
        for key,protein_list in datasets.items():

            for prot_counter,protein in enumerate(protein_list):


                pdb_list=[]
                for file in os.listdir(dir_home + protein):

                    if file.endswith(".pdb"):
                        pdb_list.append([dir_home + protein,file])



               # pdb_list=[]
               # for file in os.listdir(dir_home + protein):
                   # if file.endswith(".pdb"):
                        #pdb_list.append(args+'/'+file)
                     #   pdb_list.append([dir_home + protein,file])


                valuefile='/home/joakim/Downloads/CnM.featuresNPqDNZ'#'D:\Downloads\CnM.featuresNPqDNZ'
                value_df=pandas.read_csv(valuefile,delim_whitespace=1)



                for counter, filepath_pdb in enumerate(pdb_list):
                    #print(counter)
                    #print(pdb_list[counter][1])
                    filename_pdb=os.path.join(pdb_list[counter][0],pdb_list[counter][1])
                    print filename_pdb
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


                   # if 'dist_matrix' in dir():
                    dist_matrix.append(np.array(layers['dist']))
                        #dist_matrix=np.append(dist_matrix,layers['dist'])
                        #dist_matrix=np.concatenate((dist_matrix, layers['dist']))
                    #else:
                       # dist_matrix=np.array(layers['dist'])
                    print len(dist_matrix)

        dist_matrix=np.array(dist_matrix)    
        score=np.array(score)
        print('dist_matrix: ', dist_matrix.shape)
        print('score: ',score.shape)

    if '-s' in args and '-r' in args:
        print('Do not use save(-s) and read(-r) at the same time')
        sys.exit()
    elif '-s' in args:
        print('-Saving files...')
        np.savez_compressed('score', score)
        np.savez_compressed('dist_matrix', dist_matrix)
        #np.savez_compressed('dist_matrix', dist_matrix)
        print('-Files saved as "dist_matrix.npz","score.npz"')

    elif '-r' in args:
        print('-Loading files...')
        score = np.load('score.npz')
        dist_matrix = np.load('dist_matrix.npz')
        
        for key,array in score.items():
            score=score[key]
        for key,array in dist_matrix.items():
            dist_matrix=dist_matrix[key]
            
        print('-Files loaded')
    
    
    print('dist_matrix: ', dist_matrix.shape)
    print('score: ',score.shape)
    
    print type(dist_matrix)
    print dist_matrix[1].shape
    print type(dist_matrix[1][1])
    print type(dist_matrix[1][1][1])
    
    dist_matrix=[]
    dist_matrix=(np.zeros((56,600,3)))
    dist_matrix=np.array(dist_matrix)

    deep_learning(dist_matrix,score)
    
    
if __name__ == '__main__':
    main()

