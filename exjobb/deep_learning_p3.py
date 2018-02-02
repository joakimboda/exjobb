
# coding: utf-8

# Importing the pdb structure

# from biopandas.pdb import PandasPdb
# import sys
# import re
# import matplotlib.pyplot as plt
# from matplotlib import cm
# import math
# import numpy as np
# 
# 
# def find_structure_params(structure):
#     
#     structure_data={}
#     
#     for index, row in structure_atoms.iterrows():
#         print  row.keys#row['atom_name']
#         
# 
# 
# 

# In[47]:


import os
import Bio
from Bio.PDB import *
import sys
import re
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import numpy as np

import cPickle as pickle
from mpl_toolkits.mplot3d import Axes3D


#import pandas as pd #pandas är en python modul for att hantera s.k DataFrames.
 
        
def find_structure_params(structure):
    
    #structure_data.append([Chain_id,Seg_id, residue_name, atom_name, atom_coord_vector])
    structure_data={}
    res_key = ''
    
    for model in structure:
        for chain in model:
            
            #for atom in structure.get_atoms():  If I want all atoms in a structure, depends if I want the residue
            for residue in chain:
                if residue.get_full_id()[3][0]==' ':  #If I want to remove HOH etc (hetero-atoms) use this!
                    for atom in residue:
                        key = chain.get_id() + str(atom.serial_number)
                        structure_data[key]=([chain.get_id(),residue.get_resname(),atom.get_name(),atom.get_vector()])
                        
                        #Giving the C-terminal O its own name for when making the atom layers
                        if atom.get_name()=='OXT':
                            structure_data[res_key][2]=structure_data[res_key][2]+'Cterm'
                        res_key = key
   # for i in structure_data:
       # if i[0]=='L':
     #       print i


    return (structure_data)




# How to find different atoms for the 11 layers

# In[48]:


def make_atom_layers(structure_data):
    
    layer1_check = ['CYS:SG','MET:SD','MSE:SE']
    layer2_check = ['ASN:ND2','GLN:NE2'] #Backbone N
    layer3_check = ['HIS:ND1','HIS:NE1','TRP:NE1']
    layer4_check = ['ARG:NE','ARG:NH1','ARG:NH2','ARG:NH3']
    layer5_check = ['LYS:NZ']
    layer6_check = ['ASN:OD1','GLN:OE1'] #Backbone O?
    layer7_check = ['SER:OG','THR:OG1','TYR:OH']
    layer8_check = ['ASP:OD1','ASP:OD2','ASP:OD3','GLU:OE1','GLU:OE2','GLU:OE3']
    layer9_check = ['ARG:CZ','ASN:CG','ASP:CG','GLN:CD','GLU:CD'] #Backbone C?
    layer10_check = ['HIS:CG','HIS:CD2','HIS:CE1','PHE:CG','PHE:CD1','PHE:CD2','PHE:CD3','PHE:CE1','PHE:CE2','PHE:CE3','PHE:CZ','TRP:CG','TRP:CD1','TRP:CD2','TRP:CD3','TRP:CD3','TRP:CE1','TRP:CE2','TRP:CE3','TRP:CZ1','TRP:CZ2','TRP:CZ3','TRP:CH2','TYR:CG','TYR:CD1','TYR:CD2','TYR:CD3','TYR:CE1','TYR:CE2','TYR:CE3','TYR:CZ']
    layer11_check = ['ALA:CB','ARG:CB','ARG:CG','ARG:CD','ASN:CB','ASP:CB','CYS:CB','GLN:CB','GLN:CG','GLU:CB','GLU:CG','HIS:CB','ILE:CB','ILE:CG1','ILE:CG2','ILE:CG3','ILE:CD1','LEU:CB','LEU:CG','LEU:CD1','LEU:CD2','LEU:CD3','LYS:CB','LYS:CG','LYS:CD','LYS:CE','MET:CB','MET:CG','MET:CE','MSE:CB','MSE:CG','MSE:CE','PHE:CB','PRO:CB','PRO:CG','PRO:CD','SER:CB','THR:CB','THR:CG2','TRP:CB','TYR:CB','VAL:CB','VAL:CG1','VAL:CG2','VAL:CG3'] #Backbone CA?
    
    layers={}
    
    for key, atom in structure_data.items():
        if atom[2] == 'N':
            layers.setdefault('2', []).append(atom)
        elif atom[2] == 'O':
            layers.setdefault('6', []).append(atom)
        elif atom[2] == 'C':
            layers.setdefault('9', []).append(atom)
        elif atom[2] == 'CA':
            layers.setdefault('11', []).append(atom)
        elif atom[2] == 'OXT':
            layers.setdefault('8', []).append(atom)
        elif atom[2] == 'OCterm': #C-terminal O
            layers.setdefault('8', []).append(atom)
        else:
            atom_type = atom[1] + ':' + atom[2]
            if atom_type in layer1_check:
                layers.setdefault('1', []).append(atom)
            elif atom_type in layer2_check:
                layers.setdefault('2', []).append(atom)
            elif atom_type in layer3_check:
                layers.setdefault('3', []).append(atom)
            elif atom_type in layer4_check:
                layers.setdefault('4', []).append(atom)
            elif atom_type in layer5_check:
                layers.setdefault('5', []).append(atom)
            elif atom_type in layer6_check:
                layers.setdefault('6', []).append(atom)
            elif atom_type in layer7_check:
                layers.setdefault('7', []).append(atom)
            elif atom_type in layer8_check:
                layers.setdefault('8', []).append(atom)
            elif atom_type in layer9_check:
                layers.setdefault('9', []).append(atom)
            elif atom_type in layer10_check:
                layers.setdefault('10', []).append(atom)
            elif atom_type in layer11_check:
                layers.setdefault('11', []).append(atom)
            else:
                d=1# print str(key) + str(atom[1:3]) + ' Has not been assigned to any layer' #layers.setdefault('12', []).append(atom)
    
    #This is to check if anything is not sorted into a layer            
    #for key,atom in layers.items():
        #if key=='12':
           # for x in atom:
             #   print x

    return(layers)


# In[49]:


def find_midpoint(structure_data):
    
    allvectors = Vector([0,0,0])
	
    i=0     
    for key, atom in structure_data.items():
        allvectors=allvectors+atom[3]
        i=i+1
    midpoint=Vector([allvectors[0]/i,allvectors[1]/i,allvectors[2]/i])
    return(midpoint)    


# In[50]:


def normalize_in_origo(midpoint,structure_data):

    for key,atom in structure_data.items():
        newvector = atom[3] - midpoint
        atom_new=[atom[0],atom[1],atom[2],newvector]
        structure_data[key]=atom_new

    return (structure_data)


# In[51]:


def make_density_maps(layers,density_maps,counter):
    
        for g in range(1, 12):
            density_maps[counter].append(np.zeros((120,120,120)))
            #density_maps[counter]=np.append(np.zeros((120,120,120)))
        density_maps[counter]=np.array(density_maps[counter])    
        for key,atom in layers.items():
            for pos in atom:
                x=int(np.around(pos[3][0], decimals=0))+60
                y=int(np.around(pos[3][1], decimals=0))+60
                z=int(np.around(pos[3][2], decimals=0))+60
                
                for xx in range(-2,3):
                    for yy in range(-2,3):
                        for zz in range(-2,3):
                            r=float(max(abs(xx),abs(yy),abs(zz)))
                            density=math.exp(-((r**2)/2))
                            density_maps[counter][(int(key)-1)][x+xx][y+yy][z+zz]=(density_maps[counter][(int(key)-1)][x+xx][y+yy][z+zz]+density)

        #This is only for plotting the data in python, otherwise a numpy array is made over all layers
        x_values=[]
        y_values=[]
        z_values=[]
        density_values = []
        
        for x in range(0,120):
            for y in range(0,120):
                for z in range(0,120):
                    if density_maps[counter][1][x][y][z]>0.0:
                        x_values.append(x)
                        y_values.append(y)
                        z_values.append(z)
                        density_values.append(density_maps[counter][1][x][y][z])
        
        return (density_maps,x_values,y_values,z_values,density_values)               


# In[52]:


def plot_4d(x_values,y_values,z_values,density_value):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # choose your colormap
    cmap = plt.cm.jet

    # get a Nx4 array of RGBA corresponding to zs
    # cmap expects values between 0 and 1
    density_values = np.array(density_value) # if z_list is type `list`
    colors = cmap(density_value)

    # set the alpha values according to i_list
    # must satisfy 0 <= i <= 1
    density_value = np.array(density_value)
    colors[:,-1] = density_values / density_values.max()

    # then plot
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_values, y_values, z_values, c=colors)
    plt.show()



# In[53]:


def make_mrcfile(density_maps):
    import mrcfile
    
    density_map1=density_maps[7].astype(np.float32)

    #mrc=mrcfile.open('tmp.mrc')
    #print mrc.data
    with mrcfile.new('tmp.mrc',overwrite=True) as mrc:
        mrc.set_data(density_map1)


# In[54]:


def deep_learning(density_maps):
    from keras.models import Sequential
    from keras.layers.convolutional import Conv3D
    from keras.layers import Conv3D, MaxPooling3D,Activation,Reshape
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


    #seq.add(Conv3D(filters=11, kernel_size=(3,3,3), strides=(1,1,1), activation='relu',padding='valid', data_format='channels_first', input_shape=(11,120, 120, 120)))

    #seq.add(MaxPooling3D(pool_size=(3,3,3),strides=(2,2,2),data_format='channels_first'))

    
    #seq.add(Conv3D(filters=16, kernel_size=(3,3,3), strides=(1,1,1),padding='same', data_format='channels_first'))

    #seq.add(BatchNormalization())
    
    #seq.add(Activation('relu'))

    #seq.add(MaxPooling3D(pool_size=(3,3,3),strides=(2,2,2)))

    
    #seq.add(Conv3D(filters=32, kernel_size=(3,3,3), strides=(1,1,1),padding='same', data_format='channels_first'))

    #seq.add(BatchNormalization())

    #seq.add(Activation('relu'))

    #seq.add(Conv3D(filters=32, kernel_size=(3,3,3), strides=(1,1,1),padding='same', data_format='channels_first'))

    #seq.add(BatchNormalization())

    #seq.add(Activation('relu'))
   
    #seq.add(MaxPooling3D(pool_size=(3,3,3),strides=(2,2,2)))


    #seq.add(Conv3D(filters=64, kernel_size=(3,3,3), strides=(1,1,1),padding='same', data_format='channels_first'))

    #seq.add(BatchNormalization())

    #seq.add(Activation('relu'))

    #seq.add(Conv3D(filters=128, kernel_size=(3,3,3), strides=(1,1,1),padding='same', data_format='channels_first'))

    #seq.add(BatchNormalization())

    #seq.add(Activation('relu'))

    #seq.add(Conv3D(filters=128, kernel_size=(3,3,3), strides=(1,1,1),padding='same', data_format='channels_first'))

    #seq.add(BatchNormalization())

    #seq.add(Activation('relu'))

    #seq.add(Conv3D(filters=256, kernel_size=(3,3,3), strides=(1,1,1),padding='same', data_format='channels_first'))

    #seq.add(BatchNormalization())

    #seq.add(Activation('relu'))

    #seq.add(MaxPooling3D(pool_size=(3,3,3),strides=(2,2,2)))

    
    #seq.add(Reshape((-1,)))

    #seq.add(Activation('linear'))

    #seq.add(Activation('relu'))

    #seq.add(Activation('linear'))

    #seq.add(Activation('relu'))

    #seq.add(Activation('linear'))


    
    def mean_pred(y_true, y_pred):
        return K.mean(y_pred)
    


    adam = Adam(lr=0.0003, decay=0.01)
    seq.compile(loss='binary_crossentropy',
	      optimizer=adam,
              metrics=['accuracy', mean_pred])
    class_y = np.random.random((10, 11, 118,118,118))
	
    seq.fit(density_maps,class_y,
          epochs=20,
          batch_size=9)
    
    print('Training done')
    


# In[55]:




# In[56]:


def main():
    
    
    #Dir path
    args = sys.argv[1:]

    
    files='/home/joakim/Downloads/models/' #str(args[0])

    pdb_list=[]
    for file in os.listdir(files):
        if file.endswith(".pdb"):
            pdb_list.append(files+'/'+file)

    

    density_maps=[]    
    
    if not '-r' in args:
    	for counter, filename_pdb in enumerate(pdb_list):
        	print(filename_pdb)
		percent= np.around((np.float(counter+1) / len(pdb_list)*100), decimals=3)
		print(str((counter+1)) +'/'+ str(len(pdb_list)) + ', ' + str(percent)+'%')
        	density_maps.append([])
      
        	#filename_pdb = '/home/joakim/Downloads/5eh6.pdb'#'/home/joakim/Downloads/2HIY_A.pdb' #'/home/joakim/Downloads/D1A2K-a0a-merged.pdb'
        	try: 
        	    PDBobj = PDBParser()
        	    structure = PDBobj.get_structure(filename_pdb, filename_pdb)

        	except IOError: 
        	    print('IO Error', filename_pdb)       
        	    while 'true':
        	        input1=raw_input("Error parsing PDB file! Continue(y/n)...")
        	        if input1.lower()=='y':
        	            break
        	        elif input1.lower()=='n':
        	            sys.exit()
       #	 structure=PandasPdb().read_pdb(filename_pdb)
       #	 structure_atoms = structure.df['ATOM']


        	#df=pd.read_excel('cross_val_sets.xls')
        	#print df


        	structure_data = find_structure_params(structure)

        	midpoint = find_midpoint(structure_data) #Avrundning gör att det blir lite konstigt,!! Kolla upp

        	structure_data = normalize_in_origo(midpoint,structure_data)

        	layers = make_atom_layers(structure_data)

        	#for key,values in layers.items():
        	    #print key, len(values)

        	#Create 11 density maps (zeros)
        	(density_maps,x_values,y_values,z_values,density_values) = make_density_maps(layers,density_maps,counter)

    
      ##	#  make_mrcfile(density_maps)

      ##	#  plot_4d(x_values,y_values,z_values,density_values)

	protein_maps=np.array(density_maps)
    
#For saving and loading the array as an compressed npz file
    if '-s' in args and '-r' in args:
	print('Do not use save(-s) and read(-r) at the same time')
	sys.exit()
    elif '-s' in args:
	print('-Saving file...')
	np.savez_compressed('density_maps', protein_maps)
	print('-File saved as "density_maps.npz"')

    elif '-r' in args:
	print('-Loading file...')
	protein_maps = np.load('density_maps.npz')
    	for key,array in protein_maps.items():
    		protein_maps=protein_maps[key]
	print('-File loaded')
	
    print(type(protein_maps))
    print(type(protein_maps[0]))
    print(type(protein_maps[0][0]))
    print(type(protein_maps[0][0][0]))
    print(type(protein_maps[0][0][0][0]))
    print(type(protein_maps[0][0][0][0][0]))
		
    deep_learning(protein_maps)
        
            

if __name__ == '__main__':
  main()

