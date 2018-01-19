
# coding: utf-8

# Importing the pdb structure

# In[1]:


import Bio
from Bio.PDB import *
import sys
import re
import matplotlib.pyplot as plt
from matplotlib import cm
import math
 
         
def find_structure_params(structure):
    
    #structure_data.append([Chain_id,Seg_id, residue_name, atom_name, atom_coord_vector])
    structure_data=[]
    
    for model in structure:
        for chain in model:
            #for atom in structure.get_atoms():  If I want all atoms in a structure, depends if I want the residue
            for residue in chain:
              #  if residue.get_full_id()[3][0]==' ':  #If I want to remove HOH etc (hetero-atoms) use this!
                    for atom in residue:
                        structure_data.append([chain.get_id(),residue.get_id()[1],residue.get_resname(),atom.get_name(),atom.get_vector()])

   # for i in structure_data:
       # if i[0]=='L':
     #       print i


    return (structure_data)




# How to find different atoms for the 11 layers

# In[2]:


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
    layer10_check = ['HIS:CG','HIS:CD2','HIS:CE1','PHE:CG','PHE:CD1','PHE:CD2','PHE:CD3','PHE:CE1','PHE:CE2','PHE:CE3''PHE:CZ','TRP:CG','TRP:CD1','TRP:CD2','TRP:CD3','TRP:CD3','TRP:CE1','TRP:CE2','TRP:CE3','TRP:CZ1','TRP:CZ2','TRP:CZ3','TRP:CH2','TYR:CG','TYR:CD1','TYR:CD2','TYR:CD3','TYR:CE1','TYR:CE2','TYR:CE3','TYR:CZ']
    layer11_check = ['ALA:CB','ARG:CB','ARG:CG','ARG:CD','ASN:CB','ASP:CB','CYS:CB','GLN:CB','GLN:CG','GLU:CB','GLU:CG','HIS:CB','ILE:CB','ILE:CG1','ILE:CG2','ILE:CG3','ILE:CD1','LEU:CB','LEU:CG','LEU:CD1','LEU:CD2','LEU:CD3','LYS:CB','LYS:CG','LYS:CD','LYS:CE','MET:CB','MET:CG','MET:CE','MSE:CB','MSE:CG','MSE:CE','PHE:CB','PRO:CB','PRO:CG','PRO:CD','SER:CB','THR:CB','THR:CG2','TRP:CB','TYR:CB','VAL:CB','VAL:CG1','VAL:CG2','VAL:CG3'] #Backbone CA?
    layers=[[],[],[],[],[],[],[],[],[],[],[]]
    
    for atom in structure_data:
        if atom[3] == 'N':
            layers[1].append(atom)
        elif atom[3] == 'O':
            layers[5].append(atom)
        elif atom[3] == 'C':
            layers[8].append(atom)
        else:
            atom_type = atom[2] + ':' + atom[3]
            if atom_type in layer1_check:
                layers[0].append(atom)
            elif atom_type in layer2_check:
                layers[1].append(atom)
            elif atom_type in layer3_check:
                layers[2].append(atom)
            elif atom_type in layer4_check:
                layers[3].append(atom)
            elif atom_type in layer5_check:
                layers[4].append(atom)
            elif atom_type in layer6_check:
                layers[5].append(atom)
            elif atom_type in layer7_check:
                layers[6].append(atom)
            elif atom_type in layer8_check:
                layers[7].append(atom)
            elif atom_type in layer9_check:
                layers[8].append(atom)
            elif atom_type in layer10_check:
                layers[9].append(atom)
            elif atom_type in layer11_check:
                layers[10].append(atom)

    return(layers)


# In[3]:


def main():
    filename_pdb = '/home/joakim/Downloads/5bk1.pdb'
    
    PDBobj = PDBParser()
    structure = PDBobj.get_structure(filename_pdb, filename_pdb)
    
    structure_data = find_structure_params(structure)
    layers = make_atom_layers(structure_data)

    
if __name__ == '__main__':
  main()

