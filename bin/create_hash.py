#!/usr/bin/env python

import pickle

'''
This python script reads the tsv file containing the substrate data and
generates a python file that has a indexable dictionary of the substrate
data.

---Created by Vasu Pranav Sai Iddamsetty--- 
'''

def create_hash():
    
    #initialize the input file and the output file

    in_file = open('ChEBI_IDs_Jimmy_1.txt','rb+')
    
    #data = in_file.read()
    

    out_file = open('substrate_dict.py','w+')

    
    #initialize a list to which the dictionary entries will be appended
    substrate_info ={}
    
    #skip the first line, which contains the header information
    next(in_file)

    #from each line, create a dicitonary entry
    for line in in_file:
        
        print(line)
        
        data = line.split('\t')

        substrate_info[data[0]] = data[1:]



    #write the contents the output file and close
    pickle.dump(substrate_info,out_file)

    out_file.close()




if __name__ == "__main__":

    create_hash()


