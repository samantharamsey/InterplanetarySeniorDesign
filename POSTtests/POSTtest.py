# -*- coding: utf-8 -*-
'''
Created on Mon Feb 24 11:45:09 2020

@author: sam
'''
import pandas as pd


if __name__ == '__main__':
    
    # load in the data from POST
    filepath1 = r'C:\Users\saman\OneDrive\Desktop\POSTtests'
    filename1 = r'\tin7'
    test_data = pd.read_csv(filepath1 + filename1 + r'.hdf', delim_whitespace = True)
    # create a new DF with only the last lines
    comp_data = pd.DataFrame([])
    comp_data = comp_data.append(test_data[-1:])
    
    # load in the final orbit data
    filepath2 = r'C:\Users\saman\OneDrive\Desktop\InterplanetarySeniorDesign'
    filename2 = r'\3D_potato'
    calc_data = pd.read_hdf(filepath2 + filename2 + r'.hdf')
    
    # want to compare each line of comp_data to every line of calc_data
    # to see which comp_data cases will result in any satisfactory orbit
    listy = []
    print("hey")
    #print(comp_data)
    for i in range(1, calc_data.shape[0]):
        print(calc_data[i])
        #if comp_data == line:
            #print("found")
            #listy.append(comp_data)
    
            
    