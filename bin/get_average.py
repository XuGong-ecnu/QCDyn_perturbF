#!/usr/bin/env python

"""
Created on Tue. Dec. 21 2021 @NYU-SH
Updated on Tue. Dec. 21 2021 @NYU-SH

@author: Zhubin Hu
"""
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description='Get average result of input .csv files and save to a csv file')
parser.add_argument('input_files',nargs='*',help='specify the names of input csv files')
parser.add_argument('-p', '--prefix',help='specify the prefix of the csv file to save average result, default is \"average\"',type=str,default='average')
args = parser.parse_args()

num_files = len(args.input_files) # number of input csv files
print('The number of input csv files is', num_files)
for i in range(num_files):
    file = args.input_files[i] # filename
    if (file[-4:] != '.csv'):
        print('ERROR: Input file:', file, 'is not a csv file.')
        sys.exit()
    if (i == 0): # the first csv file
        df = pd.read_csv(file, sep=',')
    else: # accumulate another files
        df += pd.read_csv(file, sep=',')
df /= num_files # get average dataframe
df.to_csv(args.prefix + '.csv', index=False) # save the averaged results to current dir
print('The average result of them has been saved to csv file:', args.prefix + '.csv')