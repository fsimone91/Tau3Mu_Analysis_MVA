#system
import sys, os
import re
import gc

#ROOT
import ROOT
from root_pandas import read_root

import pandas as pd
import numpy as np

#data visualization
import seaborn as sn
import matplotlib.pyplot as plt

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("inputHF", type=str, help="Specify path of text file - HF")
parser.add_argument("inputW", type=str, help="Specify path of text file - W")
args = parser.parse_args()

df_HF = pd.read_csv(args.inputHF, sep='\t') 
df_HF.columns = ['run', 'lumi', 'event']
print(df_HF)
df_W = pd.read_csv(args.inputW, sep=' ') 
df_W.columns = ['run', 'lumi', 'event']
print(df_W)

#common events
df_common = pd.merge(df_HF, df_W, on=['run', 'lumi', 'event'])
print('Common')
print(df_common)
