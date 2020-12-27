#!/usr/bin/env python

import os
import argparse
import ROOT
from ROOT import TFile
import math
import array

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",help="input file. [Default: %(default)s] ", action="store", default = 'CMS_T3MSignal_13TeV_Combined.txt')
    parser.add_argument("-y", "--year",help="year. [Default: %(default)s] ", action="store", default = '')
    parser.add_argument("-c", "--combine",help="combine 2017 and 2018 input files. [Default: %(default)s] ", action="store_true", default = 'false')
    args = parser.parse_args()
    
    if (args.combine):
        command_comb = "combineCards.py CMS_T3MSignal_13TeV_Combined_2017.txt CMS_T3MSignal_13TeV_Combined_2018.txt > CMS_T3MSignal_13TeV_Combined.txt"
    print command_comb
    os.system(command_comb)

    inputfile = ''
    if (args.year=='2017' or args.year=='2018'):
        inputfile = 'CMS_T3MSignal_13TeV_Combined_'+args.year+'.txt'
    else:
        inputfile = args.input_file
    command ="combineTool.py -M T2W -o workspace_"+args.year+".root -i "+inputfile+";"
    command +="combineTool.py -M AsymptoticLimits  --run blind  -d workspace_"+args.year+".root -m 1.777 --cl 0.90;"


    print command
    os.system(command)
