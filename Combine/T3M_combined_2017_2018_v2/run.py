#!/usr/bin/env python

import os
import argparse
import ROOT
from ROOT import TFile
import math
import array


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file",help="input file. [Default: %(default)s] ", action="store", default = 'input_histograms.root')
    parser.add_argument("-y", "--year",help="input file. [Default: %(default)s] ", action="store", default = '2018')
    args = parser.parse_args()

    command="root -b -q 'createDataCards.cxx(\" " + args.input_file+ "\", \"" + args.year + "\" )';"
    command +="cd datacards; combineCards.py CMS_T3MSignal_13TeV_A1_"+args.year+".txt CMS_T3MSignal_13TeV_A2_"+args.year+".txt CMS_T3MSignal_13TeV_B1_"+args.year+".txt CMS_T3MSignal_13TeV_B2_"+args.year+".txt  CMS_T3MSignal_13TeV_C1_"+args.year+".txt CMS_T3MSignal_13TeV_C2_"+args.year+".txt > CMS_T3MSignal_13TeV_Combined_"+args.year+".txt;"


    print command
    os.system(command)
