import sys
import os
import csv
import string
import datetime

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("dataset", type=str, choices=['data', 'MC'], help="Specify if data or Monte Carlo")
parser.add_argument("anatype", type=str, choices=['tau3mu', 'control'], help="Specify analysis type")
# Optional Arguments
parser.add_argument("--run", type=str, default='', choices=[         '2016B', '2016C', '2016D', '2016E', '2016F', '2016G', '2016Hv2', '2016Hv3', 
                                                                     '2017B', '2017C', '2017D', '2017E', '2017F', 
                                                            '2018A', '2018B', '2018C', '2018D'
                                                            ], help="run in data")
parser.add_argument("--n", type=int, default=255, help="number of .root files per job")
parser.add_argument("--MCprocess", type=str, default='', choices=['2016Ds', '2016B0', '2016Bp',
                                                                  '2017Ds', '2017B0', '2017Bp',
                                                                  '2018Ds', '2018B0', '2018Bp',
                                                                  '2017DsPhiPi', '2018DsPhiPi', 'MiniBias'], help="process in Monte Carlo")
parser.add_argument("--reco", type=str, default='rereco', choices=['UL', 'rereco'], help="Specify reco campaign (rereco / UL)")
args = parser.parse_args()

#prepare output filename  and option string
if args.dataset == 'data':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+'" "'+args.run+'"'
else:
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.MCprocess+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+'" "'+args.MCprocess+'"'

startTime = datetime.datetime.now().strftime("%Y%m%d_%H%M")

# Create target Directory if don't exist
output_folder = '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/'
code_folder = os.path.dirname(os.path.abspath(__file__))

if not os.path.exists(output_folder+startTime):
    os.mkdir(output_folder+startTime)
    print('Directory '+output_folder+startTime+' created\n')
else:    
    print('Directory '+output_folder+startTime+' already exists\n')

if args.anatype == 'tau3mu' and args.reco == 'rereco':
   if args.dataset == 'data' and args.run == '2016B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run201B_AOD_v1/200217_093036'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run201B_AOD_v0/200104_124224'
   if args.dataset == 'data' and args.run == '2016C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016C_AOD_v1/200214_182431'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016C_AOD_v0/200104_130806'
   if args.dataset == 'data' and args.run == '2016D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016D_AOD_v1/200214_183128'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016D_AOD_v0/200104_130853'
   if args.dataset == 'data' and args.run == '2016E':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016E_AOD_v1/200214_183658'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016E_AOD_v0/200104_132747'
   if args.dataset == 'data' and args.run == '2016F':
      path = ''
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016F_AOD_v0/200104_133637'
   if args.dataset == 'data' and args.run == '2016G':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016G_AOD_v1bis/200214_184229'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016G_AOD_v4_IsoBS/200210_130222'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016G_AOD_v2/200204_160258'
   if args.dataset == 'data' and args.run == '2016Hv2':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv2_v1/200214_184620'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv2_v0/200104_145634'
   if args.dataset == 'data' and args.run == '2016Hv3':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv3_v1/200214_184703'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv3_v0/200104_151137'

   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_2017eraB_rereco31Mar2018_Mini_v0/201101_105048/'
      #path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017B_Mini_v10/201012_124331/'#Refitted tracks and photons
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_2017eraC_rereco31Mar2018_Mini_v0/201101_105221/'
      #path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017C_Mini_v10/201012_124349/' #Refitted tracks and photons
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_2017eraD_rereco31Mar2018_Mini_v0/201101_105248/'
      #path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017D_Mini_v10/201012_124407/' #Refitted tracks and photons
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_2017eraE_rereco31Mar2018_Mini_v0/201101_105314/'
      #path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017E_Mini_v10/201012_124426/' #Refitted tracks and photons
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_2017eraF_rereco31Mar2018_Mini_v0/201101_105340/'
      #path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017F_Mini_v10/201012_125134/' #Refitted tracks and photons

   if args.dataset == 'data' and args.run == '2018A':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_Run2018A_Mini_v0/201012_123526/'  #Refitted tracks
   if args.dataset == 'data' and args.run == '2018B':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_Run2018B_Mini_v0/201012_123546/'   #Refitted tracks
   if args.dataset == 'data' and args.run == '2018C':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_Run2018C_Mini_v0/201012_123602/'   #Refitted tracks
   if args.dataset == 'data' and args.run == '2018D':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_Run2018D_Mini_v0/201012_123618/'   #Refitted tracks
      #path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018D_AOD_v6/200908_165157/'

if args.anatype == 'tau3mu' and args.reco == 'UL':
   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_UL2017_Run2017B_Mini_v2/200715_142459/'
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_UL2017_Run2017C_Mini_v2/200715_142513/'
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_UL2017_Run2017D_Mini_v2/200715_142527/'
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_UL2017_Run2017E_Mini_v2/200715_142542/'
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_UL2017_Run2017F_Mini_v2/200715_142556/'

if args.anatype == 'control' and args.reco == 'rereco':
   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_rereco2017_Run2017B_Mini_v2/200710_142902/' #miniAOD rereco v2
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_rereco2017_Run2017C_Mini_v2/200710_142950/' #miniAOD rereco v2
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_rereco2017_Run2017D_Mini_v2/200710_143054/' #miniAOD rereco v2
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_rereco2017_Run2017E_Mini_v2/200710_143151/' #miniAOD rereco v2
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_rereco2017_Run2017F_Mini_v2/200710_143226/' #miniAOD rereco v2
   if args.dataset == 'data' and args.run == '2018A':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_Run2018A_Mini_v0/200930_160631/' #miniAOD 2018 rereco
   if args.dataset == 'data' and args.run == '2018B':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_Run2018B_Mini_v0/200930_160651/' #miniAOD 2018 rereco
   if args.dataset == 'data' and args.run == '2018C':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_Run2018C_Mini_v0/200930_160703/' #miniAOD 2018 rereco
   if args.dataset == 'data' and args.run == '2018D':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_Run2018D_Mini_v0/200930_160716/' #miniAOD 2018 rereco

if args.anatype == 'control' and args.reco == 'UL':
   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_UL2017_Run2017B_Mini_v0/200708_184129/'
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_UL2017_Run2017C_Mini_v0/200708_184148/'
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_UL2017_Run2017D_Mini_v0/200708_184206'
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_UL2017_Run2017E_Mini_v0/200708_184224/'
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimPhiPi_UL2017_Run2017F_Mini_v0/200708_184245/'

if args.dataset == 'MC' and args.MCprocess == '2017Ds' and args.reco == 'rereco':
   path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_ext_November2020/SkimTau3Mu_MC2017_DsTau3Mu_Mini_privJian_v0/201202_183157/' #additional sample by Jian
   #path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimTau3Mu_MC2017_DsTau3Mu_Mini_v6/200921_171432/' #reffitted mu tracks + photons
   #path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimTau3Mu_MC2017_DsTau3Mu_Mini_noHLT_v0/200921_172425/' #!no HLT filter applied
if args.dataset == 'MC' and args.MCprocess == '2017Ds' and args.reco == 'UL':
   path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimTau3Mu_UL2017_DsTau3Mu_Mini_v2/201010_173409/' #UL refitted mu tracks + photons

if args.dataset == 'MC' and args.MCprocess == '2017B0':
   path = '/lustre/cms/store/user/fsimone/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimTau3Mu_MC2017_BdTau3Mu_Mini_v6/201012_133409/' #reffitted mu tracks + photons

if args.dataset == 'MC' and args.MCprocess == '2017Bp':
   path = '/lustre/cms/store/user/fsimone/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimTau3Mu_MC2017_BuTau3Mu_Mini_v6/201012_133556/' #reffitted mu tracks + photons

if args.dataset == 'MC' and args.MCprocess == '2017DsPhiPi' and args.reco == 'rereco':
   path = '/lustre/cms/store/user/fsimone/DsToPhiPi_ToMuMu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimPhiPi_MC2017_DsPhiPi_Mini_v1/200707_151613/' #rereco MC DsPhiPi PV fixed
if args.dataset == 'MC' and args.MCprocess == '2017DsPhiPi' and args.reco == 'UL':
   path = '/lustre/cms/store/user/fsimone/DsToPhiPi_ToMuMu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimPhiPi_UL2017_DsPhiPi_Mini_v1/200707_142012/' #UL MC DsPhiPi PV fixed

if args.dataset == 'MC' and args.MCprocess == '2018DsPhiPi' and args.reco == 'rereco':
   path = '/lustre/cms/store/user/fsimone/DsToPhiPi_ToMuMu_MuFilter_TuneCP5_13TeV-pythia8/SkimPhiPi_MC2018_DsPhiPi_Mini_v0/200930_155602/' #rereco MC DsPhiPi 2018

if args.dataset == 'MC' and args.MCprocess == '2016Ds':
      path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignal_2016_v2/200215_162405'
if args.dataset == 'MC' and args.MCprocess == '2016B0':
      path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignal_2016_v2/200215_162625/'
if args.dataset == 'MC' and args.MCprocess == '2016Bp':
      path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignal_2016_v2/200215_162024/'

if args.dataset == 'MC' and args.MCprocess == '2018Ds':
      #path = '/lustre/cms/store/user/fsimone/DsToTau_TauTo3Mu_November2020/SkimTau3Mu_MC2018_DsTau3Mu_Mini_privJian_v0/201202_183823/' #additional sample by Jian
      path = '/lustre/cms/store/user/fsimone/DsToTau_TauTo3Mu_March2020/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_March2020_v11/200921_171906' #cov matrices + refitted tracks
      #path = '/lustre/cms/store/user/fsimone/DsToTau_TauMass1p87/SkimTau3Mu_MC2018_DsTau3Mu_TauMass1p87_Mini_privJian_v0/201209_101437/0000/' #shifted tau mass 1.87
      #path = '/lustre/cms/store/user/fsimone/DsToTau_TauMass1p67/SkimTau3Mu_MC2018_DsTau3Mu_TauMass1p67_Mini_privJian_v0/201208_091752/0000/' #shifted tau mass 1.67
      #path = '/lustre/cms/store/user/fsimone/DsToTau_TauTo3Mu_March2020/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_March2020_v10/200908_170120/' #updated ntuplizer
if args.dataset == 'MC' and args.MCprocess == '2018B0':
      path = '/lustre/cms/store/user/fsimone/B0ToTau_TauTo3Mu/SkimTau3Mu_B0ToTauTo3Mu_2018_wangjian_v0/201012_135921/'  #cov matrices + refitted tracks
if args.dataset == 'MC' and args.MCprocess == '2018Bp':
      path = '/lustre/cms/store/user/fsimone/BpToTau_TauTo3Mu_November2020/SkimTau3Mu_MC2018_BpTau3Mu_Mini_privJian_v0/201202_184603/' #additional sample by Jian
      #path = '/lustre/cms/store/user/fsimone/BuTau3Mu/SkimTau3Mu_BuTau3Mu_2018_bjoshi_v0/201012_135903/' #cov matrices + refitted tracks

#generating the list of all .root files in given directory and subdirectories
fileList = []
for r, d, f in os.walk(path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            fileList.append(os.path.join(r, file))

#prepare final script
final_script_path = output_folder+startTime+"/submit_analysis_"+startTime+".sh"
final_script = open(final_script_path, "w")
final_script.write("#!/bin/bash\n")
final_script.write("chmod 777 -R "+output_folder+startTime+"/*\n")
final_script.write("chmod 777 -R "+code_folder+"/*\n")
final_script.write("cd "+output_folder+startTime+"\n")

#loop to generate one .cpp+executable+batch system conf file for each group of "n" files
n_chunk = len(fileList)//args.n
print('using ntuples in '+path)
print('Number of files is {0:2d}'.format(len(fileList)))
print('Number of jobs is {0:2d}'.format(n_chunk+1))
for file_index in range(n_chunk+1):
      chunk = '' 
      for idx, l in enumerate(fileList):
         if idx < args.n*(file_index+1) and idx >= args.n*file_index:
             l = l.rstrip()
             l = '        chain->AddFile("{}");\n'.format(l)
             chunk = chunk + l

      #analysis.cpp template
      with open("templates/Analysis_template.cpp", "r") as in_file:
          buf = in_file.readlines()

      cpp_filename = "Analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_chunk"+str(file_index)+".cpp"
      with open(cpp_filename, "w") as out_file:
          for lb in buf:
              if lb == '        //AddFile_'+args.dataset+args.MCprocess+'_'+args.anatype+'\n':
                  #write group of files
                  out_file.write(chunk)
              elif lb == '        //OutFile_'+args.dataset+args.MCprocess+'_'+args.anatype+'\n':
                  #write output file name
                  out_file.write('        fileout = "'+out_filename+str(file_index)+'.root";\n')
              else: out_file.write(lb)

      #executable template
      with open("templates/launch_analysis_template.job", "r") as launch_infile:
          buf2 = launch_infile.readlines()

      launch_filename = "launch_analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_folder+startTime+"/"+launch_filename, "w") as launch_outfile:
          for lb2 in buf2:
              if lb2 == "#compile\n":
                  launch_outfile.write("cd "+output_folder+startTime+"\n")
                  launch_outfile.write("g++ -I $ROOTSYS/include "+code_folder+"/"+cpp_filename+" `root-config --glibs` `root-config --libs` `root-config --cflags` -lTMVA -L $ROOTSYS/lib -o executable"+str(file_index)+"\n")
              elif lb2 == "#execute\n":
                  launch_outfile.write('./executable'+str(file_index)+option_string+'\n')
              else: launch_outfile.write(lb2)

      #myCondor template
      with open("templates/my_HTCondor_template.job", "r") as myCondor_infile:
          buf3 = myCondor_infile.readlines()

      condor_filename = "my_HTCondor_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_folder+startTime+"/"+condor_filename, "w") as myCondor_outfile:
          for lb3 in buf3:
              if lb3 == "Executable = launch_analysis_template.job\n":
                  myCondor_outfile.write("Executable = "+launch_filename+"\n")
              else: myCondor_outfile.write(lb3)

      #add lines to final script
      final_script.write("echo condor_submit "+condor_filename+" -name ettore\n")
      final_script.write("condor_submit "+condor_filename+" -name ettore\n")

final_script.write("cd "+code_folder+"\n")
final_script.close()
print("to submit analysis to batch system do:\nsource "+final_script_path+"\n")
