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
parser.add_argument("anatype", type=str, choices=['muonid', 'prova'], help="Specify analysis type")
# Optional Arguments
parser.add_argument("--run", type=str, default='', choices=[
                                                            '2018A', '2018B', '2018C', '2018D'
                                                            ], help="run in data")
parser.add_argument("--n", type=int, default=255, help="number of .root files per job")
parser.add_argument("--massregion", type=str, default='', choices=['sgn', 'bkg'], help="Specify invariant mass range")
parser.add_argument("--MCprocess", type=str, default='', choices=['2018BdToKK',                 '2018BdToPiPi', '2018BdToKPi', '2018BdToKPi_2',
                                                                  '2018BsToKK', '2018BsToKK_2', '2018BsToPiPi', '2018BsToKPi',
                                                                  '2018Ds', '2018B0', '2018Bp',
                                                                  '2017BdToKK',                 '2017BdToPiPi', '2017BdToKPi', '2017BdToKPi_2',
                                                                  '2017BsToKK', '2017BsToKK_2', '2017BsToPiPi', '2017BsToKPi',
                                                                  '2017Ds', '2017B0', '2017Bp'
                                                                  ], help="process in Monte Carlo")
args = parser.parse_args()

#prepare output filename  and option string
if args.massregion:
   args.massregion = '_'+args.massregion
if args.dataset == 'data':
   out_filename = 'AnalysedTree_'+args.dataset+args.massregion+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_muonid","")+args.massregion+'" "'+args.run+'"'
else:
   out_filename = 'AnalysedTree_'+args.dataset+args.massregion+'_'+args.MCprocess+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_muonid","")+args.massregion+'" "'+args.MCprocess+'"'

startTime = datetime.datetime.now().strftime("%Y%m%d_%H%M")

# Create target Directory if don't exist
output_folder = '/lustre/cms/store/user/fsimone/MuonID/Analysis/'
code_folder = os.path.dirname(os.path.abspath(__file__))

if not os.path.exists(output_folder+startTime):
    os.mkdir(output_folder+startTime)
    print('Directory '+output_folder+startTime+' created\n')
else:    
    print('Directory '+output_folder+startTime+' already exists\n')

if args.anatype == 'muonid':
   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017B_Mini_v5/200704_142154/' #PV fixed
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017C_Mini_v6/200704_142222/' #PV fixed
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017D_Mini_v5/200704_142427/' #PV fixed
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017E_Mini_v5/200704_142250/' #PV fixed
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/fsimone/DoubleMuonLowMass/SkimTau3Mu_rereco2017_Run2017F_Mini_v5/200704_142318/' #PV fixed

   if args.dataset == 'data' and args.run == '2018A':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_v5-bis/200505_080158/'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_v4/200217_120328'
   if args.dataset == 'data' and args.run == '2018B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018B_AOD_v5/200429_184137/' #new variables
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018B_AOD_v4/200217_120246'
   if args.dataset == 'data' and args.run == '2018C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018C_AOD_v5/200429_184043/' #new variables
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018C_AOD_v4/200217_120028'
   if args.dataset == 'data' and args.run == '2018D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018D_AOD_v5/200429_183916/' #new variables
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018D_AOD_v4/200217_112055'

if args.dataset == 'MC' and args.MCprocess == '2018BdToKK':
      path = '/lustre/cms/store/user/fsimone/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKK_SoftQCDnonD_2018_v0/210212_180433/'
      #path = '/lustre/cms/store/user/rosma/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKK_SoftQCDnonD_2018_CMSSW_10_2_1_v2/200414_092653/'
if args.dataset == 'MC' and args.MCprocess == '2018BdToPiPi':
      path = '/lustre/cms/store/user/fsimone/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToPiPi_SoftQCDnonD_2018_v0/210212_180354/'
      #path = '/lustre/cms/store/user/rosma/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToPiPi_2018_CMSSW_10_2_1_v2/200414_092641/'
if args.dataset == 'MC' and args.MCprocess == '2018BdToKPi':
      path = '/lustre/cms/store/user/fsimone/BdToKPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKPi_2_SoftQCDnonD_2018_v0/210212_180420/'
      #path = '/lustre/cms/store/user/rosma/BdToKPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKPi_SoftQCDnonD_2018_CMSSW_10_2_1_v2/200414_092615/'
if args.dataset == 'MC' and args.MCprocess == '2018BdToKPi_2':
      path = '/lustre/cms/store/user/fsimone/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKPi_SoftQCDnonD_2018_v0/210212_180407/'
      #path = '/lustre/cms/store/user/rosma/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKPi_2018_CMSSW_10_2_1_v2/200414_092628/'
if args.dataset == 'MC' and args.MCprocess == '2018BsToKK':
      path = '/lustre/cms/store/user/fsimone/BsToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToKK_2_SoftQCDnonD_2018_v0/210212_180342/'
      #path = '/lustre/cms/store/user/rosma/BsToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToKK_2018_CMSSW_10_2_1_v2/200414_092715/'
if args.dataset == 'MC' and args.MCprocess == '2018BsToKK_2':
      path = '/lustre/cms/store/user/fsimone/BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToKK_SoftQCDnonD_2018_v0/210212_180329/'
      #path = '/lustre/cms/store/user/rosma/BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToKK_2_2018_CMSSW_10_2_1_v2/200414_092725/'
if args.dataset == 'MC' and args.MCprocess == '2018BsToPiPi':
      path = '/lustre/cms/store/user/fsimone/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToPiPi_SoftQCDnonD_2018_v0/210212_180316/'
      #path = '/lustre/cms/store/user/rosma/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToPiPi_2018_CMSSW_10_2_1_v2/200414_092705/'
if args.dataset == 'MC' and args.MCprocess == '2018BsToKPi':
      path = ''


if args.dataset == 'MC' and args.MCprocess == '2017BdToKK':
      path = '/lustre/cms/store/user/fsimone/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKK_SoftQCDnonD_2017_v1/210213_080022'
      #path = '/lustre/cms/store/user/fsimone/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKK_SoftQCDnonD_2017_CMSSW_10_2_18/200922_172051/'
if args.dataset == 'MC' and args.MCprocess == '2017BdToPiPi':
      path = '/lustre/cms/store/user/fsimone/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToPiPi_SoftQCDnonD_2017_v1/210213_075942'
      #path = '/lustre/cms/store/user/fsimone/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToPiPi_SoftQCDnonD_2017_CMSSW_10_2_18/200922_172021/'
if args.dataset == 'MC' and args.MCprocess == '2017BdToKPi':
      #path = ''
      path = '/lustre/cms/store/user/fsimone/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKPi_SoftQCDnonD_2017_CMSSW_10_2_18/200922_172031/'
if args.dataset == 'MC' and args.MCprocess == '2017BdToKPi_2':
      path = '/lustre/cms/store/user/fsimone/BdToKPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKPi_2_SoftQCDnonD_2017_v1/210213_080009'
      #path = '/lustre/cms/store/user/fsimone/BdToKPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BdToKPi_2_SoftQCDnonD_2017_CMSSW_10_2_18/200922_172040/'
if args.dataset == 'MC' and args.MCprocess == '2017BsToKK':
      path = '/lustre/cms/store/user/fsimone/BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToKK_SoftQCDnonD_2017_v1/210213_075917'
if args.dataset == 'MC' and args.MCprocess == '2017BsToKK_2':
      path = '/lustre/cms/store/user/fsimone/BsToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToKK_2_SoftQCDnonD_2017_v1/210213_075929'
      #path = '/lustre/cms/store/user/fsimone/BsToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToKK_2_SoftQCDnonD_2017_CMSSW_10_2_18/200922_172012/'
if args.dataset == 'MC' and args.MCprocess == '2017BsToPiPi':
      path = '/lustre/cms/store/user/fsimone/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToPiPi_SoftQCDnonD_2017_v1/210213_075904'
      #path = '/lustre/cms/store/user/fsimone/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/SkimMuID_BsToPiPi_SoftQCDnonD_2017_CMSSW_10_2_18/200922_171950/'
if args.dataset == 'MC' and args.MCprocess == '2017BsToKPi':
      path = ''


if args.dataset == 'MC' and args.MCprocess == '2018Ds':
      path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimTau3Mu_Summer20UL18_DsTau3Mu_Mini_v2/210131_223011/' #lowered skim mu pT
      #path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu_March2020/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_March2020_v9/'

if args.dataset == 'MC' and args.MCprocess == '2018B0':
      path = '/lustre/cms/store/user/fsimone/BdToTau_To3Mu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimTau3Mu_Summer20UL18_BdTau3Mu_Mini_v2/210131_222948/' #lowered skim mu pT
      #path = '/lustre/cms/store/user/rosma/B0ToTau_TauTo3Mu/SkimTau3Mu_B0ToTauTo3Mu_2018_CMSSW_10_2_1_v9/'

if args.dataset == 'MC' and args.MCprocess == '2018Bp':
      path = '/lustre/cms/store/user/fsimone/BuToTau_To3Mu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimTau3Mu_Summer20UL18_BuTau3Mu_Mini_v3/210303_144003/' #genP_GrandMotherPdgId
      #path = '/lustre/cms/store/user/rosma/BuTau3Mu/SkimTau3Mu_BuToTauTo3Mu_2018_CMSSW_10_2_1_v8/200217_122234'


if args.dataset == 'MC' and args.MCprocess == '2017Ds':
      path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimTau3Mu_Summer20UL17_DsTau3Mu_ModFilter_Mini_v0/210127_112243/'
      #path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimTau3Mu_MC2017_DsTau3Mu_Mini_v3/200704_142853/' #PV fixed

if args.dataset == 'MC' and args.MCprocess == '2017B0':
      path = '/lustre/cms/store/user/fsimone/BdToTau_To3Mu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimTau3Mu_Summer20UL17_BdTau3Mu_ModFilter_Mini_v0/210127_112223/'
      #path = '/lustre/cms/store/user/fsimone/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimTau3Mu_MC2017_BdTau3Mu_Mini_v3/200704_142949/' #PV fixed

if args.dataset == 'MC' and args.MCprocess == '2017Bp':
      path = '/lustre/cms/store/user/fsimone/BuToTau_To3Mu_MuFilter_TuneCP5_13TeV-pythia8-evtgen/SkimTau3Mu_Summer20UL17_BuTau3Mu_ModFilter_Mini_v0/210127_112234/'
      #path = '/lustre/cms/store/user/fsimone/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/SkimTau3Mu_MC2017_BuTau3Mu_Mini_v3/200704_142921/' #PV fixed

process_label = args.MCprocess
if args.MCprocess == '2018BdToPiPi' or args.MCprocess == '2018BdToKK' or args.MCprocess == '2018BdToKPi' or args.MCprocess == '2018BdToKPi_2':
    process_label = '2018Bkg'
if args.MCprocess == '2018BsToPiPi' or args.MCprocess == '2018BsToKK' or args.MCprocess == '2018BsToKK_2' or args.MCprocess == '2018BsToKPi':
    process_label = '2018Bkg'

if args.MCprocess == '2017BdToPiPi' or args.MCprocess == '2017BdToKK' or args.MCprocess == '2017BdToKPi' or args.MCprocess == '2017BdToKPi_2':
    process_label = '2017Bkg'
if args.MCprocess == '2017BsToPiPi' or args.MCprocess == '2017BsToKK' or args.MCprocess == '2017BsToKK_2' or args.MCprocess == '2017BsToKPi':
    process_label = '2017Bkg'

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

      cpp_filename = "Analysis_"+args.dataset+args.massregion+"_"+args.run+args.MCprocess+"_"+args.anatype+"_chunk"+str(file_index)+".cpp"
      with open(cpp_filename, "w") as out_file:
          for lb in buf:
              if lb == '        //AddFile_'+args.dataset+process_label+'_'+args.anatype+'\n':
                  #write group of files
                  out_file.write(chunk)
              elif lb == '        //OutFile_'+args.dataset+process_label+'_'+args.anatype+'\n':
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
