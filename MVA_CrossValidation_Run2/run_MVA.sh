#!/bin/bash
unset dataset_name
dataset_name="training_run2_CV_ntrees500depth3beta0p5cuts50_5folds"
echo "Folder will be $dataset_name"

echo "root -l -b MVA_code_run2.cpp\(\"A\"\)"
root -l -b MVA_code_run2.cpp\(\"A\"\) >& MVA_training_CV_A.log
mv BDTG_fold*.root "$dataset_name"_A/. 
mv MVA_training_CV_A.log "$dataset_name"_A/. 

echo "root -l -b MVA_code_run2.cpp\(\"B\"\)"
root -l -b MVA_code_run2.cpp\(\"B\"\) >& MVA_training_CV_B.log
mv BDTG_fold*.root "$dataset_name"_B/. 
mv MVA_training_CV_B.log "$dataset_name"_B/. 

echo "root -l -b MVA_code_run2.cpp\(\"C\"\)"
root -l -b MVA_code_run2.cpp\(\"C\"\) >& MVA_training_CV_C.log
mv BDTG_fold*.root "$dataset_name"_C/. 
mv MVA_training_CV_C.log "$dataset_name"_C/. 

#echo "root -l -b MVA_code_run2.cpp\(\"ABC\"\)"
#root -l -b MVA_code_run2.cpp\(\"ABC\"\)
