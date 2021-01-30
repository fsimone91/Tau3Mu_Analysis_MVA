#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "BDT_optimal_cut_v2.h"


void create_eventList() 
{
    TString cat_label[] = {"A", "B", "C"};

    //Open input MiniTree
    TString fout_tree_path = TMVA_inputpath+"outputTree.root";
    TChain *ttree = new TChain("outputTree");
    ttree->Add(fout_tree_path); 
    cout<<"Reading final tree from file: "<<fout_tree_path<<endl;

    //Book output text file
    std::ofstream out;
    TString fileName = TMVA_inputpath+"_eventList.txt";
    out.open(fileName);
    out << "run\tlumisection\tevt\n";
   
    double evt, run, lumi;
    ttree->SetBranchAddress("evt",&evt); 
    ttree->SetBranchAddress("run",&run); 
    ttree->SetBranchAddress("lumi",&lumi);
 
    //Loop on categories A, B, C
    for(auto i = 0; i<3; i++){
        TString category = cat_label[i];
        cout<<"Category "<<category<<endl;

        //Make sure TChain points to firts event
        ttree->LoadTree(0);

        //set optimal categorisation
        BDTcut BDTcutvalues = Get_BDT_cut_v2(category, false);   

        TString c = ""; 
        if (category=="A") c = "0";
        if (category=="B") c = "1";
        if (category=="C") c = "2";

        //data
        TString BDTa = to_string(BDTcutvalues.a);
        TString BDTb = to_string(BDTcutvalues.b);

        ttree->Draw(">>ydata", "weight*(isMC==0  && category =="+c+" && bdt >= "+BDTb+")");
        TEventList *elist = (TEventList*)gDirectory->Get("ydata");
        ttree->SetEventList(elist);
        //signal
        //ttree->Draw("tripletMass>>ysignal", "weight*(isMC>0  && category =="+c+" && bdt >= "+BDTb+")");


        for(int i=0; i<elist->GetN(); i++){
            ttree->GetEntry(elist->GetEntry(i));
            out << fixed << setprecision(0) << run  << "\t";
            out << fixed << setprecision(0) << lumi << "\t";
            out << fixed << setprecision(0) << evt  << "\n";
        }
        ttree->SetEventList(0);
    }
    out.close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fileName<<"\n\n"<<endl;
    return 0;
}
