#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <TFile.h>
#include <TString.h>
#include "TTree.h"
#include <TChain.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <TMVA/Factory.h>

////////////////////////////////////////////////////////////////////////////////
//                   Andy Perch, UCL PhD Course, Nov2012                      //
////////////////////////////////////////////////////////////////////////////////
//                   wc_trainTMVA_electronMuonCCQE.cc                         //
////////////////////////////////////////////////////////////////////////////////
//  Program to train Fisher dicriminant and multilayer perceptron with a      //
//  single hidden layer tests to distinguish signal and background events.    //
//  Adapted from http://www.pp.rhul.ac.uk/~cowan/stat/root/tmva/ for Glen     // 
//  Cowan's statistical analysis of data problem set 6 to include the MLP.    //
////////////////////////////////////////////////////////////////////////////////
//   Adapted for WCSimAnalysis output Tree June 2017                          //
//   S.Germani UCL (s.germani@ucl.ac.uk)                                      //
////////////////////////////////////////////////////////////////////////////////
//   Adapted when restructuring the PID July 2017                             //
//   J.Tingey UCL (j.tingey.16@ucl.ac.uk)                                     //
////////////////////////////////////////////////////////////////////////////////

void wc_trainTMVA_electronCCQEvsNC(TString sigFile, TString bkgNCFile, TString PIDDir)
{
	// Fill a TChain with the signal file
    TChain * sig = new TChain("PIDTree_ann");
    sig->Add(sigFile.Data());
    int nSig = sig->GetEntries();

    // Fill a TChain with the background file
    TChain * bkgNumuNC = new TChain("PIDTree_ann");
    bkgNumuNC->Add(bkgNCFile.Data());
    int nBkgNumuNC = bkgNumuNC->GetEntries();

    // Set the weight directory and open the output file
	std::string weightsDir(PIDDir);
	weightsDir = weightsDir + "weights";
	(TMVA::gConfig().GetIONames()).fWeightFileDir = weightsDir;
	//std::string outputName(PIDDir);
    TFile* outputFile = TFile::Open( (weightsDir + "/tmva_nueCCQE_vs_numuNC.root").c_str(), "RECREATE" );

    // Create the TMVA factory object
    TMVA::Factory* factory = new TMVA::Factory("tmva_nueCCQE_vs_numuNC", outputFile, "");

    // global event weights (see below for setting event-wise weights)
    double sigWeight = 1.0;
    double bkgNumuCCWeight = 1.0;
    double bkgNumuNCWeight = 1.0;
    //factory->SetInputTrees(sig, bkg, sigWeight, bkgWeight);
    factory->AddSignalTree(sig, sigWeight);
    factory->AddBackgroundTree(bkgNumuNC, bkgNumuNCWeight);

    // Define the input variables that shall be used for the MVA training
    // (the variables used in the expression must exist in the original TTree).
    //factory->AddVariable("deltaEMu",'F');

    factory->AddVariable("nHits",'F');
    factory->AddVariable("deltaCharge2LnL",'F');
    factory->AddVariable("deltaTime2LnL",'F');
    factory->AddVariable("deltaCharge2LnLOverNHits",'F');
    factory->AddVariable("fracHitsOutsideRing_mu",'F');
    factory->AddVariable("fracHitsInRing_mu",'F');
    factory->AddVariable("fracHitsInRingHole_mu",'F');
    factory->AddVariable("fracHitsOutsideRing_el",'F');
    factory->AddVariable("fracHitsInRing_el",'F');
    factory->AddVariable("fracHitsInRingHole_el",'F');
    factory->AddVariable("fracPredQOutsideRing_mu",'F');
    factory->AddVariable("fracPredQOutsideRing_el",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_mu",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_el",'F');
    factory->AddVariable("recoEOverQ_mu",'F');
    factory->AddVariable("recoEOverQ_el",'F');
    factory->AddVariable("dirX_el",'F');
    factory->AddVariable("fracHitsDownstream",'F');

    factory->AddVariable("nRings",'F');
    factory->AddVariable("firstRingHeight",'F');
    factory->AddVariable("lastRingHeight",'F');
    
    // Apply additional cuts on the signal and background sample. Can be any variable within the input tree
    // for example: TCut mycut = "abs(x)<0.5 && abs(y-0.5)<1";
    // We apply the preselection and additional for this neural network a cut on the other one to remove muon-like events...
    TCut mycut = "preselected == 1 && annNueCCQEvsNumuCCQE>0.8";

    // Use half of the events for training, half for testing
    // TString splitOpt = "NSigTrain=0:NBkgTrain=0:NSigTest=0:NBkgTest=0";
    // this version no longer works in 5.27 -- replace  7.12.10 GDC
    // factory->PrepareTrainingAndTestTree(mycut, splitOpt);
    factory->PrepareTrainingAndTestTree(mycut, mycut, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

    // Book MVA methods (see TMVA manual).  
    //factory->BookMethod(TMVA::Types::kCuts, "Cuts", "H:!V");
    //factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher");
    factory->BookMethod(TMVA::Types::kMLP, "MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5" );
    //factory->BookMethod(TMVA::Types::kKNN, "KNN", "H:!V:CreateMVAPdfs:nKNN=80");
    //factory->BookMethod(TMVA::Types::kBDT, "BDT", "H:!V");
    // Train, test and evaluate all methods

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();    

    // Save the output and finish up

    outputFile->Close();
    std::cout << "==> wrote root file TMVA.root" << std::endl;
    std::cout << "==> TMVAnalysis is done!" << std::endl; 

    delete factory;
    return;
}
