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

void trainCosmicPID()
{
    const char* sigFile = "numu_cc_combined.root";
    const char* bkgNumuCCFile = "numu_cr_combined.root";
    const char* PIDDir = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/pid/cosmicPid/";

    // Fill a TChain with the signal file
    TChain * sig = new TChain("PIDTree_ann");
    sig->Add(sigFile);
    int nSig = sig->GetEntries();

    // Fill a TChain with the background file
    TChain * bkgNumuCC = new TChain("PIDTree_ann");
    bkgNumuCC->Add(bkgNumuCCFile);
    int nBkgNumuCC = bkgNumuCC->GetEntries();

    // Set the weight directory and open the output file
	std::string weightsDir(PIDDir);
	weightsDir = weightsDir + "weights";
	(TMVA::gConfig().GetIONames()).fWeightFileDir = weightsDir;
	//std::string outputName(PIDDir);
    TFile* outputFile = TFile::Open( (weightsDir + "/tmva_cosmic.root").c_str(), "RECREATE" );

    // Create the TMVA factory object
    TMVA::Factory* factory = new TMVA::Factory("tmva_cosmic", outputFile, "");

    // global event weights (see below for setting event-wise weights)
    double sigWeight = 1.0;
    double bkgNumuCCWeight = 1.0;
    double bkgNumuNCWeight = 1.0;
    //factory->SetInputTrees(sig, bkg, sigWeight, bkgWeight);
    factory->AddSignalTree(sig, sigWeight);
    factory->AddBackgroundTree(bkgNumuCC, bkgNumuCCWeight);

    // Define the input variables that shall be used for the MVA training
    // (the variables used in the expression must exist in the original TTree).
    //factory->AddVariable("deltaEMu",'F');

    factory->AddVariable("fNVetoHits",'F');
    factory->AddVariable("fNHits",'F');
    factory->AddVariable("fFracHitsUpstream",'F');
    factory->AddVariable("fFracHitsDownstream",'F');
    //factory->AddVariable("fFracHitsInBottom",'F');
    //factory->AddVariable("fFracHitsInTop",'F');
    factory->AddVariable("fFracHitsAboveMid",'F');
    factory->AddVariable("fFracHitsBelowMid",'F');
    factory->AddVariable("fFracNHitsInRing",'F');
    factory->AddVariable("fFracNHitsOutsideRing",'F');
    factory->AddVariable("fFracNHitsInRingHole",'F');
    factory->AddVariable("fEnergy",'F');
    factory->AddVariable("fVtxX",'F');
    factory->AddVariable("fVtxY",'F');
    factory->AddVariable("fVtxZ",'F');
    factory->AddVariable("fVtxRho",'F');
    factory->AddVariable("fDirX",'F');
    factory->AddVariable("fDirY",'F');
    factory->AddVariable("fDirZ",'F');

    // Apply additional cuts on the signal and background sample. Can be any variable within the input tree
    // for example: TCut mycut = "abs(x)<0.5 && abs(y-0.5)<1";
    // We apply the preselection...

    TCut mycut = "fEscapes == 0 && fNVetoHits < 50 && fNHits > 50 && fEnergy > 550 && fEnergy < 4950 && fVtxZ > -450 && fVtxZ < 450 && fVtxRho < 1200";
    //TCut mycut;

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
