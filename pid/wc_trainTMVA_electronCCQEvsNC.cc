#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TString.h>
#include <TMVA/Factory.h>
#include <TChain.h>

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
//   Adapdet for WCSimAnalysi output Tree June 2017                           //
//   S.Germani UCL (s.germani@ucl.ac.uk)                                      //
////////////////////////////////////////////////////////////////////////////////


void wc_trainTMVA_electronCCQEvsNC(TString sigFile, TString bkgNCFile)
{
    // Fill a TChain with the signal files (electrons fitted as electrons)
    TChain * sig = new TChain("PIDTree_ann");
    sig->Add(sigFile.Data());
    int nSig = sig->GetEntries();

    TChain * bkgNumuNC = new TChain("PIDTree_ann");
    bkgNumuNC->Add(bkgNCFile.Data());
    int nBkgNumuNC = bkgNumuNC->GetEntries();

    // Create ouput file, factory object and open the input file
    TFile* outputFile = TFile::Open( "TMVA_electronMuon_signalVsNC.root", "RECREATE" );
    TMVA::Factory* factory = new TMVA::Factory("tmvaTest_nueCCQE_vs_numuNC", outputFile, "");


    // global event weights (see below for setting event-wise weights)
    double sigWeight = 1.0;
    double bkgNumuCCWeight = 1.0;
    double bkgNumuNCWeight = 1.0;
    //factory->SetInputTrees(sig, bkg, sigWeight, bkgWeight);
    factory->AddSignalTree(sig, sigWeight);
    factory->AddBackgroundTree(bkgNumuNC, bkgNumuNCWeight);

    // Define the input variables that shall be used for the MVA training
    // (the variables used in the expression must exist in the original TTree).
    /* Third try */
    factory->AddVariable("deltaCharge2LnL",'F');
    factory->AddVariable("deltaTime2LnL",'F');
    factory->AddVariable("fracQOutsideRing_mu",'F');

    // Added 
    factory->AddVariable("fracQInRing_mu",'F');
    factory->AddVariable("fracQOutsideRing_el",'F');
    factory->AddVariable("fracQInRing_el",'F');

    factory->AddVariable("fracPredQOutsideRing_mu",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_mu",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_el",'F');
    factory->AddVariable("nHits",'F');
    factory->AddVariable("recoEOverQ_el",'F');
    factory->AddVariable("recoEOverQ_mu",'F');                          
    factory->AddVariable("chargeLikelihoodRatio_el := charge2LnL_el/nHits",'F');
    factory->AddVariable("chargeLikelihoodRatio_mu := charge2LnL_mu/nHits",'F');
    ///    factory->AddVariable("chargeLikelihoodRatio_el",'F');
    ///    factory->AddVariable("chargeLikelihoodRatio_mu",'F');
    
    // Apply additional cuts on the signal and background sample.
    // for example: TCut mycut = "abs(x)<0.5 && abs(y-0.5)<1";
    TCut mycut = "nHits > 50";

    // Use half of the events for training, half for testing
    // TString splitOpt = "NSigTrain=0:NBkgTrain=0:NSigTest=0:NBkgTest=0";
    // this version no longer works in 5.27 -- replace  7.12.10 GDC
    // factory->PrepareTrainingAndTestTree(mycut, splitOpt);
    factory->PrepareTrainingAndTestTree(mycut, mycut, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

    // Book MVA methods (see TMVA manual).  
    factory->BookMethod(TMVA::Types::kCuts, "Cuts", "H:!V");   
    factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher");   
    factory->BookMethod(TMVA::Types::kMLP, "MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5" );
    factory->BookMethod(TMVA::Types::kKNN, "KNN", "H:!V:CreateMVAPdfs:nKNN=80");   
    factory->BookMethod(TMVA::Types::kBDT, "BDT", "H:!V");   
    // Train, test and evaluate all methods

    factory->TrainAllMethods();
    factory->TestAllMethods();
    // Following line used to work, causes crash with ROOT 5.20.00, should
    // be fixed in root 5.21, see https://savannah.cern.ch/bugs/?40468
    factory->EvaluateAllMethods();    

    // Save the output and finish up

    outputFile->Close();
    std::cout << "==> wrote root file TMVA.root" << std::endl;
    std::cout << "==> TMVAnalysis is done!" << std::endl; 

    delete factory;
    return;
}

