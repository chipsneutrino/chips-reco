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

//int main()
//{
    //wc_trainTMVA_electronMuonCCQE();
//}

void wc_trainTMVA_electronMuonCCQE( TString sigFile, TString bkgNumuCCFile)
{

    // Fill a TChain with the signal files (electrons fitted as electrons)
    TChain * sig = new TChain("PIDTree_ann");
    sig->Add(sigFile.Data());
    int nSig = sig->GetEntries();

    TChain * bkgNumuCC = new TChain("PIDTree_ann");
    bkgNumuCC->Add(bkgNumuCCFile.Data());
    int nBkgNumuCC = bkgNumuCC->GetEntries();

    // Create ouput file, factory object and open the input file
    TFile* outputFile = TFile::Open( "TMVA_electronMuon_signalVsNumuCCQE.root", "RECREATE" );
    TMVA::Factory* factory = new TMVA::Factory("tmvaTest_nueCCQE_vs_numuCCQE", outputFile, "");


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
    /* First try */
    /*
    factory->AddVariable("deltaCharge2LnL", 'F');
    factory->AddVariable("deltaTime2LnL", 'F');
    factory->AddVariable("insideFrac := nHitsOutsideRing_mu/nHits", 'F');
    factory->AddVariable("chargeLnLRatio_mu := charge2LnL_mu/nHits", 'F');
    factory->AddVariable("chargeLnLRatio_el := charge2LnL_el/nHits", 'F');
    */

    /* Second try */
    /*
    factory->AddVariable("charge2LnL_mu",'F');
    factory->AddVariable("charge2LnL_el",'F');
    factory->AddVariable("time2LnL_mu",'F');
    factory->AddVariable("time2LnL_el",'F');
    factory->AddVariable("fracQOutsideRing_mu",'F');
    factory->AddVariable("fracPredQOutsideRing_mu",'F');
    factory->AddVariable("predictedChargeRatio := predictedChargeOverTotalCharge_mu/predictedChargeOverTotalCharge_el");
    factory->AddVariable("nHits",'F');
    factory->AddVariable("recoE_el/totalQ",'F');
    factory->AddVariable("recoE_mu/totalQ",'F');
    */

    /* Third try */
    factory->AddVariable("deltaCharge2LnL",'F');
    factory->AddVariable("deltaTime2LnL",'F');
    factory->AddVariable("fracQOutsideRing_mu",'F');
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

    // Add these just so I get a correlation matrix...
    // factory->AddVariable("charge2LnL_el", 'F');
    // factory->AddVariable("charge2LnL_mu", 'F');
    // factory->AddVariable("endRho_el");
    // factory->AddVariable("endRho_mu");
    // factory->AddVariable("veto");
    // factory->AddVariable("escapes_el");
    // factory->AddVariable("escapes_mu");
    // factory->AddVariable("fracHitsInRingHole_el");
    // factory->AddVariable("fracHitsInRingHole_mu");
    // factory->AddVariable("fracHitsInRing_el");
    // factory->AddVariable("fracHitsInRing_mu");
    // factory->AddVariable("fracPredQInRingHole_el");
    // factory->AddVariable("fracPredQInRingHole_mu");
    // factory->AddVariable("fracPredQInRing_el");
    // factory->AddVariable("fracPredQInRing_mu");
    // factory->AddVariable("fracPredQOutsideRing_el");
    // factory->AddVariable("fracPredQOutsideRing_mu");
    // factory->AddVariable("fracQInRing_el");
    // factory->AddVariable("fracQInRing_mu");
    // factory->AddVariable("fracQOutsideRing_el");
    // factory->AddVariable("fracQOutsideRing_mu");
    // factory->AddVariable("nHits", 'I');
    // factory->AddVariable("nHitsInRing_el");
    // factory->AddVariable("nHitsInRing_mu");
    // factory->AddVariable("predictedChargeOverTotalCharge_el");
    // factory->AddVariable("predictedChargeOverTotalCharge_mu");
    // factory->AddVariable("predictedChargeRatio := predictedChargeOverTotalCharge_mu/predictedChargeOverTotalCharge_el");
    // factory->AddVariable("recoEOverCharge_el := recoE_el/totalQ");
    // factory->AddVariable("recoEOverCharge_mu := recoE_mu/totalQ");
    // factory->AddVariable("time2LnL_el");
    // factory->AddVariable("time2LnL_mu");
    // factory->AddVariable("totalPredQ_el");
    // factory->AddVariable("totalPredQ_mu");
    // factory->AddVariable("totalQInRing_el");
    // factory->AddVariable("totalQInRing_mu");
    // factory->AddVariable("charge2LnLRatio_el := charge2LnL_el/nHits");
    // factory->AddVariable("charge2LnLRatio_mu := charge2LnL_mu/nHits");







    //factory->AddVariable("recoE_mu", 'F');
    //factory->AddVariable("recoE_el", 'F');
    //factory->AddVariable("time2LnL_mu", 'F');
    //factory->AddVariable("recoE_mu", 'F');
    //factory->AddVariable("charge2LnL_el", 'F');
    //factory->AddVariable("time2LnL_el", 'F');
    //factory->AddVariable("recoE_el", 'F');
    //factory->AddVariable("totalQ", 'F');
    //factory->AddVariable("nHits", 'I');

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
