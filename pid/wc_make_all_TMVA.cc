#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TString.h>
#include <TMVA/Factory.h>
#include <TChain.h>
#include "wc_trainTMVA_electronMuonCCQE.hh"
#include "wc_trainTMVA_electronCCQEvsNC.hh"

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
int main(int argc, char *argv[])
{
  
  TString opt;
  TString sigFile;
  TString bkgFile;

  if(argc>1) opt = argv[1];
  if(argc!=4  &&  opt!="--cc"  &&  opt!="--nc" ) {
    std::cout << "   "                                                               << std::endl;
    std::cout << "usage: wc_make_all_TMVA --bkgOption  signal_file  background_file" << std::endl;
    std::cout << "   "                                                               << std::endl;
    std::cout << "   Background Options: "                                           << std::endl;
    std::cout << "   --cc                Charged Current Backgrount (Numu)  "        << std::endl;
    std::cout << "   --nc                Neutral Current Backgrount         "        << std::endl;
    std::cout << "   "                                                               << std::endl;
    std::cout << "   example: wc_make_all_TMVA --cc  nuecc.root  numucc.root"        << std::endl;
    std::cout << "   "                                                               << std::endl;

  } else if(opt=="--cc") {

    std::cout << "Training Nue CC vs Numu PID ..." << std::endl;
    std::cout << " signal file     : " << argv[2]  << std::endl;
    std::cout << " background file : " << argv[3]  << std::endl;

    sigFile = argv[2];
    bkgFile = argv[3];

    wc_trainTMVA_electronMuonCCQE(sigFile, bkgFile);

  } else if(opt=="--nc") {

    std::cout << "Training Nue CC vs CC PID ..."  << std::endl;
    std::cout << " signal file     : " << argv[2] << std::endl;
    std::cout << " background file : " << argv[3] << std::endl;

    sigFile = argv[2];
    bkgFile = argv[3];

    wc_trainTMVA_electronCCQEvsNC(sigFile, bkgFile);

  }// if argc
  
} // main
