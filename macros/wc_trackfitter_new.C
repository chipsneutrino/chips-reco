#include "WCSimFitterInterface.hh"
#include "TString.h"
#include <iostream>
#include <gperftools/profiler.h>


void PrintHelp();

void wc_trackfitter_new(const char * infile = ""){

  // Path to WCSim ROOT file
  // =======================
  TString filename(infile);
  if(filename.EqualTo("") )
  {
    filename = TString("localfile.root");
  }

  // Load Data
  // =========
  WCSimInterface::LoadData(filename.Data());
  
  WCSimFitterInterface::Instance()->SetInputFileName(filename.Data()); // For inclusion in the fitterPlots file
  WCSimFitterInterface::Instance()->SetNumTracks(1);
  WCSimFitterInterface::Instance()->SetTrackType(0, "MuonLike");

  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxX", -1200, 1200, 0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxY", -1200, 1200, 0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxZ", -900, 900, 0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxT", -900, 900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirTh", 0, TMath::Pi(), 0.5*TMath::Pi(), false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirPhi", -TMath::Pi(), TMath::Pi(), 0.25*TMath::Pi(), false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kEnergy", 800, 3000, 1500, false);


  // Plot best-fit results
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxX");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxY");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxZ");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxT");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kDirTh");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kDirPhi");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kEnergy");

  // Plot reco - true values
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxX");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxY");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxZ");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxT");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kDirTh");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kDirPhi");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kEnergy");

  // Sweep out one variable
  WCSimFitterInterface::Instance()->SetNumSurfaceBins(50);
  WCSimFitterInterface::Instance()->Make1DSurface(0, "kEnergy");
  WCSimFitterInterface::Instance()->SetMakeFits(kTRUE);
  WCSimFitterInterface::Instance()->SetMakeSurfaces(kFALSE);



  WCSimFitterInterface::Instance()->Print();
  WCSimFitterInterface::Instance()->SetNumEventsToFit(0);
  WCSimFitterInterface::Instance()->SetFirstEventToFit(0);
  ProfilerStart("profileFitter");
  WCSimFitterInterface::Instance()->Run();
  ProfilerStop();
  
  std::cout << "Done!" << std::endl;

}

void PrintHelp(){
	std::cout << "Usage instructions for fitterTest" << std::endl;
  std::cout << "This program is supposed to run the fitter over a single file so it can be profiled" << std::endl;
  std::cout << "It basically wraps around the macros/wc_trackfitter_new.C macro - normally you just want that" << std::endl;
	std::cout << "\t-h Displays the help message" << std::endl;
	std::cout << "\t-f <filename> Supply an input file" << std::endl;
}
