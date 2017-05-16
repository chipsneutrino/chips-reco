void LoadLibs(){

  //Load  Root related Libraries 
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit.so");

  //WCSim And WCSimAnalysis Libraries
  gSystem->Load("$WCSIMHOME/libWCSimRoot.so");
  gSystem->Load("$WCSIMANAHOME/lib/libWCSimAnalysis.so");
  

  //Add Include Paths
  gSystem->AddIncludePath("-I$WCSIMHOME/include/");
  gSystem->AddIncludePath("-I$WCSIMANAHOME/include/");

  gInterpreter->AddIncludePath("$WCSIMHOME/include/");
  gInterpreter->AddIncludePath("$WCSIMANAHOME/include/");

  // Load wc_readTMVA macro and compile if necessary
  gROOT->ProcessLine(".L wc_readTMVA.C+");

}
