void wc_tabulate_integrals(const char * filename="", bool includeS = false)  
{
  // Path to WCSim ROOT file
  // =======================
  TString filenameStr(filename);
  
  
  gApplication->ProcessLine(".except");

  // Load libraries
  // ==============
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
  TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
  gSystem->Load(libWCSimRoot.Data());
  gSystem->Load(libWCSimAnalysis.Data());

  TrackType::Type type = TrackType::MuonLike;
  if( includeS)
  {
    WCSimIntegralLookupMaker wcilm(type,
                                 227, 0.0, 3405.0,
                                 400, -1.0, 1.0
                                 );
    wcilm.Run(filename);
  }
  else
  {
    WCSimIntegralLookupMaker3D wcilm3D(type, 
                                       227, 0.0, 3405.0,
                                       200, -1, 1.
                                       );
    wcilm3D.Run(filename);
  }
}
 
