{

  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");
  gSystem->Load("../WCSim/libWCSimRoot.so");
  gSystem->Load("./lib/libWCSimAnalysis.so");

  WCSimEventWriter* writer = new WCSimEventWriter();

  writer->AddFile("../mcSamples/numi_numu_NC_100kT_10pC_Leigh.root");
  // writer->AddFile("/lbne/wcsim_root_files/DUSEL_100kton_10inch_40perCent/muon_vtx_20001_wcsim.root");
  // writer->AddFile("/lbne/wcsim_root_files/DUSEL_150ktonMailbox_10inch_30perCent/muon_vtx_20001mailbox.root");
  // writer->AddFile("/r01/lbne/water/wcsim_root_files/LBNE_100kton_10inch_15perCent/beamevents0.root");
  // writer->AddFile("/r01/lbne/water/wcsim_root_files/DUSEL_100kton_10inch_40perCent/muon_vtx_12001_wcsim.root");
  // writer->AddFile("/r01/lbne/water/wcsim_root_files_cosmic/DUSEL_100kton_10inch_15perCent/cosmic_DUSEL_4850ft_fix_long_100kton.001.vec*.root");
  //writer->AddFile("/r01/lbne/water/wcsim_root_files/DUSEL_100kton_10inch_15perCent/electron_vtx_12001.root");
  // writer->AddFile("/r01/lbne/water/wcsim_root_files_cosmic/DUSEL_200kton/cosmic_DUSEL_4850ft_surface_200kton.019.root");
  //writer->AddFile("/lbne/software2010/WCSim_588/wcsim.cosmic.19.10.root");

  //writer->WriteOutGeometry();
  //writer->WriteOutTracks();
  writer->WriteOutDigits();

  writer->Run(20);


}
