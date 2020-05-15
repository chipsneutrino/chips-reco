R__LOAD_LIBRARY(libGeom.so)
R__LOAD_LIBRARY(libEve.so)
R__LOAD_LIBRARY(libMinuit.so)
R__LOAD_LIBRARY(libEG.so)
R__LOAD_LIBRARY(libGui.so)
R__LOAD_LIBRARY(libSpectrum.so)
R__LOAD_LIBRARY(libWCSimRoot.so)
R__LOAD_LIBRARY(libWCSimAnalysisRoot.so)

void wc_trackfitter_cosmic(const char *infile = "", int start = 0, int fit = 1)
{
	// Setup the path to the input file
	TString filename(infile);
	if (filename.CompareTo("") == 0)
	{
		filename = getenv("CHIPSRECO");
		filename += "/config/example/example_sim_output.root";
	}

	// Get a fitter interface...
	WCSimFitterInterface myFitter;

	// Set the input file, and say if you are modifying a WCSimAnalysis output file
	// This loads the data from the corresponding WCSim file appropriately
	myFitter.SetInputFileName(filename.Data(), false);

	// Set up the fit configuration you want to run...
	WCSimFitterConfig *config_cosmic = new WCSimFitterConfig();
	config_cosmic->SetNumTracks(2);
	config_cosmic->SetTrackType(0, "ElectronLike");
	config_cosmic->SetTrackType(1, "MuonLike");
	config_cosmic->SetIsCosmicFit(true);

	config_cosmic->SetParameter(0, "kVtxX", -1200, 1200, 0.0, 10.0, false);
	config_cosmic->SetParameter(0, "kVtxY", -1200, 1200, 0.0, 10.0, false);
	config_cosmic->SetParameter(0, "kVtxZ", -950, 950, 0.0, 10.0, false);
	config_cosmic->SetParameter(0, "kVtxT", -100, 10100, 5000, 1.0, false);
	config_cosmic->SetParameter(0, "kDirTh", 0.0, TMath::Pi(), 0.5 * TMath::Pi(), 0.01, false);
	config_cosmic->SetParameter(0, "kDirPhi", -1.0 * TMath::Pi(), 1.0 * TMath::Pi(), 0, 0.02, false);
	config_cosmic->SetParameter(0, "kEnergy", 500, 5000, 1500, 250.0, false);

	config_cosmic->SetParameter(0, "kVtxX", -1200, 1200, 0.0, 10.0, false);
	config_cosmic->SetParameter(0, "kVtxY", -1200, 1200, 0.0, 10.0, false);
	config_cosmic->SetParameter(0, "kVtxZ", -950, 950, 0.0, 10.0, false);
	config_cosmic->SetParameter(0, "kVtxT", -100, 10100, 5000, 1.0, false);
	config_cosmic->SetParameter(0, "kDirTh", 0.0, TMath::Pi(), 0.5 * TMath::Pi(), 0.01, false);
	config_cosmic->SetParameter(0, "kDirPhi", -1.0 * TMath::Pi(), 1.0 * TMath::Pi(), 0, 0.02, false);
	config_cosmic->SetParameter(1, "kEnergy", 500, 5000, 850, 250.0, false);

	config_cosmic->SetNumEventsToFit(fit);
	config_cosmic->SetFirstEventToFit(start);

	// Make some changes to the default seeding parameters
	WCSimParameters::Instance()->SetSlicerClusterDistance(250);
	WCSimParameters::Instance()->SetSlicerMinSize(25);
	WCSimParameters::Instance()->SetSlicerChargeCut(0.0);
	WCSimParameters::Instance()->SetSlicerTimeCut(100);
	WCSimParameters::Instance()->SetIterateSlicing(true);

	WCSimParameters::Instance()->SetVetoClusterDistance(500);
	WCSimParameters::Instance()->SetVetoMinSize(5);
	WCSimParameters::Instance()->SetVetoMinChargeCut(2.0);
	WCSimParameters::Instance()->SetVetoMaxChargeCut(30.0);
	WCSimParameters::Instance()->SetVetoTimeCut(200);

	// Add the fitter configuration and run the fits...
	myFitter.AddFitterConfig(config_cosmic);
	myFitter.SetMakeFits(kTRUE);
	myFitter.Run();

	delete config_cosmic;

	std::cout << "Done!" << std::endl;
}
