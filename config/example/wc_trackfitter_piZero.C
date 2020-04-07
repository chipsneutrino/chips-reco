R__LOAD_LIBRARY(libWCSimRoot.so)
R__LOAD_LIBRARY(libWCSimAnalysisRoot.so)
R__LOAD_LIBRARY(libGeom.so)
R__LOAD_LIBRARY(libEve.so)
R__LOAD_LIBRARY(libMinuit.so)

void wc_trackfitter_piZero(const char *infile = "", int start = 0, int fit = 100)
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
	WCSimFitterConfig *config_piZero = new WCSimFitterConfig();
	config_piZero->SetNumTracks(2);
	config_piZero->SetTrackType(0, "PhotonLike");
	config_piZero->SetTrackType(1, "PhotonLike");
	config_piZero->SetIsPiZeroFit(true);
	config_piZero->SetForcePiZeroMass(false);

	config_piZero->SetParameter(0, "kVtxX", -1200, 1200, 0.0, 10.0, false);
	config_piZero->SetParameter(0, "kVtxY", -1200, 1200, 0.0, 10.0, false);
	config_piZero->SetParameter(0, "kVtxZ", -600, 600, 0.0, 10.0, false);
	config_piZero->SetParameter(0, "kVtxT", -100, 10100, 5000, 1.0, false);
	config_piZero->SetParameter(0, "kDirTh", 0.0, TMath::Pi(), 0.5 * TMath::Pi(), 0.01, false);
	config_piZero->SetParameter(0, "kDirPhi", -1.0 * TMath::Pi(), 1.0 * TMath::Pi(), 0, 0.02, false);
	config_piZero->SetParameter(0, "kEnergy", 500, 5000, 600, 250.0, false);
	config_piZero->SetParameter(0, "kConversionDistance", 0, 75, 50, 5.0, false);

	config_piZero->SetParameter(1, "kVtxX", -1200, 1200, 0.0, 10.0, false);
	config_piZero->SetParameter(1, "kVtxY", -1200, 1200, 0.0, 10.0, false);
	config_piZero->SetParameter(1, "kVtxZ", -600, 600, 0.0, 10.0, false);
	config_piZero->SetParameter(1, "kVtxT", -100, 10100, 5000, 1.0, false);
	config_piZero->SetParameter(1, "kDirTh", 0.0, TMath::Pi(), 0.5 * TMath::Pi(), 0.01, false);
	config_piZero->SetParameter(1, "kDirPhi", -1.0 * TMath::Pi(), 1.0 * TMath::Pi(), 0, 0.02, false);
	config_piZero->SetParameter(1, "kEnergy", 500, 5000, 600, 250.0, false);
	config_piZero->SetParameter(1, "kConversionDistance", 0, 75, 50, 5.0, false);

	config_piZero->JoinParametersTogether(0, 1, "kVtxX");
	config_piZero->JoinParametersTogether(0, 1, "kVtxY");
	config_piZero->JoinParametersTogether(0, 1, "kVtxZ");
	config_piZero->JoinParametersTogether(0, 1, "kVtxT");

	config_piZero->SetNumEventsToFit(fit);
	config_piZero->SetFirstEventToFit(start);

	// Add the fitter configuration and run the fits...
	myFitter.AddFitterConfig(config_piZero);
	myFitter.SetMakeFits(kTRUE);
	myFitter.Run();

	delete config_piZero;

	std::cout << "Done!" << std::endl;
}
