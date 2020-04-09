R__LOAD_LIBRARY(libWCSimRoot.so)
R__LOAD_LIBRARY(libWCSimAnalysisRoot.so)
R__LOAD_LIBRARY(libGeom.so)
R__LOAD_LIBRARY(libEve.so)
R__LOAD_LIBRARY(libMinuit.so)

void wc_trackfitter_mu(const char *infile = "", int start = 0, int fit = 1)
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
	WCSimFitterConfig *config_mu = new WCSimFitterConfig();
	config_mu->SetNumTracks(1);
	config_mu->SetTrackType(0, "MuonLike");

	config_mu->SetParameter(0, "kVtxX", -1200, 1200, 0.0, 10.0, false);
	config_mu->SetParameter(0, "kVtxY", -1200, 1200, 0.0, 10.0, false);
	config_mu->SetParameter(0, "kVtxZ", -600, 600, 0.0, 10.0, false);
	config_mu->SetParameter(0, "kVtxT", -100, 10100, 5000, 1.0, false);
	config_mu->SetParameter(0, "kDirTh", 0.0, TMath::Pi(), 0.5 * TMath::Pi(), 0.01, false);
	config_mu->SetParameter(0, "kDirPhi", -1.0 * TMath::Pi(), 1.0 * TMath::Pi(), 0, 0.02, false);
	config_mu->SetParameter(0, "kEnergy", 500, 5000, 1000, 250.0, false);
	config_mu->SetParameter(0, "kConversionDistance", -50, 50, 0, 5.0, true);
	config_mu->SetNumEventsToFit(fit);
	config_mu->SetFirstEventToFit(start);

	// Add the fitter configuration and run the fits...
	myFitter.AddFitterConfig(config_mu);
	myFitter.SetMakeFits(kTRUE);
	myFitter.Run();

	delete config_mu;

	std::cout << "Done!" << std::endl;
}
