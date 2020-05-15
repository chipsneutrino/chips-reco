R__LOAD_LIBRARY(libGeom.so)
R__LOAD_LIBRARY(libEve.so)
R__LOAD_LIBRARY(libMinuit.so)
R__LOAD_LIBRARY(libEG.so)
R__LOAD_LIBRARY(libGui.so)
R__LOAD_LIBRARY(libSpectrum.so)
R__LOAD_LIBRARY(libWCSimRoot.so)
R__LOAD_LIBRARY(libWCSimAnalysisRoot.so)

void wc_trackfitter_el_plots(const char *infile = "", int start = 0, int fit = 1)
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
	WCSimFitterConfig *config_el = new WCSimFitterConfig();
	config_el->SetNumTracks(1);
	config_el->SetTrackType(0, "ElectronLike");

	config_el->SetParameter(0, "kVtxX", -1200, 1200, 0.0, 10.0, false);
	config_el->SetParameter(0, "kVtxY", -1200, 1200, 0.0, 10.0, false);
	config_el->SetParameter(0, "kVtxZ", -600, 600, 0.0, 10.0, false);
	config_el->SetParameter(0, "kVtxT", -100, 10100, 5000, 1.0, false);
	config_el->SetParameter(0, "kDirTh", 0.0, TMath::Pi(), 0.5 * TMath::Pi(), 0.01, false);
	config_el->SetParameter(0, "kDirPhi", -1.0 * TMath::Pi(), 1.0 * TMath::Pi(), 0, 0.02, false);
	config_el->SetParameter(0, "kEnergy", 500, 5000, 1000, 250.0, false);
	config_el->SetParameter(0, "kConversionDistance", -50, 50, 0, 5.0, true);
	config_el->SetNumEventsToFit(fit);
	config_el->SetFirstEventToFit(start);

	// NOTE: At the moment the plots only work for one fit, so if multiple ones are created
	// the last one will overwrite the plots from the others!
	// Set up the fitter plots you want to be created...
	WCSimFitterPlots *plots = new WCSimFitterPlots();

	// Plot best-fit results
	plots->SetPlotForEachEvent("kVtxX", true);
	plots->SetPlotForEachEvent("kVtxY", true);
	plots->SetPlotForEachEvent("kVtxZ", true);
	plots->SetPlotForEachEvent("kVtxT", true);
	plots->SetPlotForEachEvent("kDirTh", true);
	plots->SetPlotForEachEvent("kDirPhi", true);
	plots->SetPlotForEachEvent("kEnergy", true);

	// Plot reco - true values
	plots->SetPlotRecoMinusTrue("kVtxX", true);
	plots->SetPlotRecoMinusTrue("kVtxY", true);
	plots->SetPlotRecoMinusTrue("kVtxZ", true);
	plots->SetPlotRecoMinusTrue("kVtxT", true);
	plots->SetPlotRecoMinusTrue("kDirTh", true);
	plots->SetPlotRecoMinusTrue("kDirPhi", true);
	plots->SetPlotRecoMinusTrue("kEnergy", true);

	// Sweep out one variable
	plots->SetNumSurfaceBins(20);
	plots->Make1DSurface("kVtxX", true, 0);
	plots->Make1DSurface("kVtxY", true, 0);
	plots->Make1DSurface("kVtxZ", true, 0);
	plots->Make1DSurface("kVtxT", true, 0);
	plots->Make1DSurface("kDirTh", true, 0);
	plots->Make1DSurface("kDirPhi", true, 0);
	plots->Make1DSurface("kEnergy", true, 0);

	// Sweep out two variables
	//plots->Make2DSurface("kVtxX", "kVtxY", true, 0, 0);
	//plots->Make2DSurface("kVtxX", "kEnergy", true, 0, 0);
	//plots->Make2DSurface("kVtxT", "kEnergy", true, 0, 0);

	// Add the fitter configuration and plots, then run the fits...
	myFitter.AddFitterConfig(config_el);
	myFitter.AddFitterPlots(plots);
	myFitter.SetMakeFits(kTRUE);
	myFitter.Run();

	delete config_el;
	delete plots;

	std::cout << "Done!" << std::endl;
}
