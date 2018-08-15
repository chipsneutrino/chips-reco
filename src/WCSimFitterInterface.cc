/*
 * WCSimFitterInterface.cc
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#include "WCSimFitterInterface.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimOutputTree.hh"
#include "WCSimPiZeroFitter.hh"
#include "WCSimCosmicFitter.hh"
#include <TString.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <cassert>
#include <stdlib.h>

ClassImp (WCSimFitterInterface)

WCSimFitterInterface::WCSimFitterInterface() :
		fNumFits(0), fModifyInputFile(false), fFileName(""), fFitter(0x0), fPiZeroFitter(0x0), fCosmicFitter(0x0), fOutputTree(
				0x0), fMakeFits(true) {

	// Clear the fitterConfig vector
	fFitterConfigs.clear();
	fOutputName = "";

	fMakePlots = false;
	fPlotsName = "";
	fFitterPlots = 0x0;

	// TODO Merge Fitter plots into a branch of the outputTree
	fOutputTree = new WCSimOutputTree();
}

WCSimFitterInterface::~WCSimFitterInterface() {
	// TODO Auto-generated destructor stub
	if (fFitter != NULL) {
		delete fFitter;
	}
	if (fPiZeroFitter != NULL) {
		delete fPiZeroFitter;
	}
	if (fCosmicFitter != NULL) {
		delete fCosmicFitter;
	}
	if (fOutputTree != NULL) {
		delete fOutputTree;
	}

	fFitterConfigs.clear();
}

void WCSimFitterInterface::SetInputFileName(const char * inputfile, bool modifyFile) {
	fModifyInputFile = modifyFile;
	char * fullPath = realpath(inputfile, NULL);
	fFileName = TString(fullPath);
	free(fullPath);

	// Load the data from the corresponding WCSim file
	LoadWCSimData();
}

void WCSimFitterInterface::SetMakeFits(bool makeFits) {
	fMakeFits = makeFits;
}

void WCSimFitterInterface::AddFitterConfig(WCSimFitterConfig * config) {
	std::cout << "Adding a fitter configuration to run order..." << std::endl;
	fNumFits += 1;
	fFitterConfigs.push_back(config);
}

void WCSimFitterInterface::AddFitterPlots(WCSimFitterPlots * plots) {
	fMakePlots = true;
	fFitterPlots = plots;
}

void WCSimFitterInterface::InitFitter(WCSimFitterConfig * config) {
	// This should reset a fitter in case it has been used already
	// It should then load the new configuration
	std::cout << "Initialise fitter..." << std::endl;
	config->Print();

	// Delete the old fitter and make a new one for each fit type...
	if (config->GetIsCosmicFit()) {
		if (fCosmicFitter != 0x0) {
			delete fCosmicFitter;
			fCosmicFitter = 0x0;
		}
		fCosmicFitter = new WCSimCosmicFitter(config);
		fCosmicFitter->SetOutputTree(fOutputTree);
		if (fMakePlots) {
			fCosmicFitter->SetFitterPlots(fFitterPlots);
		}
	} else if (config->GetIsPiZeroFit()) {
		if (fPiZeroFitter != 0x0) {
			delete fPiZeroFitter;
			fPiZeroFitter = 0x0;
		}
		fPiZeroFitter = new WCSimPiZeroFitter(config);
		fPiZeroFitter->SetOutputTree(fOutputTree);
		if (fMakePlots) {
			fPiZeroFitter->SetFitterPlots(fFitterPlots);
		}
	} else {
		if (fFitter != 0x0) {
			delete fFitter;
			fFitter = 0x0;
		}
		fFitter = new WCSimLikelihoodFitter(config);
		fFitter->SetOutputTree(fOutputTree);
		if (fMakePlots) {
			fFitter->SetFitterPlots(fFitterPlots);
		}
	}
}

void WCSimFitterInterface::LoadWCSimData() {
	std::cout << " *** WCSimFitterInterface::LoadWCSimData() *** " << std::endl;
	TString wcsimFile = fFileName;
	if (fModifyInputFile) {
		TFile * lookupFile = new TFile(fFileName.Data(), "READ");

		std::cout << "Checking the existing file..." << std::endl;
		assert(lookupFile->GetListOfKeys()->Contains("fResultsTree"));

		TTree * tree = (TTree*) (lookupFile->Get("fResultsTree"));

		//Set up the EventInfo...
		EventHeader *eventHeader = new EventHeader();
		TBranch *b_eh = tree->GetBranch("EventHeader");
		b_eh->SetAddress(&eventHeader);
		b_eh->GetEntry(0);

		wcsimFile = eventHeader->GetInputFile();

		eventHeader->Clear();
		lookupFile->Close();
	}

	WCSimInterface::LoadData(wcsimFile.Data());
}

void WCSimFitterInterface::InitFitterPlots(TString outputName, WCSimFitterConfig * config) {
	std::cout << " *** WCSimFitterInterface::InitFitterPlots() *** " << std::endl;
	fFitterPlots->SetSaveFileName(outputName.Data());
	fFitterPlots->MakeSaveFile();
	fFitterPlots->Print();
	fFitterPlots->MakeHistograms(config);
}

void WCSimFitterInterface::Run() {
	std::cout << " *** WCSimFitterInterface::Run() *** " << std::endl;
	// I want to loop through the fitter configurations and run the fits
	for (int fitNum = 0; fitNum < fNumFits; fitNum++) {
		// 1) Initialise the outputTree and fitterPlots
		if (fitNum == 0) {
			fOutputName = fOutputTree->InitOutputTree(fFileName, fModifyInputFile, fFitterConfigs[fitNum]);
			if (fMakePlots) {
				InitFitterPlots(fOutputName.ReplaceAll("tree.root", "plots.root"), fFitterConfigs[fitNum]);
			}
		} else {
			fOutputName = fOutputTree->InitOutputTree(fOutputName, fModifyInputFile, fFitterConfigs[fitNum]);
		}

		// 2) Initialise the fitter
		InitFitter (fFitterConfigs[fitNum]);

		// 3) Run the fit and fill the outputTree and fitterPlots
		if (fMakeFits) {
			if (fFitterConfigs[fitNum]->GetIsPiZeroFit()) {
				std::cout << " *** Running PiZero Fit *** " << std::endl;
				fPiZeroFitter->RunFits();
			} else if (fFitterConfigs[fitNum]->GetIsCosmicFit()) {
				std::cout << " *** Running Cosmic Fit *** " << std::endl;
				fCosmicFitter->RunFits();
			} else {
				std::cout << " *** Running Fit *** " << std::endl;
				fFitter->RunFits();
			}
		}

		// 4) If making fitterPlots fill the surfaces and write/save the plots to file
		if (fMakePlots) {
			fFitter->RunLikelihoodPlots();
			fFitter->RunSurfaces();
			fFitterPlots->SaveProfiles();
			fFitterPlots->SavePlots();
		}

		// 5) Write the tree and save the file...
		fOutputTree->SaveTree();

		// 6) Need to set the fModifyInputFile to true so the outputTree modifies the file for next fit
		fModifyInputFile = true;

		// 7) Delete and get a new output tree
		if (fOutputTree != NULL) {
			delete fOutputTree;
			fOutputTree = new WCSimOutputTree();
		}
	}
}

