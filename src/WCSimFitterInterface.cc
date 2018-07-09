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
		fNumFits(0), fModifyInputFile(false), fFileName(""), fFitter(0x0), fPiZeroFitter(
				0x0), fCosmicFitter(0x0), fOutputTree(0x0), fMakeFits(true) {

	// Clear the fitterConfig vector
	fFitterConfigs.clear();
	fOutputName = "";

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

void WCSimFitterInterface::SetInputFileName(const char * inputfile,
		bool modifyFile) {
	fModifyInputFile = modifyFile;
	char * fullPath = realpath(inputfile, NULL);
	fFileName = TString(fullPath);
	free(fullPath);
}

void WCSimFitterInterface::SetMakeFits(bool makeFits) {
	fMakeFits = makeFits;
}

void WCSimFitterInterface::AddFitterConfig(WCSimFitterConfig * config) {
	std::cout << "Adding a fitter configuration to run order..." << std::endl;
	fNumFits += 1;
	fFitterConfigs.push_back(config);
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
	} else if (config->GetIsPiZeroFit()) {
		if (fPiZeroFitter != 0x0) {
			delete fPiZeroFitter;
			fPiZeroFitter = 0x0;
		}
		fPiZeroFitter = new WCSimPiZeroFitter(config);
		fPiZeroFitter->SetOutputTree(fOutputTree);
	} else {
		if (fFitter != 0x0) {
			delete fFitter;
			fFitter = 0x0;
		}
		fFitter = new WCSimLikelihoodFitter(config);
		fFitter->SetOutputTree(fOutputTree);
	}
}

void WCSimFitterInterface::Run() {
	std::cout << " *** WCSimFitterInterface::Run() *** " << std::endl;
	// I want to loop through the fitter configurations and run the fits
	for (int fitNum = 0; fitNum < fNumFits; fitNum++) {
		// 1) First we initialise the outputTree for this fit...
		if (fitNum == 0) {
			fOutputName = fOutputTree->InitOutputTree(fFileName,
					fModifyInputFile, fFitterConfigs[fitNum]);
		} else {
			fOutputName = fOutputTree->InitOutputTree(fOutputName,
					fModifyInputFile, fFitterConfigs[fitNum]);
		}

		// 2) Initialise the fitter...
		InitFitter (fFitterConfigs[fitNum]);

		// 3) Run the fit and fill the outputTree...
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

		// 4) Write the tree and save the file...
		fOutputTree->SaveTree();

		// ?) Need to set the fModifyInputFile to true so the outputTree modifies the file for next fit
		fModifyInputFile = true;

		if (fOutputTree != NULL) {
			delete fOutputTree;
			fOutputTree = new WCSimOutputTree();
		}

	}

}
