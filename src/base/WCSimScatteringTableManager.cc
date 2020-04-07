#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>

#include "WCSimScatteringTableManager.hh"

#include "TFile.h"
#include "THnSparse.h"
#include "TString.h"
#include "TAxis.h"

WCSimScatteringTableManager* WCSimScatteringTableManager::fManager = 0x0;

WCSimScatteringTableManager::WCSimScatteringTableManager() {

	fScatteringFile = 0x0;
	fScatteringTableTop = 0x0;
	fScatteringTableBarrel = 0x0;
	fScatteringTableBottom = 0x0;

	this->LoadScatteringTables();
}

WCSimScatteringTableManager::~WCSimScatteringTableManager() {

	if (fScatteringFile != 0x0)
		delete fScatteringFile;
	if (fScatteringTableTop != 0x0)
		delete fScatteringTableTop;
	if (fScatteringTableBarrel != 0x0)
		delete fScatteringTableBarrel;
	if (fScatteringTableBottom != 0x0)
		delete fScatteringTableBottom;

}

WCSimScatteringTableManager* WCSimScatteringTableManager::Instance() {

	// Create an instance if we don't have one.
	if (fManager == 0x0) {
		fManager = new WCSimScatteringTableManager();
	}

	return fManager;
}

float WCSimScatteringTableManager::GetScatteringValue(std::vector<float> &values, int pmtLocation) {
	// Since we need an array to query the THnSparse, convert now.
	if (values.size() != 6) {
		std::cerr << "Bin array must have 6 entries" << std::endl;
		assert(0);
	}
	int binArray[6];
	for (unsigned int i = 0; i < values.size(); ++i) {
		binArray[i] = CalculateBin(i, values[i], pmtLocation);
	}
	return this->GetScatteringValue(binArray, pmtLocation);
}

float WCSimScatteringTableManager::GetScatteringValue(int *bin, int pmtLocation) {
	// Check that the pmtLocation makes sense.
	if (pmtLocation < 0 || pmtLocation > 2) {
		std::cerr << "PMT location must be 0 (top), 1 (barrel) or 2 (bottom)" << std::endl;
		return -999.9;
	}

	float scatteringValue = -999.9;
	if (pmtLocation == 0) {
		scatteringValue = fScatteringTableTop->GetBinContent(bin);
	} else if (pmtLocation == 1) {
		scatteringValue = fScatteringTableBarrel->GetBinContent(bin);
	} else {
		scatteringValue = fScatteringTableBottom->GetBinContent(bin);
	}

	return scatteringValue;
}

void WCSimScatteringTableManager::LoadScatteringTables() {

	// Todo: Should we have a geometry name based naming system?
	// - This would allow multiple scattering tables to exist and the code
	//   will choose the correct one based on the geometry settings.
	if (fScatteringFile != 0x0) {
		delete fScatteringFile;
		fScatteringFile = 0x0;
	}
	TString profileFileName = TString(getenv("CHIPSRECO"));
	profileFileName.Append("/config/scatteringTable.root");
	fScatteringFile = new TFile(profileFileName.Data());
	if (!fScatteringFile) {
		std::cerr << "Can't find the scattering table file." << std::endl;
		assert(0);
	}

	// Look for the three histograms we need.
	fScatteringTableTop = (THnSparseF*) fScatteringFile->Get("hScatteringTableTop");
	fScatteringTableBarrel = (THnSparseF*) fScatteringFile->Get("hScatteringTableBarrel");
	fScatteringTableBottom = (THnSparseF*) fScatteringFile->Get("hScatteringTableBottom");
	if (!fScatteringTableTop || !fScatteringTableBottom || !fScatteringTableBarrel) {
		std::cerr << "One or more scattering table histograms missing." << std::endl;
		assert(0);
	}

}

int WCSimScatteringTableManager::CalculateBin(int dimension, float value, int pmtLocation) {

	int bin = -1;

	THnSparseF* thisTable = 0x0;
	if (pmtLocation == 0)
		thisTable = fScatteringTableTop;
	if (pmtLocation == 1)
		thisTable = fScatteringTableBarrel;
	if (pmtLocation == 2)
		thisTable = fScatteringTableBottom;

	TAxis* axis = thisTable->GetAxis(dimension);
	bin = axis->FindBin(value);

	return bin;

}

