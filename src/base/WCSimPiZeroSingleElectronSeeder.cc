/*
 * WCSimPiZeroSingleElectronSeeder.cc
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#include "WCSimPiZeroSingleElectronSeeder.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimLikelihoodTrackFactory.hh"

#ifndef REFLEX_DICTIONARY
ClassImp (WCSimPiZeroSingleElectronSeeder)
#endif

WCSimPiZeroSingleElectronSeeder::WCSimPiZeroSingleElectronSeeder(WCSimFitterConfig * config) :
		WCSimPiZeroSeeder(config) {
	fSingleElectronSeed = 0x0;
	fMinus2LnL = 0;
	fMadeSeeds = false;
}

WCSimPiZeroSingleElectronSeeder::~WCSimPiZeroSingleElectronSeeder() {
	if (fSingleElectronSeed != 0x0) {
		delete fSingleElectronSeed;
	}
}

WCSimLikelihoodTrackBase* WCSimPiZeroSingleElectronSeeder::GetSeed() {
	if (!fMadeSeeds) {
		MakeSeeds();
	}
	return fSingleElectronSeed;
}

void WCSimPiZeroSingleElectronSeeder::MakeSeeds() {
	if (fMadeSeeds) {
		return;
	}
	WCSimFitterConfig config;
	WCSimFitterConfig * backupConfig = fFitterConfig;
	config.SetNumTracks(1);

	config.SetTrackType(0, "ElectronLike");
	config.SetParameter(0, "kVtxX", fFitterConfig->GetParMin(0, "kVtxX"), fFitterConfig->GetParMax(0, "kVtxX"),
			fFitterConfig->GetParStart(0, "kVtxX"), fFitterConfig->GetParStep(0, "kVtxX"), false);
	config.SetParameter(0, "kVtxY", fFitterConfig->GetParMin(0, "kVtxY"), fFitterConfig->GetParMax(0, "kVtxY"),
			fFitterConfig->GetParStart(0, "kVtxY"), fFitterConfig->GetParStep(0, "kVtxY"), false);
	config.SetParameter(0, "kVtxZ", fFitterConfig->GetParMin(0, "kVtxZ"), fFitterConfig->GetParMax(0, "kVtxZ"),
			fFitterConfig->GetParStart(0, "kVtxZ"), fFitterConfig->GetParStep(0, "kVtxZ"), false);
	config.SetParameter(0, "kVtxT", fFitterConfig->GetParMin(0, "kVtxT"), fFitterConfig->GetParMax(0, "kVtxT"),
			fFitterConfig->GetParStart(0, "kVtxT"), fFitterConfig->GetParStep(0, "kVtxT"), false);
	config.SetParameter(0, "kDirTh", fFitterConfig->GetParMin(0, "kDirTh"), fFitterConfig->GetParMax(0, "kDirTh"),
			fFitterConfig->GetParStart(0, "kDirTh"), fFitterConfig->GetParStep(0, "kDirTh"), false);
	config.SetParameter(0, "kDirPhi", fFitterConfig->GetParMin(0, "kDirPhi"), fFitterConfig->GetParMax(0, "kDirPhi"),
			fFitterConfig->GetParStart(0, "kDirPhi"), fFitterConfig->GetParStep(0, "kDirPhi"), false);
	config.SetParameter(0, "kEnergy", fFitterConfig->GetParMin(0, "kEnergy"), fFitterConfig->GetParMax(0, "kEnergy"),
			0.5 * (fFitterConfig->GetParMin(0, "kEnergy") + fFitterConfig->GetParMax(0, "kEnergy")),
			fFitterConfig->GetParStep(0, "kEnergy"), false);
	config.SetParameter(0, "kConversionDistance", fFitterConfig->GetParMin(0, "kConversionDistance"),
			fFitterConfig->GetParMax(0, "kConversionDistance"), 0, fFitterConfig->GetParStep(0, "kConversionDistance"),
			true);

	fFitterConfig = &config;
	fFitterTrackParMap = WCSimFitterTrackParMap(&config);
	fFitterTrackParMap.Set();

	std::cout << "Seed single electron" << std::endl;
	SeedEvent();
	std::cout << "Fit single electron track" << std::endl;
	FitEventNumber (fEvent);
	std::cout << "Done single electron fit" << std::endl;

	if (fSingleElectronSeed != 0x0) {
		delete fSingleElectronSeed;
		fSingleElectronSeed = 0x0;
	}

	fSingleElectronSeed = WCSimLikelihoodTrackFactory::MakeTrack(fFitterTrackParMap.GetTrackType(0),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kVtxX)),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kVtxY)),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kVtxZ)),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kVtxT)),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kDirTh)),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kDirPhi)),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kEnergy)),
			fFitterTrackParMap.GetCurrentValue(std::make_pair(0, FitterParameterType::kConversionDistance)));

	fFitterConfig = backupConfig;
	fFitterTrackParMap = WCSimFitterTrackParMap(fFitterConfig);
	fFitterTrackParMap.Set();
	fMadeSeeds = true;
	return;
}

double WCSimPiZeroSingleElectronSeeder::GetMinus2LnL() {
	if (!fMadeSeeds) {
		MakeSeeds();
	}
	return fMinus2LnL;
}

double WCSimPiZeroSingleElectronSeeder::FitSingleTrack() {
	std::cout << "In WCSimPiZeroSingleElectronSeeder::FitSingleTrack" << std::endl;
	double minus2LnL = FitAndGetLikelihood("Simplex");
	std::cout << "minus2LnL for single electron = " << minus2LnL << std::endl;
	return minus2LnL;
}

void WCSimPiZeroSingleElectronSeeder::SetEvent(const int &iEvent) {
	if (fSingleElectronSeed != 0x0) {
		delete fSingleElectronSeed;
		fSingleElectronSeed = 0x0;
	}
	WCSimPiZeroSeeder::SetEvent(iEvent);
}
