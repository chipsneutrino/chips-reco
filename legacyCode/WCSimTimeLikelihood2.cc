/*
 * WCSimTimeLikelihood2.cc
 *
 *  Created on: 8 Jul 2015
 *      Author: andy
 */

#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimTimePredictor.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPMTManager.hh"
#include "WCSimPMTConfig.hh"
#include "WCSimTimeLikelihood2.hh"
#include "TString.h"
#include "TMath.h"
#include <cmath>
#include <map>

#ifndef REFLEX_DICTIONARY
ClassImp (WCSimTimeLikelihood2)
#endif
WCSimTimeLikelihood2::WCSimTimeLikelihood2() {
	// TODO Auto-generated constructor stub
	fTimePredictor = 0x0;
	fLikelihoodDigitArray = 0x0;
	fEmissionProfileManager = 0x0;
	fPMTManager = new WCSimPMTManager();
}

WCSimTimeLikelihood2::WCSimTimeLikelihood2(WCSimLikelihoodDigitArray * myDigitArray,
		WCSimEmissionProfileManager * myEmissionProfileManager) {
	fLikelihoodDigitArray = myDigitArray;
	fEmissionProfileManager = myEmissionProfileManager;
	fTimePredictor = new WCSimTimePredictor(myDigitArray, fEmissionProfileManager);
	fPMTManager = new WCSimPMTManager();
	return;
}

WCSimTimeLikelihood2::~WCSimTimeLikelihood2() {
	// TODO Auto-generated destructor stub
	if (fTimePredictor != 0x0) {
		delete fTimePredictor;
	}
	if (fPMTManager != 0x0) {
		delete fPMTManager;
	}
	// I don't new the likelihooddigitarray so I don't delete it
}

void WCSimTimeLikelihood2::SetTracks(std::vector<WCSimLikelihoodTrackBase*>& myTracks) {
	//std::cout << "Setting tracks" << std::endl;
	fTracks = myTracks;
	fTimePredictor->SetTracks(fTracks);
	fAllPreds.clear();
	//std::cout << "Set tracks - there are " << fTracks.size() << " of them" << std::endl;
}

void WCSimTimeLikelihood2::ResetTracks() {
	//std::cout << "Reset tracks!" << std::endl;
	fTracks.clear();
	fAllPreds.clear();
	fTimePredictor->ClearTracks();
}

void WCSimTimeLikelihood2::ClearTracks() {
	ResetTracks();
}

Double_t WCSimTimeLikelihood2::Calc2LnL(const unsigned int& iDigit) {
	double minusTwoLnL = 0.0;

	if (fAllPreds.size() == 0) {
		fAllPreds = fTimePredictor->GetAllPredictedTimes();
	}
	WCSimLikelihoodDigit * myDigit = fLikelihoodDigitArray->GetDigit(iDigit);
	if (IsGoodDigit(myDigit)) {
		double predictedFirstHitTime = fAllPreds.at(iDigit);
		double timeResolution = GetPMTTimeResolution(myDigit);
		double actualHitTime = myDigit->GetT();

		minusTwoLnL += GetGaussianMinusTwoLnL(actualHitTime, predictedFirstHitTime, timeResolution);
	}
	return minusTwoLnL;
}

Double_t WCSimTimeLikelihood2::Calc2LnL() {
	double minusTwoLnL = 0.0;
	//std::cout << "Getting all predicted times" << std::endl;
	if (fAllPreds.size() == 0) {
		fAllPreds = fTimePredictor->GetAllPredictedTimes();
	}
	//std::cout << "Got them" << std::endl;
	for (int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit) {
		WCSimLikelihoodDigit * myDigit = fLikelihoodDigitArray->GetDigit(iDigit);
		if (myDigit->GetQ() <= 0.0) {
			continue;
		}
		double timeResolution = GetPMTTimeResolution(myDigit);
		double actualHitTime = myDigit->GetT();
		minusTwoLnL += GetGaussianMinusTwoLnL(actualHitTime, fAllPreds.at(iDigit), timeResolution);

	}
	return minusTwoLnL;
}

std::vector<Double_t> WCSimTimeLikelihood2::GetAllPredictedTimes() {
	if (fAllPreds.size() == 0) {
		fAllPreds = fTimePredictor->GetAllPredictedTimes();
	}
	return fAllPreds;
}

double WCSimTimeLikelihood2::GetPMTTimeResolution(const unsigned int& iDigit) {
	WCSimLikelihoodDigit * myDigit = fLikelihoodDigitArray->GetDigit(iDigit);
	return GetPMTTimeResolution(myDigit);
}

double WCSimTimeLikelihood2::GetPMTTimeResolution(WCSimLikelihoodDigit * myDigit) {
	TString pmtName = myDigit->GetPMTName();
	double timeConstant = 0.0;
	if (fPMTTimeConstantMap.find(pmtName) == fPMTTimeConstantMap.end()) {
		WCSimPMTConfig config = fPMTManager->GetPMTByName(std::string(pmtName.Data()));
		timeConstant = config.GetTimeConstant();
		fPMTTimeConstantMap[pmtName] = timeConstant;
	}
	return 0.33 + sqrt(fPMTTimeConstantMap[pmtName] / myDigit->GetQ()); // This is what WCSim does...
}

double WCSimTimeLikelihood2::GetGaussianMinusTwoLnL(const double& x, const double& mean, const double& sigma) {
	// If the PMT is hit but we don't expect it to be, assume we have scattered light:
	double scatteringProb = 0.005; // 0.5% chance
	double minus2LnL = 0.0;

	if (mean <= 0 && x > 0) // PMT hit when it wasn't expected to: assume scattering
			{
		minus2LnL = -2.0 * TMath::Log(scatteringProb);
	}

	if (mean > 0 && x >= 0) {
		// Calculate the likelihood:

		double tmp = 2 * TMath::Log(sqrt(2 * TMath::Pi()) * sigma) + ((x - mean) * (x - mean) / (sigma * sigma));
		if (tmp < 25) {
			minus2LnL = tmp;
		} else {
			minus2LnL = 25;
		}
	}
	//std::cout << mean << "  " << x << "  " << minus2LnL << std::endl;
	return minus2LnL;
}

bool WCSimTimeLikelihood2::IsGoodDigit(WCSimLikelihoodDigit * myDigit) {
	return (myDigit->GetQ() > 0);
}
