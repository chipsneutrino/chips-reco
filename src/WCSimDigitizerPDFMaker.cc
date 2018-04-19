/*
 * WCSimDigitizerPDFMaker.cpp
 *
 *  Created on: 5 Jan 2015
 *      Author: andy
 */

#include "WCSimDigitizerPDFMaker.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TMath.h"

#include <algorithm>
#include <iostream>
#include <cmath>

WCSimDigitizerPDFMaker::WCSimDigitizerPDFMaker() {

	// Set up the random number generator
	fRandom = new TRandom3();
	fRandom->SetSeed(0);

	fNumZero = 0;
	fMu = 0;
	fNumThrows = 1000000;
	fType = 0; // default is sk1pe for now
	totPoisson = true;

	fProbHisto = 0x0;
	fDebug = 0x0;
	SetBinning();

	fSK1pePMT = new WCSimSK1pePMT();
	fCHIPSPMT = new WCSimCHIPSPMT();
	fTOTPMT = new WCSimTOTPMT();
}

WCSimDigitizerPDFMaker::~WCSimDigitizerPDFMaker() {
	if(fProbHisto != 0x0) { delete fProbHisto; }
	delete fRandom;
}

void WCSimDigitizerPDFMaker::SetBinning(int bins, double min, double max) {
	fNumChargeBins = bins;
	fChargeMin = min;
	fChargeMax = max;
}

int WCSimDigitizerPDFMaker::GetBins() const {
	return fNumChargeBins;
}

double WCSimDigitizerPDFMaker::GetMin() const {
	return fChargeMin;
}

double WCSimDigitizerPDFMaker::GetMax() const {
	return fChargeMax;
}

void WCSimDigitizerPDFMaker::SetMu(double mu) {
	fMu = mu;
}

double WCSimDigitizerPDFMaker::GetMu() const {
	return fMu;
}

void WCSimDigitizerPDFMaker::SetNumThrows(int num) {
	fNumThrows = num;
}

int WCSimDigitizerPDFMaker::GetNumThrows() const {
	return fNumThrows;
}

void WCSimDigitizerPDFMaker::SetType(int type) {
	fType = type;
}

int WCSimDigitizerPDFMaker::GetType() const {
	return fType;
}

void WCSimDigitizerPDFMaker::SetTotPoisson(bool set) {
	totPoisson = set;
}

bool WCSimDigitizerPDFMaker::GetTotPoisson() const {
	return totPoisson;
}

void WCSimDigitizerPDFMaker::Run() {
	if(fType != 0 && fType != 1 && fType != 2){
		std::cerr << "Type not recognised -> " << fType << std::endl;
		return;
	}
	MakeHisto();
	//std::cout << "WCSimDigitizerPDFMaker::Run() LoopDigitizing..." << std::endl;
	LoopDigitize();
  	FillEmptyBins();
	//NormHistogram();
	//lnHist();
	SaveHistogram();
	//std::cout << fNumZero << std::endl;
}

int WCSimDigitizerPDFMaker::ThrowPoisson() {

	int myRandom = fRandom->Poisson(fMu);
	fDebug->Fill(myRandom);
	return myRandom;
}

void WCSimDigitizerPDFMaker::LoopDigitize() {
	for(int iThrow = 0; iThrow < fNumThrows; ++iThrow)
	{
		Digitize();
	}
}

void WCSimDigitizerPDFMaker::Digitize() {

	// Throw from a poisson with the mean at fMu
	int nPhotons;
	double peSmeared = 0.0;

	if(fType == 0){
		nPhotons = this->ThrowPoisson();
		peSmeared = fSK1pePMT->CalculateCharge(nPhotons);
	}
	else if (fType == 1){
		nPhotons = this->ThrowPoisson();
		peSmeared = fCHIPSPMT->CalculateCharge(nPhotons, 0, 0);
	} // TIME SPREAD???
	else if (fType == 2 && totPoisson){
		nPhotons = this->ThrowPoisson();
		std::string PMTName = "88mm";
		//std::string PMTName = "R6091";
		peSmeared = fTOTPMT->CalculateCharge(nPhotons, PMTName);
	}
	else if (fType == 2 && !totPoisson){
		std::string PMTName = "88mm";
		//std::string PMTName = "R6091";
		peSmeared = fTOTPMT->CalculateCharge(fMu, PMTName);
	}
	else{
		std::cerr << "Type not recognised -> " << fType << std::endl;
		return;
	}

	if(peSmeared == 0.0){
		fNumZero ++;
	}

	fProbHisto->Fill(fMu, peSmeared);
}

void WCSimDigitizerPDFMaker::FillEmptyBins()
{
	for(int iBinX = 1; iBinX <= fProbHisto->GetNbinsX(); ++iBinX)
	{
		for(int iBinY = 1; iBinY <= fProbHisto->GetNbinsY(); ++iBinY)
		{
			if( fProbHisto->GetBinContent(iBinX, iBinY) == 0 )
			{
				fProbHisto->SetBinContent(iBinX, iBinY, 1.0);
			}
		}
	}
}

void WCSimDigitizerPDFMaker::NormHistogram() {
	fProbHisto->Scale(1.0/fNumThrows);
	return;
}

void WCSimDigitizerPDFMaker::lnHist() {
	for(int iBinX = 1; iBinX <= fProbHisto->GetNbinsX(); ++iBinX)
	{
		for(int iBinY = 1; iBinY <= fProbHisto->GetNbinsY(); ++iBinY)
		{
			double content = fProbHisto->GetBinContent(iBinX, iBinY);
			content = -2 * TMath::Log(content);
			fProbHisto->SetBinContent(iBinX, iBinY, content);
		}
	}
}

void WCSimDigitizerPDFMaker::SaveHistogram() {
	fProbHisto->SetName("digiPDF");
	fProbHisto->SetDirectory(0);
	fDebug->SetDirectory(0);

	TString fileName;
	if(fType == 0){ fileName = TString::Format("sk1pe/digitizerLikelihood_%f.root", fMu); }
	else if (fType == 1){ fileName = TString::Format("pmtSim/digitizerLikelihood_%f.root", fMu); } // TIME SPREAD???
	else if (fType == 2){ fileName = TString::Format("tot/digitizerLikelihood_%f.root", fMu); }
	else{
		std::cerr << "Type not recognised -> " << fType << std::endl;
		return;
	}

	TFile * f = new TFile(fileName.Data(), "RECREATE");
	fProbHisto->Write();
	fDebug->Write();
	f->Close();
	delete f;
}

void WCSimDigitizerPDFMaker::MakeHisto()
{
	double binWidth = 0.5 * ( fChargeMax - fChargeMin ) / fNumChargeBins;
	//std::cout << "binWidth = " << binWidth << std::endl;

	// Make the actual likelihood map histogram
	fProbHisto = new TH2D("digiPDF","digiPDF",
						   fNumChargeBins, fChargeMin-binWidth, fChargeMax-binWidth,
						   fNumChargeBins, fChargeMin-binWidth, fChargeMax-binWidth);
	fProbHisto->GetXaxis()->SetTitle("Predicted mean number of photons, #mu");
	fProbHisto->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbHisto->GetZaxis()->SetTitle("P(digitized | mean photons)");

	// Make the debig histogram
	fDebug = new TH1D("fDebug","Poisson picks",(int)(ceil(fChargeMax) - floor(fChargeMin)), fChargeMin - binWidth, fChargeMax - binWidth);
	fDebug->GetXaxis()->SetTitle("Q");
	fDebug->GetYaxis()->SetTitle("Events");

	return;
}
