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

WCSimDigitizerPDFMaker::WCSimDigitizerPDFMaker()
{
	// Set up the random number generator
	fRandom = new TRandom3();
	fRandom->SetSeed(0);

	fMu = 0;
	fNumThrows = 1000000;
	fType = 0; // default is sk1pe for now

	fProbRawHisto = 0x0;
	fProbPoissonHisto = 0x0;
	fPoisson = 0x0;

	fProbRawHisto_digiNorm = 0x0;
	fProbPoissonHisto_digiNorm = 0x0;
	fProbRawHisto_ln = 0x0;
	fProbPoissonHisto_ln = 0x0;
	fProbRawHisto_digiNorm_ln = 0x0;
	fProbPoissonHisto_digiNorm_ln = 0x0;

	fSK1pePMT = new WCSimSK1pePMT();
	fCHIPSPMT = new WCSimCHIPSPMT();
	fTOTPMT = new WCSimTOTPMT();
}

WCSimDigitizerPDFMaker::~WCSimDigitizerPDFMaker()
{
	if (fProbRawHisto != 0x0) delete fProbRawHisto;
	if (fProbPoissonHisto != 0x0) delete fProbPoissonHisto;
	if (fPoisson != 0x0) delete fPoisson;

	if (fProbRawHisto_digiNorm != 0x0) delete fProbRawHisto_digiNorm;
	if (fProbPoissonHisto_digiNorm != 0x0) delete fProbPoissonHisto_digiNorm;
	if (fProbRawHisto_ln != 0x0) delete fProbRawHisto_ln;
	if (fProbPoissonHisto_ln != 0x0) delete fProbPoissonHisto_ln;
	if (fProbRawHisto_digiNorm_ln != 0x0) delete fProbRawHisto_digiNorm_ln;
	if (fProbPoissonHisto_digiNorm_ln != 0x0) delete fProbPoissonHisto_digiNorm_ln;

	delete fRandom;
}

void WCSimDigitizerPDFMaker::SetBinning(int bins, int min, int max)
{
	fNumChargeBins = bins;
	fChargeMin = min;
	fChargeMax = max;
}

int WCSimDigitizerPDFMaker::GetBins() const
{
	return fNumChargeBins;
}

int WCSimDigitizerPDFMaker::GetMin() const
{
	return fChargeMin;
}

int WCSimDigitizerPDFMaker::GetMax() const
{
	return fChargeMax;
}

void WCSimDigitizerPDFMaker::SetNumThrows(int num)
{
	fNumThrows = num;
}

int WCSimDigitizerPDFMaker::GetNumThrows() const
{
	return fNumThrows;
}

void WCSimDigitizerPDFMaker::SetType(int type)
{
	fType = type;
}

int WCSimDigitizerPDFMaker::GetType() const
{
	return fType;
}

void WCSimDigitizerPDFMaker::Run()
{
	if (fType != 0 && fType != 1 && fType != 2)
	{
		std::cerr << "Type not recognised -> " << fType << std::endl;
		return;
	}
	MakeHisto();
	LoopDigitize();
	FillEmptyBins();
	NormHists();
	CreateDigiNormHists();
	CreateLnHists();
	SaveHistogram();
}

int WCSimDigitizerPDFMaker::ThrowPoisson()
{
	int myRandom = fRandom->Poisson(fMu);
	fPoisson->Fill(myRandom);
	return myRandom;
}

void WCSimDigitizerPDFMaker::LoopDigitize()
{
	for (int bin = 0; bin < fNumChargeBins; bin++)
	{
		fMu = (double)fChargeMin + (double)(fChargeMax - fChargeMin) / (double)(fNumChargeBins)*(double)bin;
		if (bin % 50 == 0)
			std::cout << "Bin: " << bin << " = " << fMu << std::endl;
		for (int iThrow = 0; iThrow < fNumThrows; ++iThrow)
			Digitize();
	}
}

void WCSimDigitizerPDFMaker::Digitize()
{
	// Throw from a poisson with the mean at fMu
	int poissonPhotons = this->ThrowPoisson();

	double rawSmeared = 0.0;
	double poissonSmeared = 0.0;

	if (fType == 0)
	{
		rawSmeared = fSK1pePMT->CalculateCharge(fMu);
		poissonSmeared = fSK1pePMT->CalculateCharge(poissonPhotons);
	}
	else if (fType == 1)
	{
		rawSmeared = fCHIPSPMT->CalculateCharge(fMu, 0, 0);
		poissonSmeared = fCHIPSPMT->CalculateCharge(poissonPhotons, 0, 0);
	}
	else if (fType == 2)
	{
		std::string TOTPMT = "88mm"; // "R6091"
		rawSmeared = fTOTPMT->CalculateCharge(fMu, TOTPMT);
		poissonSmeared = fTOTPMT->CalculateCharge(poissonPhotons, TOTPMT);
	}

	fProbRawHisto->Fill(fMu, rawSmeared);
	fProbPoissonHisto->Fill(fMu, poissonSmeared);
}

void WCSimDigitizerPDFMaker::FillEmptyBins()
{
	for (int iBinX = 1; iBinX <= fProbRawHisto->GetNbinsX(); ++iBinX)
	{
		for (int iBinY = 1; iBinY <= fProbRawHisto->GetNbinsY(); ++iBinY)
		{
			if (fProbRawHisto->GetBinContent(iBinX, iBinY) == 0)
			{
				fProbRawHisto->SetBinContent(iBinX, iBinY, 1.0);
			}

			if (fProbPoissonHisto->GetBinContent(iBinX, iBinY) == 0)
			{
				fProbPoissonHisto->SetBinContent(iBinX, iBinY, 1.0);
			}
		}
	}
}

void WCSimDigitizerPDFMaker::NormHists()
{
	fProbRawHisto->Scale(1.0 / fNumThrows);
	fProbPoissonHisto->Scale(1.0 / fNumThrows);
	return;
}

void WCSimDigitizerPDFMaker::CreateDigiNormHists()
{
	for (int iBinY = 1; iBinY <= fProbRawHisto->GetNbinsY(); ++iBinY)
	{
		double raw_totInRow = 0.0;
		double pois_totInRow = 0.0;
		for (int iBinX = 1; iBinX <= fProbRawHisto->GetNbinsX(); ++iBinX)
		{
			raw_totInRow += fProbRawHisto->GetBinContent(iBinX, iBinY);
			pois_totInRow += fProbPoissonHisto->GetBinContent(iBinX, iBinY);
		}

		for (int iBinX = 1; iBinX <= fProbRawHisto->GetNbinsX(); ++iBinX)
		{
			double raw_tempBinContent = fProbRawHisto->GetBinContent(iBinX, iBinY);
			fProbRawHisto_digiNorm->SetBinContent(iBinX, iBinY, (raw_tempBinContent / raw_totInRow));

			double pois_tempBinContent = fProbPoissonHisto->GetBinContent(iBinX, iBinY);
			fProbPoissonHisto_digiNorm->SetBinContent(iBinX, iBinY, (pois_tempBinContent / pois_totInRow));
		}
	}
}

void WCSimDigitizerPDFMaker::CreateLnHists()
{
	for (int iBinX = 1; iBinX <= fProbRawHisto->GetNbinsX(); ++iBinX)
	{
		for (int iBinY = 1; iBinY <= fProbRawHisto->GetNbinsY(); ++iBinY)
		{
			double raw_content = fProbRawHisto->GetBinContent(iBinX, iBinY);
			raw_content = -2 * TMath::Log(raw_content);
			fProbRawHisto_ln->SetBinContent(iBinX, iBinY, raw_content);

			double pois_content = fProbPoissonHisto->GetBinContent(iBinX, iBinY);
			pois_content = -2 * TMath::Log(pois_content);
			fProbPoissonHisto_ln->SetBinContent(iBinX, iBinY, pois_content);

			double raw_content_digi_norm = fProbRawHisto_digiNorm->GetBinContent(iBinX, iBinY);
			raw_content_digi_norm = -2 * TMath::Log(raw_content_digi_norm);
			fProbRawHisto_digiNorm_ln->SetBinContent(iBinX, iBinY, raw_content_digi_norm);

			double pois_content_digi_norm = fProbPoissonHisto_digiNorm->GetBinContent(iBinX, iBinY);
			pois_content_digi_norm = -2 * TMath::Log(pois_content_digi_norm);
			fProbPoissonHisto_digiNorm_ln->SetBinContent(iBinX, iBinY, pois_content_digi_norm);
		}
	}
}

void WCSimDigitizerPDFMaker::SaveHistogram()
{
	fProbRawHisto->SetDirectory(0);
	fProbPoissonHisto->SetDirectory(0);
	fPoisson->SetDirectory(0);

	fProbRawHisto_digiNorm->SetDirectory(0);
	fProbPoissonHisto_digiNorm->SetDirectory(0);
	fProbRawHisto_ln->SetDirectory(0);
	fProbPoissonHisto_ln->SetDirectory(0);
	fProbRawHisto_digiNorm_ln->SetDirectory(0);
	fProbPoissonHisto_digiNorm_ln->SetDirectory(0);

	TString fileName;
	if (fType == 0)
	{
		fileName = "sk1pe_digitizerLikelihood.root";
	}
	else if (fType == 1)
	{
		fileName = "pmtSim_digitizerLikelihood.root";
	}
	else if (fType == 2)
	{
		fileName = "tot_digitizerLikelihood.root";
	}

	TFile *f = new TFile(fileName.Data(), "RECREATE");

	fProbRawHisto->Write();
	fProbPoissonHisto->Write();
	fPoisson->Write();

	fProbRawHisto_digiNorm->Write();
	fProbPoissonHisto_digiNorm->Write();
	fProbRawHisto_ln->Write();
	fProbPoissonHisto_ln->Write();
	fProbRawHisto_digiNorm_ln->Write();
	fProbPoissonHisto_digiNorm_ln->Write();

	f->Close();
	delete f;
}

void WCSimDigitizerPDFMaker::MakeHisto()
{
	double binWidth = 0.5 * (fChargeMax - fChargeMin) / fNumChargeBins;

	// Make the raw likelihood map histogram
	fProbRawHisto = new TH2D("rawDigiPDF", "rawDigiPDF", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
							 fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbRawHisto->GetXaxis()->SetTitle("Raw Number of Photons");
	fProbRawHisto->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbRawHisto->GetZaxis()->SetTitle("P(digitized | num photons)");

	// Make the poisson likelihood map histogram
	fProbPoissonHisto = new TH2D("poissonDigiPDF", "poissonDigiPDF", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
								 fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbPoissonHisto->GetXaxis()->SetTitle("Predicted mean number of photons, #mu");
	fProbPoissonHisto->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbPoissonHisto->GetZaxis()->SetTitle("P(digitized | mean photons)");

	// Make the poisson debug histogram
	fPoisson = new TH1D("Poisson", "Poisson picks", (int)(ceil(fChargeMax) - floor(fChargeMin)), fChargeMin - binWidth,
						fChargeMax - binWidth);
	fPoisson->GetXaxis()->SetTitle("Q");
	fPoisson->GetYaxis()->SetTitle("Events");

	// Make the raw likelihood map histogram normalised horizontally
	fProbRawHisto_digiNorm = new TH2D("rawDigiPDF_digiNorm", "rawDigiPDF_digiNorm", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
							 		  fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbRawHisto_digiNorm->GetXaxis()->SetTitle("Raw Number of Photons");
	fProbRawHisto_digiNorm->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbRawHisto_digiNorm->GetZaxis()->SetTitle("P(digitized | num photons)");

	// Make the poisson likelihood map histogram normalised horizontally
	fProbPoissonHisto_digiNorm = new TH2D("poissonDigiPDF_digiNorm", "poissonDigiPDF_digiNorm", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
								 		  fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbPoissonHisto_digiNorm->GetXaxis()->SetTitle("Predicted mean number of photons, #mu");
	fProbPoissonHisto_digiNorm->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbPoissonHisto_digiNorm->GetZaxis()->SetTitle("P(digitized | mean photons)");

	// Make the raw likelihood map histogram logged
	fProbRawHisto_ln = new TH2D("rawDigiPDF_ln", "rawDigiPDF_ln", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
							    fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbRawHisto_ln->GetXaxis()->SetTitle("Raw Number of Photons");
	fProbRawHisto_ln->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbRawHisto_ln->GetZaxis()->SetTitle("-2ln(P(digitized | num photons))");

	// Make the poisson likelihood map histogram logged
	fProbPoissonHisto_ln = new TH2D("poissonDigiPDF_ln", "poissonDigiPDF_ln", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
								    fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbPoissonHisto_ln->GetXaxis()->SetTitle("Predicted mean number of photons, #mu");
	fProbPoissonHisto_ln->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbPoissonHisto_ln->GetZaxis()->SetTitle("-2ln(P(digitized | mean photons))");

	// Make the raw likelihood map histogram normalised horizontally logged
	fProbRawHisto_digiNorm_ln = new TH2D("rawDigiPDF_digiNorm_ln", "rawDigiPDF_digiNorm_ln", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
							 		     fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbRawHisto_digiNorm_ln->GetXaxis()->SetTitle("Raw Number of Photons");
	fProbRawHisto_digiNorm_ln->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbRawHisto_digiNorm_ln->GetZaxis()->SetTitle("-2ln(P(digitized | num photons))");

	// Make the poisson likelihood map histogram normalised horizontally logged
	fProbPoissonHisto_digiNorm_ln = new TH2D("poissonDigiPDF_digiNorm_ln", "poissonDigiPDF_digiNorm_ln", fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth,
								 		     fNumChargeBins, fChargeMin - binWidth, fChargeMax - binWidth);
	fProbPoissonHisto_digiNorm_ln->GetXaxis()->SetTitle("Predicted mean number of photons, #mu");
	fProbPoissonHisto_digiNorm_ln->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbPoissonHisto_digiNorm_ln->GetZaxis()->SetTitle("-2ln(P(digitized | mean photons))");

	return;
}
