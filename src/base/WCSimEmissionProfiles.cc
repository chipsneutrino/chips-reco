/*
 * WCSimEmissionProfiles.cxx
 *
 *  Created on: 13 Nov 2014
 *      Author: andy
 */
#include "TString.h"

#include "WCSimEmissionProfiles.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimTrackParameterEnums.hh"
#include <TVector3.h>
#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCanvas.h>
#include <cassert>
#include <iostream>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimEmissionProfiles)
#endif

WCSimEmissionProfiles::WCSimEmissionProfiles()
{

	fDebug = kFALSE;

	fProfileFileName = TString("");

	fProfileFile = 0x0;

	fBinningHistogram = 0x0;
	fProfileTree = 0x0;

	fRho = 0x0;
	fGCoarse = 0x0;
	fGFine = 0x0;

	fLastPercentile = 0;
	fPercentileTrackLength = 0;

	fType = TrackType::Unknown;
	fEnergy = -999.9;

	fTimeCosThetaMax = 1.0;
	fTimeCosThetaMin = -1.0;
	fSForTime = 0x0;
	fSCosThetaForTime = 0x0;
}

WCSimEmissionProfiles::~WCSimEmissionProfiles()
{
	if (fProfileFile != 0x0)
	{
		// This also deletes fBinningHistogram
		fProfileFile->Close();
		delete fProfileFile;
		fProfileFile = 0x0;
		fProfileTree = 0x0;
		fRho = 0x0;
		fGFine = 0x0;
		fGCoarse = 0x0;
		fSCosThetaForTime = 0x0;
		fSForTime = 0x0;
	}
}

WCSimEmissionProfiles::WCSimEmissionProfiles(const TrackType::Type &type, const double &energy)
{

	fProfileFileName = TString("");

	fProfileFile = 0x0;

	fBinningHistogram = 0x0;
	fProfileTree = 0x0;

	fRho = 0x0;
	fGCoarse = 0x0;
	fGFine = 0x0;

	fSCosThetaForTime = 0x0;
	fSForTime = 0x0;
	fTimeCosThetaMin = -1.0;
	fTimeCosThetaMax = 1.0;

	LoadFile(type, energy);
}

Double_t WCSimEmissionProfiles::GetRho(const Double_t &s)
{
	return fRho->GetBinContent(fRho->GetXaxis()->FindBin(s));
}

Double_t WCSimEmissionProfiles::GetRhoWidth(const Double_t &s)
{
	return fRho->GetBinWidth(fRho->GetXaxis()->FindBin(s));
}

Double_t WCSimEmissionProfiles::GetG(const Double_t &s, const Double_t &cosTheta)
{
	TH2F *hist = fGFine;
	if (cosTheta < fGFine->GetXaxis()->GetXmin())
	{
		hist = fGCoarse;
	}
	return hist->GetBinContent(hist->GetXaxis()->FindBin(cosTheta), hist->GetYaxis()->FindBin(s));
}

Double_t WCSimEmissionProfiles::GetGWidth(const Double_t &cosTheta)
{
	TH2F *hist = fGFine;
	if (cosTheta < fGFine->GetXaxis()->GetXmin())
	{
		hist = fGCoarse;
	}
	return hist->GetXaxis()->GetBinWidth(hist->GetXaxis()->FindBin(cosTheta));
}

Double_t WCSimEmissionProfiles::GetIntegral(WCSimLikelihoodTrackBase *myTrack, EmissionProfile_t::Type type,
											WCSimLikelihoodDigit *myDigit, Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth)
{
	if (type == EmissionProfile_t::kRho)
	{
		return GetRhoIntegral(myTrack, sPower, 0, cutoffS, multiplyByWidth);
	}
	else if (type == EmissionProfile_t::kRhoTimesG)
	{
		return GetRhoGIntegral(myTrack, myDigit, sPower, cutoffS, multiplyByWidth);
	}
	assert(0);
	return 0;
}

Double_t WCSimEmissionProfiles::GetRhoIntegral(WCSimLikelihoodTrackBase *myTrack, Int_t sPower, Double_t startS,
											   Double_t endS, Bool_t multiplyByWidth)
{
	std::vector<Int_t> sVec(1, sPower);
	return GetRhoIntegrals(myTrack, sVec, startS, endS, multiplyByWidth).at(0);
}

std::vector<Double_t> WCSimEmissionProfiles::GetRhoIntegrals(WCSimLikelihoodTrackBase *myTrack,
															 std::vector<Int_t> sPowers, Double_t startS, Double_t endS, Bool_t multiplyByWidth)
{
	// The integrals we'll return
	std::vector<Double_t> integrals(sPowers.size(), 0.0);

	// Need to subtract off the conversion distance (it's zero unless we have a photon track)
	startS -= myTrack->GetConversionDistance();
	endS -= myTrack->GetConversionDistance();

	// Now find the right bins:
	Int_t startBin = fRho->GetXaxis()->FindBin(startS);
	Int_t endBin = fRho->GetXaxis()->FindBin(endS);
	if (endBin < 1)
	{
		return integrals;
	} // We don't get to any emission, so can just return here
	if (startBin < 1)
	{
		startBin = 1;
	} // Are we starting off the end of the histogram? Let's not.

	//	std::cout << fRho->GetBinContent(endBin) << std::endl;;
	for (Int_t iBin = startBin; iBin < endBin; ++iBin)
	{
		Double_t binWidth = multiplyByWidth ? fRho->GetXaxis()->GetBinWidth(iBin) : 1;
		for (UInt_t iPower = 0; iPower < sPowers.size(); ++iPower)
		{
			integrals.at(iPower) += binWidth * fRho->GetBinContent(iBin) * pow(fRho->GetBinCenter(iBin), sPowers.at(iPower));
		}
	}
	return integrals;
}

Double_t WCSimEmissionProfiles::GetRhoGIntegral(WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit,
												Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth)
{
	std::vector<Int_t> sVec(1, sPower);
	return GetRhoGIntegrals(myTrack, myDigit, sVec, cutoffS, multiplyByWidth).at(0);
}

std::vector<Double_t> WCSimEmissionProfiles::GetRhoGIntegrals(WCSimLikelihoodTrackBase *myTrack,
															  WCSimLikelihoodDigit *myDigit, std::vector<Int_t> sPowers, Double_t cutoffS, Bool_t multiplyByWidth)
{

	// std::cout << "WCSimEmissionProfiles::GetRhoGIntegrals()" << std::endl;
	std::vector<Double_t> integrals(3, 0.0);

	// If the track has a conversion distance we need to subtract it off the cutoff to
	// work out how far through the emission profile we need to integrate
	Int_t startBin = fRho->GetXaxis()->FindBin(0.);
	double emissionCutoff = myTrack->GetConversionDistance();
	if (emissionCutoff <= 0.0)
	{
		return integrals;
	}

	Int_t endBin = fRho->GetXaxis()->FindBin(emissionCutoff);

	TVector3 pmtPos = myDigit->GetPos();
	TVector3 vtxPos = myTrack->GetVtx();
	TVector3 vtxDir = myTrack->GetDir();

	// std::cout << "Integrating from " << 0.0 << " to " << cutoffS << std::endl;

	for (Int_t iBin = startBin; iBin < endBin; ++iBin)
	{
		Double_t s = fRho->GetBinCenter(iBin);
		if (s < myTrack->GetConversionDistance())
		{
			continue;
		}

		Double_t emissionDistance = s - myTrack->GetConversionDistance();

		TVector3 toPMT = pmtPos - myTrack->GetPropagatedPos(s);
		Double_t cosTheta = TMath::Cos(vtxDir.Angle(toPMT));

		Double_t binWidthS = multiplyByWidth ? GetRhoWidth(emissionDistance) : 1;

		//    if( myDigit->GetTubeId() == 3364)
		//    {
		//      std::cout << "s = " << s << "   toPMT = "  << toPMT.Mag() << "   costheta = " << cosTheta << std::endl;
		//    }
		double adding = 0;
		for (UInt_t iPower = 0; iPower < sPowers.size(); ++iPower)
		{
			// std::cout << "Binwidth   = " << binWidth << std::endl
			//           << "fRho = " << fRho->GetBinContent(iBin) << std::endl
			//           << "fG         = " << fG->GetBinContent(fG->GetXaxis()->FindBin(cosTheta), iBin) << std::endl;

			if (fRho->GetBinContent(iBin) == 0)
			{
				continue;
			}
			adding = binWidthS * fRho->GetBinContent(iBin) * GetG(emissionDistance, cosTheta);
			integrals.at(iPower) += adding * pow(emissionDistance, sPowers.at(iPower));
			// std::cout << "Now integrals.at(" << iPower << ") = " << integrals.at(iPower) <<  " sPower = " << sPowers.at(iPower) << std::endl;
		}
		//
		//    if( myDigit->GetTubeId() == 3364)
		//    {
		//      std::cout << "adding " << adding << "   now integrals[0] = " << integrals.at(0) << "  integrals[1] = " << integrals.at(1) << "  integrals[2] = " << integrals.at(2) << std::endl;
		//    }
		//
	}
	return integrals;
}

UInt_t WCSimEmissionProfiles::GetTreeEntry(Double_t energy) const
{
	// std::cout << "Binning histogram is " << fBinningHistogram << std::endl;
	int iBin = fBinningHistogram->GetXaxis()->FindBin(energy) - 1;
	if (iBin < 0)
	{
		std::cerr << "Warning, array bin number is less than zero - returning 0 instead - bin was " << iBin
				  << std::endl;
		iBin = 0;
	}

	return iBin;
}

UInt_t WCSimEmissionProfiles::GetHistBin(Double_t energy) const
{
	return GetTreeEntry(energy) + 1;
}

Double_t WCSimEmissionProfiles::GetFractionThroughBin(Double_t energy) const
{
	UInt_t iBin = GetHistBin(energy);
	// std::cout << "Hist bin = " << iBin << std::endl;
	// std::cout << "Its low edge = " << fBinningHistogram->GetXaxis()->GetBinLowEdge(iBin) << std::endl;
	// std::cout << "And its width = " << fBinningHistogram->GetXaxis()->GetBinWidth(iBin);
	Double_t fraction = (energy - fBinningHistogram->GetXaxis()->GetBinLowEdge(iBin)) / (fBinningHistogram->GetXaxis()->GetBinWidth(iBin));
	return fraction;
}

void WCSimEmissionProfiles::LoadFile(const TrackType::Type &type, const double &energy)
{
	std::cout << " *** WCSimEmissionProfiles::LoadFile - Loading profile" << std::endl;
	fProfileFileName = TString(getenv("CHIPSRECO"));
	switch (type)
	{
	case TrackType::PhotonLike:
		// Fall through to the electron track
	case TrackType::ElectronLike:
		fProfileFileName.Append("/config/emissionProfilesElectrons.root");
		break;
	case TrackType::MuonLike:
		fProfileFileName.Append("/config/emissionProfilesMuons.root");
		break;
	default:
		std::cerr << "Error: unknown track type in WCSimLikelihoodTuner::LoadEmissionProfiles" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Open the files and ge the profile tree
	std::cout << "Profile file name = " << fProfileFileName << std::endl;
	fProfileFile = new TFile(fProfileFileName.Data(), "READ");
	fProfileFile->GetObject("fProfileTree", fProfileTree);
	fProfileFile->GetObject("hWhichHisto", fBinningHistogram);

	// Set the tree branches
	std::cout << "Set branch for rho" << std::endl;
	fProfileTree->SetBranchAddress("hRho", &fRho);
	std::cout << "Set branch for GFine" << std::endl;
	fProfileTree->SetBranchAddress("hGFine", &fGFine);
	std::cout << "Set branch for GCoarse" << std::endl;
	fProfileTree->SetBranchAddress("hGCoarse", &fGCoarse);
	std::cout << "Set branch for SCosThetaForTime" << std::endl;
	fProfileTree->SetBranchAddress("hSCosThetaForTime", &fSCosThetaForTime);
	std::cout << "Set branch for SForTime" << std::endl;
	fProfileTree->SetBranchAddress("hSForTime", &fSForTime);
	std::cout << "Set branch for timeCosThetaMin" << std::endl;
	fProfileTree->SetBranchAddress("timeCosThetaMin", &fTimeCosThetaMin);
	std::cout << "Set branch for timeCosThetaMax" << std::endl;
	fProfileTree->SetBranchAddress("timeCosThetaMax", &fTimeCosThetaMax);

	fType = type;
	SetEnergy(energy); // This calls GetEntry

	if (fDebug)
	{
		SaveProfiles(type, energy);
	}

	fStoppingDistance = -999.9;
	return;
}

void WCSimEmissionProfiles::SetEnergy(const double &energy)
{
	if (energy != fEnergy)
	{
		fEnergy = energy;
		int bin = fBinningHistogram->GetXaxis()->FindBin(energy) - 1;
		assert(bin <= fProfileTree->GetEntries());
		fProfileTree->GetEntry(bin);
		fStoppingDistance = -999.9;
		fLastPercentile = 0;
		fPercentileTrackLength = 0;
	}

	return;
}

Double_t WCSimEmissionProfiles::GetCriticalDistance(TrackType::Type type, Double_t energy) const
{
	assert(type == TrackType::MuonLike || type == TrackType::ElectronLike);
	switch (type)
	{
	case TrackType::MuonLike:
		return -398 + 0.756 * energy - 8e-5 * energy * energy;
		break;
	default:
		assert(0);
	}
	return 0;
}

std::vector<Double_t> WCSimEmissionProfiles::GetProfileEnergies() const
{
	// std::cout << "fBinningHistogram = " << fBinningHistogram << std::endl;
	// std::cout << "bins = "; std::cout << fBinningHistogram->GetNbinsX() << std::endl;
	std::vector<Double_t> energies;

	if (fBinningHistogram == 0x0)
	{
		std::cerr << "Error: binning histogram not yet loaded" << std::endl;
		assert(fBinningHistogram != 0x0);
	}

	int nBins = fBinningHistogram->GetNbinsX();
	for (int i = 1; i <= nBins; ++i)
	{
		energies.push_back(fBinningHistogram->GetXaxis()->GetBinLowEdge(i));
	}
	return energies;
}

Double_t WCSimEmissionProfiles::GetLightFlux(const TrackType::Type &type, const double &energy) const
{
	assert(type == TrackType::MuonLike || type == TrackType::ElectronLike || type == TrackType::PhotonLike);
	Double_t flux = 0.0;

	// There was an overall normalization of 2pi that used to be neglected
	// because the solid angle subtended by a cone is 2pi(1 - cosTheta)
	// I've absorbed this into the SolidAngle now

	Double_t nPhotons = 0.0;

	switch (type)
	{
	case TrackType::MuonLike:
		// From a linear fit to nPhotons as a function of energy - it's fairly
		// linear for muons and very linear for electrons

		// Fitted on 12th May 2016
		nPhotons = -30849 + 350.289 * energy;
		if (nPhotons > 0)
		{
			flux = nPhotons;
		}
		break;
	case TrackType::PhotonLike:
		// Fall through to ElectronLike
	case TrackType::ElectronLike:
		nPhotons = 198.152 + 383.238 * energy;
		if (nPhotons > 0)
		{
			flux = nPhotons;
		}
		break;
	default:
		assert(0);
		break;
	}
	// std::cout << "numPhotons = " << nPhotons << std::endl;
	return flux;
}

Double_t WCSimEmissionProfiles::GetTrackLengthForPercentile(Double_t percentile)
{
	if (fLastPercentile == percentile)
	{
		return fPercentileTrackLength;
	}
	Double_t runningTotal = 0.0;
	Int_t iBin = 0;
	while (iBin < fRho->GetNbinsX() && runningTotal < percentile)
	{
		iBin++;
		runningTotal += fRho->GetBinContent(iBin) * fRho->GetXaxis()->GetBinWidth(iBin);
	}
	fLastPercentile = percentile;
	fPercentileTrackLength = fRho->GetXaxis()->GetBinCenter(iBin);
	return fPercentileTrackLength;
}

void WCSimEmissionProfiles::SaveProfiles(const TrackType::Type &type, const double &energy)
{
	return;
	TCanvas *canRho = new TCanvas("canRho", "canRho", 800, 600);
	fRho->Draw();
	fRho->SetLineWidth(2);
	fRho->SetLineColor(kBlue);
	fRho->GetXaxis()->SetTitle("s/cm");
	fRho->GetYaxis()->SetTitle("rho(s)");
	fRho->GetXaxis()->CenterTitle();
	fRho->GetYaxis()->CenterTitle();
	TString rhoTitle("");
	rhoTitle.Form("rho_%d_", (int)energy);
	rhoTitle += TrackType::AsString(type);
	fRho->SetTitle(rhoTitle.Data());

	canRho->SaveAs((rhoTitle + TString(".png")).Data());
	canRho->SaveAs((rhoTitle + TString(".C")).Data());

	TCanvas *canG = new TCanvas("canG", "canG", 800, 600);
	fGCoarse->Draw("COLZ");
	fGCoarse->GetXaxis()->SetTitle("cos#theta");
	fGCoarse->GetYaxis()->SetTitle("s (cm)");
	fGCoarse->GetXaxis()->CenterTitle();
	fGCoarse->GetYaxis()->CenterTitle();
	TString gTitle("");
	gTitle.Form("gCoarse_%d_", (int)energy);
	gTitle += TrackType::AsString(type);
	fGCoarse->SetTitle(gTitle.Data());

	canG->SaveAs((gTitle + TString(".png")).Data());
	canG->SaveAs((gTitle + TString(".C")).Data());

	fGFine->Draw("COLZ");
	fGFine->GetXaxis()->SetTitle("cos#theta");
	fGFine->GetYaxis()->SetTitle("s (cm)");
	fGFine->GetXaxis()->CenterTitle();
	fGFine->GetYaxis()->CenterTitle();
	gTitle.Form("gFine_%d_", (int)energy);
	gTitle += TrackType::AsString(type);
	fGFine->SetTitle(gTitle.Data());

	canG->SaveAs((gTitle + TString(".png")).Data());
	canG->SaveAs((gTitle + TString(".C")).Data());
	delete canRho;
	delete canG;

	return;
}

Double_t WCSimEmissionProfiles::GetStoppingDistance()
{
	if (fStoppingDistance == -999.9)
	{
		Int_t bin = fRho->FindLastBinAbove(0);
		fStoppingDistance = fRho->GetBinCenter(bin);
	}
	return fStoppingDistance;
}

TH1F *WCSimEmissionProfiles::GetRho()
{
	return fRho;
}

TH2F *WCSimEmissionProfiles::GetGCoarse()
{
	return GetG().first;
}

TH2F *WCSimEmissionProfiles::GetGFine()
{
	return GetG().second;
}

std::pair<TH2F *, TH2F *> WCSimEmissionProfiles::GetG()
{
	return std::make_pair(fGCoarse, fGFine);
}

TH1F *WCSimEmissionProfiles::GetEnergyHist()
{
	return fBinningHistogram;
}

TH2F *WCSimEmissionProfiles::GetSCosThetaForTime()
{
	return fSCosThetaForTime;
}

TH1F *WCSimEmissionProfiles::GetSForTime()
{
	return fSForTime;
}

double WCSimEmissionProfiles::GetTimeCosThetaMin()
{
	return fTimeCosThetaMin;
}

double WCSimEmissionProfiles::GetTimeCosThetaMax()
{
	return fTimeCosThetaMax;
}
