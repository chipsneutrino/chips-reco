/*
 * WCSimEmissionProfiles.cxx
 *
 *  Created on: 13 Nov 2014
 *      Author: andy
 */
#include "TString.h"

#include "WCSimEmissionProfiles.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodDigit.hh"
#include <TVector3.h>
#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <cassert>
#include <iostream>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimEmissionProfiles)
#endif

WCSimEmissionProfiles::WCSimEmissionProfiles() {

  fDebug = kFALSE;
  
	fLastEnergy = -999.9;
	fLastType = WCSimLikelihoodTrack::Unknown;

	fProfileFileName = TString("");

	fProfileFile = 0x0;

	fRhoArray = 0x0;
	fGCoarseArray = 0x0;
	fGFineArray = 0x0;
	fBinningHistogram = 0x0;

	fRhoInterp = 0x0;
	fRhoInterpLo = 0x0;
	fRhoInterpHi = 0x0;
	fGCoarse = 0x0;
	fGFine = 0x0;
}

WCSimEmissionProfiles::~WCSimEmissionProfiles() {
	if(fProfileFile != 0x0)
	{
		// This also deletes fBinningHisto
		fProfileFile->Close();
		delete fProfileFile;
		fProfileFile = 0x0;
	}
	if(fRhoArray != 0x0)
	{
		delete fRhoArray;
		fRhoArray = 0x0;
	}
	if(fGCoarseArray != 0x0)
	{
		delete fGCoarseArray;
		fGCoarseArray = 0x0;
	}
	if(fGFineArray != 0x0)
	{
		delete fGFineArray;
		fGFineArray = 0x0;
	}

	if( fRhoInterp != 0x0 )
	{
		delete fRhoInterp;
		fRhoInterp = 0x0;
	}

	if(	fGCoarse != 0x0 )
	{
		delete fGCoarse;
		fGCoarse = 0x0;
	}

	if(	fGFine != 0x0 )
	{
		delete fGFine;
		fGFine = 0x0;
	}

	if(fRhoInterpLo != 0x0)
	{
		delete fRhoInterpLo;
		fRhoInterpLo = 0x0;
	}

	if(fRhoInterpHi != 0x0)
	{
		delete fRhoInterpHi;
		fRhoInterpHi = 0x0;
	}

}

WCSimEmissionProfiles::WCSimEmissionProfiles(WCSimLikelihoodTrack* myTrack) {

	fProfileFileName = TString("");

	fProfileFile = 0x0;

	fRhoArray = 0x0;
	fGCoarseArray = 0x0;
	fGFineArray = 0x0;
	fBinningHistogram = 0x0;

	fRhoInterp = 0x0;
	fGCoarse = 0x0;
	fGFine= 0x0;

	fRhoInterp = 0x0;
	fRhoInterpLo = 0x0;
	fRhoInterpHi = 0x0;

	SetTrack(myTrack);

}

void WCSimEmissionProfiles::SetTrack(WCSimLikelihoodTrack* myTrack) {
	if( myTrack->GetE() != fLastEnergy || myTrack->GetType() != fLastType)
	{
		LoadFile(myTrack);
	}
	return;
}

Double_t WCSimEmissionProfiles::GetRho(Double_t s) {
	return fRhoInterp->GetBinContent(fRhoInterp->GetXaxis()->FindBin(s));
}

Double_t WCSimEmissionProfiles::GetG(Double_t s, Double_t cosTheta) {
	TH2F * hist = fGFine;
	if(cosTheta < fGFine->GetXaxis()->GetXmin())
	{
		hist = fGCoarse;
	}
	return hist->GetBinContent(hist->GetXaxis()->FindBin(cosTheta), hist->GetYaxis()->FindBin(s));
}


Double_t WCSimEmissionProfiles::GetIntegral(EmissionProfile_t::Type type, WCSimLikelihoodTrack* myTrack,
		WCSimLikelihoodDigit * myDigit, Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth) {
	if( myTrack->GetType() != fLastType || myTrack->GetE() != fLastEnergy )
	{
		SetTrack(myTrack);
	}
	if( type == EmissionProfile_t::kRho )
	{
		return GetRhoIntegral(sPower, 0, cutoffS, multiplyByWidth);
	}
	else if( type == EmissionProfile_t::kRhoTimesG )
	{
		return GetRhoGIntegral(myTrack, myDigit, sPower, cutoffS, multiplyByWidth);
	}
	assert(0);
	return 0;
}

Double_t WCSimEmissionProfiles::GetRhoIntegral(Int_t sPower, Double_t energy,
		Double_t startS, Double_t endS, Bool_t multiplyByWidth) {
  std::vector<Int_t> sVec(1, sPower);
  return GetRhoIntegrals(sVec, energy, startS, endS, multiplyByWidth).at(0);
}

std::vector<Double_t> WCSimEmissionProfiles::GetRhoIntegrals(std::vector<Int_t> sPowers, Double_t energy,
		Double_t startS, Double_t endS, Bool_t multiplyByWidth) {
	Int_t startBin = fRhoInterp->GetXaxis()->FindBin(startS);
	Int_t endBin = fRhoInterp->GetXaxis()->FindBin(endS);
//	std::cout << fRhoInterp->GetBinContent(endBin) << std::endl;;
	std::vector<Double_t> integrals(sPowers.size(), 0.0);
	for(Int_t iBin = startBin; iBin < endBin; ++iBin)
	{
		Double_t binWidth = multiplyByWidth ? fRhoInterp->GetXaxis()->GetBinWidth(iBin) : 1;
    for(UInt_t iPower = 0; iPower < sPowers.size(); ++iPower)
    {
     	integrals.at(iPower) += binWidth * fRhoInterp->GetBinContent(iBin) * pow(fRhoInterp->GetBinCenter(iBin), sPowers.at(iPower));
	  }
  }
	return integrals;
}


Double_t WCSimEmissionProfiles::GetRhoGIntegral(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit,
		Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth) {
  std::vector<Int_t> sVec(1, sPower);
  return GetRhoGIntegrals(myTrack, myDigit, sVec, cutoffS, multiplyByWidth).at(0);
}

std::vector<Double_t> WCSimEmissionProfiles::GetRhoGIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit,
		std::vector<Int_t> sPowers, Double_t cutoffS, Bool_t multiplyByWidth) {

	// std::cout << "WCSimEmissionProfiles::GetRhoGIntegrals()" << std::endl;
	std::vector<Double_t> integrals(3, 0.0);
   
	Int_t startBin = fRhoInterp->GetXaxis()->FindBin(0.);
	Int_t endBin = fRhoInterp->GetXaxis()->FindBin(cutoffS);

    TVector3 pmtPos = myDigit->GetPos();
    TVector3 vtxPos = myTrack->GetVtx();
    TVector3 vtxDir = myTrack->GetDir();

   // std::cout << "Integrating from " << 0.0 << " to " << cutoffS << std::endl;

	for(Int_t iBin = startBin; iBin < endBin; ++iBin)
	{
		Double_t s = fRhoInterp->GetBinCenter(iBin);
		TVector3 toPMT = pmtPos - myTrack->GetPropagatedPos(s);
		Double_t cosTheta = TMath::Cos(vtxDir.Angle( toPMT ));

		Double_t binWidth = multiplyByWidth ? fRhoInterp->GetXaxis()->GetBinWidth(iBin) : 1;
    
		for(UInt_t iPower = 0 ; iPower < sPowers.size() ; ++iPower)
		{
			// std::cout << "Binwidth   = " << binWidth << std::endl
			//           << "fRhoInterp = " << fRhoInterp->GetBinContent(iBin) << std::endl
			//           << "fG         = " << fG->GetBinContent(fG->GetXaxis()->FindBin(cosTheta), iBin) << std::endl;

			if(fRhoInterp->GetBinContent(iBin) == 0) { continue; }
			integrals.at(iPower) += binWidth
									* fRhoInterp->GetBinContent(iBin)
									* GetG(s, cosTheta)
									* pow(fRhoInterp->GetBinCenter(iBin), sPowers.at(iPower));
			// std::cout << "Now integrals.at(" << iPower << ") = " << integrals.at(iPower) <<  " sPower = " << sPowers.at(iPower) << std::endl;
		}
	}
	return integrals;

}

UInt_t WCSimEmissionProfiles::GetArrayBin(Double_t energy) const{
	// std::cout << "Binning histogram is " << fBinningHistogram << std::endl;
	int iBin = fBinningHistogram->GetXaxis()->FindBin(energy) - 1;
	if(iBin < 0)
	{
		std::cerr << "Warning, array bin number is less than zero - returning 0 instead - bin was " << iBin << std::endl;
		iBin = 0;
	}
  
	return iBin;
}

UInt_t WCSimEmissionProfiles::GetHistBin(Double_t energy) const{
	return GetArrayBin(energy)+1;
}

Double_t WCSimEmissionProfiles::GetFractionThroughBin(Double_t energy) const{
	UInt_t iBin = GetHistBin(energy);
	// std::cout << "Hist bin = " << iBin << std::endl;
	// std::cout << "Its low edge = " << fBinningHistogram->GetXaxis()->GetBinLowEdge(iBin) << std::endl;
	// std::cout << "And its width = " << fBinningHistogram->GetXaxis()->GetBinWidth(iBin);
	Double_t fraction =   (energy - fBinningHistogram->GetXaxis()->GetBinLowEdge(iBin))
						/ (fBinningHistogram->GetXaxis()->GetBinWidth(iBin));
	return fraction;
}

void WCSimEmissionProfiles::LoadFile(WCSimLikelihoodTrack* myTrack) {
	std::cout << " *** WCSimEmissionProfiles::LoadFile - Loading profile" << std::endl;
	// std::cout << "Track type is " << WCSimLikelihoodTrack::TrackTypeToString(myTrack->GetType()) << std::endl;
	if( myTrack->GetType() != fLastType )
	{

		fProfileFileName = TString(getenv("WCSIMANAHOME"));
		fProfileFileName = "";
		switch( myTrack->GetType() )
		{
			case WCSimLikelihoodTrack::ElectronLike:
				fProfileFileName.Append("/home/ajperch/CHIPS/WCSim_fixGeneration/rootfiles/sum/emissionProfilesElectronsCoarseFine.root");
				break;
			case WCSimLikelihoodTrack::MuonLike:
				fProfileFileName.Append("/home/ajperch/CHIPS/WCSim_fixGeneration/rootfiles/sum/emissionProfilesMuonsCoarseFine.root");
				break;
			default:
        myTrack->Print();
				std::cerr << "Error: unknown track type in WCSimLikelihoodTuner::LoadEmissionProfiles" << std::endl;
				exit(EXIT_FAILURE);
		}


		if( fProfileFile != NULL )
		{
			delete fProfileFile;
			delete fRhoArray;
			delete fGCoarseArray;
			delete fGFineArray;
		}

		fProfileFile      = new TFile(fProfileFileName.Data(),"READ");
		fProfileFile->ls();
		fLastType         = myTrack->GetType();
		fRhoArray         = (TObjArray *) fProfileFile->Get("histArray");
    std::cout << fGFineArray << std::endl;
    TObjArray * arr = new TObjArray();
    fProfileFile->GetObject("angHistFineArray", arr);
    std::cout << arr << std::endl;
		fGFineArray 	  = (TObjArray *) fProfileFile->Get("angHistFineArray");
		fGCoarseArray 	  = (TObjArray *) fProfileFile->Get("angHistCoarseArray");
		fBinningHistogram = (TH1F*) fProfileFile->Get("hWhichHisto");

		fRhoArray->SetOwner(kTRUE);
		fGCoarseArray->SetOwner(kTRUE);
		fGFineArray->SetOwner(kTRUE);
		std::cout << "Loaded array" << std::endl;
	}

	fLastEnergy = myTrack->GetE();
	fLastType = myTrack->GetType();
	InterpolateProfiles(myTrack);

	if(fDebug) { SaveProfiles(); }
	return;

}

Double_t WCSimEmissionProfiles::GetCriticalDistance(
		WCSimLikelihoodTrack::TrackType type, Double_t energy) const {
	assert(type == WCSimLikelihoodTrack::MuonLike || type == WCSimLikelihoodTrack::ElectronLike);
	switch(type)
	{
	case WCSimLikelihoodTrack::MuonLike:
		return -398 + 0.756*energy - 8e-5 * energy * energy;
    break;
  default:
    assert(0);
	}
  return 0;
}

void WCSimEmissionProfiles::InterpolateProfiles(WCSimLikelihoodTrack* myTrack) {
	InterpolateRho(myTrack);
	InterpolateG(myTrack);
}

void WCSimEmissionProfiles::InterpolateRho(WCSimLikelihoodTrack* myTrack) {
	if(fRhoInterp != 0x0)
	{
    // std::cout << fRhoInterp << std::endl;
    // std::cout << fRhoInterp->GetName() << std::endl;
		delete fRhoInterp;
		fRhoInterp = 0x0;
	}

  int bin = GetArrayBin(myTrack->GetE());
	std::cerr << "Warning: not doing any interpolation, just getting the nearest profile - bin = " << bin << std::endl;
	fRhoInterp = (TH1F*)(fRhoArray->At(GetArrayBin(myTrack->GetE()))->Clone());
	fRhoInterp->SetDirectory(0);
	// std::cout << "Energy = " << myTrack->GetE() << "  integral = " << fRhoInterp->Integral("W");
	return;


	// There are two kinds of interpolation that need to happen here:
	// 1) Between emission profiles - need to smoothly transition from one profile to another
  //	  as energy increases, such that if energy exactly matches the emission profile binning
	//    we get that exact plot, and otherwise we combine the two plots either side
	// 2) Between bins of s - emission profiles are binned in s, so when we work out how much
	//    of each plot we take, we need to make sure there aren't discrete steps when the energy
	//    increases by just enough to take the track length into the next bin


	// Get the emission profiles for the energy bin below and above our current energy
	Double_t energy = myTrack->GetE();
  // std::cout << "Getting profiles from array bin " << GetArrayBin(energy) << std::endl;
	TH1F * profileLo = (TH1F*)fRhoArray->At(GetArrayBin(energy));
	TH1F * profileHi = (TH1F*)fRhoArray->At(GetArrayBin(energy)+1);

	// Get the energies these two profiles correspond to
	Double_t upEnergy = fBinningHistogram->GetXaxis()->GetBinUpEdge(GetArrayBin(energy)+1);
	Double_t downEnergy = fBinningHistogram->GetXaxis()->GetBinLowEdge(GetArrayBin(energy)+1);
	// std::cout << "upEnergy = " << upEnergy << std::endl;
  // std::cout << "downEnergy = " << downEnergy << std::endl;

  // What fraction of the way between the two energies is the energy we're considering?
	Double_t fracHi = GetFractionThroughBin(energy);
  if(fracHi < 0 || fracHi > 1.0){
    // std::cout << "energy = " << energy << "   fracHi = " << fracHi << std::endl;
    assert(!(fracHi < 0 || fracHi > 1.0));
  }

  // Now make some histograms to add together in these fractions
  fRhoInterpLo = (TH1F*)profileLo->Clone(); // Clone to steal the binning
  fRhoInterpHi = (TH1F*)profileHi->Clone();
  fRhoInterpLo->Reset();
  fRhoInterpHi->Reset();
  fRhoInterpLo->SetDirectory(0);
  fRhoInterpHi->SetDirectory(0);
  
	// This will be our final interpolated histogram
	fRhoInterp = new TH1F("fRho","fRho", profileLo->GetNbinsX(), profileLo->GetXaxis()->GetXmin(), profileLo->GetXaxis()->GetXmax());
  fRhoInterp->SetDirectory(0);

  // Get the distances to the shoulder for each of the three energies
  Double_t hiDist = GetCriticalDistance(fLastType, upEnergy);
  Double_t loDist = GetCriticalDistance(fLastType, downEnergy);
  Double_t truDist = GetCriticalDistance(fLastType, energy);
  // std::cout << "hiDist  = " << hiDist << std::endl;
  // std::cout << "loDist  = " << loDist << std::endl;
  // std::cout << "truDist  = " << truDist << std::endl;

  // Work out which bins to steal from the higher energy profile
  double distanceToSteal = truDist - loDist;
  double distanceToSkip = hiDist - distanceToSteal;
  double exactNumBins = distanceToSteal / (fRhoInterpHi->GetXaxis()->GetBinWidth(1));
  double exactStartBin = distanceToSkip / (fRhoInterpHi->GetXaxis()->GetBinWidth(1));
  // e.g. we need to fill in 3.2 bins' worth of distance, starting from bin 110.5
  // so bins to steal for roudning down:
  int nBinsToStealLo = (int)floor(exactNumBins);
  int nBinsToStealHi = nBinsToStealLo+1;
  int firstBinToStealLo = (int)ceil(exactStartBin);
  int firstBinToStealHi = firstBinToStealLo - 1;
  double scaleStolenBinsLo = 1.0 - (exactNumBins - (int)exactNumBins); // Does the ratio between high and low x bins
  double normStolenBinsLo = (1.0/profileHi->GetBinContent(firstBinToStealLo+nBinsToStealLo)) * profileLo->GetBinContent(1);
  double normStolenBinsHi = (1.0/profileHi->GetBinContent(firstBinToStealHi+nBinsToStealHi)) * profileLo->GetBinContent(1);
  // std::cout << "Scale stolen bins " << scaleStolenBinsLo << std::endl;
  // std::cout << "Norm stolen bins low: " << normStolenBinsLo << std::endl;
  // std::cout << "Norm stolen bins high: " << normStolenBinsHi << std::endl;



  // Work out how far to push to high energy histogram across
  double distanceToPush = hiDist - truDist;
  double exactNumBinsToPush = distanceToPush / (fRhoInterpHi->GetXaxis()->GetBinWidth(1));
  int firstBinForHiHistoLo = (int)floor(exactNumBinsToPush);
  int firstBinForHiHistoHi = firstBinForHiHistoLo+1;
  double scalePushHistoLo = 1.0 - (exactNumBinsToPush - (int)exactNumBinsToPush);


  // There's probably a much more efficient way of doing this that doesn't involve multiple loops, but
  // for now I want it to be readable...

  ////////////////////////////////////
  // First fill the low histogram:
  for(int iBin = 1; iBin <= nBinsToStealLo ; ++iBin)
  {
  	double content = fRhoInterpLo->GetBinContent(iBin);
  	fRhoInterpLo->SetBinContent(iBin, content + scaleStolenBinsLo * normStolenBinsLo * profileHi->GetBinContent(firstBinToStealLo+iBin-1));
  }
  for( int iBin = 1; iBin <= nBinsToStealHi; ++iBin)
  {
  	double content = fRhoInterpLo->GetBinContent(iBin);
  	fRhoInterpLo->SetBinContent(iBin, content + (1-scaleStolenBinsLo) * normStolenBinsHi * profileHi->GetBinContent(firstBinToStealHi+iBin-1));
  }


  // Now fill the shifted part
  for(int iBin = nBinsToStealLo+1; iBin < fRhoInterpLo->GetNbinsX(); ++iBin)
  {
  	double content = fRhoInterpLo->GetBinContent(iBin);
  	fRhoInterpLo->SetBinContent(iBin, content + scaleStolenBinsLo * profileLo->GetBinContent(iBin-nBinsToStealLo));
  }
  for(int iBin = nBinsToStealHi+1; iBin < fRhoInterpLo->GetNbinsX(); ++iBin)
  {
  	double content = fRhoInterpLo->GetBinContent(iBin);
  	fRhoInterpLo->SetBinContent(iBin, content + (1-scaleStolenBinsLo) * profileLo->GetBinContent(iBin - nBinsToStealHi));
  }
  ////////////////////////////////////

  //////////////////////////////////////
  // Now we'll fill the high energy one
  for(int iBin = 1; iBin < fRhoInterpHi->GetNbinsX(); ++iBin)
  {
  	double content = fRhoInterpHi->GetBinContent(iBin);
  	if(firstBinForHiHistoLo+iBin <= profileHi->GetNbinsX())
  	{
  		fRhoInterpHi->SetBinContent( iBin, content + scalePushHistoLo * profileHi->GetBinContent(firstBinForHiHistoLo + iBin));
  	}
  	content = fRhoInterpHi->GetBinContent(iBin);
  	if(firstBinForHiHistoHi+iBin <= profileLo->GetNbinsX())
  	{
  		fRhoInterpHi->SetBinContent( iBin, content + (1-scalePushHistoLo) * profileHi->GetBinContent(firstBinForHiHistoHi + iBin));
  	}
  }
  /////////////////////////////////////






  /*

  for(int i = 1 ; i <= nBinsToAverageLargeShift; ++i)
  {
    double hiContent = profileHi->GetBinContent(profileHi->GetXaxis()->FindBin(hiDist - truDist) + i - 1);
    if( hiContent > 0 )
    {
    	double toAddSmallShift = (profileLo->GetBinContent(i) / hiContent) / nBinsToAverageSmallShift;
    	double toAddLargeShift = (profileLo->GetBinContent(i) / hiContent) / nBinsToAverageLargeShift;
      scaleAddedLargeShift += toAddSmallShift;
      if(i <= nBinsToAverageSmallShift)
      {
      	scaleAddedSmallShift += toAddLargeShift;
      }
    }
  }

	for(Int_t iBin = 1; iBin <= fRhoInterp->GetNbinsX(); ++iBin)
	{
		// Double_t binContent = 0.0;
    Double_t binContentLo = 0.0;
    Double_t binContentHi = 0.0;
    /////

    // This is shifting the low histogram to the right
    if(iBin >= binForLowHistLo)
    {
      binContentLo += binFracForLoHist * profileLo->GetBinContent(iBin - binForLowHist);
    }

    if(iBin >= binForLowHistHi)
    {
    	binContentLo += (1.0 - binFracForLoHist) * profileLo->GetBinContent(iBin - binForLowHist);
    }

    // This is filling in the gap where we've shifted the low histogram
    // TODO: Need to double this up for low and high bins
    if(iBin < binForLowHistLo && iBin + nBinsToAverageSmallShift - 1 < fRhoInterp->GetXaxis()->GetNbins())
    {
      binContentLo = binFracForLoHist * profileHi->GetBinContent(iBin + nBinsToAverageSmallShift) * scaleAddedSmallShift;
    }
    if(iBin < binForLowHistLo && iBin + nBinsToAverageLargeShift - 1 < fRhoInterp->GetXaxis()->GetNbins())
    {
      binContentLo = (1-binFracForLoHist) * profileHi->GetBinContent(iBin + nBinsToAverageLargeShift) * scaleAddedLargeShift;
    }



    // This is shifting the high histogram to the left
    if( iBin + binForHiHistHi - 1 <= fRhoInterp->GetXaxis()->GetNbins() )
    {
      binContentHi += binFracForHiHist * profileHi->GetBinContent(iBin + binForHiHistHi - 1);
      if( iBin + binForHiHistLo - 1 <= fRhoInterp->GetXaxis()->GetNbins() )
      {
      	binContentHi += (1-binFracForHiHist) * profileHi->GetBinContent(iBin + binForHiHistLo - 1);
      }

    }

 */


    /////
    /*
		if(iBin + upShiftBin -1 <= fRhoInterp->GetNbinsX())
		{
			binContent += upScale * profileHi->GetBinContent(iBin + upShiftBin - 1 ) * fracHi;
      binContentHi = upScale * profileHi->GetBinContent(iBin + upShiftBin - 1) * fracHi;

		}

		if(iBin + toAddFrom <= toAddTo)
		{
			binContent += (1-fracHi) * profileLo->GetBinContent(iBin - toAddFrom + 1);
      binContentLo = (1-fracHi) * profileLo->GetBinContent(iBin - toAddFrom + 1);

		}
		else
		{
			binContent += (1-fracHi) * scaleAdded * profileHi->GetBinContent(iBin+toAddFrom);
      binContentLo = (1-fracHi) * scaleAdded * profileHi->GetBinContent(iBin+toAddFrom);
		}
  }
    */

  TString title = TString::Format("Energy = %.01f", energy);
  fRhoInterp->SetTitle(title.Data());
  fRhoInterpLo->SetLineColor(kAzure);
  fRhoInterpHi->SetLineColor(kMagenta);
  fRhoInterpLo->SetLineWidth(2);
  fRhoInterpHi->SetLineWidth(2);


  // TString canName1 = TString::Format("can_%01f_before.png", energy);
  // TCanvas * can1 = new TCanvas("can","can",700,500);
  // fRhoInterpLo->Draw();
  // fRhoInterpHi->Draw("SAME");
  // can1->SaveAs(canName1);
  // delete can1;

  fRhoInterpLo->Scale(1.0/fRhoInterpLo->Integral("width"));
  fRhoInterpHi->Scale(1.0/fRhoInterpHi->Integral("width"));

  fRhoInterp->Reset();
  fRhoInterp->Add(fRhoInterpLo, fRhoInterpHi, 1.-fracHi, fracHi);

  // TString canName = TString::Format("can_%.01f_rho.png",energy);

  // TCanvas * can = new TCanvas("can","can",700,500);
  // fRhoInterp->Draw();
  // fRhoInterpLo->Draw("SAME");
  // fRhoInterpHi->Draw("SAME");
  // can->SaveAs(canName.Data());
  // delete can;

  return;
}

void WCSimEmissionProfiles::InterpolateG(WCSimLikelihoodTrack* myTrack) {

	if(fGCoarse != 0x0)
	{
		delete fGCoarse;
		fGCoarse = 0x0;
	}

	if(fGFine != 0x0)
	{
		delete fGFine;
		fGFine = 0x0;
	}

  int bin = GetArrayBin(myTrack->GetE());
	std::cerr << "Warning: not doing any interpolation, just getting the nearest profile - bin = " << bin << std::endl;
	fGCoarse = (TH2F*)(fGCoarseArray->At(GetArrayBin(myTrack->GetE()))->Clone());
	fGCoarse->SetDirectory(0);
	fGFine = (TH2F*)(fGFineArray->At(GetArrayBin(myTrack->GetE()))->Clone());
	fGFine->SetDirectory(0);
	return;

/*
    TH2F * hGLo = (TH2F*)fGArray->At(GetArrayBin(energy));
    TH2F * hGHi = (TH2F*)fGArray->At(GetArrayBin(energy)+1);

    Double_t fracHi =  (myTrack->GetE() - fBinningHistogram->GetXaxis()->GetBinLowEdge(GetArrayBin(energy)+1))
    		 	            /(fBinningHistogram->GetXaxis()->GetBinWidth(GetArrayBin(energy)+1));
    Double_t fracLo = (1.0 - fracHi);

    if( fG != 0x0 ){ std::cout << "fG = " << fG << std::endl; std::cout << "Name = " << fG->GetName() << std::endl; delete fG; fG = 0x0;}
    fG = new TH2F("fG","fG", hGLo->GetNbinsX(), hGLo->GetXaxis()->GetXmin(), hGLo->GetXaxis()->GetXmax(),
    						 hGLo->GetNbinsY(), hGLo->GetYaxis()->GetXmin(), hGLo->GetYaxis()->GetXmax());
    fG->Add(hGLo, hGHi, fracLo, fracHi);
    fG->SetDirectory(0);

    return;
*/
}



Double_t WCSimEmissionProfiles::GetLightFlux(WCSimLikelihoodTrack * myTrack) const
{
  //std::cout << "*** WCSimChargeLikelihood::GetLightFlux() *** Getting the light flux as a function of the track KE" << std::endl;
  Double_t energy      = myTrack->GetE();
  WCSimLikelihoodTrack::TrackType type = myTrack->GetType();
  return this->GetLightFlux(type, energy);

}

std::vector<Double_t> WCSimEmissionProfiles::GetProfileEnergies() const
{
  // std::cout << "fBinningHistogram = " << fBinningHistogram << std::endl;
  // std::cout << "bins = "; std::cout << fBinningHistogram->GetNbinsX() << std::endl;
	std::vector<Double_t> energies;
	
  if(fBinningHistogram == 0x0)
  {
    std::cerr << "Error: binning histogram not yet loaded" << std::endl;
    assert(fBinningHistogram != 0x0);
  }

  int nBins = fBinningHistogram->GetNbinsX();
	for(int i = 1; i <= nBins; ++i)
	{
		energies.push_back(fBinningHistogram->GetXaxis()->GetBinLowEdge(i));
	}
	return energies;
}

Double_t WCSimEmissionProfiles::GetLightFlux(
		WCSimLikelihoodTrack::TrackType type, Double_t energy) const {
	assert(type == WCSimLikelihoodTrack::MuonLike || type == WCSimLikelihoodTrack::ElectronLike);
	Double_t flux = 0.0;
	Double_t factorFromG = 1/(4.0*TMath::Pi()); // We want solid angle fraction

	Double_t factorFromSolidAngle = 1; //2*TMath::Pi();
	// There was an overall normalization of 2pi that used to be neglected
	// because the solid angle subtended by a cone is 2pi(1 - cosTheta)
	// I've absorbed this into the SolidAngle now

  Double_t nPhotons = 0.0;

  double fudge = 1.0;
	switch(type)
	{
	case WCSimLikelihoodTrack::MuonLike:
		// From a linear fit to nPhotons as a function of energy (it's very linear)
	  // New fit on 10/March/15
	  nPhotons = 135.247 * energy -14823.6;

	  // Things get sketchy at really low energies
	  if(nPhotons > 0){ flux = factorFromG * factorFromSolidAngle * nPhotons; }
      break;
    case WCSimLikelihoodTrack::ElectronLike:
      nPhotons = 116.272 + 123.399 * energy;  // New fit on 10/March/15
      if(nPhotons > 0){ flux = factorFromG * factorFromSolidAngle * nPhotons; }
      break;
    default:
      assert(0);
      break;
	}
  flux = flux * fudge;
	return flux;
}

Double_t WCSimEmissionProfiles::GetTrackLengthForPercentile(Double_t percentile)
{
  Double_t runningTotal = 0.0;
  Int_t iBin = 0;
  while(iBin < fRhoInterp->GetNbinsX() && runningTotal < 0.75)
  {
      iBin++;
      runningTotal += fRhoInterp->GetBinContent(iBin) * fRhoInterp->GetXaxis()->GetBinWidth(iBin);
  }
  return fRhoInterp->GetXaxis()->GetBinCenter(iBin);

}

void WCSimEmissionProfiles::SaveProfiles()
{
  return;
  TCanvas * canRho = new TCanvas("canRho", "canRho", 800, 600);
  fRhoInterp->Draw();
  fRhoInterp->SetLineWidth(2);
  fRhoInterp->SetLineColor(kBlue);
  fRhoInterp->GetXaxis()->SetTitle("s/cm");
  fRhoInterp->GetYaxis()->SetTitle("rho(s)");
  fRhoInterp->GetXaxis()->CenterTitle();
  fRhoInterp->GetYaxis()->CenterTitle();
  TString rhoTitle("");
  rhoTitle.Form("rho_%d_", (int)fLastEnergy);
  rhoTitle += WCSimLikelihoodTrack::TrackTypeToString(fLastType);
  fRhoInterp->SetTitle(rhoTitle.Data());

  canRho->SaveAs( (rhoTitle + TString(".png")).Data());
  canRho->SaveAs( (rhoTitle + TString(".C")).Data());
  
  TCanvas * canRhoLo = new TCanvas("canRhoLo", "canRhoLo", 800, 600);
  fRhoInterpLo->Draw();
  fRhoInterpLo->SetLineWidth(2);
  fRhoInterpLo->SetLineColor(kRed);
  fRhoInterpLo->GetXaxis()->SetTitle("s/cm");
  fRhoInterpLo->GetYaxis()->SetTitle("rho(s)");
  fRhoInterpLo->GetXaxis()->CenterTitle();
  fRhoInterpLo->GetYaxis()->CenterTitle();
  TString rhoLoTitle("");
  rhoLoTitle.Form("rhoLo_%d_", (int)fLastEnergy);
  rhoLoTitle += WCSimLikelihoodTrack::TrackTypeToString(fLastType);
  fRhoInterpLo->SetTitle(rhoLoTitle.Data());
  canRhoLo->SaveAs( (rhoLoTitle + TString(".png")).Data());
  canRhoLo->SaveAs( (rhoLoTitle + TString(".C")).Data());
  
  TCanvas * canRhoHi = new TCanvas("canRhoHi", "canRhoHi", 800, 600);
  fRhoInterpHi->Draw();
  fRhoInterpHi->SetLineWidth(2);
  fRhoInterpHi->SetLineColor(kRed);
  fRhoInterpHi->GetXaxis()->SetTitle("s/cm");
  fRhoInterpHi->GetYaxis()->SetTitle("rho(s)");
  fRhoInterpHi->GetXaxis()->CenterTitle();
  fRhoInterpHi->GetYaxis()->CenterTitle();
  TString rhoHiTitle("");
  rhoHiTitle.Form("rhoHi_%d_", (int)fLastEnergy);
  rhoHiTitle += WCSimLikelihoodTrack::TrackTypeToString(fLastType);
  fRhoInterpHi->SetTitle(rhoHiTitle.Data());
  canRhoHi->SaveAs( (rhoHiTitle + TString(".png")).Data());
  canRhoHi->SaveAs( (rhoHiTitle + TString(".C")).Data());



  TCanvas * canG = new TCanvas("canG", "canG", 800, 600);
  fGCoarse->Draw("COLZ");
  fGCoarse->GetXaxis()->SetTitle("cos#theta");
  fGCoarse->GetYaxis()->SetTitle("s (cm)");
  fGCoarse->GetXaxis()->CenterTitle();
  fGCoarse->GetYaxis()->CenterTitle();
  TString gTitle("");
  gTitle.Form("gCoarse_%d_", (int)fLastEnergy);
  gTitle += WCSimLikelihoodTrack::TrackTypeToString(fLastType);
  fGCoarse->SetTitle(gTitle.Data());

  canG->SaveAs( (gTitle + TString(".png")).Data());
  canG->SaveAs( (gTitle + TString(".C")).Data());


  fGFine->Draw("COLZ");
  fGFine->GetXaxis()->SetTitle("cos#theta");
  fGFine->GetYaxis()->SetTitle("s (cm)");
  fGFine->GetXaxis()->CenterTitle();
  fGFine->GetYaxis()->CenterTitle();
  gTitle.Form("gFine_%d_", (int)fLastEnergy);
  gTitle += WCSimLikelihoodTrack::TrackTypeToString(fLastType);
  fGFine->SetTitle(gTitle.Data());

  canG->SaveAs( (gTitle + TString(".png")).Data());
  canG->SaveAs( (gTitle + TString(".C")).Data());
  delete canRho;
  delete canG;
  delete canRhoHi;
  delete canRhoLo;

  return;

}

Double_t WCSimEmissionProfiles::GetStoppingDistance(WCSimLikelihoodTrack * track)
{
	SetTrack(track);
	Int_t bin = fRhoInterp->FindLastBinAbove(0);
	return (fRhoInterp->GetBinCenter(bin));
}

TH1F * WCSimEmissionProfiles::GetRho(WCSimLikelihoodTrack::TrackType particle, double energy)
{
	WCSimLikelihoodTrack tempTrack;
	tempTrack.SetType(particle);
	tempTrack.SetE(energy);
	SetTrack(&tempTrack);
	return fRhoInterp;
}

TH2F * WCSimEmissionProfiles::GetGCoarse(WCSimLikelihoodTrack::TrackType particle, double energy)
{
	return GetG(particle, energy).first;
}

TH2F * WCSimEmissionProfiles::GetGFine(WCSimLikelihoodTrack::TrackType particle, double energy)
{
	return GetG(particle, energy).second;
}

std::pair<TH2F *, TH2F* > WCSimEmissionProfiles::GetG(WCSimLikelihoodTrack::TrackType particle, double energy)
{
	WCSimLikelihoodTrack tempTrack;
	tempTrack.SetType(particle);
	tempTrack.SetE(energy);
	SetTrack(&tempTrack);
	return std::make_pair(fGCoarse, fGFine);
}

TH1F * WCSimEmissionProfiles::GetEnergyHist(WCSimLikelihoodTrack::TrackType particle)
{
	WCSimLikelihoodTrack tempTrack;
	tempTrack.SetType(particle);
	tempTrack.SetE(1000.0);
	SetTrack(&tempTrack);
	return fBinningHistogram;
}
