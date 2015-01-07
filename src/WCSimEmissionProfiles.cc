/*
 * WCSimEmissionProfiles.cxx
 *
 *  Created on: 13 Nov 2014
 *      Author: andy
 */

#include "WCSimEmissionProfiles.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodDigit.hh"
#include <TVector3.h>
#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TH1D.h>
#include <TH2D.h>
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
	fGArray = 0x0;
	fBinningHistogram = 0x0;

	fRhoInterp = 0x0;
	fG = 0x0;

}

WCSimEmissionProfiles::~WCSimEmissionProfiles() {
	if(fProfileFile != 0x0)
	{
    fProfileFile->Close();
    delete fProfileFile;
    fProfileFile = 0x0;
	}

	if( fRhoInterp != 0x0 )
	{
		delete fRhoInterp;
		fRhoInterp = 0x0;
	}

	if(	fG != 0x0 )
	{
		delete fG;
		fG = 0x0;
	}
}

WCSimEmissionProfiles::WCSimEmissionProfiles(WCSimLikelihoodTrack* myTrack) {

	fProfileFileName = TString("");

	fProfileFile = 0x0;

	fRhoArray = 0x0;
	fGArray = 0x0;
	fBinningHistogram = 0x0;

	fRhoInterp = 0x0;
	fG = 0x0;

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
	return fG->GetBinContent(fG->GetXaxis()->FindBin(cosTheta), fG->GetYaxis()->FindBin(s));
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
      //          << "fRhoInterp = " << fRhoInterp->GetBinContent(iBin) << std::endl 
      //          << "fG         = " << fG->GetBinContent(fG->GetXaxis()->FindBin(cosTheta), iBin) << std::endl;

      integrals.at(iPower) += binWidth
			  	                   * fRhoInterp->GetBinContent(iBin)
				                     * fG->GetBinContent(fG->GetXaxis()->FindBin(cosTheta), iBin)
				                     * pow(fRhoInterp->GetBinCenter(iBin), sPowers.at(iPower));
      // std::cout << "Now integrals.at(" << iPower << ") = " << integrals.at(iPower) << std::endl;
	  }
  }
	return integrals;

}

UInt_t WCSimEmissionProfiles::GetArrayBin(Double_t energy) const{
	int iBin = fBinningHistogram->GetXaxis()->FindBin(energy) - 1;
  if(iBin < 0) { std::cerr << "Warning, array bin number is less than zero - returning 0 instead" << std::endl;}
  return iBin; 
}

UInt_t WCSimEmissionProfiles::GetHistBin(Double_t energy) const{
  return GetArrayBin(energy)+1; 
}

Double_t WCSimEmissionProfiles::GetFractionThroughBin(Double_t energy) const{
	UInt_t iBin = GetHistBin(energy);
	Double_t fraction =   (energy - fBinningHistogram->GetXaxis()->GetBinLowEdge(iBin))
						/ (fBinningHistogram->GetXaxis()->GetBinWidth(iBin));
	return fraction;
}

void WCSimEmissionProfiles::LoadFile(WCSimLikelihoodTrack* myTrack) {
	std::cout << " *** WCSimEmissionProfiles::LoadFile - Loading profile" << std::endl;
	if( myTrack->GetType() != fLastType )
	{

		fProfileFileName = TString(getenv("WCSIMANAHOME"));
		switch( myTrack->GetType() )
		{
			case WCSimLikelihoodTrack::ElectronLike:
				fProfileFileName.Append("/config/emissionProfilesElectron.root");
				break;
			case WCSimLikelihoodTrack::MuonLike:
				fProfileFileName.Append("/config/emissionProfiles.root");
				break;
			default:
				std::cerr << "Error: unknown track type in WCSimLikelihoodTuner::LoadEmissionProfiles" << std::endl;
				exit(EXIT_FAILURE);
		}


		if( fProfileFile != NULL )
		{
			delete fProfileFile;
			delete fRhoArray;
			delete fGArray;
		}

		fProfileFile      = new TFile(fProfileFileName.Data(),"READ");
		fLastType         = myTrack->GetType();
		fRhoArray         = (TObjArray *) fProfileFile->Get("histArray");
		fGArray 		  = (TObjArray *) fProfileFile->Get("angHistArray");
		fBinningHistogram = (TH1D*) fProfileFile->Get("hWhichHisto");

		fRhoArray->SetOwner(kTRUE);
		fGArray->SetOwner(kTRUE);
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
    std::cout << fRhoInterp << std::endl;
    std::cout << fRhoInterp->GetName() << std::endl;
		delete fRhoInterp;
		fRhoInterp = 0x0;
	}

	Double_t energy = myTrack->GetE();
  std::cout << "Getting profiles from array bin " << GetArrayBin(energy) << std::endl;
	TH1D * profileLo = (TH1D*)fRhoArray->At(GetArrayBin(energy));
	TH1D * profileHi = (TH1D*)fRhoArray->At(GetArrayBin(energy)+1);

  fRhoInterpLo = (TH1D*)profileLo->Clone();
  fRhoInterpHi = (TH1D*)profileHi->Clone();
  fRhoInterpLo->Reset();
  fRhoInterpHi->Reset();
  fRhoInterpLo->SetDirectory(0);
  fRhoInterpHi->SetDirectory(0);
  
	Double_t fracHi = GetFractionThroughBin(energy);
  if(fracHi < 0 || fracHi > 1.0){ 
    // std::cout << "energy = " << energy << "   fracHi = " << fracHi << std::endl;
    assert(!(fracHi < 0 || fracHi > 1.0));
  }


	// Chop some distance off the front of the higher energy histogram and rescale:
	Double_t upEnergy = fBinningHistogram->GetXaxis()->GetBinUpEdge(GetArrayBin(energy)+1);
  /*
	Double_t upShift = GetCriticalDistance(fLastType, upEnergy) - GetCriticalDistance(fLastType, energy);
  if( !(upShift >= 0 && upShift <= fBinningHistogram->GetXaxis()->GetBinWidth(GetArrayBin(energy)+1)))
  {
    std::cout << "Ruh-roh, Raggy - upShift = " << upShift << " and bin width = " << fBinningHistogram->GetXaxis()->GetBinWidth(GetArrayBin(energy)+1) << std::endl;
    assert(upShift >= 0 && upShift <= fBinningHistogram->GetXaxis()->GetBinWidth(GetArrayBin(energy)+1));
  }
	Int_t upShiftBin = profileHi->GetXaxis()->FindBin(upShift);
  Double_t upScale = GetLightFlux(fLastType, upEnergy) / GetLightFlux(fLastType, energy);
  std::cout << "Energy: " << energy << "   upshift = " << upShift << "    bins = " << upShiftBin << "   scale = " << upScale << std::endl;
  */

	// Add some distance on to the lower energy histogram and rescale
  std::cout << "upEnergy = " << upEnergy << std::endl;
	Double_t downEnergy = fBinningHistogram->GetXaxis()->GetBinLowEdge(GetArrayBin(energy)+1);
  std::cout << "downEnergy = " << downEnergy << std::endl;
	/*
  Int_t toAddTo   = profileHi->GetXaxis()->FindBin(GetCriticalDistance(fLastType, upEnergy) - GetCriticalDistance(fLastType, energy));
	*/

	// Now add the two histograms together fractionally so that we smoothly transition from the profile
	// in the lower energy bin to the profile in the higher one
	fRhoInterp = new TH1D("fRho","fRho", profileLo->GetNbinsX(), profileLo->GetXaxis()->GetXmin(), profileLo->GetXaxis()->GetXmax());
  fRhoInterp->SetDirectory(0);

  // Get the distances to the shoulder for each of the three energies
  Double_t hiDist = GetCriticalDistance(fLastType, upEnergy);
  Double_t loDist = GetCriticalDistance(fLastType, downEnergy);
  Double_t truDist = GetCriticalDistance(fLastType, energy);
  // std::cout << "hiDist  = " << hiDist << std::endl;
  // std::cout << "loDist  = " << loDist << std::endl;
  // std::cout << "truDist  = " << truDist << std::endl;

  // Convert the high and low distances into bin numbers
  Int_t binForLowHist = fRhoInterp->GetXaxis()->FindBin(truDist - loDist);
  Int_t binForHiHist = fRhoInterp->GetXaxis()->FindBin(hiDist - truDist);
  // std::cout << "binForLowHist = " << binForLowHist << std::endl;
  // std::cout << "binForHiHist  = " << binForHiHist << std::endl;
  
  Double_t scaleAdded = 0;
  Int_t nBinsToAverage = (profileHi->GetXaxis()->FindBin(hiDist - loDist) - profileHi->GetXaxis()->FindBin(hiDist - truDist));

  for(int i = 1 ; i <= nBinsToAverage; ++i)
  {
    double hiContent = profileHi->GetBinContent(profileHi->GetXaxis()->FindBin(hiDist - truDist) + i - 1);
    if( hiContent > 0 )
    {
      scaleAdded += (profileLo->GetBinContent(i) / hiContent) / nBinsToAverage;
    }
  }

	for(Int_t iBin = 1; iBin <= fRhoInterp->GetNbinsX(); ++iBin)
	{
		// Double_t binContent = 0.0;
    Double_t binContentLo = 0.0;
    Double_t binContentHi = 0.0;
    /////
    if(iBin >= binForLowHist)
    {
      binContentLo = profileLo->GetBinContent(iBin - binForLowHist + 1);
    }

    if(iBin < binForLowHist && iBin + binForHiHist - 1 < fRhoInterp->GetXaxis()->GetNbins())
    {
      binContentLo = profileHi->GetBinContent(iBin + binForHiHist - 1) * scaleAdded;
    }

    if( iBin + binForHiHist - 1 <= fRhoInterp->GetXaxis()->GetNbins() )
    {
      binContentHi = profileHi->GetBinContent(iBin + binForHiHist - 1);
    }
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
    */

    fRhoInterpLo->SetBinContent(iBin, binContentLo);
    fRhoInterpHi->SetBinContent(iBin, binContentHi);
	}
  fRhoInterpLo->Scale(1.0/fRhoInterpLo->Integral());
  fRhoInterpHi->Scale(1.0/fRhoInterpHi->Integral());

  fRhoInterp->Reset();
  fRhoInterp->Add(fRhoInterpLo, fRhoInterpHi, 1.-fracHi, fracHi);
  
  // TCanvas * can = new TCanvas("can","can",700,500);
  // fRhoInterp->Draw();
  // fRhoInterpLo->SetLineColor(kAzure);
  // fRhoInterpHi->SetLineColor(kMagenta);
  // fRhoInterpLo->SetLineWidth(2);
  // fRhoInterpHi->SetLineWidth(2);
  // fRhoInterpLo->Draw("SAME");
  // fRhoInterpHi->Draw("SAME");
  // can->SaveAs("can.png");
  // delete can;
  return;
}

void WCSimEmissionProfiles::InterpolateG(WCSimLikelihoodTrack* myTrack) {
	Double_t energy = myTrack->GetE();
    TH2D * hGLo = (TH2D*)fGArray->At(GetArrayBin(energy));
    TH2D * hGHi = (TH2D*)fGArray->At(GetArrayBin(energy)+1);

    Double_t fracHi =  (myTrack->GetE() - fBinningHistogram->GetXaxis()->GetBinLowEdge(GetArrayBin(energy)+1))
    		 	            /(fBinningHistogram->GetXaxis()->GetBinWidth(GetArrayBin(energy)+1));
    Double_t fracLo = (1.0 - fracHi);

    if( fG != 0x0 ){ std::cout << "fG = " << fG << std::endl; std::cout << "Name = " << fG->GetName() << std::endl; delete fG; fG = 0x0;}
    fG = new TH2D("fG","fG", hGLo->GetNbinsX(), hGLo->GetXaxis()->GetXmin(), hGLo->GetXaxis()->GetXmax(),
    						 hGLo->GetNbinsY(), hGLo->GetYaxis()->GetXmin(), hGLo->GetYaxis()->GetXmax());
    fG->Add(hGLo, hGHi, fracLo, fracHi);
    fG->SetDirectory(0);

    return;

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
	switch(type)
	{
		case WCSimLikelihoodTrack::MuonLike:
			// From a linear fit to nPhotons as a function of energy (it's very linear)
			nPhotons = 246.564 * energy - 20700;
			// Things get sketchy at really low energies
			if(nPhotons > 0){ flux = factorFromG * factorFromSolidAngle * nPhotons; }
      break;
    default:
      assert(0);
      break;
	}

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
  fG->Draw("COLZ");
  fG->GetXaxis()->SetTitle("cos#theta");
  fG->GetYaxis()->SetTitle("s (cm)");
  fG->GetXaxis()->CenterTitle();
  fG->GetYaxis()->CenterTitle();
  TString gTitle("");
  gTitle.Form("g_%d_", (int)fLastEnergy);
  gTitle += WCSimLikelihoodTrack::TrackTypeToString(fLastType);
  fG->SetTitle(gTitle.Data());

  canG->SaveAs( (gTitle + TString(".png")).Data());
  canG->SaveAs( (gTitle + TString(".C")).Data());

  delete canRho;
  delete canG;
  delete canRhoHi;
  delete canRhoLo;

  return;

}


Double_t WCSimEmissionProfiles::GetLightFlux(WCSimLikelihoodTrack * myTrack)
{
  //std::cout << "*** WCSimChargeLikelihood::GetLightFlux() *** Getting the light flux as a function of the track KE" << std::endl;
  Double_t energy      = myTrack->GetE();
  Double_t factorFromG = 1/(4.0*TMath::Pi()); // We want solid angle fraction

  Double_t factorFromSolidAngle = 1; //2*TMath::Pi();
  // There was an overall normalization of 2pi that used to be neglected
  // because the solid angle subtended by a cone is 2pi(1 - cosTheta)
  // I've absorbed this into the SolidAngle now

  // Double_t nPhotons = 341.726*fEnergy;    // <-- EARLIER PARAMETRISATION
  // Double_t nPhotons = 36.9427* - 3356.71; // <-- EARLIER PARAMETRISATION

  // From a linear fit to nPhotons as a function of energy (it's very linear)
  Double_t nPhotons = 17.91 + 39.61*energy;

  // Things get sketchy at really low energies
  if(nPhotons < 0) return 0;

  return factorFromG * factorFromSolidAngle * nPhotons;

}
