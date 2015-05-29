/*
 * WCSimIntegralLookupMaker.cc
 *
 *  Created on: 9 Mar 2015
 *      Author: ajperch
 */

#include <stdexcept>
#include "WCSimIntegralLookupMaker.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimIntegralLookupMaker)
#endif

WCSimIntegralLookupMaker::WCSimIntegralLookupMaker( WCSimLikelihoodTrack::TrackType particle,
													int nR0Bins, double R0Min, double R0Max,
													int nCosTh0Bins, double cosTh0Min, double cosTh0Max )
{
	std::cout << " *** WCSimIntegralLookupMaker for track type = " << WCSimLikelihoodTrack::TrackTypeToString(particle) << " *** " << std::endl;
	fRhoInt = 0x0;
	fRhoSInt = 0x0;
	fRhoSSInt = 0x0;
	fRhoGInt = 0x0;
	fRhoGSInt = 0x0;
	fRhoGSSInt = 0x0;
	fCutoffS = 0x0;
	fType = particle;

	TH1F * binningHist = fEmissionProfiles.GetEnergyHist(particle);
	TH1F * hRhoTemp = fEmissionProfiles.GetRho(particle, binningHist->GetXaxis()->GetXmin());

	TCanvas * can1=  new TCanvas("can1","can1",800,600);
	binningHist->Draw();
	can1->SaveAs("can1.C");
	can1->SaveAs("can1.root");

	int eBins = binningHist->GetXaxis()->GetNbins();
	double eMax = binningHist->GetXaxis()->GetBinLowEdge(eBins) + binningHist->GetBinWidth(eBins-1);

	SetBins(binningHist->GetXaxis()->GetNbins(), binningHist->GetXaxis()->GetXmin(), eMax,
			hRhoTemp->GetXaxis()->GetNbins(), hRhoTemp->GetXaxis()->GetXmin(), hRhoTemp->GetXaxis()->GetXmax(),
			nR0Bins, R0Min, R0Max, nCosTh0Bins, cosTh0Min, cosTh0Max);
}

WCSimIntegralLookupMaker::WCSimIntegralLookupMaker()
{
	fRhoInt = 0x0;
	fRhoSInt = 0x0;
	fRhoSSInt = 0x0;
	fRhoGInt = 0x0;
	fRhoGSInt = 0x0;
	fRhoGSSInt = 0x0;
	fCutoffS = 0x0;


}

WCSimIntegralLookupMaker::~WCSimIntegralLookupMaker() {
	// TODO Auto-generated destructor stub
}


void WCSimIntegralLookupMaker::SetBins(const int &nEBins, const double &eMin, const double &eMax,
								  const int &nSBins, const double &sMin, const double &sMax,
								  const int &nR0Bins, const double &R0Min, const double &R0Max,
								  const int &nCosTh0Bins, const double &cosTh0Min, const double &cosTh0Max) {

	if(fRhoInt != 0x0){ delete fRhoInt; }
	if(fRhoSInt != 0x0){ delete fRhoSInt; }
	if(fRhoSSInt != 0x0){ delete fRhoSSInt; }

	if(fRhoGInt != 0x0){delete fRhoGInt;}
	if(fRhoGSInt != 0x0){delete fRhoGSInt;}
	if(fRhoGSSInt != 0x0){delete fRhoGSSInt;}

	if(fCutoffS != 0x0){ delete fCutoffS; }

	fNEBins = nEBins;
	fEMin = eMin;
	fEMax = eMax;

	fNSBins = nSBins;
	fSMin = sMin;
	fSMax = sMax;

	fNR0Bins = nR0Bins;
	fR0Min = R0Min;
	fR0Max = R0Max;

	fNCosTh0Bins = nCosTh0Bins;
	fCosTh0Min = cosTh0Min;
	fCosTh0Max = cosTh0Max;

	int binsRho[2] = {fNEBins, fNSBins};
	double minRho[2] = {fEMin, fSMin};
	double maxRho[2] = {fEMax, fSMax};
	fRhoInt = new THnSparseF("fRhoInt","Integral of #rho(s)",2, binsRho, minRho, maxRho);
	fRhoSInt = new THnSparseF("fRhoSInt","Integral of #rho(s) #times s",2, binsRho, minRho, maxRho);
	fRhoSSInt = new THnSparseF("fRhoSSInt","Integral of #rho(s) #times s^{2}",2, binsRho, minRho, maxRho);

	int binsG[4] = {fNEBins, fNSBins, fNR0Bins, fNCosTh0Bins};
	double minG[4] = {fEMin, fSMin, fR0Min, fCosTh0Min};
	double maxG[4] = {fEMax, fSMax, fR0Max, fCosTh0Max};
	fRhoGInt = new THnSparseF("fRhoGInt","Integral of #rho(s) #times g(s, cos#theta)", 4, binsG, minG, maxG);
	fRhoGSInt = new THnSparseF("fRhoGSInt","Integral of #rho(s) #times g(s, cos#theta) #times s", 4, binsG, minG, maxG);
	fRhoGSSInt = new THnSparseF("fRhoGSSInt","Integral of #rho(s) #times g(s, cos#theta) #times s^{2}", 4, binsG, minG, maxG);


	fCutoffS = new TH1F("fCutoffS","Distance after which there is no more emission", fNEBins, fEMin, fEMax);
	fIntegrals.SetHists(fRhoInt, fRhoSInt, fRhoSSInt, fRhoGInt, fRhoGSInt, fRhoGSSInt, fCutoffS);
}

void WCSimIntegralLookupMaker::MakeLookupTables() {
	MakeCutoffS();
	MakeRhoTables();
	MakeRhoGTables();
}

void WCSimIntegralLookupMaker::MakeCutoffS() {
	for( int iEBin = 1; iEBin <= fCutoffS->GetXaxis()->GetNbins(); ++iEBin)
	{
		std::cout << "Filling cutoff histogram" << std::endl;
		TH1F * hRho = fEmissionProfiles.GetRho(fType, fCutoffS->GetXaxis()->GetBinLowEdge(iEBin));
		int bin = hRho->FindLastBinAbove(0);
		std::cout << "Energy = " << fCutoffS->GetXaxis()->GetBinLowEdge(iEBin) << " cutoff = " << hRho->GetXaxis()->GetBinCenter(bin) << std::endl;
		fCutoffS->SetBinContent(iEBin, hRho->GetBinCenter(bin));
	}
}

void WCSimIntegralLookupMaker::MakeRhoTables() {
	std::cout << "*** Making tables for rho *** " << std::endl;
	int binsFilled = 0;
    TH1F * binningHisto = fEmissionProfiles.GetEnergyHist(fType);

    for( int iEBin = 1; iEBin <= binningHisto->GetXaxis()->GetNbins(); ++iEBin)
    {
      std::cout << "Energy bin = " << iEBin << "/" << binningHisto->GetXaxis()->GetNbins() << std::endl;
      Double_t runningTotal[3] = {0.0, 0.0, 0.0};
	  Double_t EandS[2] = {binningHisto->GetXaxis()->GetBinLowEdge(iEBin), 0.0}; // E and s for each bin - THnSparseF needs an array
	  TH1F * hRho = fEmissionProfiles.GetRho(fType, EandS[0]);
	  double cutoffS = fCutoffS->GetBinContent(iEBin);

	  for(int iSBin = 1; iSBin <= hRho->GetXaxis()->GetNbins(); ++iSBin)
	  {
		  EandS[1] = hRho->GetXaxis()->GetBinCenter(iSBin);
		  if(EandS[1] > cutoffS){
			  std::cout << "Exceeded the cutoff at s = " << EandS[1] << std::endl;
			  break;
		  }
		  double toAdd = hRho->GetBinContent(iSBin) * hRho->GetBinWidth(iSBin);
		  if( toAdd > 0.0 )
		  {
			  runningTotal[0] += toAdd;
			  runningTotal[1] += toAdd * EandS[1];
			  runningTotal[2] += toAdd * EandS[1] * EandS[1];

			  ++binsFilled;
			  fIntegrals.GetRhoInt()->Fill(EandS, runningTotal[0]);
			  fIntegrals.GetRhoSInt()->Fill(EandS, runningTotal[1]);
			  fIntegrals.GetRhoSSInt()->Fill(EandS, runningTotal[2]);
		  }
	  }
   }
   std::cout << "Rho bins filled=  " << binsFilled << std::endl;
}

void WCSimIntegralLookupMaker::MakeRhoGTables() {
	std::cout << " *** MakeRhoGTables *** " << std::endl;
    TH1F * binningHisto = fEmissionProfiles.GetEnergyHist(fType);

    TAxis axisR0(fNR0Bins, fR0Min, fR0Max);
    TAxis axisCosTh0(fNCosTh0Bins, fCosTh0Min, fCosTh0Max);

    int binsFilled = 0;
	for( int iEBin = 1; iEBin <= binningHisto->GetXaxis()->GetNbins(); ++iEBin)
	{
		//if(iEBin != binningHisto->GetXaxis()->GetNbins() ) { continue; }
		std::cout << "Energy bin = " << iEBin << "/" << binningHisto->GetXaxis()->GetNbins() << std::endl;
		Double_t vars[4] = { binningHisto->GetXaxis()->GetBinLowEdge(iEBin), 0.0, 0.0, 0.0 };
		// Entries are E, s, R0, cosTh0
		TH1F * hRho = fEmissionProfiles.GetRho(fType, vars[0]);
		std::pair<TH2F *, TH2F *> hG = fEmissionProfiles.GetG(fType, vars[0]);
		double cutoffS = fCutoffS->GetBinContent(iEBin);

		for( Int_t iR0Bin = 1; iR0Bin <= fNR0Bins; ++iR0Bin)
		{
			Double_t R0 = axisR0.GetBinCenter(iR0Bin);
			vars[2] = R0;

			for( Int_t iCosTh0Bin = 1; iCosTh0Bin <= fNCosTh0Bins; ++iCosTh0Bin)
			{
				Double_t cosTh0 = axisCosTh0.GetBinCenter(iCosTh0Bin);
				vars[3] = cosTh0;

				// Three entries are powers of s = 0,1,2
			    Double_t runningTotal[3] = {0.0, 0.0, 0.0};

				for( Int_t iSBin = 0; iSBin < fNSBins; ++iSBin)
				{
					// Work out the co-ordinates needed to look up the emission profiles at this step
					Double_t s = hRho->GetBinCenter(iSBin);
					if(s > cutoffS){
						// std::cout << "Breaking at s = " << s << std::endl;
						break;
					}
					vars[1] = s;

					Double_t sWidth = hRho->GetBinWidth(iSBin); // +1 is array numbering vs. bin numbering
					Double_t cosTh = (R0*cosTh0 - s)/TMath::Sqrt(R0*R0 + s*s - 2*R0*s*cosTh0);

          // Work out cos(theta) as a function of R0, s and cosTh0
          // I've broken down the fraction into part to save computations but
          // this is just applying cosine rule twice an rearranging
          double twoS = 2 * s;
          double sMinuR0CosTh0 = s - R0*cosTh0;
          double R0SqMinusSSq = (R0+s)*(R0-s);
          cosTh = (twoS * sMinuR0CosTh0) / ( R0SqMinusSSq + twoS*sMinuR0CosTh0);

					// Now get the values of the emission profiles
					int thetaBin = 0;
					double toAdd = 0.0;
					if( cosTh < hG.second->GetXaxis()->GetXmin())
					{
						thetaBin = hG.first->GetXaxis()->FindBin(cosTh);
						toAdd = hRho->GetBinContent(iSBin) * hG.first->GetBinContent(thetaBin, iSBin) * sWidth;
					}
					else
					{
						thetaBin = hG.second->GetXaxis()->FindBin(cosTh);
						toAdd = hRho->GetBinContent(iSBin) * hG.second->GetBinContent(thetaBin, iSBin) * sWidth;
					}

					if( toAdd > 0.0 )
					{
						runningTotal[0] += toAdd;
						runningTotal[1] += toAdd * vars[1];
						runningTotal[2] += toAdd * vars[1] * vars[1];
						// std::cout << "vars = " << vars[0] << " " << vars[1] << " " << vars[2] << " " << vars[3] << " bin = " << fIntegrals.GetRhoGInt()->GetBin(vars) << std::endl;

						fIntegrals.GetRhoGInt()->Fill(vars, runningTotal[0]);
						fIntegrals.GetRhoGSInt()->Fill(vars, runningTotal[1]);
						fIntegrals.GetRhoGSSInt()->Fill(vars, runningTotal[2]);
						++binsFilled;
					}
				}
			}
		}
	}
	std::cout << "g bins filled = " << binsFilled << std::endl;
	std::cout << "Fraction of bins = " << fIntegrals.GetRhoGInt()->GetSparseFractionBins() << " and of memory " << fIntegrals.GetRhoGInt()->GetSparseFractionMem() << std::endl;
	std::cout << fIntegrals.GetRhoGInt()->GetNbins() << std::endl << std::endl;
}

void WCSimIntegralLookupMaker::Run(TString fileName) {
	MakeLookupTables();
	SaveLookupTables(fileName);
}

void WCSimIntegralLookupMaker::SaveLookupTables(TString fileName) {
	TFile * toSave = new TFile(fileName.Data(), "RECREATE");
    if (!toSave)
    {
        throw std::runtime_error("Could not open file");
    }
    if( toSave->IsZombie())
    {
    	throw std::runtime_error("Integral lookup TFile is a zombie");
    }

    fIntegrals.GetCutoffS()->Write();
    fIntegrals.GetRhoInt()->Write();
    fIntegrals.GetRhoSInt()->Write();
    fIntegrals.GetRhoSSInt()->Write();
    fIntegrals.GetRhoGInt()->Write();
    fIntegrals.GetRhoGSInt()->Write();
    fIntegrals.GetRhoGSSInt()->Write();
    toSave->Close();


}
