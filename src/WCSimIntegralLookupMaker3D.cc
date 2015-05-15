/*
 * WCSimIntegralLookupMaker3D.cc
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#include <stdexcept>
#include "WCSimIntegralLookupMaker3D.hh"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimIntegralLookupMaker3D)
#endif


WCSimIntegralLookupMaker3D::WCSimIntegralLookupMaker3D( WCSimLikelihoodTrack::TrackType particle,
													int nR0Bins, double R0Min, double R0Max,
													int nCosTh0Bins, double cosTh0Min, double cosTh0Max )
: WCSimIntegralLookupMaker(particle, nR0Bins, R0Min, R0Max, nCosTh0Bins, cosTh0Min, cosTh0Max)
{
	std::cout << " *** WCSimIntegralLookupMaker3D for track type = " << WCSimLikelihoodTrack::TrackTypeToString(particle) << " *** " << std::endl;
	fRhoInt = 0x0;
	fRhoSInt = 0x0;
	fRhoSSInt = 0x0;
	fRhoGInt = 0x0;
	fRhoGSInt = 0x0;
	fRhoGSSInt = 0x0;
	fType = particle;

	TH1D * binningHist = fEmissionProfiles.GetEnergyHist(particle);
	TH1D * hRhoTemp = fEmissionProfiles.GetRho(particle, binningHist->GetXaxis()->GetXmin());

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

WCSimIntegralLookupMaker3D::~WCSimIntegralLookupMaker3D() {
	// TODO Auto-generated destructor stub
}


void WCSimIntegralLookupMaker3D::SetBins(const int &nEBins, const double &eMin, const double &eMax,
								  const int &nSBins, const double &sMin, const double &sMax,
								  const int &nR0Bins, const double &R0Min, const double &R0Max,
								  const int &nCosTh0Bins, const double &cosTh0Min, const double &cosTh0Max) {

	if(fRhoInt != 0x0){ delete fRhoInt; }
	if(fRhoSInt != 0x0){ delete fRhoSInt; }
	if(fRhoSSInt != 0x0){ delete fRhoSSInt; }

	if(fRhoGInt != 0x0){delete fRhoGInt;}
	if(fRhoGSInt != 0x0){delete fRhoGSInt;}
	if(fRhoGSSInt != 0x0){delete fRhoGSSInt;}

	fNEBins = nEBins;
	fEMin = eMin;
	fEMax = eMax;

	fNR0Bins = nR0Bins;
	fR0Min = R0Min;
	fR0Max = R0Max;

	fNCosTh0Bins = nCosTh0Bins;
	fCosTh0Min = cosTh0Min;
	fCosTh0Max = cosTh0Max;

	int binsRho[1] = {fNEBins};
	double minRho[1] = {fEMin};
	double maxRho[1] = {fEMax};
	fRhoInt = new THnSparseF("fRhoInt","Integral of #rho(s)",1, binsRho, minRho, maxRho);
	fRhoSInt = new THnSparseF("fRhoSInt","Integral of #rho(s) #times s",2, binsRho, minRho, maxRho);
	fRhoSSInt = new THnSparseF("fRhoSSInt","Integral of #rho(s) #times s^{2}",2, binsRho, minRho, maxRho);

	int binsG[3] = {fNEBins, fNR0Bins, fNCosTh0Bins};
	double minG[3] = {fEMin, fR0Min, fCosTh0Min};
	double maxG[3] = {fEMax, fR0Max, fCosTh0Max};
	fRhoGInt = new THnSparseF("fRhoGInt","Integral of #rho(s) #times g(s, cos#theta)", 3, binsG, minG, maxG);
	fRhoGSInt = new THnSparseF("fRhoGSInt","Integral of #rho(s) #times g(s, cos#theta) #times s", 3, binsG, minG, maxG);
	fRhoGSSInt = new THnSparseF("fRhoGSSInt","Integral of #rho(s) #times g(s, cos#theta) #times s^{2}", 3, binsG, minG, maxG);


	fIntegrals.SetHists(fRhoInt, fRhoSInt, fRhoSSInt, fRhoGInt, fRhoGSInt, fRhoGSSInt, 0x0);
}

void WCSimIntegralLookupMaker3D::MakeLookupTables() {
	MakeCutoffS();
	MakeRhoTables();
	MakeRhoGTables();
}

void WCSimIntegralLookupMaker3D::MakeRhoTables() {
	std::cout << "*** Making tables for rho *** " << std::endl;
	int binsFilled = 0;
    TH1D * binningHisto = fEmissionProfiles.GetEnergyHist(fType);

    for( int iEBin = 1; iEBin <= binningHisto->GetXaxis()->GetNbins(); ++iEBin)
    {
      std::cout << "Energy bin = " << iEBin << "/" << binningHisto->GetXaxis()->GetNbins() << std::endl;
      Double_t runningTotal[3] = {0.0, 0.0, 0.0};
	  Double_t ELowEdge[1] = {binningHisto->GetXaxis()->GetBinLowEdge(iEBin)}; // E and s for each bin - THnSparseF needs an array
	  TH1D * hRho = fEmissionProfiles.GetRho(fType, ELowEdge[0]);
	  for(int iSBin = 1; iSBin <= hRho->GetXaxis()->GetNbins(); ++iSBin)
	  {
		  double s = hRho->GetXaxis()->GetBinCenter(iSBin);
		  double toAdd = hRho->GetBinContent(iSBin) * hRho->GetBinWidth(iSBin);
		  if( toAdd > 0.0 )
		  {
			  runningTotal[0] += toAdd;
			  runningTotal[1] += toAdd * s;
			  runningTotal[2] += toAdd * s * s;
		  }
		  else
		  {
			  break;
		  }
	  }
	  ++binsFilled;
	  fIntegrals.GetRhoInt()->Fill(ELowEdge, runningTotal[0]);
	  fIntegrals.GetRhoSInt()->Fill(ELowEdge, runningTotal[1]);
	  fIntegrals.GetRhoSSInt()->Fill(ELowEdge, runningTotal[2]);
   }
   std::cout << "Rho bins filled=  " << binsFilled << std::endl;
}

void WCSimIntegralLookupMaker3D::MakeRhoGTables() {
	std::cout << " *** MakeRhoGTables *** " << std::endl;
    TH1D * binningHisto = fEmissionProfiles.GetEnergyHist(fType);

    TAxis axisR0(fNR0Bins, fR0Min, fR0Max);
    TAxis axisCosTh0(fNCosTh0Bins, fCosTh0Min, fCosTh0Max);

    int binsFilled = 0;
	for( int iEBin = 1; iEBin <= binningHisto->GetXaxis()->GetNbins(); ++iEBin)
	{
		//if(iEBin != binningHisto->GetXaxis()->GetNbins() ) { continue; }
		std::cout << "Energy bin = " << iEBin << "/" << binningHisto->GetXaxis()->GetNbins() << std::endl;
		Double_t vars[3] = { binningHisto->GetXaxis()->GetBinLowEdge(iEBin), 0.0, 0.0 };
		// Entries are E, s, R0, cosTh0
		TH1D * hRho = fEmissionProfiles.GetRho(fType, vars[0]);
		TH2D * hG = fEmissionProfiles.GetG(fType, vars[0]);

		for( Int_t iR0Bin = 1; iR0Bin <= fNR0Bins; ++iR0Bin)
		{
			Double_t R0 = axisR0.GetBinCenter(iR0Bin);
			vars[1] = R0;

			for( Int_t iCosTh0Bin = 1; iCosTh0Bin <= fNCosTh0Bins; ++iCosTh0Bin)
			{
				Double_t cosTh0 = axisCosTh0.GetBinCenter(iCosTh0Bin);
				vars[2] = cosTh0;

				// Three entries are powers of s = 0,1,2
			    Double_t runningTotal[3] = {0.0, 0.0, 0.0};

				for( Int_t iSBin = 0; iSBin < hRho->GetNbinsX(); ++iSBin)
				{
					// Work out the co-ordinates needed to look up the emission profiles at this step
					Double_t s = hRho->GetBinCenter(iSBin);

					Double_t sWidth = hRho->GetBinWidth(iSBin); // +1 is array numbering vs. bin numbering
					Double_t cosTh = (R0*cosTh0 - s)/TMath::Sqrt(R0*R0 + s*s - 2*R0*s*cosTh0);

					// Now get the values of the emission profiles
					double toAdd = hRho->GetBinContent(iSBin) * hG->GetBinContent(hG->GetXaxis()->FindBin(cosTh), iSBin) * sWidth;
					if( toAdd > 0.0 )
					{
						runningTotal[0] += toAdd;
						runningTotal[1] += toAdd * s;
						runningTotal[2] += toAdd * s * s;
						// std::cout << "vars = " << vars[0] << " " << vars[1] << " " << vars[2] << " " << vars[3] << " bin = " << fIntegrals.GetRhoGInt()->GetBin(vars) << std::endl;
					}
					else if(hRho->GetBinContent(iSBin) > 0) {break;}
				}
				fIntegrals.GetRhoGInt()->Fill(vars, runningTotal[0]);
				fIntegrals.GetRhoGSInt()->Fill(vars, runningTotal[1]);
				fIntegrals.GetRhoGSSInt()->Fill(vars, runningTotal[2]);
				if(runningTotal[0] > 0) { // std::cout << "E = " << vars[0] << "  R0 = " << vars[1] << "  cosTh0 = " << vars[2] << "  rhogint = " << runningTotal[0] << std::endl;
					++binsFilled;
				}
			}
		}
	}
	std::cout << "g bins filled = " << binsFilled << std::endl;
	std::cout << "Fraction of bins = " << fIntegrals.GetRhoGInt()->GetSparseFractionBins() << " and of memory " << fIntegrals.GetRhoGInt()->GetSparseFractionMem() << std::endl;
	std::cout << fIntegrals.GetRhoGInt()->GetNbins() << std::endl << std::endl;
}

void WCSimIntegralLookupMaker3D::Run(TString fileName) {
	MakeLookupTables();
	SaveLookupTables(fileName);
}

void WCSimIntegralLookupMaker3D::SaveLookupTables(TString fileName) {
	TFile * toSave = new TFile(fileName.Data(), "RECREATE");
    if (!toSave)
    {
        throw std::runtime_error("Could not open file");
    }
    if( toSave->IsZombie())
    {
    	throw std::runtime_error("Integral lookup TFile is a zombie");
    }

    fIntegrals.GetRhoInt()->Write();
    fIntegrals.GetRhoSInt()->Write();
    fIntegrals.GetRhoSSInt()->Write();
    fIntegrals.GetRhoGInt()->Write();
    fIntegrals.GetRhoGSInt()->Write();
    fIntegrals.GetRhoGSSInt()->Write();
    toSave->Close();


}
