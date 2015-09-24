/*
 * WCSimIntegralLookupMaker3D.cc
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#include <stdexcept>
#include "WCSimIntegralLookupMaker3D.hh"
#include "WCSimTrackParameterEnums.hh"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimIntegralLookupMaker3D)
#endif


WCSimIntegralLookupMaker3D::WCSimIntegralLookupMaker3D( TrackType::Type particle,
													int nR0Bins, double R0Min, double R0Max,
													int nCosTh0Bins, double cosTh0Min, double cosTh0Max )
{
	std::cout << " *** WCSimIntegralLookupMaker3D for track type = " << TrackType::AsString(particle) << " *** " << std::endl;
	fRhoInt = 0x0;
	fRhoSInt = 0x0;
	fRhoSSInt = 0x0;
	fRhoGInt = 0x0;
	fRhoGSInt = 0x0;
	fRhoGSSInt = 0x0;
	fType = particle;

	TH1F * binningHist = fEmissionProfileManager.GetEnergyHist(particle);
	TH1F * hRhoTemp = fEmissionProfileManager.GetRho(particle, binningHist->GetXaxis()->GetXmin());

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

	fRhoInt = new TH1F("fRhoInt","Integral of #rho(s)", fNEBins, fEMin, fEMax);
	fRhoSInt = new TH1F("fRhoSInt","Integral of #rho(s) #times s", fNEBins, fEMin, fEMax);
	fRhoSSInt = new TH1F("fRhoSSInt","Integral of #rho(s) #times s^{2}", fNEBins, fEMin, fEMax);

	fRhoGInt = new TH3F("fRhoGInt","Integral of #rho(s) #times g(s, cos#theta)", fNEBins, fEMin, fEMax, fNR0Bins, fR0Min, fR0Max, fNCosTh0Bins, fCosTh0Min, fCosTh0Max);
	fRhoGSInt = new TH3F("fRhoGSInt","Integral of #rho(s) #times g(s, cos#theta) #times s", fNEBins, fEMin, fEMax, fNR0Bins, fR0Min, fR0Max, fNCosTh0Bins, fCosTh0Min, fCosTh0Max);
	fRhoGSSInt = new TH3F("fRhoGSSInt","Integral of #rho(s) #times g(s, cos#theta) #times s^{2}", fNEBins, fEMin, fEMax, fNR0Bins, fR0Min, fR0Max, fNCosTh0Bins, fCosTh0Min, fCosTh0Max);

	fIntegrals.SetHists(fRhoInt, fRhoSInt, fRhoSSInt, fRhoGInt, fRhoGSInt, fRhoGSSInt, 0x0);
}

void WCSimIntegralLookupMaker3D::MakeLookupTables() {
	MakeRhoTables();
	MakeRhoGTables();
}

void WCSimIntegralLookupMaker3D::MakeRhoTables() {
	std::cout << "*** Making tables for rho *** " << std::endl;
	int binsFilled = 0;
    TH1F * binningHisto = fEmissionProfileManager.GetEnergyHist(fType);

    for( int iEBin = 1; iEBin <= binningHisto->GetXaxis()->GetNbins(); ++iEBin)
    {
      std::cout << "Energy bin = " << iEBin << "/" << binningHisto->GetXaxis()->GetNbins() << std::endl;
      Double_t runningTotal[3] = {0.0, 0.0, 0.0};
	  Double_t ELowEdge = binningHisto->GetXaxis()->GetBinLowEdge(iEBin); // E and s for each bin - THnSparseF needs an array
	  TH1F * hRho = fEmissionProfileManager.GetRho(fType, ELowEdge);
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
   // std::cout << "Rho bins filled=  " << binsFilled << std::endl;
}

void WCSimIntegralLookupMaker3D::MakeRhoGTables() {
	std::cout << " *** MakeRhoGTables *** " << std::endl;
    TH1F * binningHisto = fEmissionProfileManager.GetEnergyHist(fType);

    TAxis axisR0(fNR0Bins, fR0Min, fR0Max);
    TAxis axisCosTh0(fNCosTh0Bins, fCosTh0Min, fCosTh0Max);

    int binsFilled = 0;
	for( int iEBin = 1; iEBin <= binningHisto->GetXaxis()->GetNbins(); ++iEBin)
	{
		//if(iEBin != binningHisto->GetXaxis()->GetNbins() ) { continue; }
		std::cout << "Energy bin = " << iEBin << "/" << binningHisto->GetXaxis()->GetNbins() << std::endl;
    // std::cout << fIntegrals.GetRhoGSSInt()->GetBinContent(1613039);
		Double_t vars[3] = { binningHisto->GetXaxis()->GetBinLowEdge(iEBin), 0.0, 0.0 };
		// Entries are E, R0, cosTh0
		TH1F * hRho = fEmissionProfileManager.GetRho(fType, vars[0]);
		std::pair<TH2F*, TH2F*> hG = fEmissionProfileManager.GetG(fType, vars[0]);

		for( Int_t iR0Bin = 1; iR0Bin <= fNR0Bins; ++iR0Bin)
		{
			Double_t R0 = axisR0.GetBinLowEdge(iR0Bin);
			vars[1] = R0;

			for( Int_t iCosTh0Bin = 1; iCosTh0Bin <= fNCosTh0Bins; ++iCosTh0Bin)
			{
				Double_t cosTh0 = axisCosTh0.GetBinLowEdge(iCosTh0Bin);
				vars[2] = cosTh0;

				// Three entries are powers of s = 0,1,2
			    Double_t runningTotal[3] = {0.0, 0.0, 0.0};

				for( Int_t iSBin = 1; iSBin <= hRho->GetNbinsX(); ++iSBin)
				{
					// Work out the co-ordinates needed to look up the emission profiles at this step
					Double_t s = hRho->GetBinCenter(iSBin);

					// Work out cos(theta) as a function of R0, s and cosTh0
					// I've broken down the fraction into part to save computations but
					// this is just applying cosine rule twice an rearranging
					double R0CosTh0 = R0 * cosTh0;
					double cosTh = (R0CosTh0 - s) / TMath::Sqrt(s*s + R0*R0 - 2*s*R0CosTh0);

					double sWidth = hRho->GetBinWidth(iSBin); // +1 is array numbering vs. bin numbering
          double thetaWidth = 1; // We've made the emission profiles so that all we need are the bin contents
          /*if(thetaWidth > 5e-4)
          {
            std::cout << "next " << cosThNext << "   current " << cosTh << std::endl;
            std::cout << "sWidth = " << sWidth << "  thetaWidth = " << thetaWidth << std::endl;
          }*/
					// Now get the values of the emission profiles
					int thetaBin = 0;
					double toAdd = 0.0;
					if( cosTh < hG.second->GetXaxis()->GetXmin())
					{
						thetaBin = hG.first->GetXaxis()->FindBin(cosTh);
            // thetaWidth = hG.first->GetXaxis()->GetBinWidth(thetaBin);
						toAdd = hRho->GetBinContent(iSBin) * hG.first->GetBinContent(thetaBin, iSBin) * sWidth * thetaWidth;
					}
					else
					{
						thetaBin = hG.second->GetXaxis()->FindBin(cosTh);
            // thetaWidth = hG.second->GetXaxis()->GetBinWidth(thetaBin);
						toAdd = hRho->GetBinContent(iSBin) * hG.second->GetBinContent(thetaBin, iSBin) * sWidth * thetaWidth;
					}

					if( toAdd > 0.0 )
					{
						runningTotal[0] += toAdd;
						runningTotal[1] += toAdd * s;
						runningTotal[2] += toAdd * s * s;
					}
					//else if(hRho->GetBinContent(iSBin) == 0 && iSBin > 10) {break;} // Stop adding if we've reached the end of the profile
//          if( vars[0] == 1000 && 1285 <= vars[1] && 1290 >= vars[1] && 0.784 <= vars[2] && 0.786 > vars[2] )
//          {
// //            std::cout << "s = " << s << "  cosTh = " << cosTh << " - adding " << toAdd << "  integral[0] = " << runningTotal[0] << "  integral[1] = " << runningTotal[1] << "   integral[2] = " << runningTotal[2] << std::endl;
//					std::cout << "vars = " << vars[0] << " " << vars[1] << " " << vars[2] << " bin = " << fIntegrals.GetRhoGInt()->GetBin(vars, false) << std::endl;
//          }
				}
				fIntegrals.GetRhoGInt()->Fill(vars[0], vars[1], vars[2], runningTotal[0]);
				fIntegrals.GetRhoGSInt()->Fill(vars[0], vars[1], vars[2], runningTotal[1]);
				fIntegrals.GetRhoGSSInt()->Fill(vars[0], vars[1], vars[2], runningTotal[2]);
				if(runningTotal[0] > 0) {
		//			std::cout << "E = " << vars[0] << "  R0 = " << vars[1] << "  cosTh0 = " << vars[2] << "  rhogint = " << runningTotal[0] << std::endl;
	  			++binsFilled;
				}
			}
		}
	}
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
