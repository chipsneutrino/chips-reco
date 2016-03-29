/*
 * WCSimIntegralLookupMaker3D.cc
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */


#include <algorithm>
#include <stdexcept>
#include "WCSimIntegralLookupMaker3D.hh"
#include "WCSimTrackParameterEnums.hh"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TString.h>
#include <TSpline.h>
#include <TTree.h>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimIntegralLookupMaker3D)
#endif


WCSimIntegralLookupMaker3D::WCSimIntegralLookupMaker3D( TrackType::Type particle,
													int nR0Bins, double R0Min, double R0Max,
													int nCosTh0Bins, double cosTh0Min, double cosTh0Max )
{
	std::cout << " *** WCSimIntegralLookupMaker3D for track type = " << TrackType::AsString(particle) << " *** " << std::endl;
	fType = particle;

	TH1F * binningHist = fEmissionProfileManager.GetEnergyHist(particle);
	TH1F * hRhoTemp = fEmissionProfileManager.GetRho(particle, binningHist->GetXaxis()->GetXmin());

	int eBins   = binningHist->GetXaxis()->GetNbins();
	double eMax = binningHist->GetXaxis()->GetBinLowEdge(eBins) + binningHist->GetBinWidth(eBins-1);

    std::cout << "About to set bins" << std::endl;
    fIntegrals.Verify(); 
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
	fEMin   = eMin;
	fEMax   = eMax;

	fNR0Bins = nR0Bins;
	fR0Min   = R0Min;
	fR0Max   = R0Max;

	fNCosTh0Bins = nCosTh0Bins;
	fCosTh0Min   = cosTh0Min;
	fCosTh0Max   = cosTh0Max;

    fIntegrals = WCSimIntegralLookupHistArray(fR0Max, fR0Min, fNR0Bins, fCosTh0Max, fCosTh0Min, fNCosTh0Bins);
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
      std::cout << "\rEnergy bin = " << iEBin << "/" << binningHisto->GetXaxis()->GetNbins() << std::flush;
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
      
      // This is to avoid having graphs full of zeros - we don't 
      // need to was memory/disk space by splining those, we can just make the 
      // lookup return zero
      if(runningTotal[0] > 0 || fIntegrals.GetRhoInt()->GetN() > 0)
      {
	    fIntegrals.GetRhoInt()->SetPoint(fIntegrals.GetRhoInt()->GetN(), ELowEdge, runningTotal[0]);
	    fIntegrals.GetRhoSInt()->SetPoint(fIntegrals.GetRhoSInt()->GetN(), ELowEdge, runningTotal[1]);
	    fIntegrals.GetRhoSSInt()->SetPoint(fIntegrals.GetRhoSSInt()->GetN(), ELowEdge, runningTotal[2]);
      }
   }
   std::cout << std::endl;
   // std::cout << "Rho bins filled=  " << binsFilled << std::endl;
}

void WCSimIntegralLookupMaker3D::MakeRhoGTables() {
	std::cout << " *** MakeRhoGTables *** " << std::endl;
    TH1F * binningHisto = fEmissionProfileManager.GetEnergyHist(fType);

    TAxis axisR0(fNR0Bins, fR0Min, fR0Max);
    TAxis axisCosTh0(fNCosTh0Bins, fCosTh0Min, fCosTh0Max);

    int binsFilled = 0;

    std::vector<TH1F*> rhoHistVec;
    std::vector<TH2F*> gCoarseHistVec;
    std::vector<TH2F*> gFineHistVec;

    for(int iEBin = 1; iEBin <= binningHisto->GetNbinsX(); ++iEBin)
    {
        double energy = binningHisto->GetXaxis()->GetBinLowEdge(iEBin);
        rhoHistVec.push_back(static_cast<TH1F*>(fEmissionProfileManager.GetRho(fType, energy)->Clone()));
        gCoarseHistVec.push_back(static_cast<TH2F*>(fEmissionProfileManager.GetG(fType, energy).first->Clone()));
        gFineHistVec.push_back(static_cast<TH2F*>(fEmissionProfileManager.GetG(fType, energy).second->Clone()));
        rhoHistVec.back()->SetDirectory(0);
        gCoarseHistVec.back()->SetDirectory(0);
        gFineHistVec.back()->SetDirectory(0);

    }


    int nBins = fNR0Bins * fNCosTh0Bins;
    int lastPercent = -1;
	for( Int_t iR0Bin = 1; iR0Bin <= fNR0Bins; ++iR0Bin)
	{
		Double_t R0 = axisR0.GetBinLowEdge(iR0Bin);

		for( Int_t iCosTh0Bin = 1; iCosTh0Bin <= fNCosTh0Bins; ++iCosTh0Bin)
		{
            int percent = (iCosTh0Bin-1) * (iR0Bin-1) * 100 / nBins;
            if(percent > lastPercent)
            {
                std::cout << "\rProcessed " << percent << "%" << std::flush;
                lastPercent = percent;
            }
			Double_t cosTh0 = axisCosTh0.GetBinLowEdge(iCosTh0Bin);

			// Three entries are powers of s = 0,1,2
			for( Int_t iEBin = 1; iEBin <= binningHisto->GetNbinsX(); ++iEBin)
			{
				Double_t runningTotal[3] = {0.0, 0.0, 0.0};
				Double_t energy = binningHisto->GetXaxis()->GetBinLowEdge(iEBin);
				TH1F * hRho = rhoHistVec.at(iEBin-1);
				std::pair<TH2F*, TH2F*> hG = std::make_pair(gCoarseHistVec.at(iEBin-1), gFineHistVec.at(iEBin-1));
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


					// Now get the values of the emission profiles
					int thetaBin = 0;
					double toAdd = 0.0;

					if( cosTh < hG.second->GetXaxis()->GetXmin())
					{
						thetaBin = hG.first->GetXaxis()->FindBin(cosTh);
						toAdd = hRho->GetBinContent(iSBin) * hG.first->GetBinContent(thetaBin, iSBin) * sWidth * thetaWidth;
					}
					else
					{
						thetaBin = hG.second->GetXaxis()->FindBin(cosTh);
						toAdd = hRho->GetBinContent(iSBin) * hG.second->GetBinContent(thetaBin, iSBin) * sWidth * thetaWidth;
					}

					if( toAdd > 0.0 )
					{
						runningTotal[0] += toAdd;
						runningTotal[1] += toAdd * s;
						runningTotal[2] += toAdd * s * s;
					}
				}

				TGraph * rhoGInt = fIntegrals.GetRhoGInt(R0, cosTh0);
				TGraph * rhoGSInt = fIntegrals.GetRhoGSInt(R0, cosTh0);
				TGraph * rhoGSSInt = fIntegrals.GetRhoGSSInt(R0, cosTh0);

				// This is to avoid having graphs full of zeros - we don't 
                // need to was memory/disk space by splining those, we can just make the 
                // lookup return zero
                if( runningTotal[0] > 0 || rhoGInt->GetN() > 0)
                {
                    rhoGInt->SetPoint(rhoGInt->GetN(), energy, runningTotal[0]);
				    rhoGSInt->SetPoint(rhoGSInt->GetN(), energy, runningTotal[1]);
				    rhoGSSInt->SetPoint(rhoGSSInt->GetN(), energy, runningTotal[2]);
                }
				if(runningTotal[0] > 0)
				{
					//std::cout << "E = " << vars[0] << "  R0 = " << vars[1] << "  cosTh0 = " << vars[2] << "  rhogint = " << runningTotal[0] << std::endl;
					++binsFilled;
				}
			}
		}
	}
    std::cout << std::endl;

    for(int iEBin = 1; iEBin <= binningHisto->GetNbinsX(); ++iEBin)
    {
        delete rhoHistVec.at(iEBin-1);
        delete gCoarseHistVec.at(iEBin-1);
        delete gFineHistVec.at(iEBin-1);
    }
    rhoHistVec.clear();
    gCoarseHistVec.clear();
    gFineHistVec.clear();

}

void WCSimIntegralLookupMaker3D::Run(TString fileName) {
	MakeLookupTables();
    SmoothLookupTables();
	MakeSplines();
	SaveLookupTables(fileName);
}

void WCSimIntegralLookupMaker3D::MakeRhoGSplines() {
	std::cout << "*** MakeRhoGSplines ***" << std::endl;
    int nBins = fNR0Bins * fNCosTh0Bins;
    int lastPercent = -1;
	for(int iR0Bin = 0; iR0Bin < fNR0Bins; ++iR0Bin)
	{
		double R0  = fR0Min + iR0Bin * (fR0Max - fR0Min)/fNR0Bins;
		for(int iCosTh0Bin = 0; iCosTh0Bin < fNCosTh0Bins; ++iCosTh0Bin)
		{
            int bin = iCosTh0Bin + (iR0Bin * fNCosTh0Bins);
            int percent = bin * 100 / nBins;
            if(percent > lastPercent)
            {
                std::cout << "\rProcessed " << percent << "%" << std::flush;
                lastPercent = percent;
            }
			double cosTh0 = fCosTh0Min + iCosTh0Bin * (fCosTh0Max - fCosTh0Min)/fNCosTh0Bins;

			TString name = TString::Format("rhoGInt_%f_%f", R0, cosTh0);
			TGraph * gr  = fIntegrals.GetRhoGInt(R0, cosTh0);
            TSpline3 * rhoGSpline = NULL;
            if(gr->GetN() > 0)
            {
                std::cout << "Making spline" << std::endl;
			    rhoGSpline = new TSpline3(name.Data(), gr);
                rhoGSpline->SetName(name.Data());
            }
			fIntegrals.SetRhoGSpline(R0, cosTh0, rhoGSpline);


			name = TString::Format("rhoGSInt_%f_%f", R0, cosTh0);
            TGraph * grS  =  fIntegrals.GetRhoGSInt(R0, cosTh0);
            if(R0 == 1700 && cosTh0 == 0.96){ std::cout << "RhoGSInt = " << grS << std::endl; }
            TSpline3 * rhoGSSpline = NULL;
            if(grS->GetN() > 0)
            {
			    rhoGSSpline = new TSpline3(name.Data(), grS);
                rhoGSSpline->SetName(name.Data());
            }
			fIntegrals.SetRhoGSSpline(R0, cosTh0, rhoGSSpline);

			name		 	   = TString::Format("rhoGSSInt_%f_%f", R0, cosTh0);
			TGraph * grSS      = fIntegrals.GetRhoGSSInt(R0, cosTh0);
            if(R0 == 1700 && cosTh0 == 0.96){ std::cout << "RhoGSSInt = " << grSS << std::endl; }
            TSpline3 * rhoGSSSpline = NULL;
            if(grSS->GetN() > 0)
            {
			    rhoGSSSpline = new TSpline3(name.Data(), grSS);
                rhoGSSSpline->SetName(name.Data());
            }
			fIntegrals.SetRhoGSSSpline(R0, cosTh0, rhoGSSSpline);
		}
	}
    std::cout << std::endl;
}

void WCSimIntegralLookupMaker3D::MakeRhoSplines() {
	std::cout << "*** MakeRhoSplines ***" << std::endl;
	TSpline3 rhoSSpline("fRhoSSpline", fIntegrals.GetRhoSInt());
	TSpline3 rhoSSSpline("fRhoSSSpline", fIntegrals.GetRhoSSInt());
	fIntegrals.SetRhoSSpline(&rhoSSpline);
	fIntegrals.SetRhoSSSpline(&rhoSSSpline);
	return;
}

void WCSimIntegralLookupMaker3D::MakeSplines() {
	MakeRhoSplines();
	MakeRhoGSplines();
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

    fIntegrals.Verify();
    toSave->cd();
    fIntegrals.Write();
    toSave->Close();


}

void WCSimIntegralLookupMaker3D::SmoothLookupTables()
{
    // SmoothRhoTables();
    SmoothRhoGTables();
}

void WCSimIntegralLookupMaker3D::SmoothRhoTables()
{
    // Replace the points in the rho graph with smoothed ones
	TGraph * gRho   = fIntegrals.GetRhoInt();
    ReplaceGraphPoints(gRho, SmoothGraph(gRho));

    // Replace the points in the rhoS graph with smoothed ones
	TGraph * gRhoS  = fIntegrals.GetRhoSInt();
    ReplaceGraphPoints(gRhoS, SmoothGraph(gRhoS));

    // Replace the points in the rhoSS graph with smoothed ones
	TGraph * gRhoSS = fIntegrals.GetRhoSSInt();
    ReplaceGraphPoints(gRhoSS, SmoothGraph(gRhoSS));

    return;
}

void WCSimIntegralLookupMaker3D::ReplaceGraphPoints(TGraph * original, const TGraph& replacement)
{
    if(original->GetN() == 0){ return; }
    // Clear the original graph
    original->Set(0);

    // Copy across the points from the replacement graph
    for(int i = 0; i < replacement.GetN(); ++i)
    {
        original->SetPoint(original->GetN(), replacement.GetX()[i], replacement.GetY()[i]);
    }
    return;
}

void WCSimIntegralLookupMaker3D::SmoothRhoGTables()
{
    //TGraph * grTemp = fIntegrals.GetRhoGInt(1700,0.72);
    //TCanvas * canTmp = new TCanvas("canTmp","",800,600);
    //grTemp->Draw("ALP");
    //canTmp->SaveAs("atStart_1700.000000_0.720000.png");
    //delete canTmp;


    // Make some TAxis objects to get R0 and cosTh0 values easily
    TAxis axisR0(fNR0Bins, fR0Min, fR0Max);
    TAxis axisCosTh0(fNCosTh0Bins, fCosTh0Min, fCosTh0Max);

    int lastPercent = -1; // For keeping track of the % processed
    int nBins = fNR0Bins * fNCosTh0Bins; // Total # of bins for tracking % processed

    // Loop over R0 and cosTh0 and smooth each set of three integral graphs
	for( Int_t iR0Bin = 1; iR0Bin <= fNR0Bins; ++iR0Bin)
	{
		Double_t R0 = axisR0.GetBinLowEdge(iR0Bin);

		for( Int_t iCosTh0Bin = 1; iCosTh0Bin <= fNCosTh0Bins; ++iCosTh0Bin)
		{
			Double_t cosTh0 = axisCosTh0.GetBinLowEdge(iCosTh0Bin);

            // Print out the percentage completed
            int percent = (iCosTh0Bin-1) * (iR0Bin-1) * 100 / nBins;
            if(percent > lastPercent)
            {
                std::cout << "\rSmoothing " << percent << "%" << std::flush;
                lastPercent = percent;
            }

            // Smooth the RhoG integral graph
			TGraph * gr  = fIntegrals.GetRhoGInt(R0, cosTh0);
            ReplaceGraphPoints(gr, SmoothGraph(gr, 3));

            // Smooth the RhoGS integral graph
			TGraph * grS  = fIntegrals.GetRhoGInt(R0, cosTh0);
            grS->Draw("ALP");
            ReplaceGraphPoints(grS, SmoothGraph(grS, 3));
            
            // Smooth the RhoGSS integral graph
			TGraph * grSS  = fIntegrals.GetRhoGInt(R0, cosTh0);
            ReplaceGraphPoints(grSS, SmoothGraph(grSS, 3));
        }
    }

    return;
}
    
TGraph WCSimIntegralLookupMaker3D::SmoothGraph(TGraph * graph, int nTimes)
{
    
    // This is the graph we're going to return
    TGraph grSmooth;
    grSmooth.SetName(TString::Format("%s_smoothed%d", graph->GetName(), nTimes).Data());
    grSmooth.SetTitle(grSmooth.GetName());

    // Converting it into a histogram relies on there being at least 2 points
    // because we set the upper bin edge equal to the width of the bin below
    if(graph->GetN() > 1)
    {

        // We need to contruct a histogram out of the graph points
        // So first we'll get the x-coordinates and make them into an
        // array of ascending bin edges
        std::vector<double> xVec;

        // Get the graphs x-coordinates 
        double x, y;
        

        // Copy these n coordinates into our n+1 sized array
        for(int i = 0; i < graph->GetN(); ++i)
        {
            graph->GetPoint(i, x, y);
            xVec.push_back(x);
        }

        // Ensure they're ordered, then set the last point to have the same bin width as the second-last bin
        std::sort(xVec.begin(), xVec.end());
        xVec.push_back(2*xVec.at(xVec.size()-1) - xVec.at(xVec.size()-2));

        // Now make the histogram and fill it with the y-values of the graph
        TH1D * hist = new TH1D("hist","",graph->GetN(), &xVec[0]);
        for(int i = 0; i < graph->GetN(); ++i)
        {
            graph->GetPoint(i, x, y);
            hist->Fill(x, y);
        }

        // Smooth the histogram
        hist->Smooth(nTimes);

        // Now fill the graph to return with the smoothed points
        for(int i = 0; i < hist->GetNbinsX(); ++i)
        {
            grSmooth.SetPoint(i, xVec.at(i), hist->GetBinContent(i+1));
        }

        //TCanvas * can = new TCanvas(TString::Format("can_%s", graph->GetName()).Data(),"", 1280, 800);
        //TGraph * graphClone = (TGraph*)(graph->Clone());
        //graphClone->SetMarkerColor(kBlue);
        //graphClone->SetMarkerSize(0.7);
        //graphClone->SetLineColor(kBlue);
        //hist->SetLineColor(kRed);
        //hist->Draw();
        //graphClone->Draw("LP SAME");
        //grSmooth.SetMarkerColor(kMagenta);
        //grSmooth.SetMarkerSize(0.5);
        //grSmooth.SetLineColor(kMagenta);
        //grSmooth.Draw("LP SAME");

        //can->SaveAs(TString::Format("%s.png", can->GetName()).Data());
        //can->SaveAs(TString::Format("%s.root", can->GetName()).Data());
        //delete can;
        //delete graphClone;

        // Tidy up
        delete hist;
    }
    else
    {
        for(int i = 0; i < graph->GetN(); ++i)
        {
            grSmooth.SetPoint(i, graph->GetX()[i], graph->GetY()[i]);
        }
    }
    // And we're done
    return grSmooth;
}
