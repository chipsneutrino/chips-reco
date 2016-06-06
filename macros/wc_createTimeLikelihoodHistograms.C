/*
 * wc_createTimeLikelihoodHistograms.C
 *
 *  Created on: 15 Dec 2015
 *      Author: ajperch
 */
#include "TSystem.h"
#include <cassert>
#include <utility>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TTree.h>
#include <TLine.h>
#include <TFile.h>
#include <iostream>

void SaveCosThetaProjection(TH1F * proj);
void SaveSCosThetaForTime(TH2F * hist);
std::pair<double, double> GetPercentiles(TH1F * hist, double percentile);

void wc_createTimeLikelihoodHistograms(const char * inputName)
{
	  // Load libraries
	  // ==============
	  gApplication->ProcessLine(".except");
	  gSystem->Load("libGeom");
	  gSystem->Load("libEve");
	  gSystem->Load("libMinuit");
	  gApplication->ProcessLine("#include <vector>");

	  TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
	  TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
	  gSystem->Load(libWCSimRoot.Data());
	  gSystem->Load(libWCSimAnalysis.Data());

	  TFile profileFile(inputName,"UPDATE");
	  if( profileFile.IsZombie())
	  {
		  std::cerr << "Error: couldn't open file " << inputName << std::endl;
		  return;
	  }

	  // Read the saved histograms of where photons are emitted from the summary files
	  TH2F * sThetaCoarseHist = 0x0;
	  TH2F * sThetaFineHist = 0x0;
	  TH1F * sHist = 0x0;

	  profileFile.GetObject("fS", sHist);
	  profileFile.GetObject("fSCosThetaCoarse", sThetaCoarseHist);
	  profileFile.GetObject("fSCosThetaFine", sThetaFineHist);
	  std::cout << "fS = " << sHist << " and sThetaCoarseHist = " << sThetaCoarseHist << " and sThetaFineHist = " << sThetaFineHist << std::endl;

	  // I assume this is constant in order to rebin.  If you start from WCSim rather than from (cosTheta, s) histograms
	  // for photons then you can get around this
	  double sBinWidth = sHist->GetXaxis()->GetBinWidth(1);

	  // Make a single histogram of cosTheta copying the binning of the coarse and fine histograms
	  std::vector<double> thetaBins;
	  for(int iBin = 1; iBin <= sThetaCoarseHist->GetNbinsX(); ++iBin)
	  {
		  thetaBins.push_back(sThetaCoarseHist->GetXaxis()->GetBinLowEdge(iBin));
	  }
	  for(int iBin = 1; iBin <= sThetaFineHist->GetNbinsX(); ++iBin)
	  {
		  thetaBins.push_back(sThetaFineHist->GetXaxis()->GetBinLowEdge(iBin));
	  }
	  thetaBins.push_back(sThetaFineHist->GetXaxis()->GetBinUpEdge(sThetaFineHist->GetNbinsX()));
	  TH1F * cosThetaProjection = new TH1F("cosThetaProjection",";cos#theta;s (cm)", thetaBins.size()-1, &(thetaBins[0]));

	  // Fill the cosTheta projection histogram:
	  for(int iSBin = 1 ; iSBin <= sThetaCoarseHist->GetNbinsY(); ++iSBin)
	  {
		  for( int jThBin = 1; jThBin <= sThetaCoarseHist->GetNbinsX(); ++jThBin)
		  {
			  cosThetaProjection->Fill(sThetaCoarseHist->GetXaxis()->GetBinLowEdge(jThBin),
					  	  	  	  	   sThetaCoarseHist->GetBinContent(jThBin, iSBin) / sThetaCoarseHist->GetXaxis()->GetBinWidth(jThBin));
		  }
		  for( int jThBin = 1; jThBin <= sThetaFineHist->GetNbinsX(); ++jThBin)
		  {
			  cosThetaProjection->Fill(sThetaFineHist->GetXaxis()->GetBinLowEdge(jThBin),
					  	  	  	  	   sThetaFineHist->GetBinContent(jThBin, iSBin) / sThetaFineHist->GetXaxis()->GetBinWidth(jThBin));
		  }
	  }
	  SaveCosThetaProjection(cosThetaProjection);


	  // Save the angles at which we need to cut off to get rid of the last 2.5% of tails
	  // on either side of the peak
	  double percentile = 2.5;
	  std::pair<double, double> cutoffs = GetPercentiles(cosThetaProjection, percentile);
	  TTree * angleCutoffTree = new TTree("angleCutoffTree","angleCutoffTree");
	  angleCutoffTree->Branch("timePercentile",&percentile);
	  angleCutoffTree->Branch("timeCosThetaMin", &cutoffs.first);
	  angleCutoffTree->Branch("timeCosThetaMax",&cutoffs.second);
	  angleCutoffTree->Fill();
	  // angleCutoffTree->Write("",TObject::kOverwrite);
    profileFile.WriteTObject(angleCutoffTree, "angleCutoffTree", "Overwrite");

	  // Now make some coarsely-binned plots of s vs cosTheta
	  double binWidth = sBinWidth * (floor(25.0/sBinWidth) + (((25.0 / sBinWidth) - floor(25.0 / sBinWidth)) != 0));
	  int nBins = (int)(sHist->GetXaxis()->GetXmax() / binWidth) + ((25.0/sBinWidth) - floor(25.0/sBinWidth) != 0);
	  double sMax = nBins * binWidth;
	  TH2F * sCosThetaForTime = new TH2F("hSCosThetaForTime","Fraction of photons;cos#theta;s (cm)",
			  	  	  	  	  	  	  	 thetaBins.size()-1, &(thetaBins[0]),
										 nBins, 0, sMax);

	  // Fill the coarse s-binned histogram with raw event numbers
	  TH1F * sForTime = new TH1F("hSForTime","Fraction of photons;s (cm)", nBins, 0, sMax);
	  for(int iSBin = 1; iSBin <= sThetaCoarseHist->GetNbinsY(); ++iSBin)
	  {
		  for( int jThBin = 1; jThBin <= sThetaCoarseHist->GetNbinsX(); ++jThBin)
		  {
			  sForTime->Fill(sThetaCoarseHist->GetYaxis()->GetBinLowEdge(iSBin), sThetaCoarseHist->GetBinContent(jThBin, iSBin));
		  }
		  for( int jThBin = 1; jThBin <= sThetaFineHist->GetNbinsX(); ++jThBin)
		  {
			  sForTime->Fill(sThetaFineHist->GetYaxis()->GetBinLowEdge(iSBin), sThetaFineHist->GetBinContent(jThBin, iSBin));
		  }
	  }

	  // Now normalise the s plot correctly: we want sum(bin contents * bin width) = 1
	  double sIntegral = sForTime->Integral();
	  for(int sBin = 1; sBin <= sForTime->GetNbinsX(); ++sBin)
	  {
		  sForTime->SetBinContent(sBin, sForTime->GetBinContent(sBin) / (sIntegral * sForTime->GetXaxis()->GetBinWidth(sBin)));
	  }

	  // Draw and save the s plot
	  TCanvas * can1 = new TCanvas("can1","",800,600);
	  sForTime->Draw();
	  // can1->SaveAs("sForTime.png");
	  // can1->SaveAs("sForTime.C");
	  delete can1;


	  // Now some  reshuffling to basically copy sThetaCoarse and Find into the new more-coarsely s-binned histogram
	  for(int iSBin = 1 ; iSBin <= sThetaCoarseHist->GetNbinsY(); ++iSBin)
	  {
		  double x, y;
		  y = sThetaCoarseHist->GetYaxis()->GetBinLowEdge(iSBin);
		  for( int jThBin = 1; jThBin <= sThetaCoarseHist->GetNbinsX(); ++jThBin)
		  {
			  x = sThetaCoarseHist->GetXaxis()->GetBinLowEdge(jThBin);
			  sCosThetaForTime->Fill(x, y,
					  	  	  	  	 sThetaCoarseHist->GetBinContent(jThBin, iSBin));///sThetaCoarseHist->GetXaxis()->GetBinWidth(jThBin));
		  }
		  for( int jThBin = 1; jThBin <= sThetaFineHist->GetNbinsX(); ++jThBin)
		  {
			  x = sThetaFineHist->GetXaxis()->GetBinLowEdge(jThBin);
			  sCosThetaForTime->Fill(x, y,
					  	  	    	 sThetaFineHist->GetBinContent(jThBin, iSBin));// /sThetaFineHist->GetXaxis()->GetBinWidth(jThBin));
		  }
	  }

	  for(int iSBin = 1; iSBin <= sCosThetaForTime->GetNbinsY(); ++iSBin )
	  {
		  double nPhotons = 0.0;
		  for( int jThBin = 1; jThBin <= sCosThetaForTime->GetNbinsX(); ++jThBin)
		  {
			  nPhotons += sCosThetaForTime->GetBinContent(jThBin, iSBin);
		  }
		  if(nPhotons > 0)
		  {
			  double normFactor = 1.0 / nPhotons;
			  for(int jThBin = 1; jThBin <= sCosThetaForTime->GetNbinsX(); ++jThBin)
			  {
				  sCosThetaForTime->SetBinContent(jThBin, iSBin,
						  	  	  	  	  	  	  sCosThetaForTime->GetBinContent(jThBin, iSBin) * normFactor / sCosThetaForTime->GetXaxis()->GetBinWidth(jThBin));
			  }
		  }
	  }

	  SaveSCosThetaForTime(sCosThetaForTime);

	  // sCosThetaForTime->Write("",TObject::kOverwrite);
    profileFile.WriteTObject(sCosThetaForTime, "hSCosThetaForTime", "Overwrite");
    // sForTime->Write("",TObject::kOverwrite);
    profileFile.WriteTObject(sForTime, "hSForTime", "Overwrite");
	  profileFile.Close();
}

void SaveCosThetaProjection(TH1F * proj)
{
	TCanvas * can = new TCanvas("projCan","",800,600);
	proj->Draw();
	proj->GetXaxis()->CenterTitle();
	proj->GetYaxis()->CenterTitle();
	proj->SetLineColor(kAzure-5);
	proj->SetLineWidth(2);
	proj->SetTitle("Projection of number of photons in cos#theta, normalised to 1");
	proj->Scale(1.0/proj->Integral());

	std::pair<double, double> xPercentiles = GetPercentiles(proj, 2.5);
	double maximum = proj->GetMaximum();

	TLine lineMin(xPercentiles.first, 0, xPercentiles.first, maximum);
	TLine lineMax(xPercentiles.second, 0, xPercentiles.second, maximum);
	lineMin.SetLineColor(kGray);
	lineMax.SetLineColor(kGray);
	lineMin.SetLineWidth(2);
	lineMax.SetLineWidth(2);
	lineMin.SetLineStyle(kDashed);
	lineMax.SetLineStyle(kDashed);

	lineMin.Draw();
	lineMax.Draw();

	// can->SaveAs("cosThetaProjection.png");
	// can->SaveAs("cosThetaProjection.C");
	// can->SaveAs("cosThetaProjection.root");

	TH1F * integral = (TH1F*)proj->Clone();
	integral->Reset();
	double runningTot = 0.0;
	for(int iBin = 1; iBin < proj->GetNbinsX(); ++iBin)
	{
		runningTot += proj->GetBinContent(iBin);
		integral->SetBinContent(iBin, runningTot);
	}
	integral->SetLineColor(kAzure-5);
	integral->Draw();
	// can->SaveAs("cosThetaProjIntegral.png");

	delete can;
}

std::pair<double, double> GetPercentiles(TH1F * hist, double percentile)
{
	double integral = hist->Integral("width");
	double runningTotal = 0.0;
	double min = 0.0, max = 0.0;
	for(int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin)
	{
		runningTotal += hist->GetBinContent(iBin) * hist->GetBinWidth(iBin);
		if( runningTotal > percentile / 100.0 * integral && min == 0 )
		{
			min = hist->GetXaxis()->GetBinLowEdge(iBin);
			std::cout << "Min = " << min << " runningTotal = " << runningTotal << std::endl;

		}
		if( runningTotal > (100 - percentile)/100.0 * integral && max == 0.0 )
		{
			max = hist->GetXaxis()->GetBinUpEdge(iBin);
			std::cout << "Max = " << max << " runningTotal = " << runningTotal << std::endl;
			break;
		}
	}
	return std::make_pair(min, max);
}

void SaveSCosThetaForTime(TH2F * hist)
{
	TCanvas * can = new TCanvas("projCan","",800,600);
	hist->Draw("COLZ");
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->CenterTitle();

	// can->SaveAs("sCosThetaForTime.png");
	// can->SaveAs("sCosThetaForTime.C");
	// can->SaveAs("sCosThetaForTime.root");
	delete can;
}
