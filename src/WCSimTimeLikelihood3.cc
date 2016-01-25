/*
 * WCSimTimeLikelihood3.cc
 *
 *  Created on: 8 Jul 2015
 *      Author: andy
 */


#include "WCSimAnalysisConfig.hh"
#include "WCSimEmissionProfileManager.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimPMTConfig.hh"
#include "WCSimPMTManager.hh"
#include "WCSimRootGeom.hh"
#include "WCSimTimeLikelihood3.hh"
#include "WCSimTimePredictor.hh"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TString.h"
#include <cmath>
#include <map>

#include "WCSimInterface.hh"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimTimeLikelihood3)
#endif
double ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons)
{
    double prob = 0.0;
    if(nPhotons == 0 || sigma == 0) { return 0; }
    if(nPhotons < 25)
    {
	    prob =   nPhotons * TMath::Sqrt2() / (pow(2, nPhotons) * sigma * sqrt(M_PI))
		 	  * pow(1.0 - TMath::Erf((t - mean)/(TMath::Sqrt2()*sigma)), nPhotons-1)
		   	  * exp(-(t - mean)*(t - mean) / (2 * sigma * sigma));
    }
    else // Use logarithms to take care of the two large n powers almost cancelling
    {
        double log =  TMath::Log(nPhotons * TMath::Sqrt2() / (sigma * TMath::Sqrt(TMath::Pi())))
                    - nPhotons * TMath::Log(2)
                    + (nPhotons - 1) * TMath::Log(1.0 - TMath::Erf((t - mean) / (2 * sigma * sigma)))
                    - (t - mean) * (t - mean) / (2 * sigma * sigma);
        prob = TMath::Exp(log);

    }
	return prob;
}

double ConvertMeanTimeToMeanFirstTime(double * x, double * par)
{
	double nPhotons = par[0];
	double mean = par[1];
	double sigma = par[2];
	double t = x[0];
	return ConvertMeanTimeToMeanFirstTime(t, mean, sigma, nPhotons);
}

double IntegrateMeanTimeToMeanFirstTime( const double &t, const double &mean, const double &sigma, const double &nPhotons)
{
  double integral(0.0);
  if(nPhotons <= 25)
  {
    integral = 1.0 - 1.0/(pow(2,nPhotons)) * TMath::Power((1.0 - TMath::Erf((t - mean)/(TMath::Sqrt2() * sigma))), nPhotons);
  }
  else
  {
    double exponent = nPhotons * (TMath::Log(1 - TMath::Erf((t - mean) / (TMath::Sqrt2() * sigma))) - TMath::Log(2.0));
    integral = 1.0 - TMath::Exp(exponent);
  }
  return integral;
}

double IntegrateMeanTimeToMeanFirstTime(double * x, double * par)
{
  double &nPhotons = (par[0]);
  double &mean = (par[1]);
  double &sigma = (par[2]);
  double &t = (x[0]);

  return IntegrateMeanTimeToMeanFirstTime(t, mean, sigma, nPhotons);
}

double MinimumOfArrivalAndResolution(const double &t, const double &meanArr, const double rmsArr, const double pmtTime, const double &pmtRes, const double &nPhotons)
{
  double arrivalPDF = ConvertMeanTimeToMeanFirstTime(t, meanArr, rmsArr, nPhotons);
  double resolutionPDF = 1.0 / (pmtRes * TMath::Sqrt(2*TMath::Pi())) * TMath::Exp( - (t - pmtTime)*(t - pmtTime) / (2 * pmtRes * pmtRes));
  if( TMath::IsNaN(arrivalPDF) || TMath::IsNaN(resolutionPDF))
  {

    std::cout << "Arrival is " << arrivalPDF << " and resolution PDF = " << resolutionPDF << std::endl;
    std::cout << "t = " << t << " and array mean = " << meanArr << " and array RMS = " << rmsArr << " and PMT hit at " << pmtTime << " and PMT res = " << pmtRes << " and nPhotons = " << nPhotons << std::endl;
    assert(0);
  }
    
  if( resolutionPDF < arrivalPDF)
  {
    return resolutionPDF;
  }
  return arrivalPDF;
}

double MinimumOfArrivalAndResolution( double * x, double * par )
{
  double &nPhotons = (par[0]);
  double &meanArr = (par[1]);
  double &rmsArr = (par[2]);
  double &pmtTime = (par[3]);
  double &pmtRes = (par[4]);
  double &t = (x[0]);

  return MinimumOfArrivalAndResolution(t, meanArr, rmsArr, pmtTime, pmtRes, nPhotons);
}

WCSimTimeLikelihood3::WCSimTimeLikelihood3() : fSpeedOfParticle(1.0 * TMath::C() * 1e-7){
	// TODO Auto-generated constructor stub
	fLikelihoodDigitArray = 0x0;
  fEmissionProfileManager = 0x0;
	fPMTManager = new WCSimPMTManager();
}

WCSimTimeLikelihood3::WCSimTimeLikelihood3( WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager)
{
  fLikelihoodDigitArray = myDigitArray;
  fEmissionProfileManager = myEmissionProfileManager;
  fPMTManager = new WCSimPMTManager();
	fDistVsPred = new TGraphErrors();
	fDistVsPred->SetName("fDistVsPred");
	fDistVsProb = new TGraph();
	fDistVsProb->SetName("fDistVsProb");
	fDistVsSpreadProb = new TGraph();
	fDistVsSpreadProb->SetName("fDistVsSpreadProb");
	fTMinusTPredAll = new TGraphAsymmErrors();
	fTMinusTPredSource = new TGraphAsymmErrors();
	fTMinusTPredAll->SetName("fTMinusTPredAll");
	fTMinusTPredSource->SetName("fTMinusTPredSource");
  fTMinusTPredVsQ = new TGraphErrors();
  fTMinusTPredVsQ->SetName("fTMinusTPredVsQ");
  TString histName = TString::Format("fTMinusTPredHist_%d", WCSimInterface::GetEventNumber());
  fTMinusTPredHist = new TH1D(histName.Data(),"All PMTs;t_{hit} - t_{pred} (ns);PMTs",500,-50,50);
  fTMinusTPredSharpHist = new TH1D("fTMinusTPredSharpHist","All PMTs with #sigma < 2.0ns;t_{hit} - t_{pred} (ns);PMTs",500,-50,50);
  fHitQVsRMS = new TH2D("fHitQVsRMS", ";ceil(Total number of photons);Predicted RMS arrival time",20,1,21,20,0,10);

  return;
}

WCSimTimeLikelihood3::~WCSimTimeLikelihood3() {
	// TODO Auto-generated destructor stub
	if(fPMTManager != 0x0) { delete fPMTManager; }
	// I don't new the likelihooddigitarray so I don't delete it
}

void WCSimTimeLikelihood3::SetTracks(
		std::vector<WCSimLikelihoodTrackBase*>& myTracks) {
  //std::cout << "Setting tracks" << std::endl;
	fTracks = myTracks;
  fAllPreds.clear();
  fAllPreds.resize(fLikelihoodDigitArray->GetNDigits()); 
  //std::cout << "Set tracks - there are " << fTracks.size() << " of them" << std::endl;
}

void WCSimTimeLikelihood3::ResetTracks() {
	//std::cout << "Reset tracks!" << std::endl;
  fTracks.clear();
  fAllPreds.clear();
}

void WCSimTimeLikelihood3::ClearTracks() {
  ResetTracks();
}

double WCSimTimeLikelihood3::Calc2LnL(const unsigned int& iDigit) {


	double lnL = 0.0;

	WCSimLikelihoodDigit * myDigit = fLikelihoodDigitArray->GetDigit(iDigit);
	WCSimLikelihoodTrackBase * myTrack = fTracks.at(0);
	if(IsGoodDigit(myDigit))
	{
    // std::cout << "Digit " << myDigit << " is good!" << std::endl;
		std::pair<double, double> arrivalMeanSigma = GetArrivalTimeMeanSigma(myTrack, myDigit);
    // std::cout << "Arrival time: mean = " << arrivalMeanSigma.first << " and rms = " << arrivalMeanSigma.second << " and t = " << myDigit->GetT() << " and nPhotons = " << myDigit->GetQ() << std::endl;

		double reso = GetPMTTimeResolution(myDigit);
		double pScatter = 0.01;

		double sigmaSquared = arrivalMeanSigma.second * arrivalMeanSigma.second;
		double& sigma = arrivalMeanSigma.second;
        double& mean = arrivalMeanSigma.first;
		// std::cout << "sigma = " << sigma << std::endl;
		double nPhotons = myDigit->GetQ();

		double timeLikelihood = 0.0;
		if(arrivalMeanSigma.second > 0) // Floating point arithmetic can struggle with (1-erf)^{n-} / 2^{n} for large n
		{
            TF1 overlapFunc("overlapFunc", MinimumOfArrivalAndResolution, myDigit->GetT() - 5.0 * reso, myDigit->GetT() + 5.0 * reso, 5);
            overlapFunc.SetParameters(nPhotons, mean, sigma, myDigit->GetT(), reso);
            timeLikelihood = overlapFunc.Integral(myDigit->GetT() - 5.0 * reso, myDigit->GetT() + 5.0 * reso);
        }
		
        double likelihood =   (1. - pScatter) * (timeLikelihood) + pScatter;
		if(likelihood <= 0 || likelihood > 1)
		{
		  //std::cout << "iDigt = " << iDigit << " and likelihood = " << likelihood << std::endl;
		  assert(likelihood > 0 && likelihood <= 1);
		}
			lnL += log(likelihood);
		if(TMath::IsNaN(lnL))
		{
		  std::cout << "Digit is " << iDigit << " and likelihood is " << likelihood << " because arrives at " << myDigit->GetT() << " and predict " << arrivalMeanSigma.first << " with width " << arrivalMeanSigma.second << " and res " << reso << std::endl;
		  assert(!TMath::IsNaN(lnL));
		}
		double distance = (myDigit->GetPos() - myTrack->GetVtx()).Mag();
        double sqrtPi = 1.772453851;
        double start = arrivalMeanSigma.first - 10 * arrivalMeanSigma.second;
        double end   = arrivalMeanSigma.first + 10 * arrivalMeanSigma.second;

        TF1 * fMean = new TF1("fMean", ConvertMeanTimeToMeanFirstTime, start, end, 3);
       // fMean->SetParameters(nPhotons, arrivalMeanSigma.first, arrivalMeanSigma.second);
       // fMean->SetNpx(1000);
        TF1 * fInt = new TF1("fInt", IntegrateMeanTimeToMeanFirstTime, start, end, 3);
       //  fInt->SetParameters(nPhotons, arrivalMeanSigma.first, arrivalMeanSigma.second);
       // fInt->SetNpx(1000);
       //  TF1 * fOverlap = new TF1("fOverlap", MinimumOfArrivalAndResolution, start, end, 5);
       //  fOverlap->SetParameters(nPhotons, mean, sigma, myDigit->GetT(), reso);
       //  fOverlap->SetNpx(1000);
       //  TF1 * fPMTTimeRes = new TF1("fPMTTimeRes","1.0 / ([1] * TMath::Sqrt(2 * TMath::Pi())) * TMath::Exp(-(x - [0])*(x-[0])/(2 * [1]*[1]))",start, end);
       //  fPMTTimeRes->SetParameters(myDigit->GetT(), reso);
       //  TCanvas* can11 = new TCanvas("can11","",800,600);
       //  TH1D hTemplate("hTemplate","Calculating the time likelihood;Arrival time (ns);P(t)",100, start, end);
       //  hTemplate.SetMaximum(1.2 * fMean->GetMaximum());
       //  TLegend leg1(0.55, 0.85, 0.8, 0.6);
       //  leg1.AddEntry(fMean,"PDF for arrival time","L");
       //  leg1.AddEntry(fPMTTimeRes,"PMT time measurement","L");
       //  leg1.AddEntry(fOverlap,"Overlap","FL");
        
        // hTemplate.Draw("AXIS");
        // hTemplate.GetXaxis()->CenterTitle();
        // hTemplate.GetYaxis()->CenterTitle();
        // hTemplate.SetStats(0);
        // fMean->Draw("SAME");
        // leg1.Draw();
        // fPMTTimeRes->Draw("SAME");
        // fOverlap->Draw("SAME");
        // fMean->SetLineColor(kRed);
        // fPMTTimeRes->SetLineColor(kBlue);
        // fOverlap->SetLineColor(kViolet-1);
        // fOverlap->SetFillColor(kViolet-1);
        // fOverlap->SetFillStyle(3002);
        // TString name = TString::Format("overlap_%d", myDigit->GetTubeId());
        // can11->SaveAs((name + ".png").Data());
        // can11->SaveAs((name + ".pdf").Data());
        // can11->SaveAs((name + ".C").Data());
        // can11->SaveAs((name + ".root").Data());
        //  delete can11;
        /*
        
        double predMean = fMean->GetX(fMean->Mean(start, end, fMean->GetParameters()));
        double errLow = predMean - fInt->GetX(0.32);
        double errHi = fInt->GetX(0.68) - predMean;
        //std::cout << "nPhotons = " << nPhotons << "  ArrivalMeanSigma = (" << arrivalMeanSigma.first << ", " << arrivalMeanSigma.second << ")   PredMean = " << predMean << " and errors are " << fInt->GetX(0.32) << "  " << fInt->GetX(0.68) << std::endl;
        double delta = myDigit->GetT() - predMean;
        if( delta > 50 ) { delta = 50.0; }
        fTMinusTPredSource->SetPoint(fTMinusTPredSource->GetN(), distance, delta);
        fTMinusTPredSource->SetPointError(fTMinusTPredSource->GetN()-1, 0, 0, errLow, errHi);



		fDistVsPred->SetPoint(fDistVsPred->GetN(), distance, arrivalMeanSigma.first);
		fDistVsPred->SetPointError(fDistVsPred->GetN()-1, 0, arrivalMeanSigma.second);

		fDistVsProb->SetPoint(fDistVsProb->GetN(), distance, -2.0*lnL);
		fDistVsSpreadProb->SetPoint(fDistVsSpreadProb->GetN(), distance, -2.0*lnL);
        //fTMinusTPredSource->SetPoint(fTMinusTPredSource->GetN(), distance, (myDigit->GetT()-arrivalMeanSigma.first));
		//fTMinusTPredSource->SetPointError(fTMinusTPredSource->GetN()-1, 0, arrivalMeanSigma.second);
		fTMinusTPredAll->SetPoint(fTMinusTPredAll->GetN(), distance, delta);
		fTMinusTPredAll->SetPointError(fTMinusTPredAll->GetN()-1, 0, 0, sqrt(reso * reso + errLow*errLow), sqrt(reso * reso  +errHi * errHi));
		fTMinusTPredVsQ->SetPoint(fTMinusTPredVsQ->GetN(), myDigit->GetQ(), myDigit->GetT()-arrivalMeanSigma.first);
		fTMinusTPredVsQ->SetPointError(fTMinusTPredVsQ->GetN()-1, 0, sqrt(sigmaSquared));
		fHitQVsRMS->Fill(ceil(myDigit->GetQ()), arrivalMeanSigma.second);
        fTMinusTPredHist->Fill(delta);
        if( errLow < 1.0 ){
          fTMinusTPredSharpHist->Fill(delta);
        }
        */
        delete fInt;
        delete fMean;
        // delete fOverlap;
        // delete fPMTTimeRes;

	}
	return -2.0 * lnL;
}

double WCSimTimeLikelihood3::Calc2LnL() {
	double minusTwoLnL = 0.0;
  std::cout << "Calling ::Calc2LnL" << std::endl;
	
  for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit)
	{
		minusTwoLnL += Calc2LnL(iDigit);
	}
	return minusTwoLnL;
}

double WCSimTimeLikelihood3::GetPMTTimeResolution(const unsigned int& iDigit) {
	WCSimLikelihoodDigit * myDigit = fLikelihoodDigitArray->GetDigit(iDigit);
	return GetPMTTimeResolution(myDigit);
}

double WCSimTimeLikelihood3::GetPMTTimeResolution(WCSimLikelihoodDigit * myDigit)
{
  TString pmtName = myDigit->GetPMTName();
  double timeConstant = 0.0;
  if( fPMTTimeConstantMap.find(pmtName) == fPMTTimeConstantMap.end())
  {
	  WCSimPMTConfig config = fPMTManager->GetPMTByName(std::string(pmtName.Data()));
	  timeConstant = config.GetTimeConstant();
    fPMTTimeConstantMap[pmtName] = timeConstant;
  }
  return 0.33 + sqrt(fPMTTimeConstantMap[pmtName] / myDigit->GetQ()); // This is what WCSim does...
}

bool WCSimTimeLikelihood3::IsGoodDigit(WCSimLikelihoodDigit * myDigit)
{
    return (myDigit->GetQ() > 0);
}

std::pair<double, double> WCSimTimeLikelihood3::GetArrivalTimeMeanSigma(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
  bool debug = (myDigit->GetX() > 1200 && myDigit->GetZ() > 700 && myDigit->GetZ() < 730);
  debug = false;

  TH2F * hSCosTheta = fEmissionProfileManager->GetEmissionProfile(myTrack)->GetSCosThetaForTime();
  TH1F * hS = fEmissionProfileManager->GetEmissionProfile(myTrack)->GetSForTime();
  double minCosTheta = fEmissionProfileManager->GetEmissionProfile(myTrack)->GetTimeCosThetaMin();
  double maxCosTheta = fEmissionProfileManager->GetEmissionProfile(myTrack)->GetTimeCosThetaMax();
  double n = myDigit->GetAverageRefIndex();
  // std::cout << "Loading hS and hSCosTheta " << std::endl;

  // We'll avoid doing lots of GetXaxis()->FindBin calls by calculating all the cosTheta values ahead of time
  // Then we can sort it and just loop through the array of cosTheta bins once

  // Elements are (cosTheta, sBin, dCosTheta, distance to PMT)
  double init[4] = {-999., -999., -999., -999.};
  std::vector< std::vector<double> > cosSAndDelta(hS->GetNbinsX(), std::vector<double>(init, init+4));
  
  TVector3 dir = myTrack->GetDir();
  TVector3 pmt = myDigit->GetPos();

  double dCos = -999;
  double magToPMT = 0.0;
  double cosTheta = 0.0;
  TVector3 toPMT;
  for( int iSBin = 1; iSBin <= hS->GetNbinsX(); ++iSBin )
  {
	  double s = hS->GetXaxis()->GetBinLowEdge(iSBin);
	  if( dCos == -999)
	  {
		  toPMT = pmt - myTrack->GetPropagatedPos(s);
		  magToPMT = toPMT.Mag();
		  cosTheta = dir.Dot(toPMT)/magToPMT;
	  }

	  // Update it to work out the width in cosTheta that we traverse
	  // Then remember it because it'll be the starting point for the next loop iteration
	  toPMT = pmt - myTrack->GetPropagatedPos(hS->GetXaxis()->GetBinUpEdge(iSBin));
	  magToPMT = toPMT.Mag();
	  double newCos = dir.Dot(toPMT) / magToPMT;
	  if(! (cosTheta < minCosTheta || cosTheta >= maxCosTheta) )
	  {
		  cosSAndDelta[iSBin-1][0] = cosTheta;
		  cosSAndDelta[iSBin-1][1] = iSBin;
		  cosSAndDelta[iSBin-1][2] = fabs(newCos - cosTheta);
		  cosSAndDelta[iSBin-1][3] = magToPMT;
	  }
	  cosTheta = newCos;
  }
  // std::cout << "Sorting " << cosSAndDelta.size() << " entries" << " where first is " << cosSAndDelta[0][0] << std::endl;
  std::sort(cosSAndDelta.begin(), cosSAndDelta.end());
  // std::cout << "Now first is " << cosSAndDelta[0][0] << std::endl;

  // Now do the loop over cosTheta bins and look up the emission profile histograms to work out
  // the photon arrival time at the PMT and a probability weightinh
  int lastBin = 1;
  TAxis * axis = hSCosTheta->GetXaxis();

  double weightedSum = 0.0;
  double sumOfWeights = 0.0;

  double speedOfLightInCmPerNs = TMath::C() * 1e-7;
   // std::cout << "n = " << n << " so speed of light = " << speedOfLightInCmPerNs/n << std::endl;
  if(WCSimAnalysisConfig::Instance()->GetUseCustomParticleSpeed())
  {
	  fSpeedOfParticle = WCSimAnalysisConfig::Instance()->GetCustomParticleSpeed() * speedOfLightInCmPerNs;
  }
  else
  {
	  fSpeedOfParticle = myTrack->GetPropagationSpeedFrac() * speedOfLightInCmPerNs;
  }

  if(WCSimAnalysisConfig::Instance()->GetUseCustomSpeedOfLight())
  {
    n = 1.0 / WCSimAnalysisConfig::Instance()->GetCustomSpeedOfLight();
  }

  std::vector<double> times(cosSAndDelta.size());
  std::vector<double> probs(cosSAndDelta.size());
  int index = 0;
  for(std::vector<std::vector<double> >::iterator thItr = cosSAndDelta.begin(); thItr != cosSAndDelta.end(); ++thItr)
  {
	  // Had cosTheta outside the allowed range
	  if(thItr->at(0) == -999){ continue; }

	  for(int thBin = lastBin; thBin <= hSCosTheta->GetNbinsX(); ++thBin)
	  {
		  if(axis->GetBinLowEdge(thBin) <= thItr->at(0) && thItr->at(0) < axis->GetBinUpEdge(thBin))
		  {
			  lastBin = thBin;

			  double t =   myTrack->GetT()
					     + hS->GetBinLowEdge(thItr->at(1)) / fSpeedOfParticle
					     + n * thItr->at(3) / speedOfLightInCmPerNs; // Make c into cm/ns
  //      std::cout << "Digit " << (int)myDigit->GetTubeId() << "  Track t = " << myTrack->GetT() << " and pred t  = " << t << std::endl;
			  double prob =   hS->GetBinContent(thItr->at(1))
							* hS->GetBinWidth(thItr->at(1))
					        * hSCosTheta->GetBinContent(thBin, thItr->at(1))
							* thItr->at(2);
			  weightedSum += prob * t;
			  sumOfWeights += prob;

			  times[index] = t;
			  probs[index] = prob;
			  if(debug)
			  {
				  std::cout << "s = " << hS->GetBinLowEdge(thItr->at(1)) << " and cosTheta = " << thItr->at(0)
					    	<< " and t0 = " << myTrack->GetT()
							<< " and particle time = " << hS->GetBinLowEdge(thItr->at(1)) / fSpeedOfParticle
							<< " and photon time = " << n * thItr->at(3) / speedOfLightInCmPerNs
							<< " so t = " << t << " = " << myTrack->GetT() + (hS->GetBinLowEdge(thItr->at(1)) / fSpeedOfParticle) + n * thItr->at(3) / speedOfLightInCmPerNs << std::endl
							<< "Prob = " << prob << " of which s->" << hS->GetBinContent(thItr->at(1))
							<< " cosTheta -> " << hSCosTheta->GetBinContent(thBin, thItr->at(1))
							<< " and widths are " << hS->GetBinWidth(thItr->at(1)) << " x " << thItr->at(2) << std::endl << std::endl;
			  }
			  index++;
		  }
	  }
  }



  // Finally work out the mean and rms of the predicted times
  double mean = 0;
  double rms = 0;
  //std::cout << "Sum of weights = " << sumOfWeights << " and sum of weighted times = " << weightedSum << std::endl;
  if(sumOfWeights > 0)
  {
	  mean = weightedSum / sumOfWeights;
    for(int i = 0; i < index ; ++i)
    {
	  rms += (mean - times[i])*(mean-times[i])*probs[i] / sumOfWeights;
    }
    rms = sqrt(rms);
  }
  else
  {
    // std::cerr << "Somehow I don't think PMT number " << myDigit->GetTubeId() << "could be hit" << std::endl;
  }


  // printf("So mean = %f and RMS = %f at (%.01f, %.01f, %.01f) for digit %d hit at %f", mean, rms, myDigit->GetX(), myDigit->GetY(), myDigit->GetZ(), myDigit->GetTubeId(), myDigit->GetT());
  // std::cout << "So mean = " << mean << " and RMS = " << rms << " at " << myDigit->GetX() << "," << myDigit->GetY() << ", " << myDigit->GetZ() << std::endl;
  fAllPreds[myDigit->GetTubeId()] = mean;

  if(debug && false)
  {
    TGraph * gr = new TGraph();
    gr->SetName("Arrival times");
    for(int i = 0; i < index; ++i)
    {
      gr->SetPoint(gr->GetN(), times[i], probs[i]);
    }
    TCanvas * can = new TCanvas();
    gr->Draw("AP");
    gr->GetXaxis()->SetTitle("Hit time (ns)");
    gr->GetYaxis()->SetTitle("Relative Probability");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kBlue);
    can->SaveAs(TString::Format("fDigitHitProbabilities_%d.png", myDigit->GetTubeId()).Data());
    can->SaveAs(TString::Format("fDigitHitProbabilities_%d.C", myDigit->GetTubeId()).Data());
    can->SaveAs("fDigitHitProbabilities.C");
    delete can;
    delete gr;
    // assert(0);

    
  }
  return std::make_pair(mean, rms);
}

std::vector<double> WCSimTimeLikelihood3::GetAllPredictedTimes()
{
	return fAllPreds;
}

void WCSimTimeLikelihood3::SavePlot()
{
    return;
	TCanvas * can = new TCanvas("can","",800,600);
	fDistVsPred->Draw("AP");
	fDistVsPred->GetXaxis()->SetTitle("Distance from vertex");
	fDistVsPred->GetYaxis()->SetTitle("Predicted arrival time");
	fDistVsPred->SetMarkerStyle(20);
	fDistVsPred->SetMarkerSize(0.8);
	fDistVsPred->SetMarkerColor(kRed);
	can->SaveAs("fDistVsPred.png");
	can->SaveAs("fDistVsPred.C");
	can->SaveAs("fDistVsPred.root");
	fDistVsProb->Draw("AP");
	fDistVsProb->GetXaxis()->SetTitle("Distance from vertex");
	fDistVsProb->GetYaxis()->SetTitle("-2Ln(L)");
	fDistVsProb->SetMarkerStyle(20);
	fDistVsProb->SetMarkerSize(0.8);
	fDistVsProb->SetMarkerColor(kRed);
	can->SaveAs("fDistVsProb.png");
	can->SaveAs("fDistVsProb.C");
	can->SaveAs("fDistVsProb.root");
	fDistVsSpreadProb->Draw("AP");
	fDistVsSpreadProb->GetXaxis()->SetTitle("Distance from vertex");
	fDistVsSpreadProb->GetYaxis()->SetTitle("-2Ln(L)");
	fDistVsSpreadProb->SetMarkerStyle(20);
	fDistVsSpreadProb->SetMarkerSize(0.8);
	fDistVsSpreadProb->SetMarkerColor(kRed);
	can->SaveAs("fDistVsSpreadProb.png");
	can->SaveAs("fDistVsSpreadProb.C");
	can->SaveAs("fDistVsSpreadProb.root");
  TH1D hTemplate("hTemplate",";Distance from vertex (cm);T_{hit} - T_{first, pred} (ns)", 250, 0, 2500);
  hTemplate.SetStats(0);
  hTemplate.Draw();
  hTemplate.SetMinimum(-50.0);
  hTemplate.SetMaximum(50.0);
	fTMinusTPredAll->Draw("P SAME");
	fTMinusTPredSource->Draw("|| SAME");
	fTMinusTPredAll->GetXaxis()->SetTitle("Distance from vertex");
	fTMinusTPredAll->GetYaxis()->SetTitle("T_{hit} - T_{first, pred} (ns)");
	fTMinusTPredAll->SetMarkerStyle(20);
	fTMinusTPredAll->SetMarkerSize(0.6);
	fTMinusTPredAll->SetMarkerColor(kRed);
	fTMinusTPredAll->SetLineColor(kBlack);
	fTMinusTPredSource->GetXaxis()->SetTitle("Distance from vertex");
	fTMinusTPredSource->GetYaxis()->SetTitle("T_{hit} - T_{first, pred} (ns)");
	fTMinusTPredSource->SetMarkerStyle(20);
	fTMinusTPredSource->SetMarkerSize(0.6);
	fTMinusTPredSource->SetMarkerColor(kRed+1);
	fTMinusTPredSource->SetLineColor(kRed+1);
	TString name = Form("fTMinusTPred_%.02fc", fSpeedOfParticle / (1e-7 * TMath::C()));
	can->SaveAs((name+TString(".png")).Data());
	can->SaveAs((name+TString(".C")).Data());
	can->SaveAs((name+TString(".root")).Data());
	fTMinusTPredVsQ->Draw("AP");
	fTMinusTPredVsQ->GetXaxis()->SetTitle("Hit charge (p.e.)");
	fTMinusTPredVsQ->GetYaxis()->SetTitle("t_{hit} - t_{pred} (ns)");
	fTMinusTPredVsQ->SetMarkerStyle(20);
	fTMinusTPredVsQ->SetMarkerSize(0.8);
	fTMinusTPredVsQ->SetMarkerColor(kCyan+2);
  fTMinusTPredVsQ->Draw("AP");
  can->SaveAs("fTMinusTPRedVsQ.png");
  can->SaveAs("fTMinusTPRedVsQ.C");
  can->SaveAs("fTMinusTPRedVsQ.root");
  int event = WCSimInterface::GetEventNumber();
  fHitQVsRMS->Draw("BOX");
  fHitQVsRMS->SetDirectory(0);
  TDirectory * tmpDir = gDirectory;
  TFile f(TString::Format("fHitQVsRMSFile_%d.root", event).Data(), "RECREATE");
  fHitQVsRMS->Write();
  f.Close();
  tmpDir->cd();
  can->SaveAs(TString::Format("fHitQVsRMS_%d.png", event).Data());
  can->SaveAs(TString::Format("fHitQVsRMS_%d.C", event).Data());
  can->SaveAs(TString::Format("fHitQVsRMS_%d.root", event).Data());

  fTMinusTPredHist->SetDirectory(0);
  TFile f2("fTMinusTPredHistFile.root","UPDATE");
  fTMinusTPredHist->Write();
  f2.Close();
  tmpDir->cd();
  fTMinusTPredHist->Draw();
  fTMinusTPredHist->SetLineColor(kRed+1);
  fTMinusTPredHist->SetLineWidth(2);
  fTMinusTPredHist->SetFillColor(kRed);
  fTMinusTPredHist->SetFillStyle(3001);
  can->SaveAs("fTMinusTPredHist.png");
  can->SaveAs("fTMinusTPredHist.C");
  can->SaveAs("fTMinusTPredHist.root");
  fTMinusTPredSharpHist->Draw();
  fTMinusTPredSharpHist->SetLineColor(kAzure-1);
  fTMinusTPredSharpHist->SetLineWidth(2);
  fTMinusTPredSharpHist->SetFillColor(kAzure);
  fTMinusTPredSharpHist->SetFillStyle(3001);

  can->SaveAs("fTMinusTPredSharpHist.png");
  can->SaveAs("fTMinusTPredSharpHist.C");
  can->SaveAs("fTMinusTPredSharpHist.root");
	delete can;
}

