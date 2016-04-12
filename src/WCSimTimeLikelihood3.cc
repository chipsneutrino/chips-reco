/*
 * WCSimTimeLikelihood3.cc
 *
 *  Created on: 8 Jul 2015
 *      Author: andy
 */
#include "WCSimAnalysisConfig.hh"
#include "WCSimEmissionProfileManager.hh"
#include "WCSimFastMath.hh"
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
#include <algorithm>
#include <cmath>
#include <map>

#include "WCSimInterface.hh"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TSpline.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimTimeLikelihood3)
#endif

const double WCSimTimeLikelihood3::fMaximumLnL = 25.0/2.0; // So the maximum of -2LnL = 25
const double WCSimTimeLikelihood3::fMinimumLikelihood = TMath::Exp(-fMaximumLnL); // So the maximum of -2LnL = 25
const double WCSimTimeLikelihood3::fSqrtPi = TMath::Sqrt(TMath::Pi());
const double WCSimTimeLikelihood3::fLog2 = TMath::Log(2.0);

double ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons)
{
    double prob = 0.0;
    if(nPhotons == 0 || sigma == 0) { return 0; }
    if(nPhotons < 25)
    {
	    prob =   nPhotons * TMath::Sqrt2() / (pow(2, nPhotons) * sigma * sqrt(M_PI))
		 	  * pow(1.0 - WCSimFastMath::erf((t - mean)/(TMath::Sqrt2()*sigma)), nPhotons-1)
		   	  * TMath::Exp(-(t - mean)*(t - mean) / (2 * sigma * sigma));
    }
    else // Use logarithms to take care of the two large n powers almost cancelling
    {
        double logL =  log(nPhotons * TMath::Sqrt2() / (sigma * WCSimTimeLikelihood3::fSqrtPi ))
                     - nPhotons * WCSimTimeLikelihood3::fLog2
                     + (nPhotons - 1) * log(1.0 - WCSimFastMath::erf((t - mean) / (2 * sigma * sigma)))
                     - (t - mean) * (t - mean) / (2 * sigma * sigma);
        prob = TMath::Exp(logL);

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
    integral = 1.0 - 1.0/(pow(2,nPhotons)) * pow((1.0 - WCSimFastMath::erf((t - mean)/(TMath::Sqrt2() * sigma))), nPhotons);
  }
  else
  {
    double exponent = nPhotons * (log(1 - WCSimFastMath::erf((t - mean) / (TMath::Sqrt2() * sigma))) - WCSimTimeLikelihood3::fLog2);
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
  double root2pi = TMath::Sqrt(2 * TMath::Pi());
  double arrivalPDF = ConvertMeanTimeToMeanFirstTime(t, meanArr, rmsArr, nPhotons);
  double resolutionPDF = 1.0 / (pmtRes * root2pi) * TMath::Exp( - (t - pmtTime)*(t - pmtTime) / (2 * pmtRes * pmtRes));
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
    fProbToCharge = new TGraph();
    fProbToCharge->SetName("fProbToCharge");
    fProbToCharge->SetMarkerStyle(20);
    fProbToCharge->SetMarkerSize(0.6);
    fProbToCharge->SetMarkerColor(kAzure);
	fDistVsPred = 0x0;
	fDistVsProb = 0x0;
	fDistVsSpreadProb = 0x0;
	fTMinusTPredAll = 0x0;
	fTMinusTPredSource = 0x0;
    fTMinusTPredVsQ = 0x0;
    fTMinusTPredHist = 0x0;
    fTMinusTPredSharpHist = 0x0;
    fHitQVsRMS = 0x0;
    fLastE = -999;
    fLastType = TrackType::Unknown;


	fLikelihoodDigitArray = 0x0;
    fEmissionProfileManager = 0x0;
	fPMTManager = new WCSimPMTManager();
}

WCSimTimeLikelihood3::WCSimTimeLikelihood3( WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager) :
			fSpeedOfParticle(1.0 * TMath::C() * 1e-7)
{
    fProbToCharge = new TGraph();
    fProbToCharge->SetName("fProbToCharge");
    fProbToCharge->SetMarkerStyle(20);
    fProbToCharge->SetMarkerSize(0.6);
    fProbToCharge->SetMarkerColor(kAzure);
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
  fLastE = -999;
  fLastType = TrackType::Unknown;
  return;
}

WCSimTimeLikelihood3::~WCSimTimeLikelihood3() {
	// TODO Auto-generated destructor stub
	if(fPMTManager != 0x0) { delete fPMTManager; }
    if( fTMinusTPredHist != 0x0){ delete fTMinusTPredHist; }
    if( fTMinusTPredSharpHist != NULL ){ delete fTMinusTPredSharpHist; }
    if( fHitQVsRMS != NULL ){ delete fHitQVsRMS; }
    if( fDistVsPred != NULL ){ delete fDistVsPred; }
    if( fDistVsProb != NULL ){ delete fDistVsProb; }
    if( fDistVsSpreadProb != NULL ){ delete fDistVsSpreadProb; }
    if( fTMinusTPredAll != NULL ){ delete fTMinusTPredAll; }
    if( fTMinusTPredSource != NULL ){ delete fTMinusTPredSource; }
    if( fTMinusTPredVsQ != NULL ){ delete fTMinusTPredVsQ; }
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
	double lnL = 0;  // Default to this so we return 25 as our default penalty

	WCSimLikelihoodDigit * myDigit = fLikelihoodDigitArray->GetDigit(iDigit);
	WCSimLikelihoodTrackBase * myTrack = fTracks.at(0);
	if(IsGoodDigit(myDigit))
	{
        double prob     = 1.0;
		double reso     = GetPMTTimeResolution(myDigit);
		double nPhotons = myDigit->GetQ();
		TimePrediction prediction = PredictArrivalTime(myTrack, myDigit);


        // We have four possibilities:
        // -------------------------------------------------------------------------------  
        // |                 | Predicted charge = 0         | Predicted charge > 0       |  
        // |------------------------------------------------|-----------------------------  
        // |Hit charge = 0   | (0x00) likelihood = 1        | (0x01) likelihood = exp(-Q)|  
        // |Hit charge > 0   | (0x10)likelihood = (1 - prob)| (0x11) calculate as normal |  
        // -------------------------------------------------------------------------------  
        int whichCase = 0x10 * (myDigit->GetQ() > 0) + 0x01 * (prediction.GetRMS() > 1e-6);
        // std::cout << "myDigit->GetTubeId() = " << myDigit->GetTubeId() << "whichCase = " << whichCase << std::endl;

		double timeLikelihood = fMinimumLikelihood;
        if(whichCase == 0x11) // Predict and detect a hit - use PDF overlap
        {
            // std::cout << " -- Doing overlap" << std::endl;
            float pmtLo  = myDigit->GetT() - 5.0 * reso;
            float pmtHi  = myDigit->GetT() + 5.0 * reso;
            float predLo = prediction.GetMean() - 5 * prediction.GetRMS();
            float predHi = prediction.GetMean() + 5 * prediction.GetRMS();

            // Function containing the overlap between the predicted first arrival PDF and the smeared true hit time
            TF1 overlapFunc("overlapFunc", MinimumOfArrivalAndResolution,
            				prediction.GetMean()-5*prediction.GetRMS(),
							prediction.GetMean()+5*prediction.GetRMS(), 5);
            overlapFunc.SetParameters(nPhotons,
            						  prediction.GetMean(), prediction.GetRMS(),
									  myDigit->GetT(), reso);
            overlapFunc.SetNpx(10000);
            timeLikelihood = overlapFunc.Integral(prediction.GetMean() - 5*prediction.GetRMS(), prediction.GetMean() + 5*prediction.GetRMS());
            if(timeLikelihood == 0){ whichCase = 0x10; }
        }
        else if(whichCase == 0x10) // Detect a hit but didn't predict one - use exp(-Q)
        {
            // std::cout << " -- Doing expo" << std::endl;
            timeLikelihood = TMath::Exp(-myDigit->GetQ());
            //timeLikelihood = 1.0;
        }
        else if(whichCase == 0x01) // Predict a hit but don't detect one - use predicted probability
        {   
            // std::cout << " -- Using prob" << std::endl;
            //timeLikelihood = (1.0 - prediction.GetProb());
            timeLikelihood = 1.0;
        }
        else if(whichCase == 0x00) // Didn't predict or detect a hit - carry on
        {
            // std::cout << " -- Using 1" << std::endl;
            timeLikelihood = 1.0;
        }
        else
        {
            std::cerr << "Didn't end up in a case I understand how to process, whichCase = " << whichCase << std::endl;
            assert(false);
        }

        if(timeLikelihood < 0 || timeLikelihood > 1)
		{
		  std::cout << "iDigt = " << iDigit << " and timeLikelihood = " << timeLikelihood << std::endl;
          std::cout << "timeLikelihood is " << timeLikelihood << " and is this > 1? " << (timeLikelihood > 1) << std::endl;
          std::cout << "timeLikelihood - 1 = " << timeLikelihood - 1 << std::endl;
		  assert(timeLikelihood > 0 && timeLikelihood <= 1);
		}
        else if(timeLikelihood >= 0 && timeLikelihood <= fMinimumLikelihood)
        {
            lnL = -fMaximumLnL;
        }
        else{
			lnL = log(timeLikelihood);
        }
		

		if(TMath::IsNaN(lnL) || !TMath::Finite(lnL))
		{
		  std::cout << "Digit is " << iDigit << " and timeLikelihood is " << timeLikelihood
				    << " because arrives at " << myDigit->GetT() << " and predict "
					<< prediction.GetMean() << " with width " << prediction.GetRMS()
					<< " and res " << reso << std::endl;
		  assert(!TMath::IsNaN(lnL));
          assert(TMath::Finite(lnL));
		}

        double start = prediction.GetMean() - 10 * prediction.GetRMS();
        if( start > myDigit->GetT() - reso){ start = myDigit->GetT() - reso; }
        double end   = prediction.GetMean() + 10 * prediction.GetRMS();
        if( end < myDigit->GetT() + reso){ end = myDigit->GetT() + reso; }

        /*
        TF1 * fMean = new TF1("fMean", ConvertMeanTimeToMeanFirstTime, start, end, 3);
        fMean->SetParameters(nPhotons, arrivalMeanSigma.first, arrivalMeanSigma.second);
        fMean->SetNpx(1000);
        TF1 * fInt = new TF1("fInt", IntegrateMeanTimeToMeanFirstTime, start, end, 3);
        fInt->SetParameters(nPhotons, arrivalMeanSigma.first, arrivalMeanSigma.second);
        fInt->SetNpx(1000);
        TF1 * fOverlap = new TF1("fOverlap", MinimumOfArrivalAndResolution, start, end, 5);
        fOverlap->SetParameters(nPhotons, mean, sigma, myDigit->GetT(), reso);
        fOverlap->SetNpx(1000);
        TF1 * fPMTTimeRes = new TF1("fPMTTimeRes","1.0 / ([1] * TMath::Sqrt(2 * TMath::Pi())) * TMath::Exp(-(x - [0])*(x-[0])/(2 * [1]*[1]))",start, end);
        fPMTTimeRes->SetParameters(myDigit->GetT(), reso);
        TCanvas* can11 = new TCanvas("can11","",800,600);
        TH1D hTemplate("hTemplate","Calculating the time likelihood;Arrival time (ns);P(t)",100, start, end);
        hTemplate.SetMaximum(1.2 * fMean->GetMaximum());
        fAllPreds[myDigit->GetTubeId()-1] = fMean->Mean(start, end, fMean->GetParameters());
        TLegend leg1(0.55, 0.85, 0.8, 0.6);
        leg1.AddEntry(fMean,"PDF for arrival time","L");
        leg1.AddEntry(fPMTTimeRes,TString::Format("PMT time measurement, #sigma = %.02f", reso).Data(),"L");
        leg1.AddEntry(fOverlap,TString::Format("Overlap = %.02e",timeLikelihood).Data() ,"FL");
        
        hTemplate.Draw("AXIS");
        hTemplate.GetXaxis()->CenterTitle();
        hTemplate.GetYaxis()->CenterTitle();
        hTemplate.SetStats(0);
        fMean->Draw("SAME");
        leg1.Draw();
        fPMTTimeRes->Draw("SAME");
        fOverlap->Draw("SAME");
        fMean->SetLineColor(kRed);
        fPMTTimeRes->SetLineColor(kBlue);
        fOverlap->SetLineColor(kViolet-1);
        fOverlap->SetFillColor(kViolet-1);
        fOverlap->SetFillStyle(3002);
        TString name = TString::Format("overlap_%d_%d", myDigit->GetTubeId(), (int)myTrack->GetE());
        can11->SetLogy();
        can11->SaveAs((name + ".png").Data());
        can11->SaveAs((name + ".pdf").Data());
        can11->SaveAs((name + ".C").Data());
        can11->SaveAs((name + ".root").Data());
        delete can11;
        //double predMean = fMean->GetX(fMean->Mean(start, end, fMean->GetParameters()));
        //double errLow = predMean - fInt->GetX(0.32);
        //double errHi = fInt->GetX(0.68) - predMean;
        ////std::cout << "nPhotons = " << nPhotons << "  ArrivalMeanSigma = (" << arrivalMeanSigma.first << ", " << arrivalMeanSigma.second << ")   PredMean = " << predMean << " and errors are " << fInt->GetX(0.32) << "  " << fInt->GetX(0.68) << std::endl;
        //double delta = myDigit->GetT() - predMean;
        //if( delta > 50 ) { delta = 50.0; }
        //fTMinusTPredSource->SetPoint(fTMinusTPredSource->GetN(), distance, delta);
        //fTMinusTPredSource->SetPointError(fTMinusTPredSource->GetN()-1, 0, 0, errLow, errHi);



		//fDistVsPred->SetPoint(fDistVsPred->GetN(), distance, arrivalMeanSigma.first);
		//fDistVsPred->SetPointError(fDistVsPred->GetN()-1, 0, arrivalMeanSigma.second);

		//fDistVsProb->SetPoint(fDistVsProb->GetN(), distance, -2.0*lnL);
		//fDistVsSpreadProb->SetPoint(fDistVsSpreadProb->GetN(), distance, -2.0*lnL);
        ////fTMinusTPredSource->SetPoint(fTMinusTPredSource->GetN(), distance, (myDigit->GetT()-arrivalMeanSigma.first));
		////fTMinusTPredSource->SetPointError(fTMinusTPredSource->GetN()-1, 0, arrivalMeanSigma.second);
		//fTMinusTPredAll->SetPoint(fTMinusTPredAll->GetN(), distance, delta);
		//fTMinusTPredAll->SetPointError(fTMinusTPredAll->GetN()-1, 0, 0, sqrt(reso * reso + errLow*errLow), sqrt(reso * reso  +errHi * errHi));
		//fTMinusTPredVsQ->SetPoint(fTMinusTPredVsQ->GetN(), myDigit->GetQ(), myDigit->GetT()-arrivalMeanSigma.first);
		//fTMinusTPredVsQ->SetPointError(fTMinusTPredVsQ->GetN()-1, 0, sqrt(sigmaSquared));
		//fHitQVsRMS->Fill(ceil(myDigit->GetQ()), arrivalMeanSigma.second);
        //fTMinusTPredHist->Fill(delta);
        //if( errLow < 1.0 ){
        //  fTMinusTPredSharpHist->Fill(delta);
        //}
        //delete fInt;
        //delete fMean;
        // delete fOverlap;
        // delete fPMTTimeRes;

	}
    if(myDigit->GetTubeId() > 518400){
  TCanvas * can  = new TCanvas("can","",800,600);
  fProbToCharge->Draw("AP");
  can->SaveAs("probToCharge.C");
  can->SaveAs("probToCharge.png");
  delete can;
    */
    }
    //std::cout << "TubeId " << myDigit->GetTubeId() << " Energy " << myTrack->GetE() << " -2LnL " << -2.0 * lnL << std::endl;
	return -2.0 * lnL;
}

double WCSimTimeLikelihood3::Calc2LnL() {
	double minusTwoLnL = 0.0;
	
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

void WCSimTimeLikelihood3::SetParticleSpeed(WCSimLikelihoodTrackBase * myTrack) {
  double speedOfLightInCmPerNs = TMath::C() * 1e-7;
  if(WCSimAnalysisConfig::Instance()->GetUseCustomParticleSpeed())
  {
	  fSpeedOfParticle = WCSimAnalysisConfig::Instance()->GetCustomParticleSpeed() * speedOfLightInCmPerNs;
  }
  else
  {
	  fSpeedOfParticle = myTrack->GetPropagationSpeedFrac() * speedOfLightInCmPerNs;
  }
}

double WCSimTimeLikelihood3::GetRefractiveIndex(WCSimLikelihoodDigit * myDigit) {
  double n = 1.0;
  if(WCSimAnalysisConfig::Instance()->GetUseCustomSpeedOfLight())
  {
    n = 1.0 / WCSimAnalysisConfig::Instance()->GetCustomSpeedOfLight();
  }
  else if(WCSimAnalysisConfig::Instance()->GetUseFittedSpeedOfLight())
  {
    n = 1.0 / WCSimAnalysisConfig::Instance()->GetFittedSpeedOfLight();
  }
  else
  {
    n = myDigit->GetAverageRefIndex();
  }
  return n;
}

bool WCSimTimeLikelihood3::IsGoodDigit(WCSimLikelihoodDigit * myDigit)
{
    return (myDigit->GetQ() > 0);
    return true;
}

TimePrediction WCSimTimeLikelihood3::PredictArrivalTime(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
  bool debug = false;
  if(myDigit->GetTubeId() == 5168 || myDigit->GetTubeId() == 5167)
  {
      debug = true;
  }


  // Load the emission profiles for the four emission profiles two either side of E
  // ==============================================================================
  GetFourNearestHists(myTrack);

  double minCosTheta = fSCosThetaVec[0]->GetXaxis()->GetXmin();
  double maxCosTheta = fSCosThetaVec[0]->GetXaxis()->GetXmax();
  const int nSBins = fSVec[0]->GetNbinsX();
  const int nThetaBins = fSCosThetaVec[0]->GetNbinsX();



  // Pre-calculate the values of (s,cosTheta) we need to look up from the emission profiles
  // for every distance step along the track
  // ======================================================================================
  std::vector<StepParameters> stepParameterVec(nSBins);
  TVector3 dir    = myTrack->GetDir();  // These get used a lot inside the loop
  TVector3 pmt    = myDigit->GetPos();  // So declare them here to only do it once
  double dCos     = -999;  // Width of the current cosTheta bin
  double magToPMT = 0.0;   // Distance from start of step to PMT
  double cosTheta = 0.0;   // cos of angle from start of step to PMT
  TVector3 toPMT;

  double s = fSVec[0]->GetXaxis()->GetBinCenter(1);     // Assumes all bins have the same width
  double sWidth = fSVec[0]->GetXaxis()->GetBinWidth(1); // to save on GetBinLowEdge() calls
  s = s - sWidth;   // Want bin 1 to start in the right place after adding sWidth

  for( int iSBin = 1; iSBin <= nSBins; ++iSBin )
  {
	  double s = s + sWidth;  // Take a step
	  if(dCos == -999) // i.e. the first step
	  {
		  toPMT    = pmt - myTrack->GetPropagatedPos(s);
		  magToPMT = toPMT.Mag() - myDigit->GetExposeHeight()/10.0; // Expose height is in mm
		  cosTheta = dir.Dot(toPMT)/magToPMT;
	  }

	  // Update it to work out the width in cosTheta that we traverse
	  // Then remember it because it'll be the starting point for the next loop iteration
	  toPMT = pmt - myTrack->GetPropagatedPos(s + sWidth);
	  magToPMT = toPMT.Mag() - myDigit->GetExposeHeight()/10.;
	  double newCos = dir.Dot(toPMT) / magToPMT;
	  if(! (cosTheta < minCosTheta || cosTheta >= maxCosTheta) )
	  {
		  stepParameterVec[iSBin-1] = StepParameters(cosTheta, iSBin, fabs(newCos - cosTheta), magToPMT);
	  }
	  cosTheta = newCos;
  }

  // If we sort it into cosTheta order then it's in order of ascending histogram
  // bin number and we only have to sweep through the axis once
  std::sort(stepParameterVec.begin(), stepParameterVec.end());

  // Work out the speeds for the charge particle and light
  SetParticleSpeed(myTrack);
  double n = GetRefractiveIndex(myDigit);
  double speedOfLightInCmPerNs = TMath::C() / 1e7;




  // For each of the four energies, calculate a time, RMS and weight:
  // ================================================================
  std::vector<double> meanVec(fSVec.size()); // Vectors to hold the values for each energy
  std::vector<double> rmsVec(fSVec.size());
  std::vector<double> probVec(fSVec.size());

  for(size_t iEnergy = 0; iEnergy < fSVec.size(); ++iEnergy)
  {
	  // A vector of the predicted hit times from each step, and one of each step's weight
	  std::vector<double> times(stepParameterVec.size(), 0.);
	  std::vector<double> probs(stepParameterVec.size(), 0.);
	  int index = 0;

	  // Now do the loop over cosTheta bins and look up the emission profile histograms to work out
	  // the photon arrival time at the PMT and a probability weighting
	  TAxis * axis = fSCosThetaVec[0]->GetXaxis(); // For looking up cosTheta bins
	  int lastBin = 1;
	  double weightedSum = 0.0;
	  double sumOfWeights = 0.0;

	  // Step along the track, calculating the hit time and probability for each step
	  // ============================================================================
	  for(std::vector<StepParameters>::iterator stepItr = stepParameterVec.begin(); stepItr != stepParameterVec.end(); ++stepItr)
	  {
		  // Had cosTheta outside the allowed range
          double cosTheta = stepItr->GetCosTheta();
		  if( cosTheta == -999){ continue; }

		  // Step through the theta bins keeping trackof the last bin so we only do it once
		  // for(int thBin = lastBin; thBin <= nThetaBins; ++thBin)
		  // {
			  // Does the current theta bin contain a cosTheta value we're interested in for this step?
              int thBin = fSCosThetaVec[iEnergy]->GetXaxis()->FindBin(cosTheta);

			  // if(axis->GetBinLowEdge(thBin) <= cosTheta && cosTheta < axis->GetBinUpEdge(thBin))
			  // {
				  // lastBin = thBin;

				  // Do the (cosTheta, s) profile first seeing as it's more likely to be zero
				  double prob = fSCosThetaVec[iEnergy]->GetBinContent(thBin, stepItr->GetSBin());
				  if(prob == 0)
				  {  // We can skip the other lookups if it's zero
					  continue;
				  }
                  //prob = prob * fSVec[iEnergy]->GetBinContent(stepItr->GetSBin()) 
                   //           * fSVec[iEnergy]->GetBinWidth(stepItr->GetSBin());

				  double t = 0;
				  // We only need to evaluate the arrival time if the probabilty is nonzero
				  if(prob > 0)
				  {
					  t =   myTrack->GetT()
								+ fSVec[iEnergy]->GetBinCenter(stepItr->GetSBin()) / fSpeedOfParticle
								+ n * stepItr->GetMagToPMT() / speedOfLightInCmPerNs;
				  }
				  else{
					  prob = 0;
				  }

				  // Keep track of everything for the moving average
				  weightedSum += prob * t;
				  sumOfWeights += prob;
                  // std::cout << "Digit " << myDigit->GetTubeId() << " thBin = " << thBin << " sBin = " << stepItr->GetSBin() << " prob = " << prob << " and sumOfWeights = " << sumOfWeights << "  sContent = " << fSVec[iEnergy]->GetBinContent(stepItr->GetSBin()) << " and thetaConstent = " << fSCosThetaVec[iEnergy]->GetBinContent(thBin, stepItr->GetSBin()) << std::endl;
                  if(sumOfWeights > 1){ TCanvas * can = new TCanvas("can","",800,600); fSCosThetaVec[iEnergy]->Draw(); can->SaveAs("above1.root"); assert(0); }

				  times[index] = t;
				  probs[index] = prob;

				  index++;
			  // }
		  // }
	  }


	  // Finally work out the mean and rms of the predicted times for this energy
	  // ========================================================================
	  double mean = 0;
	  double rms = 0;

	  if(sumOfWeights > 0)
	  {
		mean = weightedSum / sumOfWeights;
		for(int i = 0; i < index ; ++i)
		{
		  if(probs[i] == 0){ continue; }
		  rms += (mean - times[i])*(mean-times[i])*probs[i] / sumOfWeights;
		}
		rms = sqrt(rms);
	  }

	  meanVec[iEnergy] = mean;
	  rmsVec[iEnergy]  = rms;
	  // std::cout << myDigit->GetTubeId() << " at " << fEnergiesVec[iEnergy] << " has t=  " << mean << " and rms = " << rms << " and prob = " << sumOfWeights << std::endl;
	  probVec[iEnergy] = sumOfWeights;
  }


  // Construct a spline through each of the four values:
  // ===================================================
  double finalMean = WCSimFastMath::CatmullRomSpline(&(fEnergiesVec[0]), &(meanVec[0]), myTrack->GetE());
  double finalRMS  = WCSimFastMath::CatmullRomSpline(&(fEnergiesVec[0]), &(rmsVec[0]), myTrack->GetE());
  double finalProb = WCSimFastMath::CatmullRomSpline(&(fEnergiesVec[0]), &(probVec[0]), myTrack->GetE());
  if(finalMean < 0){ finalMean = 0; }
  if(finalRMS  < 0){ finalRMS  = 0; }
  if(finalProb < 0){ finalProb = 0; }
  fAllPreds[myDigit->GetTubeId()-1] = finalMean;


  return TimePrediction(finalMean, finalRMS, finalProb);
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

void WCSimTimeLikelihood3::GetFourNearestHists(WCSimLikelihoodTrackBase * myTrack)
{
    if(myTrack->GetType() == fLastType)
    {
        if(myTrack->GetE() == fLastE){ return; }
        if( fEmissionProfileManager->EnergiesInSameBin(fLastType, myTrack->GetE(), fLastE))
        {
            fLastE = myTrack->GetE();
            return; 
        }
    }

    fSCosThetaVec = fEmissionProfileManager->GetFourNearestSCosThetaForTimeHists(myTrack);
    fSVec         = fEmissionProfileManager->GetFourNearestSForTimeHists(myTrack);
    fEnergiesVec  = fEmissionProfileManager->GetFourNearestEnergies(myTrack);
    fLastE = myTrack->GetE();
    fLastType = myTrack->GetType();

    return;
}
