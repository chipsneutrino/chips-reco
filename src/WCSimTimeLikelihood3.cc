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
#include "TMarker.h"
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
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimTimeLikelihood3)
#endif

const double WCSimTimeLikelihood3::fMaximumLnL = 25.0/2.0; // So the maximum of -2LnL = 25
const double WCSimTimeLikelihood3::fMinimumLikelihood = TMath::Exp(-fMaximumLnL); // So the maximum of -2LnL = 25
const double WCSimTimeLikelihood3::fSqrtPi = TMath::Sqrt(TMath::Pi());
const double WCSimTimeLikelihood3::fLog2 = TMath::Log(2.0);

/**
 * Calculate the log of a Gaussian PDF normalised to unit area
 *
 * @param x The value at which to calculate the log
 * @param mean The mean of the Gaussian
 * @param sigma The width of the Gaussian
 * @return The natural log of the Gaussian PDF at this x value
 */
double LogGaussian(const double x, const double mean, const double sigma)
{
    return -log(sigma * WCSimTimeLikelihood3::fSqrtPi) - (x-mean)*(x-mean)/(2*sigma*sigma);
}

/**
 * Calculate the log of the PDF for the arrival time of the first photon
 * given the number of photons and the parameters of a Gaussian PDF describing
 * their arrival time
 *
 * @param x The value at which to calculate the log
 * @param n The number of photons, i.e. of repeated samples drawn from a Gaussian
 * @param mean The centre of the Gaussian being repeatedly sampled
 * @param rms The 1-sigma width of the Gaussian being repeatedly sampled
 * @return The log of the PDF for the earliest of n samples from the Gaussian
 */
double LogFirstArrival(const double x, const int n, const double mean, const double rms)
{
	double xmm = x-mean;
    return (   log(n * TMath::Sqrt2() / (rms * WCSimTimeLikelihood3::fSqrtPi)) 
             - n * WCSimTimeLikelihood3::fLog2 
             + (n-1)*log(TMath::Erfc((xmm)/(rms*TMath::Sqrt2())))
             - (xmm)*(xmm) / (2*rms*rms) );
}

/**
 * Calculate the difference between the log of a Gaussian and of the PDF for the first
 * arrival time, used for working out where the PMT resolution PDF and the photon arrival
 * PDF cross
 *
 * @param x The value at which to calculate the difference in logs
 * @param n The number of photons, i.e. of repeated samples drawn from a Gaussian
 * @param predMean The centre of the Gaussian repeatedly sampled to generate the first arrival time
 * @param predRMS The 1-sigma width of the Gaussian repeatedly sampled to generate the first arrival time
 * @param t The centre of the Gaussian (i.e the PMT hit time)
 * @param reso The 1-sigma width of the Gaussian (i.e. the PMT time resolution)
 * @return The difference between the log of the first arrival time PDF and the PMT resolution PDF
 */
double DiffLogGaussFirstArrival(const double x, const int n, const double predMean, const double predRMS, const double t, const double reso)
{
    return LogFirstArrival(x, n, predMean, predRMS) - LogGaussian(x, t, reso);
}

/**
 * Calculate the difference between the log of a Gaussian and of the PDF for the first
 * arrival time.  This has double-pointer arguments so we can make it into a TF1
 * @param x Pointer to the value at which to calculate the difference in logs
 * @param par An array containing: number of photons, centre of repeatedly-sampled Gaussian, RMS of same Gaussian, PMT hit time, PMT time resolution
 * @return The difference between the log of the first arrival time PDF and the PMT resolution PDF
 */
double DiffLogGaussFirstArrival(double * x, double * par)
{
    return DiffLogGaussFirstArrival(x[0], (int)par[0], par[1], par[2], par[3], par[4]);
}

/**
 * Calculate the time at which the PMT time resolution PDF and the first arrival time
 * PDF cross
 * @param n The number of photons, i.e. of repeated samples drawn from a Gaussian
 * @param predMean The centre of the Gaussian repeatedly sampled to generate the first arrival time
 * @param predRMS The 1-sigma width of the Gaussian repeatedly sampled to generate the first arrival time
 * @param t The centre of the Gaussian (i.e the PMT hit time)
 * @param reso The 1-sigma width of the Gaussian (i.e. the PMT time resolution)
 * @return The x-value at which the two PDFs are equal
 */
double FindCrossing(const int n, const double predMean, const double predRMS, const double t, const double reso)
{
    double min = t < predMean ? t : predMean;
    double max = t > predMean ? t : predMean;
    assert(max != min);
    TF1 f("func", DiffLogGaussFirstArrival, min, max, 5);
    f.SetParameters(n, predMean, predRMS, t, reso);

    // Create the function and wrap it
    ROOT::Math::WrappedTF1 wf1(f);
  
    // Create the Integrator
    ROOT::Math::BrentRootFinder brf;
  
    // Set parameters of the method
    brf.SetFunction( wf1, min, max );
    brf.Solve();

    double crossingX = brf.Root();
    return crossingX;
}


/**
 * Calculate an approximation to the aera of overlap between the PMT resolution and first arrival time PDFS
 * using the trapezium rule
 * @param n The number of photons, i.e. of repeated samples drawn from a Gaussian
 * @param predMean The centre of the Gaussian repeatedly sampled to generate the first arrival time
 * @param predRMS The 1-sigma width of the Gaussian repeatedly sampled to generate the first arrival time
 * @param t The centre of the Gaussian (i.e the PMT hit time)
 * @param reso The 1-sigma width of the Gaussian (i.e. the PMT time resolution)
 * @return Approximate area of overlap between the two PDFs
 */
double ApproximateIntegral(const int n, const double predMean, const double predRMS, const double t, const double reso)
{
	// Work out where the two PDFs intersect
    double x = FindCrossing(n, predMean, predRMS, t, reso);
    double logYCrossing = LogGaussian(x, t, reso);
    double gausExponent  = 0;
    double firstExponent =  0;

    // We want to integrate the PMT PDF to a distance reso from the crossing and the
    // arrival time PDF to 3*predRMS - so work out which side is which and do that
    bool predGreater = (predMean > t);
    if(predGreater)
    {
        gausExponent = LogGaussian(x-reso, t, reso);
        firstExponent = LogFirstArrival(x+predRMS, n, predMean, predRMS);
    }
    else
    {
        gausExponent = LogGaussian(x+reso, t, reso);
        firstExponent = LogFirstArrival(x-predRMS, n, predMean, predRMS);
    }

    // Use the trapezium rule to approximate the integral
    double middle = exp(logYCrossing);
    double gaus   = exp(gausExponent);
    double erf    = exp(firstExponent);
    double area   = 0.5*(middle+gaus)*reso + 0.5*(middle+erf)*predRMS;
    return area;
}


/**
 * Given n samples drawn from a Gaussian distribution, return the PDF for the earliest sample
 * @param t The time of the first photon, whose PDF we want to evaluate
 * @param mean The mean of the Gaussian describing the individual photon times
 * @param sigma The width of the Gaussian describing the individual photon times
 * @param nPhotons The number of photons
 * @return The PDF for the first arrival time being equal to t
 */
double ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons)
{
    double prob = 0.0;
    if(nPhotons == 0 || sigma == 0) { return 0; }
    if(nPhotons < 25)
    {
	    prob =   nPhotons * TMath::Sqrt2() / (pow(2, nPhotons) * sigma * WCSimTimeLikelihood3::fSqrtPi)
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

/**
 * Given n samples drawn from a Gaussian distribution, return the PDF for the earliest sample,
 * with double pointer arguments for making it into a TF1
 *
 * @param x Pointer to the value at which we want to evaluate the PDF
 * @param par Array containing the underlying Gaussian's mean and RMS, and number of photons
 * @return The PDF for the first arrival time being equal to x[0]
 * @param par Array containing the underlying Gaussian's mean and RMS, and number of photons
 */
double ConvertMeanTimeToMeanFirstTime(double * x, double * par)
{
	double nPhotons = par[0];
	double mean = par[1];
	double sigma = par[2];
	double t = x[0];
	return ConvertMeanTimeToMeanFirstTime(t, mean, sigma, nPhotons);
}

/**
 * Evaluate the integral of the PDF for the first arrival time
 * @param t Integral evaluated from -infinity to t
 * @param mean Mean of the underlying Gaussian
 * @param sigma Width of the underlying Gaussian
 * @param nPhotons Number of photons, i.e. number of times that Gaussian is sampled
 * @return Integral of the PDF for the first arrival time from -infinity to t
 */
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

/**
 * Evaluate the integral of the PDF for the first arrival time - with double
 * pointer arguments so it can be made into a TF1
 * @param x Integral evaluated from -infinity to x[0]
 * @param par Array containing the underlying Gaussian's mean and RMS, and number of photons
 * @return Integral of the PDF for the first arrival time from -infinity to x[0]
 */
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

            // Make sure we evaluate the TF1 at a reasonable number of points
            int npx = (int)((prediction.GetMean() - myDigit->GetT())*1000);
            if(npx < 1000){ npx = 1000; }
            if(npx > 10000000){ npx = 10000000;}
            overlapFunc.SetNpx(npx);
            timeLikelihood = overlapFunc.Integral(prediction.GetMean() - 5*prediction.GetRMS(), prediction.GetMean() + 5*prediction.GetRMS());

            // If the integral was too small there's probably some floating point complexity breaking it
            if(timeLikelihood <= fMinimumLikelihood)
            {
                // Try constructing a zoom of the overlap region by numerically solving for where
                // the two PDFs cross, using their logarithms
                double crossAt = FindCrossing(nPhotons, prediction.GetMean(), prediction.GetRMS(), myDigit->GetT(), reso);
                double start   = crossAt - 3*prediction.GetRMS();
                double end     = crossAt + 3*prediction.GetRMS();

                // Plot the overlap between +/- 3sigma around the overlap points and integrate it
                TF1 overlapFuncZoom("overlapFuncZoom", MinimumOfArrivalAndResolution, start, end, 5);
                overlapFuncZoom.SetParameters(nPhotons,
            						      prediction.GetMean(), prediction.GetRMS(),
									      myDigit->GetT(), reso);
                overlapFuncZoom.SetNpx(100000);
                timeLikelihood = overlapFuncZoom.Integral(start, end);

                
                // Even this integral was too small!  Normally the problem comes from erf^(n-1) / 2^n so we can 
                // use logs to get something sensible for lnL and then exponentiate it.  We'll do that, using
                // the trapezium rule around the crossing point
                if(timeLikelihood == 0)
                {
                    timeLikelihood = ApproximateIntegral(nPhotons, prediction.GetMean(), prediction.GetRMS(), myDigit->GetT(), reso);
                }
            }

        }
        else if(whichCase == 0x10) // Detect a hit but didn't predict one - use exp(-Q)
        {
            // std::cout << " -- Doing expo" << std::endl;
            // timeLikelihood = TMath::Exp(-myDigit->GetQ());
            timeLikelihood = 1.0;
        }
        else if(whichCase == 0x01) // Predict a hit but don't detect one - use predicted probability
        {   
            // std::cout << " -- Using prob" << std::endl;
            // timeLikelihood = (1.0 - prediction.GetProb());
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


        // Check likelihood is sensible
        // ============================
        if(timeLikelihood < 0 || timeLikelihood > 1)
		{
		  std::cout << "iDigit = " << iDigit << " and timeLikelihood = " << timeLikelihood << std::endl;
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
    }


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
  int numNonZero = 0; // Count how many nonzero probabilities we get - need 4 to make a spline
  std::vector<double> meanVec(fSVec.size()); // Vectors to hold the values for each energy
  std::vector<double> rmsVec(fSVec.size());
  std::vector<double> probVec(fSVec.size());

  std::vector<TMarker> sThetaMarkers[4];
  std::vector<TMarker> timeMarkers[4];

  for(size_t iEnergy = 0; iEnergy < fSVec.size(); ++iEnergy)
  {
      int lastColor = 30;
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

		  // Does the current theta bin contain a cosTheta value we're interested in for this step?
          int thBin = fSCosThetaVec[iEnergy]->GetXaxis()->FindBin(cosTheta);

		  // Do the (cosTheta, s) profile first seeing as it's more likely to be zero
		  double prob = fSCosThetaVec[iEnergy]->GetBinContent(thBin, stepItr->GetSBin());
		  if(prob == 0)
		  {  // We can skip the other lookups if it's zero
		      continue;
		  }

		  double t = 0;

		  // We only need to evaluate the arrival time if the probability is nonzero
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
		  weightedSum  += prob * t;
		  sumOfWeights += prob;
		  times[index] = t;
		  probs[index] = prob;

		  index++;
	  }


	  // Finally work out the mean and rms of the predicted times for this energy
	  // ========================================================================
	  double mean = 0;
	  double rms  = 0.1; // We have to start this off non-zero in case there's only one non-zero-weighted time
                         // 0.1ns seems high enough to give a smooth likelihood surface as we crank up the energy

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

      // If the mean comes out zero we'll get a really distorted spline because it
      // has to go through three nearby points and zero. So default the time to the
      // straight line time taken by a photon from the vertex, but keep the weight zero
      if(mean != 0)
      {
          numNonZero++;
      }


	  meanVec[iEnergy] = mean;
	  rmsVec[iEnergy]  = rms;
	  probVec[iEnergy] = sumOfWeights;
  }

  // Construct a spline through each of the four values:
  // ===================================================
  double finalMean = 0;
  double finalRMS  = 0;
  double finalProb = 0;

  // We need four nonzero points to do the spline, otherwise the zeros make it unnaturally curvy
  if(numNonZero == 4)
  {
     finalMean = WCSimFastMath::CatmullRomSpline(&(fEnergiesVec[0]), &(meanVec[0]), myTrack->GetE());
     finalRMS  = WCSimFastMath::CatmullRomSpline(&(fEnergiesVec[0]), &(rmsVec[0]), myTrack->GetE());
     finalProb = WCSimFastMath::CatmullRomSpline(&(fEnergiesVec[0]), &(probVec[0]), myTrack->GetE());
  }
  else
  {
     // We can't do the spline, so just get ROOT to interpolate/extrapolate using the nearest points
     TGraph gMeans;
     TGraph gRMSs;
     TGraph gProbs;
     for(size_t i = 0; i < 4; ++i)
     {
        if(meanVec[i] == 0){ continue; }
        gMeans.SetPoint(gMeans.GetN(), fEnergiesVec[i], meanVec[i]);
        gRMSs.SetPoint(gRMSs.GetN(), fEnergiesVec[i], rmsVec[i]);
        gProbs.SetPoint(gProbs.GetN(), fEnergiesVec[i], probVec[i]);
     }
     if(numNonZero == 1)
     {
         // We don't have enough points to draw a straight line - just take the values we have
         finalMean = gMeans.GetY()[0];
         finalRMS = gRMSs.GetY()[0];
         finalProb = gProbs.GetY()[0];
     }
     else
     {
        // Eval will interpolate between the two nearest points, or extrapolate, both linearly
        finalMean = gMeans.Eval(myTrack->GetE());
        finalRMS  = gRMSs.Eval(myTrack->GetE());
        finalProb = gProbs.Eval(myTrack->GetE());
     }
  }

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
