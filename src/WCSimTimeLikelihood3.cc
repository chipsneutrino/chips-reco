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
#include "Math/RootFinderAlgorithms.h"
#include "Math/RootFinder.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimTimeLikelihood3)
#endif

const double WCSimTimeLikelihood3::fMaximumLnL = 100.0/2.0; // So the maximum of -2LnL = 25
const double WCSimTimeLikelihood3::fMinimumLikelihood = TMath::Exp(-fMaximumLnL); // So the maximum of -2LnL = 25
const double WCSimTimeLikelihood3::fSqrtPi = TMath::Sqrt(TMath::Pi());
const double WCSimTimeLikelihood3::fLog2 = TMath::Log(2.0);
/**
 * Calculate the log of a Gaussian PDF normalised to unit area
 *
 * @param x The value at which to calculate the log
 * @return The natural log of the Gaussian PDF at this x value
 */
double LogFinder::LogGaussian(const double x) const
{
    return -log(fReso * TMath::Sqrt2() * WCSimTimeLikelihood3::fSqrtPi) - (x-fT)*(x-fT)/(2*fReso*fReso);
}

/**
 * Calculate the gradient of the log of a Gaussian PDF normalised to unit area
 *
 * @param x The value at which to calculate the log
 * @return The gradient of the natural log of the Gaussian PDF at this x value
 */
double LogFinder::LogGaussianGradient(const double x) const
{
    return (fT - x)/(fReso*fReso);
}

/**
 * Calculate the log of the PDF for the arrival time of the first photon
 * given the number of photons and the parameters of a Gaussian PDF describing
 * their arrival time
 *
 * @param x The value at which to calculate the log
 * @return The log of the PDF for the earliest of n samples from the Gaussian
 */
double LogFinder::LogFirstArrival(const double x) const
{
	double xmm = x-fPredMean;
    return (   log(fN * TMath::Sqrt2() / (fPredRMS * WCSimTimeLikelihood3::fSqrtPi)) 
             - fN * WCSimTimeLikelihood3::fLog2 
             + (fN-1)*log(TMath::Erfc((xmm)/(fPredRMS*TMath::Sqrt2())))
             - (xmm)*(xmm) / (2*fPredRMS*fPredRMS) );
}

/**
 * Calculate the gradient of the log of the PDF for the arrival time of the first photon
 * given the number of photons and the parameters of a Gaussian PDF describing
 * their arrival time
 *
 * @param x The value at which to calculate the gradient of the log
 * @param n The number of photons, i.e. of repeated samples drawn from a Gaussian
 * @param mean The centre of the Gaussian being repeatedly sampled
 * @param rms The 1-sigma width of the Gaussian being repeatedly sampled
 * @return The gradient of the log of the PDF for the earliest of n samples from the Gaussian
 */
double LogFinder::LogFirstArrivalGradient(const double x) const
{
	double xmm = x-fPredMean;
	return    1.0 / (fPredRMS * TMath::Sqrt2() * WCSimTimeLikelihood3::fSqrtPi)
			* 2*TMath::Exp(-xmm*xmm/(2*fPredRMS*fPredRMS))
			* ((1-fN)/TMath::Erfc(xmm/(fPredRMS*TMath::Sqrt2())))
			- (xmm)/(fPredRMS*fPredRMS);
}

/**
 * Calculate the difference between the log of a Gaussian and of the PDF for the first
 * arrival time, used for working out where the PMT resolution PDF and the photon arrival
 * PDF cross
 *
 * @param x The value at which to calculate the difference in logs
 * @return The difference between the log of the first arrival time PDF and the PMT resolution PDF
 */
double LogFinder::DiffLogGaussFirstArrival(const double x) const
{
    double gaus = LogGaussian(x);
    double arr  = LogFirstArrival(x);
    
    if(TMath::Finite(gaus - arr))
    {
        return gaus - arr;
    }
    else if(TMath::Finite(gaus))
    {
        return -arr;
    }
    else if(TMath::Finite(arr))
    {
        return gaus;
    }
    return -arr; // The arrival time probably falls off faster than the Gaussian, but all bets are off by this point
    return LogGaussian(x) - LogFirstArrival(x);
}

/**
 * Calculate the gradient of the difference between the log of a Gaussian and of the PDF for the first
 * arrival time, used for working out where the PMT resolution PDF and the photon arrival
 * PDF cross
 *
 * @param x The value at which to calculate the gradient of the difference in logs
 * @return The gradient of the difference between the log of the first arrival time PDF and the PMT resolution PDF
 */
double LogFinder::DiffLogGaussFirstArrivalGradient(const double x) const
{
    return LogGaussianGradient(x) - LogFirstArrivalGradient(x);
}



WCSimTimeLikelihood3::WCSimTimeLikelihood3() : fSpeedOfParticle(1.0 * TMath::C() * 1e-7){
	// TODO Auto-generated constructor stub
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
  fLastE = -999;
  fLastType = TrackType::Unknown;
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
	double lnL = 0;  // Default to this so we return 25 as our default penalty

	WCSimLikelihoodDigit * myDigit = fLikelihoodDigitArray->GetDigit(iDigit);
	WCSimLikelihoodTrackBase * myTrack = fTracks.at(0);
	if(IsGoodDigit(myDigit))
	{
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
        int whichCase = 0x10 * (myDigit->GetQ() > 0) + 0x01 * (prediction.GetProb() > fMinimumLikelihood);
        // std::cout << "myDigit->GetTubeId() = " << myDigit->GetTubeId() << "whichCase = " << whichCase << std::endl;

		double timeLikelihood = fMinimumLikelihood;
        if(whichCase == 0x11) // Predict and detect a hit - use PDF overlap
        {
            // Function containing the overlap between the predicted first arrival PDF and the smeared true hit time
            timeLikelihood = FindOverlap(nPhotons, prediction.GetMean(), prediction.GetRMS(), myDigit->GetT(), reso);

            // If the integral was too small there's probably some floating point complexity breaking it
            if(timeLikelihood <= fMinimumLikelihood || timeLikelihood > 1)
            {
                // Try constructing a zoom of the overlap region by numerically solving for where
                // the two PDFs cross, using their logarithms
                double start   = myDigit->GetT() - 10 * reso;
                double end     = myDigit->GetT() + 10 * reso;

                // Plot the overlap between +/- 3sigma around the overlap points and integrate it
                TF1 overlapFuncZoom("overlapFuncZoom", this, &WCSimTimeLikelihood3::MinimumOfArrivalAndResolutionTF1,
                					start, end, 5, "WCSimTimeLikelihood3", "MinimumOfArrivalAndResolutionTF1");
                overlapFuncZoom.SetParameters(nPhotons,
            						      prediction.GetMean(), prediction.GetRMS(),
									      myDigit->GetT(), reso);
                ROOT::Math::WrappedTF1 f1(overlapFuncZoom);
                ROOT::Math::Integrator ig(f1, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
                timeLikelihood = ig.Integral(prediction.GetMean() - 5*prediction.GetRMS(), prediction.GetMean() + 5*prediction.GetRMS());

                // Even this integral was too small!  Normally the problem comes from erf^(n-1) / 2^n so we can 
                // use logs to get something sensible for lnL and then exponentiate it.  We'll do that, using
                // the trapezium rule around the crossing point
                /* Causes GSL errors because it often has infinites at the either end of the t range - let's write
                 * these events off for now
                if(timeLikelihood == 0)
                {
    
                    timeLikelihood = ApproximateIntegral(nPhotons, prediction.GetMean(), prediction.GetRMS(), myDigit->GetT(), reso);
                }
                */
            }
            
            // Add on a flat distribution of scattered light times
            timeLikelihood += GetScatteredTimeLikelihood(myDigit);
        }
        else if(whichCase == 0x10) // Detect a hit but didn't predict one - use exp(-Q)
        {
            timeLikelihood = GetScatteredTimeLikelihood(myDigit);
            //timeLikelihood = 1.0;
        }
        else if(whichCase == 0x01) // Predict a hit but don't detect one - use predicted probability
        {   
            // std::cout << " -- Using prob" << std::endl;
            //timeLikelihood = (1.0 - prediction.GetProb());
            //timeLikelihood = 1.0;
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

  for(size_t iEnergy = 0; iEnergy < fSVec.size(); ++iEnergy)
  {
	  // A vector of the predicted hit times from each step, and one of each step's weight
	  std::vector<double> times(stepParameterVec.size(), 0.);
	  std::vector<double> probs(stepParameterVec.size(), 0.);
	  std::vector<double> squareUncertainties(stepParameterVec.size(), 0.);
	  int index = 0;

	  // Now do the loop over cosTheta bins and look up the emission profile histograms to work out
	  // the photon arrival time at the PMT and a probability weighting
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
          int thBin = 1 + (cosTheta - minCosTheta)/(maxCosTheta - minCosTheta) * nThetaBins;

		  // Do the (cosTheta, s) profile first seeing as it's more likely to be zero
		  double prob = fSCosThetaVec[iEnergy]->GetBinContent(thBin, stepItr->GetSBin());
		  if(prob == 0)
		  {  // We can skip the other lookups if it's zero
		      continue;
		  }

		  double t = 0;
          double uncertaintySquared = 0;

		  // We only need to evaluate the arrival time if the probability is nonzero
		  if(prob > 0)
		  {
              double particleDistance = fSVec[iEnergy]->GetBinCenter(stepItr->GetSBin());
              double photonDistance = stepItr->GetMagToPMT();
		      t =   myTrack->GetT()
		    			+ fSVec[iEnergy]->GetBinCenter(stepItr->GetSBin()) / fSpeedOfParticle
		    			+ n * stepItr->GetMagToPMT() / speedOfLightInCmPerNs;

              // Propagate a 3% uncertainty on the muon and photon speeds
              // to allow for a bit of scattering/wavelength dependence of n
              uncertaintySquared =   (0.02 * particleDistance / fSpeedOfParticle) * (0.02 * particleDistance / fSpeedOfParticle)
                                   + (0.02 * photonDistance * n / speedOfLightInCmPerNs) * (0.02 * photonDistance * n / speedOfLightInCmPerNs);
		  }
		  else{
		      prob = 0;
		  }

		  // Keep track of everything for the moving average
		  weightedSum  += prob * t;
		  sumOfWeights += prob;
		  times[index] = t;
		  probs[index] = prob;
          squareUncertainties[index] = uncertaintySquared;

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


        double speedUncertaintySq = 0.0;
        for(int i = 0; i < index; ++i)
        {
            speedUncertaintySq += probs[i] * probs[i] * squareUncertainties[i];
        }
        speedUncertaintySq = speedUncertaintySq / (sumOfWeights*sumOfWeights);


		rms = sqrt(rms + speedUncertaintySq);
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
  else if(numNonZero > 0)
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
        if(finalRMS < 0.1){ finalRMS = 0.1; } 
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

double WCSimTimeLikelihood3::FindOverlap(const double n, const double predMean, const double predRMS, const double t, const double reso)
{
	std::vector<double> crossingsX = FindCrossings(n, predMean, predRMS, t, reso);
    LogFinder finder(n, predMean, predRMS, t, reso);
    if(crossingsX.empty()){ return 0; }

	double integral = 0.0;

	size_t numCrossings = crossingsX.size();

	// The first interval is from -infinity to the first root:
	bool gaussianLower = (finder.DiffLogGaussFirstArrival(crossingsX[0]-0.1) < 0);
	if(gaussianLower)
	{
        double beginning = IntegrateGaussianFromMinusInfinity(t, reso, crossingsX[0]);
        integral += beginning;
        // std::cout << "Adding beginning: gaus= " << beginning << " total = " << integral << std::endl;
	}
	else
	{
        double beginning = IntegrateMeanTimeToMeanFirstTime(crossingsX[0], predMean, predRMS, n);
        integral += beginning;
        // std::cout << "Adding beginning: erfc= " << beginning << " total = " << integral << std::endl;
	}

	// The last interval is from the last root to +infinity
	// So we need to start from slightly past the root
	gaussianLower = (finder.DiffLogGaussFirstArrival(crossingsX.back()+0.1) < 0);

	if(gaussianLower)
	{
		double ending = IntegrateGaussianToInfinity(t, reso, crossingsX.back());
        integral += ending;
        // std::cout << "Adding ending: gaus= " << ending << " total = " << integral << std::endl;
	}
	else
	{
        double ending = 1 - IntegrateMeanTimeToMeanFirstTime(crossingsX.back(), predMean, predRMS, n);
        integral += ending;
        // std::cout << "Adding ending: erf = " << ending << " total = " << integral << std::endl;
	}

	// Now loop over the intervals in between
	// Interval i goes from crossingsX[i] to crossingsX[i+1]
	for(size_t interval = 0; interval < numCrossings - 1; ++interval)
	{
		double min = crossingsX[interval];
		double max = crossingsX[interval+1];
		double mid = 0.5 * (min+max);
		gaussianLower = (finder.DiffLogGaussFirstArrival(mid) < 0);

		if(gaussianLower)
		{
            double intervalInt = IntegrateGaussian(t, reso, min, max);
            integral+= intervalInt;
            // std::cout << "Interval " << interval << " and gauss int = " << intervalInt << " so total = " << integral << std::endl;
		}
		else
		{
            double intervalInt =   IntegrateMeanTimeToMeanFirstTime(max, predMean, predRMS, n)
						         - IntegrateMeanTimeToMeanFirstTime(min, predMean, predRMS, n);
            integral += intervalInt;
            // std::cout << "Interval " << interval << " and erf   int = " << intervalInt << " so total = " << integral << std::endl;
		}
	}
    if(integral > 1)
    {
    	integral = 1;
    }
	return integral;
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
double WCSimTimeLikelihood3::FindCrossing(const double n, const double predMean, const double predRMS, const double t, const double reso)
{
    double min = t < predMean ? t : predMean;
    double max = t > predMean ? t : predMean;
    assert(max != min);

    // RootFinder with derivative functions
    LogFinder finder(n, predMean, predRMS, t, reso);
    ROOT::Math::GradFunctor1D  f(finder);
    ROOT::Math::RootFinder rfn(ROOT::Math::RootFinder::kGSL_BRENT);
    rfn.SetFunction(f, min, max);
    int result = rfn.Solve();
    double crossingX = rfn.Root();

    return crossingX;
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
std::vector<double> WCSimTimeLikelihood3::FindCrossings(const double n, const double predMean, const double predRMS, const double t, const double reso)
{
   std::vector<double> minima;
   // Go within 5 sigma of the measured or predicted hit time, 
   // using whichever of the predicted or measured values gives the wider range
   double min = t-5*reso < predMean-5*predRMS ? t-5*reso : predMean-5*predRMS;
   double max = t+5*reso > predMean+5*predRMS ? t+5*reso : predMean + 5*predRMS;

   // RootFinder with derivative functions
   LogFinder finder(n, predMean, predRMS, t, reso);
   ROOT::Math::GradFunctor1D  f(finder);

   // Break the full range into windows of width 0.2 * rms
   double width = 5 * (predRMS > reso ? predRMS : reso);
   int nWindows = TMath::CeilNint(width / 0.2*(predRMS < reso ? predRMS : reso));
   double dx = 2 * width / nWindows;

   // Lower and upper edges of windows where the function changes sign
   std::vector<std::pair<double, double> > windowEdges; 

   double start = t - width;
   double x = start;

   double atMin = finder.Eval(min);
   double lastEdge = finder.Eval(x);

   // Look for a sign change in this first step
   if(TMath::Finite(lastEdge) && TMath::Finite(lastEdge) && atMin * lastEdge < 0 && min < x)
   { 
       windowEdges.push_back(std::make_pair(min, x));
   }
   for(int i = 1; i < nWindows; ++i)
   {
         x += dx;
         double y = finder.Eval(x);

         bool cross = (y*lastEdge <= 0);
         if(cross)
         {
            if(TMath::Finite(y) && TMath::Finite(lastEdge))
            {
				windowEdges.push_back(std::make_pair(x-dx, x));
            }
            else if(TMath::Finite(y))
            {
                for(int step = 0; step < 100; ++step)
                {
                    double newLowEdge = x - step*dx/100.0;
                    double newLast = finder.Eval(newLowEdge);
                    if(TMath::Finite(newLast) && newLast * y < 0)
                    {
				        windowEdges.push_back(std::make_pair(newLast, x));
                        break;
                    }
                }

            }
            else if(TMath::Finite(lastEdge))
            {
                for(int step = 0; step < 100; ++step)
                {
                    double newUpEdge = x - dx + step*dx/100.0;
                    double newY = finder.Eval(newUpEdge);
                    if(TMath::Finite(newY) && lastEdge * newY < 0)
                    {
				        windowEdges.push_back(std::make_pair(x-dx, newUpEdge));
                        break;
                    }
                }
            }
			lastEdge = y;
         }
   }
   double atMax = finder.Eval(max);
   if(TMath::Finite(atMax) && atMax * lastEdge < 0 && max > x){ windowEdges.push_back(std::make_pair(x, max)); }


   ROOT::Math::RootFinder rfn(ROOT::Math::RootFinder::kGSL_BRENT);
   size_t nFound = windowEdges.size();
   for(size_t iWindow = 0; iWindow < nFound; ++iWindow)
   {
     double lo = windowEdges[iWindow].first;
     double hi = windowEdges[iWindow].second;


     if(!TMath::Finite(finder.Eval(lo)) || !TMath::Finite(finder.Eval(hi)))
     {
         continue;
     }
     rfn.SetFunction(f, lo, hi);
     int result = rfn.Solve();
     double root = rfn.Root();

     if ( (lo <= root || finder.ApproxEqual(root, lo, 1e-6) ) && ( root <= hi || finder.ApproxEqual(root, hi, 1e-6)) )
     {
         minima.push_back(root);
     }
     assert(lo < hi);
   }

   if(minima.size() == 1)
   {
       return minima;
   }

   std::vector<double> unique;
   std::sort(minima.begin(), minima.end());
   if( minima.size() > 1)
   {
	   unique.push_back(minima[0]);
   }
   for(size_t i = 1; i < minima.size(); ++i)
   {

     if(!finder.ApproxEqual(minima[i], minima[i-1], 2.5e-3))
     {
        unique.push_back(minima[i]);
     }
   }
   return unique;
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
double WCSimTimeLikelihood3::ApproximateIntegral(const double n, const double predMean, const double predRMS, const double t, const double reso)
{
	// Work out where the two PDFs intersect
    LogFinder finder(n, predMean, predRMS, t, reso);
    double x = FindCrossing(n, predMean, predRMS, t, reso);
    double logYCrossing = finder.LogGaussian(x);
    double gausExponent  = 0;
    double firstExponent =  0;

    // We want to integrate the PMT PDF to a distance reso from the crossing and the
    // arrival time PDF to 3*predRMS - so work out which side is which and do that
    bool predGreater = (predMean > t);
    if(predGreater)
    {
        gausExponent = finder.LogGaussian(x-reso);
        firstExponent = finder.LogFirstArrival(x+predRMS);
    }
    else
    {
        gausExponent = finder.LogGaussian(x+reso);
        firstExponent = finder.LogFirstArrival(x-predRMS);
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
double WCSimTimeLikelihood3::ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons)
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

double WCSimTimeLikelihood3::IntegrateGaussian(const double mean, const double sigma,
		const double from, const double to) {
	double denom = sigma*TMath::Sqrt2();
	return 0.5 * (TMath::Erf((to - mean)/denom) - TMath::Erf((from - mean)/denom));
}

double WCSimTimeLikelihood3::IntegrateGaussianFromMinusInfinity(const double mean, const double sigma,
		const double to) {
	if(to < mean)
	{
		return 0.5 - IntegrateGaussian(mean, sigma, to, mean);
	}
	return 0.5+IntegrateGaussian(mean, sigma, mean, to);
}

double WCSimTimeLikelihood3::IntegrateGaussianToInfinity(const double mean, const double sigma,
		const double from) {
	if(from > mean)
	{
		return 0.5 - IntegrateGaussian(mean, sigma, mean, from);
	}
	return 0.5 + IntegrateGaussian(mean, sigma, from, mean);

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
double WCSimTimeLikelihood3::ConvertMeanTimeToMeanFirstTime(double * x, double * par)
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
double WCSimTimeLikelihood3::IntegrateMeanTimeToMeanFirstTime( const double &t, const double &mean, const double &sigma, const double &nPhotons)
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
double WCSimTimeLikelihood3::IntegrateMeanTimeToMeanFirstTime(double * x, double * par)
{
  double &nPhotons = (par[0]);
  double &mean = (par[1]);
  double &sigma = (par[2]);
  double &t = (x[0]);

  return IntegrateMeanTimeToMeanFirstTime(t, mean, sigma, nPhotons);
}

double WCSimTimeLikelihood3::MinimumOfArrivalAndResolution(const double &t, const double &meanArr, const double rmsArr, const double pmtTime, const double &pmtRes, const double &nPhotons)
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

double WCSimTimeLikelihood3::MinimumOfArrivalAndResolutionTF1( double * x, double * par )
{
  double &nPhotons = (par[0]);
  double &meanArr = (par[1]);
  double &rmsArr = (par[2]);
  double &pmtTime = (par[3]);
  double &pmtRes = (par[4]);
  double &t = (x[0]);

  return MinimumOfArrivalAndResolution(t, meanArr, rmsArr, pmtTime, pmtRes, nPhotons);
}

double WCSimTimeLikelihood3::GetScatteredTimeLikelihood(WCSimLikelihoodDigit * myDigit)
{
    double duration = fLikelihoodDigitArray->GetDuration();
    if(duration < 20){ duration = 20; } // Put a floor of 20ns on the event length
    double flatProb = 0.01 / duration; // Assume a 1% chance of scattered light distributed uniformly in time throughout the event

    // Want the overlap of y = flatProb with the Gaussian describing the PMT resolutionA
	double reso = GetPMTTimeResolution(myDigit);
    double rhs = flatProb * reso * fSqrtPi * TMath::Sqrt2();

    if( rhs >= 1) { return fMinimumLikelihood; }
    double crossX = myDigit->GetT() - TMath::Sqrt(-2 * reso * TMath::Log(flatProb * reso * fSqrtPi * TMath::Sqrt2()));
    double prob = 2 * IntegrateGaussianFromMinusInfinity(myDigit->GetT(), reso, crossX);
    if(prob < fMinimumLikelihood){ prob = fMinimumLikelihood;}
    return prob;
}

