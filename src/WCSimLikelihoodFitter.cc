#include "WCSimAnalysisConfig.hh"
#include "WCSimChargePredictor.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimFitterTree.hh"
#include "WCSimFitterInterface.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimReco.hh"
#include "WCSimRecoSeed.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoFactory.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimTrackParameterEnums.hh"

#include "TClonesArray.h"
#include "TCollection.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <string>
#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include <algorithm>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodFitter)
#endif

bool WCSimLikelihoodFitter::RingSort(const std::pair<WCSimRecoRing*,double> &a, const std::pair<WCSimRecoRing*,double> &b){
  return (a.first)->GetHeight() > (b.first)->GetHeight();
}

/**
 * @todo Remove hardcoding of track type
 */
WCSimLikelihoodFitter::WCSimLikelihoodFitter(WCSimFitterConfig * config) : fFitterTrackParMap(config)
{
  fCalls = 0;
  fEvent = -999;
  fFitterConfig = config;
  fFitterPlots = NULL;
  fFitterTree = 0x0;
  fTotalLikelihood = NULL;
  fRootEvent = NULL;
  fLikelihoodDigitArray = NULL;
  fTrueLikelihoodTracks = NULL;

  fUseHoughFitterForSeed = WCSimAnalysisConfig::Instance()->GetUseHoughFitterForSeed();
  ResetEvent();
}

WCSimLikelihoodFitter::~WCSimLikelihoodFitter()
{
}

UInt_t WCSimLikelihoodFitter::GetNPars()
{
  // Do we know how to fit this number of tracks?
  unsigned int nPars = fFitterConfig->GetNumIndependentParameters();
  return nPars;
}



void WCSimLikelihoodFitter::Minimize2LnL()
{
}

/**
 * Wrapper: constructs the correct number of track objects
 * then works out the total likelihood
 */
Double_t WCSimLikelihoodFitter::WrapFunc(const Double_t * x)
{
  UInt_t nTracks = fFitterConfig->GetNumTracks();

  // If we've fixed some track parameters together, then our array of fit parameters won't just be of size
  // n tracks * m parameters per track, and we need to work out which entry corresponds to which parameter
  // of which track(s) - this is probably expensive, so we'll do it once and look it up in a map thereafter

  std::vector<WCSimLikelihoodTrackBase*> tracksToFit;
  for(UInt_t iTrack = 0 ; iTrack < nTracks; ++iTrack)
  {
	  TrackType::Type trackType = fFitterTrackParMap.GetTrackType(iTrack);
	  std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes(trackType);
	  TrackAndType trackParX = std::make_pair(iTrack, FitterParameterType::kVtxX);
	  TrackAndType trackParY = std::make_pair(iTrack, FitterParameterType::kVtxY);
	  TrackAndType trackParZ = std::make_pair(iTrack, FitterParameterType::kVtxZ);
	  TrackAndType trackParT = std::make_pair(iTrack, FitterParameterType::kVtxT);
	  TrackAndType trackParTh = std::make_pair(iTrack, FitterParameterType::kDirTh);
	  TrackAndType trackParPhi = std::make_pair(iTrack, FitterParameterType::kDirPhi);
	  TrackAndType trackParE = std::make_pair(iTrack, FitterParameterType::kEnergy);

	  TrackAndType trackParConversionDistance = std::make_pair(iTrack, FitterParameterType::kConversionDistance);
	  std::map<FitterParameterType::Type, double> extraPars;
	  extraPars[FitterParameterType::kConversionDistance] = x[fFitterTrackParMap.GetIndex(trackParConversionDistance)];

	  double energy = x[fFitterTrackParMap.GetIndex(trackParE)];
	  if(iTrack == 1 && GetUsePiZeroMassConstraint())
	  {
		  energy = GetPiZeroSecondTrackEnergy(x);
	  }

	  WCSimLikelihoodTrackBase * track = WCSimLikelihoodTrackFactory::MakeTrack(
			  	  	  	  	  	  fFitterTrackParMap.GetTrackType(iTrack),
			  	  	  	  	  	  x[fFitterTrackParMap.GetIndex(trackParX)],
		          						  x[fFitterTrackParMap.GetIndex(trackParY)],
		          						  x[fFitterTrackParMap.GetIndex(trackParZ)],
		          						  x[fFitterTrackParMap.GetIndex(trackParT)],
		          						  x[fFitterTrackParMap.GetIndex(trackParTh)],
		          						  x[fFitterTrackParMap.GetIndex(trackParPhi)],
		          						  energy,
		          						  extraPars
		          						  );
	  tracksToFit.push_back(track);
  }



  // if(fIsFirstCall ) // TEMP: Print it out all the time for now
  // {
  //   std::cout << "Tracks used for first call:" << std::endl;
  //   std::vector<WCSimLikelihoodTrackBase>::iterator itr = tracksToFit.begin();
  //   for( ; itr < tracksToFit.end(); ++itr)
  //   {
  //     (*itr).Print();
  //   }
  // }

  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(tracksToFit);
  minus2LnL = fTotalLikelihood->Calc2LnL();

  for(size_t iTrack = 0; iTrack < tracksToFit.size(); ++iTrack)
  {
    tracksToFit.at(iTrack)->Print();
    // std::cout << iTrack << "/" << tracksToFit.size() << std::endl;
	  delete (tracksToFit.at(iTrack));
  }
  tracksToFit.clear();
  std::cout << "Function evaluated " << ++fCalls << " times for event " << fEvent << " ... " << " -2Ln(L) = " << minus2LnL << std::endl << std::endl;
  return minus2LnL;
}

double WCSimLikelihoodFitter::WrapFuncAlongTrack(const Double_t * x)
{
  UInt_t nTracks = fFitterConfig->GetNumTracks();

  // If we've fixed some track parameters together, then our array of fit parameters won't just be of size
  // n tracks * m parameters per track, and we need to work out which entry corresponds to which parameter
  // of which track(s) - this is probably expensive, so we'll do it once and look it up in a map thereafter

  std::vector<WCSimLikelihoodTrackBase*> tracksToFit;
  TrackType::Type trackType = fFitterTrackParMap.GetTrackType(0);
  std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes(trackType);
  TrackAndType trackParX = std::make_pair(0, FitterParameterType::kVtxX);
  TrackAndType trackParY = std::make_pair(0, FitterParameterType::kVtxY);
  TrackAndType trackParZ = std::make_pair(0, FitterParameterType::kVtxZ);
  TrackAndType trackParT = std::make_pair(0, FitterParameterType::kVtxT);
  TrackAndType trackParTh = std::make_pair(0, FitterParameterType::kDirTh);
  TrackAndType trackParPhi = std::make_pair(0, FitterParameterType::kDirPhi);
  TrackAndType trackParE = std::make_pair(0, FitterParameterType::kEnergy);

  TVector3 vertex(fFitterTrackParMap.GetCurrentValue(trackParX),
		  	  	  fFitterTrackParMap.GetCurrentValue(trackParY),
		  	  	  fFitterTrackParMap.GetCurrentValue(trackParZ));
  
  std::vector<TVector3> directions;
  std::vector<double> energies;
  for(unsigned int i = 0; i < nTracks; ++i)
  {
    TVector3 directionTmp;
    directionTmp.SetMagThetaPhi(1.0, 
    							fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirTh),
    							fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirPhi)
                               );
    directions.push_back(directionTmp);
    energies.push_back(x[i+1]);
  }
  double time = fFitterTrackParMap.GetCurrentValue(trackParT) + x[0]/30.0 + x[nTracks+1];

  TVector3 direction(0,0,0);
  for(unsigned int i = 0; i < nTracks; ++i)
  {
    direction += (directions.at(i));
  }
  direction *= 1.0 / directions.size();


  TVector3 newVertex = vertex + direction * x[0];
  std::cout << "x = " << x[0] << " So vertex goes from (" << vertex.X() << ", " << vertex.Y() << ", " << vertex.Z() << ") to ("
            << newVertex.X() << ", " << newVertex.Y() << ", " << newVertex.Z() << ")" << std::endl;
  std::cout << "Time is now " << time << " and energy is " << energies.at(0) << std::endl;

  for(unsigned int i  = 0 ; i < nTracks; ++i)
  {
	  WCSimLikelihoodTrackBase * track = WCSimLikelihoodTrackFactory::MakeTrack(
		    	  	  	  	  	  fFitterTrackParMap.GetTrackType(i),
		    	  	  	  	  	  newVertex.X(),
		    	  	  	  	  	  newVertex.Y(),
		    	  	  	  	  	  newVertex.Z(),
                            time,
		    	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirTh),
		    	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirPhi),
		    	  	  	  	  	  energies.at(i)
	  	  	  	  	  	  	 );
	  tracksToFit.push_back(track);
  }

  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(tracksToFit);
  minus2LnL = fTotalLikelihood->Calc2LnL();

  // std::cout << "deleting tracks" << std::endl;
  for(size_t iTrack = 0; iTrack < tracksToFit.size(); ++iTrack)
  {
	  // std::cout << iTrack << "/" << tracksToFit.size() << std::endl;
	  delete (tracksToFit.at(iTrack));
  }
  // std::cout << "Clearing" << std::endl;
  tracksToFit.clear();
  // std::cout << "Done" << std::endl;
  std::cout << "Function evaluated " << ++fCalls << " times for event " << fEvent << " ... " << " -2Ln(L) = " << minus2LnL << std::endl << std::endl;
  return minus2LnL;
}


// A special fit function for neutral pions where we can constrain the energy and direction
// of the two photon tracks so that they have the invariant mass of the parent pi0
double WCSimLikelihoodFitter::WrapFuncPiZero(const Double_t * x)
{
  UInt_t nTracks = fFitterConfig->GetNumTracks();
  assert(nTracks == 2);

  // If we've fixed some track parameters together, then our array of fit parameters won't just be of size
  // n tracks * m parameters per track, and we need to work out which entry corresponds to which parameter
  // of which track(s) - this is probably expensive, so we'll do it once and look it up in a map thereafter

  std::vector<WCSimLikelihoodTrackBase*> tracksToFit;
  
  ////////////////////////////////////////////////////////////
  // First track: everything for this one is a free parameter
  ////////////////////////////////////////////////////////////
  TrackType::Type trackType = fFitterTrackParMap.GetTrackType(0);
  assert(trackType == TrackType::PhotonLike);

  std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes(trackType);
  TrackAndType trackParX = std::make_pair(0, FitterParameterType::kVtxX);
  TrackAndType trackParY = std::make_pair(0, FitterParameterType::kVtxY);
  TrackAndType trackParZ = std::make_pair(0, FitterParameterType::kVtxZ);
  TrackAndType trackParT = std::make_pair(0, FitterParameterType::kVtxT);
  TrackAndType trackParTh = std::make_pair(0, FitterParameterType::kDirTh);
  TrackAndType trackParPhi = std::make_pair(0, FitterParameterType::kDirPhi);
  TrackAndType trackParE = std::make_pair(0, FitterParameterType::kEnergy);

  TrackAndType trackParConversionDistance = std::make_pair(0, FitterParameterType::kConversionDistance);
  std::map<FitterParameterType::Type, double> extraPars;
  extraPars[FitterParameterType::kConversionDistance] = x[fFitterTrackParMap.GetIndex(trackParConversionDistance)];

  WCSimLikelihoodTrackBase * track = WCSimLikelihoodTrackFactory::MakeTrack(
		  	  	  	  	  	  trackType,
		  	  	  	  	  	  x[fFitterTrackParMap.GetIndex(trackParX)],
	          						  x[fFitterTrackParMap.GetIndex(trackParY)],
	          						  x[fFitterTrackParMap.GetIndex(trackParZ)],
	          						  x[fFitterTrackParMap.GetIndex(trackParT)],
	          						  x[fFitterTrackParMap.GetIndex(trackParTh)],
	          						  x[fFitterTrackParMap.GetIndex(trackParPhi)],
	          						  x[fFitterTrackParMap.GetIndex(trackParE)],
	          						  extraPars
	          						  );
	tracksToFit.push_back(track);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Second track: this one has its energy fixed so that the invariant mass of the two-track pair is that of a pi0
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  trackParX = std::make_pair(1, FitterParameterType::kVtxX);
  trackParY = std::make_pair(1, FitterParameterType::kVtxY);
  trackParZ = std::make_pair(1, FitterParameterType::kVtxZ);
  trackParT = std::make_pair(1, FitterParameterType::kVtxT);
  trackParTh = std::make_pair(1, FitterParameterType::kDirTh);
  trackParPhi = std::make_pair(1, FitterParameterType::kDirPhi);
  double track1Energy = GetPiZeroSecondTrackEnergy(x);

  std::map<FitterParameterType::Type, double> extraPars2;
  trackParConversionDistance = std::make_pair(1, FitterParameterType::kConversionDistance);
  extraPars2[FitterParameterType::kConversionDistance] = x[fFitterTrackParMap.GetIndex(trackParConversionDistance)];

  track = WCSimLikelihoodTrackFactory::MakeTrack(
		  	  	  	  	  	  	  	  	  	  	     trackType,
                                                 x[fFitterTrackParMap.GetIndex(trackParX)],
                                                 x[fFitterTrackParMap.GetIndex(trackParY)],
                                                 x[fFitterTrackParMap.GetIndex(trackParZ)],
                                                 x[fFitterTrackParMap.GetIndex(trackParT)],
                                                 x[fFitterTrackParMap.GetIndex(trackParTh)],
                                                 x[fFitterTrackParMap.GetIndex(trackParPhi)],
                                                 track1Energy,
                                                 extraPars2
  	  	  	  	  	  	  	  	  	  	  	  	);
  tracksToFit.push_back(track);


  //if(fIsFirstCall ) // TEMP: Print it out all the time for now
  //{
  //   std::cout << "Tracks used for first call:" << std::endl;
  //   std::vector<WCSimLikelihoodTrackBase*>::iterator itr = tracksToFit.begin();
  //   int i = 0;
  //   for( ; itr != tracksToFit.end(); ++itr)
  //   {
  //  	 ++i;
  //     (*itr)->Print();
  //   }
  //}

  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(tracksToFit);
  minus2LnL = fTotalLikelihood->Calc2LnL();

  // std::cout << "deleting tracks" << std::endl;
  for(size_t iTrack = 0; iTrack < tracksToFit.size(); ++iTrack)
  {
    tracksToFit.at(iTrack)->Print();
    // std::cout << iTrack << "/" << tracksToFit.size() << std::endl;
    delete (tracksToFit.at(iTrack));
  }
  // std::cout << "Clearing" << std::endl;
  tracksToFit.clear();
  // std::cout << "Done" << std::endl;
  std::cout << "Function evaluated " << ++fCalls << " times for event " << fEvent << " ... " << " -2Ln(L) = " << minus2LnL << std::endl << std::endl;
  return minus2LnL;

}

// Use the previous Hough transform-based reconstruction algorithm to seed the multi-track one
void WCSimLikelihoodFitter::SeedParams()
{
	return;
}


std::vector<WCSimLikelihoodTrackBase*> WCSimLikelihoodFitter::GetBestFit()
{
  return fBestFit;
}

Double_t WCSimLikelihoodFitter::GetMinimum()
{
  return fMinimum;
}

void WCSimLikelihoodFitter::FitEventNumber(Int_t iEvent) {
	fIsFirstCall = true;
	SetEvent(iEvent);


	if(fFitterConfig->GetMakeFits() && CanFitEvent()) {
    std::cout << "Fitting event " << iEvent << std::endl;

    if(fFitterConfig->GetIsPiZeroFit())
    {
      FitPiZeroEvent();
    }
    else
    {
  
      fCalls = 0;
  	  fFitterTrackParMap.Set();
  
      // Run seed
      SeedEvent();
  
      // Fit directions
      if(fFitterConfig->GetNumTracks() > 1)
      {
        // Leigh: Fit the time first.
        if(!fFitterTrackParMap.GetIsFixed(0,FitterParameterType::kVtxT)){
          std::cout << "LEIGH: About to fit time for two track event" << std::endl;
          FixVertex();
          FixDirection();
          FixEnergy();
          FreeTime();
          FitTime();
          FreeVertex();
          FreeDirection();
          FreeEnergy();
        }

        FixVertex();
        FreeDirection();
        FixEnergy();
        FitVertex();
        FreeVertex();
        FreeEnergy();
      }
      
      // Now fit the energy
      bool switchBack = false;
      if( WCSimAnalysisConfig::Instance()->GetUseChargeAndTime())
      {
        WCSimAnalysisConfig::Instance()->SetUseChargeOnly();
        switchBack = true;
      }
      //
      //
      std::cout << "Fitting energy only" << std::endl;
      FixVertex();
      FixDirection();
      FixTime();
      if( WCSimAnalysisConfig::Instance()->GetUseChargeAndTime())
      {
        WCSimAnalysisConfig::Instance()->SetUseChargeOnly();
        FitEnergy();
        WCSimAnalysisConfig::Instance()->SetUseChargeAndTime();
      }
      else
      {
        FitEnergy();
      }
      FreeVertex();
      FreeDirection();
      FreeTime();

      // Fix the directions and energy and move the vertex
      // along the (average) momentum vector (helps because the seed
      // often gets this wrong in the direction of the track by 
      // changing the cone angle)
      std::cout << "Fitting along track" << std::endl;
      FixDirection();
      // FixEnergy();
      if(WCSimAnalysisConfig::Instance()->GetUseTimeOnly())
      {
        // MetropolisHastingsAlongTrack(500);
        FreeTime();
        FitAlongTrack();
      }
      else if( WCSimAnalysisConfig::Instance()->GetUseChargeAndTime())
      {
        WCSimAnalysisConfig::Instance()->SetUseChargeOnly();
        FixTime();
        FitAlongTrack();
        FreeTime();
        WCSimAnalysisConfig::Instance()->SetUseChargeAndTime();
        FitAlongTrack();
      }
      else if( WCSimAnalysisConfig::Instance()->GetUseChargeOnly())
      {
    	  FixTime();
    	  FitAlongTrack();
      }
      // MetropolisHastingsAlongTrack(250);
      FreeDirection();
      FreeEnergy();
  
      // Now freely fit the vertex and direction    
      FixEnergy();
      FitVertex();
      FreeEnergy();
  
      // And finally tidy up the energy
      FixVertex();
      FixDirection();
      FreeEnergy();
      if( WCSimAnalysisConfig::Instance()->GetUseChargeAndTime())
      {
    	  WCSimAnalysisConfig::Instance()->SetUseChargeOnly();
		  FitEnergy();
		  WCSimAnalysisConfig::Instance()->SetUseChargeAndTime();
      }
      else
      {
    	  FitEnergy();
      }
      // Leave everything free at the end
      FreeVertex();
      FreeDirection();    
      FreeEnergy();
      //Fit("Simplex");
      // MetropolisHastings(2500);
      //
      if( switchBack )
      {
        std::cout << "Switching back to charge AND time" << std::endl;
        WCSimAnalysisConfig::Instance()->SetUseTimeOnly();
        FixVertex();
        FixDirection();
        FixEnergy();
        FreeTime();
        FitTime();
        WCSimAnalysisConfig::Instance()->SetUseChargeAndTime();
        FreeVertex();
        FreeDirection();
        FitVertex();

        WCSimAnalysisConfig::Instance()->SetUseChargeOnly();
        FreeEnergy();
        FixVertex();
        FixDirection();
        FixTime();
        FitEnergy();
        FreeVertex();
        FreeDirection();
        FreeEnergy();
        FreeTime();
        WCSimAnalysisConfig::Instance()->SetUseChargeAndTime();
      }
      
    }
  /*
  
  
  
  		Minimize2LnL();
  */
		fTrueLikelihoodTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
  }
  std::cout << "Fitted event number " << iEvent << std::endl;
	return;
}


void WCSimLikelihoodFitter::SetFitterPlots(WCSimFitterPlots* fitterPlots) {
  std::cout << "Setting fFitterPlots" << std::endl;
	fFitterPlots = fitterPlots;
  std::cout << "Number of surface bins" << fFitterPlots->GetNumSurfaceBins() << std::endl;
}

void WCSimLikelihoodFitter::Make1DSurface(TrackAndType trackPar) {
  std::cout << " *** WCSimLikelihoodFitter::Make1DSurface() *** " << std::endl;
 
  double min = fFitterTrackParMap.GetMinValue(trackPar);
  double max = fFitterTrackParMap.GetMaxValue(trackPar);

  std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
  Double_t * x = &(startVals[0]);

  for(unsigned int iBin = 1; iBin <= fFitterPlots->GetNumSurfaceBins(); ++iBin)
  {
	  double stepVal = min + (max - min) * iBin / (double)fFitterPlots->GetNumSurfaceBins();
	  startVals.at(fFitterTrackParMap.GetIndex(trackPar)) = stepVal;
    double lnL = WrapFunc(x);
    std::cout << "Bin = " << iBin << "   StepVal = " << stepVal << "   -2LnL = " << lnL <<std::endl;

	  fFitterPlots->Fill1DProfile(trackPar, iBin, lnL);
  }
}

void WCSimLikelihoodFitter::Make2DSurface(
		std::pair<TrackAndType, TrackAndType > trackPar)
{
  std::cout << " *** WCSimLikelihoodFitter::Make1DSurface() *** " << std::endl;
  
  unsigned int arrIndexX = fFitterTrackParMap.GetIndex(trackPar.first);
  double minX = fFitterTrackParMap.GetMinValue(trackPar.first);
  double maxX = fFitterTrackParMap.GetMaxValue(trackPar.first);

  unsigned int arrIndexY = fFitterTrackParMap.GetIndex(trackPar.second);
  double minY = fFitterTrackParMap.GetMinValue(trackPar.second);
  double maxY = fFitterTrackParMap.GetMaxValue(trackPar.second);

  std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
  Double_t * x = &(startVals[0]);


  // Now step over each variable
  for(unsigned int iXBin = 0; iXBin < fFitterPlots->GetNumSurfaceBins(); ++iXBin)
  {
	  double stepValX = minX + (maxX - minX) * iXBin / (double)fFitterPlots->GetNumSurfaceBins();
	  x[arrIndexX] = stepValX;

	  for(unsigned int iYBin = 0; iYBin < fFitterPlots->GetNumSurfaceBins(); ++iYBin)
	  {
		  double stepValY = minY + (maxY - minY) * iYBin / (double)fFitterPlots->GetNumSurfaceBins();
		  x[arrIndexY] = stepValY;
		  fFitterPlots->Fill2DProfile(trackPar, iXBin, iYBin, WrapFunc(x));
	  }
  }
  return;
}

void WCSimLikelihoodFitter::PerformEnergyGridSearch(Double_t& best2LnL,
		std::vector<Double_t>& bestEnergies)
{


}

void WCSimLikelihoodFitter::FillHitComparison() {
  fTotalLikelihood->SetTracks(fBestFit);
  fTotalLikelihood->Calc2LnL();
	std::vector<double> predictedCharges = fTotalLikelihood->GetPredictedChargeVector();
	std::vector<double> measuredCharges = fTotalLikelihood->GetMeasuredChargeVector();
	std::vector<double> best2LnLs = fTotalLikelihood->GetTotal2LnLVector();
	std::vector<double> charge2LnLs = fTotalLikelihood->GetCharge2LnLVector();
	std::vector<double> time2LnLs = fTotalLikelihood->GetTime2LnLVector();
  std::vector<double> predictedTimes = fTotalLikelihood->GetPredictedTimeVector();

	std::vector<WCSimLikelihoodTrackBase*> *correctTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
  
	fTotalLikelihood->SetTracks(*correctTracks);
	fTotalLikelihood->Calc2LnL();
	std::vector<double> correctPredictedCharges = fTotalLikelihood->GetPredictedChargeVector();
	std::vector<double> correctPredictedTimes = fTotalLikelihood->GetPredictedTimeVector();
	std::vector<double> correctCharge2LnLs = fTotalLikelihood->GetCharge2LnLVector();
	std::vector<double> correctTime2LnLs = fTotalLikelihood->GetTime2LnLVector();
	std::vector<double> correct2LnLs = fTotalLikelihood->GetTotal2LnLVector();
	fFitterTree->FillHitComparison(fEvent, fLikelihoodDigitArray, predictedCharges, correctPredictedCharges, measuredCharges, predictedTimes, correctPredictedTimes, best2LnLs, correct2LnLs, charge2LnLs, correctCharge2LnLs, time2LnLs, correctTime2LnLs);

}

WCSimLikelihoodTrackBase * WCSimLikelihoodFitter::RescaleParams(Double_t x, Double_t y, Double_t z, Double_t t,
                                                          Double_t th, Double_t phi, Double_t E, 
                                                          TrackType::Type type)
{
  // Currently do nothing...
  return WCSimLikelihoodTrackFactory::MakeTrack(type, x,y,z,t,th,phi,E);
/*
  // But could do this - scale variables in the range (0,1) up to proper track parameers:
  double pi = TMath::Pi();
  double x2,y2,z2,t2,th2,phi2,E2;
  // x2   = (1900.+x)/3800.;
  // y2   = (1900.+y)/3800.;
  // z2   = (940. +z)/1880.;
  // t2   = t;
  // th2  = th / pi;
  // phi2 = (phi+pi)/(2*pi);
  // E2   = E / 1500;

  // Dirty hack: detector dimensions hardcoded for quick testing
  x2 = 3800.0 * x - 1900.;
  y2 = 3800.0 * y - 1900.;
  z2 = 1880.0 * z - 940.;
  th2 = pi * th;
  std::cout << "theta2 = " << th2 << std::endl;
  phi2 = 2*pi*phi - pi;
  E2 = 1500*E;
  WCSimLikelihoodTrackBase myTrack(x2,y2,z2,t2,th2,phi2,E2,type);

  return myTrack;
*/
}

Int_t WCSimLikelihoodFitter::GetStatus()
{
  return fStatus;
}

void WCSimLikelihoodFitter::RunFits() {
  std::cout << "WCSimLikelihoodFitter::RunFits() " << std::endl;
	UInt_t firstEvent = fFitterConfig->GetFirstEventToFit();
	UInt_t numEventsToFit = fFitterConfig->GetNumEventsToFit();


	for(UInt_t iEvent = firstEvent; iEvent < firstEvent + numEventsToFit ; ++iEvent)
	{
		FitEventNumber(iEvent);
	    if( fMinimum > 0 )
	    {
			  FillPlots();
			  FillTree();
			  FillHitComparison();
	    }
	    else
	    {
	      fFitterTree->FillRecoFailures(iEvent);
	    }
	    fFitterTree->SaveTree();
		fFitterPlots->SavePlots();
	}
}

void WCSimLikelihoodFitter::RunSurfaces() {
	fRootEvent = WCSimInterface::Instance()->GetWCSimEvent(fFitterConfig->GetFirstEventToFit());
	fLikelihoodDigitArray = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(fFitterConfig->GetFirstEventToFit());
  fFitterTrackParMap.Set();
	std::cout << "There are " << fLikelihoodDigitArray->GetNDigits() << " digits" << std::endl;
	if(fTotalLikelihood != NULL)
	{
		delete fTotalLikelihood;
		fTotalLikelihood = NULL;
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}
	else
	{
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}

	
  std::vector<std::pair<unsigned int, FitterParameterType::Type> > all1DProfile = fFitterPlots->GetAll1DSurfaceKeys();
  std::cout << "Size of 1D surface vector = " << all1DProfile.size() << std::endl;
	for(std::vector<std::pair<unsigned int, FitterParameterType::Type> >::iterator itr1D = all1DProfile.begin() ;
		itr1D != all1DProfile.end(); ++itr1D)
	{
		Make1DSurface(*itr1D);
	}

	std::vector<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > > all2DProfile = fFitterPlots->GetAll2DSurfaceKeys();
	for(std::vector<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > >::iterator itr2D = all2DProfile.begin() ;
			itr2D != all2DProfile.end(); ++itr2D)
	{
			Make2DSurface(*itr2D);
	}
	fFitterPlots->SaveProfiles();
	return;
}

void WCSimLikelihoodFitter::ResetEvent() {
    fStatus = -999;
    fRootEvent = NULL;
    fLikelihoodDigitArray = NULL; // WCSimInterface cleans up after itself so no need to delete these here
    if(fTotalLikelihood != NULL)
    {
      delete fTotalLikelihood;
      fTotalLikelihood = NULL;
    }

    fStatus = -999;
    fTypes.clear();

    fParMap[1] = 8; // The number of fit parameters for n tracks, plus 1 for nTracks itself (fixed)
    fParMap[2] = 11;
    fMinimum  = 0;
    fIsFirstCall = true;
}

void WCSimLikelihoodFitter::FillPlots() {
	if(fFitterPlots != NULL)
	{
		fFitterPlots->FillPlots(fBestFit);
		fFitterPlots->FillRecoMinusTrue(fBestFit, fTrueLikelihoodTracks);
	}
	return;
}

void WCSimLikelihoodFitter::FillTree() {
	if(fFitterTree != NULL)
	{
		std::vector<Bool_t> bestFitEscapes;
		std::vector<Bool_t> trueTrackEscapes;
		for(unsigned int i = 0; i < fBestFit.size(); ++i)
		{
			bestFitEscapes.push_back(this->GetFitTrackEscapes(i));
		}
		for(unsigned int j = 0; j < fTrueLikelihoodTracks->size(); ++j)
		{
			trueTrackEscapes.push_back(this->GetTrueTrackEscapes(j));
		}
		fFitterTree->Fill(fEvent, fBestFit, bestFitEscapes, *fTrueLikelihoodTracks, trueTrackEscapes, fMinimum);
	}
	return;
}

void WCSimLikelihoodFitter::SetFitterTree(WCSimFitterTree* fitterTree) {
  std::cout << "Setting fFitterTree" << std::endl;
  fFitterTree = fitterTree;
}

Bool_t WCSimLikelihoodFitter::GetTrueTrackEscapes(unsigned int iTrack) const{
	WCSimLikelihoodTrackBase * track = fTrueLikelihoodTracks->at(iTrack);
  std::cout << "Seeing if true track escapes" << std::endl;
  track->Print();
	return GetTrackEscapes(track);
}

Bool_t WCSimLikelihoodFitter::GetFitTrackEscapes( unsigned int iTrack) const{
	WCSimLikelihoodTrackBase * track = (fBestFit.at(iTrack));
	return GetTrackEscapes(track);
}

Bool_t WCSimLikelihoodFitter::GetTrackEscapes(WCSimLikelihoodTrackBase * track) const{
	Double_t distanceToEdge = WCSimGeometry::Instance()->ForwardProjectionToEdge(track->GetX(), track->GetY(), track->GetZ(),
																	 	 	 	 track->GetDirX(), track->GetDirY(), track->GetDirZ());
  Double_t stoppingDistance = 0.0;
  if( track->GetType() != TrackType::Unknown )
  {
    fTotalLikelihood->GetEmissionProfileManager()->GetStoppingDistance(track);
  }
  else
  {
    std::cerr << "Warning, track was of unknown type" << std::endl;
    std::cerr << "Will assume it was contained" << std::endl; 
  }
	return (stoppingDistance > distanceToEdge);
}

void WCSimLikelihoodFitter::SeedEvent()
{
	std::cout << " *** WCSimLikelihoodFitter::SeedEvent() *** " << std::endl;

  // Run the old Hough transform reco
  WCSimRecoSeed * myReco = dynamic_cast<WCSimRecoSeed*>(WCSimRecoFactory::Instance()->MakeReco("seed")); // This calls new() - delete when done
  myReco->SetNumberOfTracks(fFitterConfig->GetNumTracks());
  WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
  std::vector<WCSimRecoEvent*> slicedEvents = myReco->RunSeed(recoEvent);

  // Make a vector of all of the available rings we have. 
  std::vector<WCSimRecoRing*> ringVec;
  std::vector<std::pair<WCSimRecoRing*,double> > otherRings;
  std::vector<double> ringTime;
  // Add the primary ring from each slice first.
  for(unsigned int e = 0; e < slicedEvents.size(); ++e){
    //ringVec.push_back(slicedEvents[e]->GetRing(0));
    //ringTime.push_back(slicedEvents[e]->GetVtxTime());
    for(int r = 0; r < slicedEvents[e]->GetNRings(); ++r){
      std::pair<WCSimRecoRing*,double> tempPair;
      tempPair = std::make_pair<WCSimRecoRing*,double>(slicedEvents[e]->GetRing(r),slicedEvents[e]->GetVtxTime());
      otherRings.push_back(tempPair);
    }
  }
  // Sort the otherRings vector by height and add to the main vector
  std::sort(otherRings.begin(),otherRings.end(),RingSort);
  for(unsigned int r = 0; r < otherRings.size(); ++r){
    ringVec.push_back(otherRings[r].first);
    ringTime.push_back(otherRings[r].second);
  }

  // Now we need to loop over all the track parameters available to the fitter
  // and set them to the corresponding seed parameter
  for(unsigned int iTrack = 0; iTrack < fFitterConfig->GetNumTracks(); ++iTrack)
  {
    double seedX, seedY, seedZ, seedT;
    double dirX, dirY, dirZ;
    double seedTheta, seedPhi;
    std::cout << "Track number " <<  iTrack << " compared to " << slicedEvents.size() << std::endl;
    if(iTrack < ringVec.size()){
      seedX = ringVec[iTrack]->GetVtxX();
      seedY = ringVec[iTrack]->GetVtxY();
      seedZ = ringVec[iTrack]->GetVtxZ();
      seedT = ringTime[iTrack];
      dirX = ringVec[iTrack]->GetDirX();
      dirY = ringVec[iTrack]->GetDirY();
      dirZ = ringVec[iTrack]->GetDirZ();
    }
    else{
      seedX = ringVec[0]->GetVtxX();
      seedY = ringVec[0]->GetVtxY();
      seedZ = ringVec[0]->GetVtxZ();
      seedT = ringTime[0];
      dirX = ringVec[0]->GetDirX();
      dirY = ringVec[0]->GetDirY();
      dirZ = ringVec[0]->GetDirZ();
    }

    seedTheta = TMath::ACos(dirZ);
    if(dirY != 0.0){ seedPhi = TMath::ATan2(dirY,dirX); }// Ensure range is -pi to pi
    else{ seedPhi = (dirX < 0.0)? 0.5*TMath::Pi() : -0.5*TMath::Pi(); }

    std::cout << "Track seed " << iTrack << ": " << seedX << ", " << seedY << ", " << seedZ << " :: "
                                                 << dirX << ", " << dirY << ", " << dirZ << ", " << seedTheta << ", " << seedPhi << std::endl; 

    // Only see the parameters if they are not requested to be fixed:
    // Also check if any of the parameters are joined with any other tracks: we should
    // let the lower-numbered track win in this case, because it has the higher Hough peak
    WCSimFitterConfig * fitterConfig = fFitterConfig;
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxX)){
      if(fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxX) >= iTrack){
        fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxX, seedX);
      }
    }
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxY)){
      if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxY) >= iTrack){
        fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxY, seedY);
      }
    }
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxZ)){
      if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxZ) >= iTrack){
        fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxZ, seedZ);
      }
    }
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxT)){
      if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxT) >= iTrack){
        fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxT, seedT);
      }
    }
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kDirTh)){
      if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kDirTh) >= iTrack){
        fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kDirTh, seedTheta);
      }
    }
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kDirPhi)){
      if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kDirPhi) >= iTrack){
        fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kDirPhi, seedPhi);
      }
    }
  }
  delete myReco;

  // Need to delete the elements of slicedEvents as they are not destroyed by WCSimRecoSlicer
  for(unsigned int v = 0; v < slicedEvents.size(); ++v){
    delete slicedEvents[v];
  }

}


void WCSimLikelihoodFitter::FixVertex()
{
    fFitterTrackParMap.FixVertex();
}

void WCSimLikelihoodFitter::FreeVertex()
{
	fFitterTrackParMap.FreeVertex();
}

void WCSimLikelihoodFitter::FixDirection()
{
  fFitterTrackParMap.FixDirection();
}

void WCSimLikelihoodFitter::FreeDirection()
{
  fFitterTrackParMap.FreeDirection();
}

void WCSimLikelihoodFitter::FixEnergy()
{
	fFitterTrackParMap.FixEnergy();
}

void WCSimLikelihoodFitter::FreeEnergy()
{
	fFitterTrackParMap.FreeEnergy();
}

void WCSimLikelihoodFitter::FixTime()
{
  fFitterTrackParMap.FixTime();
}

void WCSimLikelihoodFitter::FreeTime()
{
  fFitterTrackParMap.FreeTime();
}

void WCSimLikelihoodFitter::FitEnergyGridSearch()
{
	TrackAndType trackAndEnergy;
	if( fFitterConfig->GetNumTracks() > 2 )
	{
		std::cerr << "Error: Performing a grid search with more than two tracks will take an extremely long time" << std::endl;
		assert(fFitterConfig->GetNumTracks() <= 2);
	}

	// Get energies allowed for the first track
	// We have to make an emission profile using a temporary track set to this energy
	// Then get the energy bins for this track into a vector, and delete the track
	TrackAndType trackZeroEnergy(0, FitterParameterType::kEnergy);
	TrackType::Type trackZeroType = fFitterConfig->GetTrackType(0);

	//std::cout << "Track zero type should be " << TrackType::AsString(trackZeroType) << std::endl;
	WCSimLikelihoodTrackBase * tempTrack = WCSimLikelihoodTrackFactory::MakeTrack(trackZeroType);
	tempTrack->SetE(fFitterTrackParMap.GetMinValue(0, FitterParameterType::kEnergy));

	//std::cout << "Making empty emission profile" << std::endl;
	WCSimEmissionProfiles * tempProfileZero = new WCSimEmissionProfiles(tempTrack->GetType(), tempTrack->GetE());
	std::vector<Double_t> energyBinsZero = tempProfileZero->GetProfileEnergies();

	// Emission profiles can only be used for energies between (not equal to) the first and last bins
	assert(energyBinsZero.size() > 2);
	energyBinsZero.pop_back();
	energyBinsZero.erase(energyBinsZero.begin(), energyBinsZero.begin() + 1);

	delete tempTrack;
	delete tempProfileZero;

	// Now we have to do the same thing to get the energy bins for the second track
	TrackType::Type trackOneType = TrackType::Unknown;
	std::pair<UInt_t, FitterParameterType::Type> trackOneEnergy = std::make_pair(1,FitterParameterType::kEnergy);
	std::vector<Double_t> energyBinsOne;
	if( fFitterConfig->GetNumTracks() > 1)
	{
		trackOneType = fFitterConfig->GetTrackType(1);
		WCSimLikelihoodTrackBase * tempTrack = WCSimLikelihoodTrackFactory::MakeTrack(trackOneType);
		tempTrack->SetE(fFitterConfig->GetParStart(1, "kEnergy"));
		WCSimEmissionProfiles * tempProfileOne = new WCSimEmissionProfiles(tempTrack->GetType(), tempTrack->GetE());
		energyBinsOne = tempProfileOne->GetProfileEnergies();
		assert(energyBinsOne.size() > 2);
		energyBinsOne.pop_back();
		energyBinsOne.erase(energyBinsOne.begin(), energyBinsOne.begin() + 1);
		delete tempTrack;
		delete tempProfileOne;
	}

	// This loop assumes we have at most two tracks
	double best2LnL = -999.9;
	std::vector<double> bestEnergies;

	// Set everything to the default value

	bool isFirstLoop = true;
  for(unsigned int iBin = 0; iBin < energyBinsZero.size(); ++iBin)
	//for(unsigned int iBin = 0; iBin < energyBinsZero.size(); iBin += 4) // LEIGH
	{
		// std::cout << "Energy bin " << iBin << " is " << energyBinsZero.at(iBin) << std::endl;
    // std::cout << "Setting initial energy to " << energyBinsZero.at(iBin) << std::endl;
		fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kEnergy, energyBinsZero.at(iBin));
		// std::cout << "about to do the loop" << std::endl;

		unsigned int jBin = 0;
		do
		{
		  // std::cout << "jBin = " << jBin << std::endl;
		  // std::cout << "Current energy bin " << iBin << "," << jBin << std::endl; // LEIGH
	          // TString toPrint = Form("Current energy bin = %f", energyBinsZero.at(iBin));
			if( fFitterConfig->GetNumTracks() > 1)
			{
				fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kEnergy, energyBinsOne.at(jBin));
        			// toPrint = toPrint + Form(" and %f", energyBinsOne.at(jBin));
			}
			std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
      			for(unsigned int i = 0; i < startVals.size(); ++i)
      			{
       				// std::cout << startVals.at(i) << std::endl;
      			}
			Double_t * x = &(startVals[0]);
     			// std::cout << "GETTING TEMPMIN" << std::endl;
			double tempMin = WrapFunc(x);
      			// toPrint = toPrint + Form(" and -2LnL = %f, c.f. previous best of %f", tempMin, best2LnL);
      			// std::cout << toPrint << std::endl;
			if( tempMin < best2LnL || isFirstLoop )
			{
				best2LnL = tempMin;
				isFirstLoop = false;
				bestEnergies.clear();
				// std::cout << "Pushing back for iBin = " << iBin << std::endl;
				bestEnergies.push_back(energyBinsZero.at(iBin));
				if( fFitterConfig->GetNumTracks() == 2 )
				{
						// std::cout << "And for jBin = " << std::endl;
						bestEnergies.push_back(energyBinsOne.at(jBin));
				}
			}
			jBin++;

			//jBin+=4; // LEIGH
		} while(jBin < energyBinsOne.size());

	}
	// Update everything with the best fit
	fMinimum = best2LnL;
	for(unsigned int iTrack = 0; iTrack < fFitterConfig->GetNumTracks(); ++iTrack)
	{
		fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kEnergy, bestEnergies.at(iTrack));
	}

	std::cout << "Grid search finished" << std::endl;
	std::cout << "Best 2 Ln(L) = " << fMinimum << std::endl;
  UpdateBestFits();

	return;
}

void WCSimLikelihoodFitter::FitEnergy()
{
  Fit("Simplex"); 
}

void WCSimLikelihoodFitter::FitVertex()
{
  Fit("Simplex");
}

void WCSimLikelihoodFitter::FitTime()
{
  Fit("Simplex");
}

void WCSimLikelihoodFitter::Fit(const char * minAlgorithm)
{
	double likelihood = FitAndGetLikelihood(minAlgorithm);
	std::cout << "Best-fit -2Ln(likelihood) was " << likelihood << std::endl;
	return;
}

double WCSimLikelihoodFitter::FitAndGetLikelihood(const char * minAlgorithm)
{
	// Set up the minimizer
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", minAlgorithm);

	// Alternatively: use a different algorithm to check the minimizer works
	// ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "GSLSimAn");


	min->SetMaxFunctionCalls(500);
	min->SetMaxIterations(1);
	min->SetPrintLevel(3);
	//min->SetTolerance(0.);
	min->SetErrorDef(1.0);
	min->SetStrategy(2);
	std::cout << " Tolerance = " << min->Tolerance() << std::endl;

	// Convert nTracks into number of parameters
	const unsigned int nPars = this->GetNPars();

	// Tell the minimizer the function to minimize
	// We have to wrap it in this functor because it's a member function of this class
	ROOT::Math::Functor func(this,&WCSimLikelihoodFitter::WrapFunc, nPars);


	// Tell the minimizer the functor and variables to consider
	min->SetFunction(func);

	// Set the fixed and free variables
	Bool_t isAnythingFree = false;

	std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
	std::vector<double> minVals = fFitterTrackParMap.GetMinValues();
	std::vector<double> maxVals = fFitterTrackParMap.GetMaxValues();
	std::vector<double> steps = fFitterTrackParMap.GetSteps();
	std::vector<bool> currentlyFixed = fFitterTrackParMap.GetCurrentlyFixed();
	std::vector<std::string> names = fFitterTrackParMap.GetNames();
	for(UInt_t i = 0; i < nPars; ++i)
	{
		// Sometimes the fast seed goes outside the allowed range:
		if( currentlyFixed[i] ){
      std::cout << "Fixed: " << names[i] << "  Start: " << startVals[i] << std::endl;
			min->SetFixedVariable(i, names[i], startVals[i]);
		}
		else{
			isAnythingFree = true;
      std::cout << "Free: " << names[i] << "  Start: " << startVals[i] << "   Min: " << minVals[i] << "   Max: " << maxVals[i] << std::endl;
			min->SetLimitedVariable(i, names[i], startVals[i], steps[i], minVals[i], maxVals[i]);
		}

		if(GetUsePiZeroMassConstraint())
		{
			if(i == fFitterTrackParMap.GetIndex(1, FitterParameterType::kEnergy))
			{
				min->SetFixedVariable(i, names[i], startVals[i]);
			}
		}
	}

	// Print the parameters we're using
	for(UInt_t j = 0; j < nPars; ++j)
	{
		std::cout << j << "   " << names[j] << "   " << startVals[j] << "   " << minVals[j] << "   " << maxVals[j] << std::endl;
	}




	  // Perform the minimization
	  if( isAnythingFree )
	  {
	    min->Minimize();
	    fMinimum = min->MinValue();
	    fStatus = min->Status();
	  }
	  // This is what we ought to do:
	  const Double_t * outPar = min->X();
	  // Future tracks need to start from this best fit
	  for( unsigned int i = 0; i < nPars; ++i)
	  {
	    fFitterTrackParMap.SetCurrentValue(i, outPar[i]);

		  if(GetUsePiZeroMassConstraint())
		  {
		  	if(i == fFitterTrackParMap.GetIndex(1, FitterParameterType::kEnergy))
		  	{
		  		fFitterTrackParMap.SetCurrentValue(i, GetPiZeroSecondTrackEnergy(outPar));
		  	}
		  }
	  }

	  // Minuit gives us a const double array, but (for now) we want to do the grid search and then modify it
	  // std::cout << "Here's outPar..." << std::endl;
	  // std::cout << outPar[0] << "  " << outPar[1] << std::endl;

	  // Get and print the fit results
	  // std::cout << "Best fit track: " << std::endl;
	  // RescaleParams(outPar2[0],  outPar2[1],  outPar2[2],  outPar2[3],  outPar2[4],  outPar2[5],  outPar2[6], fType).Print();
    UpdateBestFits();
    return min->MinValue();
}


void WCSimLikelihoodFitter::FitPiZero(const char * minAlgorithm)
{
	// Set up the minimizer
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", minAlgorithm);

	// Alternatively: use a different algorithm to check the minimizer works
	// ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "GSLSimAn");


	min->SetMaxFunctionCalls(500);
	min->SetMaxIterations(1);
	min->SetPrintLevel(3);
	//min->SetTolerance(0.);
  min->SetErrorDef(1.0);
	min->SetStrategy(2);
	std::cout << " Tolerance = " << min->Tolerance() << std::endl;

	// Convert nTracks into number of parameters
	const unsigned int nPars = this->GetNPars(); // We can constrain one of the energies using the pion invariant mass

	// Tell the minimizer the function to minimize
	// We have to wrap it in this functor because it's a member function of this class
	ROOT::Math::Functor func(this,&WCSimLikelihoodFitter::WrapFuncPiZero, nPars);


	// Tell the minimizer the functor and variables to consider
	min->SetFunction(func);

	// Set the fixed and free variables
	Bool_t isAnythingFree = false;

	std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
	std::vector<double> minVals = fFitterTrackParMap.GetMinValues();
	std::vector<double> maxVals = fFitterTrackParMap.GetMaxValues();
	std::vector<double> steps = fFitterTrackParMap.GetSteps();
	std::vector<bool> currentlyFixed = fFitterTrackParMap.GetCurrentlyFixed();
	std::vector<std::string> names = fFitterTrackParMap.GetNames();
	
  UInt_t trk1EnergyIndex = (fFitterTrackParMap.GetIndex(1, FitterParameterType::kEnergy));

	for(UInt_t i = 0; i < nPars; ++i)
	{

		// Sometimes the fast seed goes outside the allowed range:
		if( currentlyFixed[i] ){ 
      // The fitter shouldn't vary track 1's energy because the pion mass constrains it
      // But we still need to set it so the indices don't get off by one from fFitterTrackParMap
      // So we'll fix it, the fit function won't use it, and we can update it manually at the end
      //
      // std::cout << "Fixed: " << names[i] << "  Start: " << startVals[i] << std::endl;
			min->SetFixedVariable(i, names[i], startVals[i]);
		}
		else{
			isAnythingFree = true;
      // std::cout << "Free: " << names[i] << "  Start: " << startVals[i] << "   Min: " << minVals[i] << "   Max: " << maxVals[i] << std::endl;
			min->SetLimitedVariable(i, names[i], startVals[i], steps[i], minVals[i], maxVals[i]);
		}
	}

	  // Print the parameters we're using
	  for(UInt_t j = 0; j < nPars; ++j)
	  {
      std::string isFixedStr = currentlyFixed[j] ? "Fixed" : "Free";
	  	std::cout << j << "   " << names[j] << "   " << isFixedStr << "   " << startVals[j] << "   " << minVals[j] << "   " << maxVals[j] << std::endl;
	  }

	  // Perform the minimization
	  if( isAnythingFree )
	  {
	    min->Minimize();
	    fMinimum = min->MinValue();
	    fStatus = min->Status();
	  }
	  // This is what we ought to do:
	  const Double_t * outPar = min->X();
	  // Future tracks need to start from this best fit
	  for( unsigned int i = 0; i < nPars; ++i)
	  {
	    fFitterTrackParMap.SetCurrentValue(i, outPar[i]);

	    if(i == trk1EnergyIndex )
	    {
	    	fFitterTrackParMap.SetCurrentValue(i, GetPiZeroSecondTrackEnergy(outPar));
	    }
	  }

	  // Minuit gives us a const double array, but (for now) we want to do the grid search and then modify it
	  // std::cout << "Here's outPar..." << std::endl;
	  // std::cout << outPar[0] << "  " << outPar[1] << std::endl;

	  // Get and print the fit results
	  // std::cout << "Best fit track: " << std::endl;
	  // RescaleParams(outPar2[0],  outPar2[1],  outPar2[2],  outPar2[3],  outPar2[4],  outPar2[5],  outPar2[6], fType).Print();
    UpdateBestFits();
	  return;
}

void WCSimLikelihoodFitter::UpdateBestFits()
{
	  // Now save the best-fit tracks
	  fBestFit.clear();
	  for(UInt_t iTrack = 0; iTrack < fFitterConfig->GetNumTracks(); ++iTrack)
	  {
		  std::cout << "Track number " << iTrack << std::endl;
		  TrackAndType trackParX = std::make_pair(iTrack, FitterParameterType::kVtxX);
		  TrackAndType trackParY = std::make_pair(iTrack, FitterParameterType::kVtxY);
		  TrackAndType trackParZ = std::make_pair(iTrack, FitterParameterType::kVtxZ);
		  TrackAndType trackParT = std::make_pair(iTrack, FitterParameterType::kVtxT);
		  TrackAndType trackParTh = std::make_pair(iTrack, FitterParameterType::kDirTh);
		  TrackAndType trackParPhi = std::make_pair(iTrack, FitterParameterType::kDirPhi);
		  TrackAndType trackParE = std::make_pair(iTrack, FitterParameterType::kEnergy);
		  TrackAndType trackParConversionDistance = std::make_pair(iTrack, FitterParameterType::kConversionDistance);
		  std::map<FitterParameterType::Type, Double_t> extraPars;
		  extraPars[FitterParameterType::kConversionDistance] = fFitterTrackParMap.GetCurrentValue(trackParConversionDistance);

		  std::cout << "Type is " << fFitterTrackParMap.GetTrackType(iTrack) << std::endl;
		  WCSimLikelihoodTrackBase * track = WCSimLikelihoodTrackFactory::MakeTrack(
				  	  fFitterTrackParMap.GetTrackType(iTrack),
				  	  fFitterTrackParMap.GetCurrentValue(trackParX),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParY),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParZ),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParT),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParTh),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParPhi),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParE),
                      extraPars);
		  std::cout << "Best-fit track number " << iTrack << std::endl;
		  track->Print();
		  fBestFit.push_back(track);
	  }
}



void WCSimLikelihoodFitter::MetropolisHastings(const int nTries)
{
	std::vector<std::vector<double> > markovChain;

	std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
	std::vector<double> minVals = fFitterTrackParMap.GetMinValues();
	std::vector<double> maxVals = fFitterTrackParMap.GetMaxValues();
	std::vector<double> steps = fFitterTrackParMap.GetSteps();
	std::vector<bool> currentlyFixed = fFitterTrackParMap.GetCurrentlyFixed();

	markovChain.push_back(startVals);
	double * startX = &(startVals[0]);
	double old2LnL = WrapFunc(startX);


	int iTry = 0;
	TRandom3 ran;
	int accepted = 0;
	while(iTry < nTries)
	{
		std::cout << "Metropolis-Hastings iteration " << iTry << std::endl;
		std::vector<double> oldParams = markovChain.at(markovChain.size() - 1);
		std::vector<double> newParams;
		for(unsigned int iParam = 0; iParam < oldParams.size(); ++iParam)
		{
			if(currentlyFixed.at(iParam))
			{
				newParams.push_back(oldParams.at(iParam));
			}
		  else
		  {
			  double newParam = ran.Gaus(oldParams.at(iParam), steps.at(iParam));
			  newParams.push_back(newParam);			
		  }
		}
		double * x = &(newParams[0]);
		double new2LnL = WrapFunc(x);

		if( new2LnL < old2LnL ) // they're -ve log likelihoods
		{
			markovChain.push_back(newParams);
			accepted++;
		}
		else
		{
			double uniform = ran.Rndm();
			if(uniform < TMath::Exp(0.5*(old2LnL - new2LnL)))
			{
				markovChain.push_back(newParams);
				accepted++;
			}
			else 
			{
				markovChain.push_back(oldParams);
			}
		}
		++iTry;
	}
	std::vector<std::vector<double> > cleaned;
	for(unsigned int i = markovChain.size()/5 + 1; i < markovChain.size(); ++i)
	{
		cleaned.push_back(markovChain.at(i));
	}
	

	for(unsigned int iParam = 0; iParam < cleaned.at(0).size(); ++iParam)
	{
	  double runningTotal = 0.0;
	  int numToAverage = 0;
		for(unsigned int iChain = 0; iChain < cleaned.size(); ++iChain)
		{
			runningTotal += markovChain.at(iChain).at(iParam);
			numToAverage++;
		}
		assert(numToAverage);
		fFitterTrackParMap.SetCurrentValue(iParam, runningTotal/numToAverage);
	}
    	UpdateBestFits();
	return;
}

void WCSimLikelihoodFitter::MetropolisHastingsAlongTrack(const int nTries)
{
	const unsigned int nTracks = fFitterConfig->GetNumTracks();
	const double cInCmPerNs = 29.9792458;  // Speed of light

	if(nTracks != 1) { return;}

	// If we don't have only one track, make sure they share a vertex
	double sharedVertex = true;
	if(nTracks > 1)
	{
		sharedVertex = true;
		for(unsigned int i = 1; i < nTracks; ++i)
		{
			if(   fFitterConfig->GetTrackIsJoinedWith(i, "kVtxX") != 0
					|| fFitterConfig->GetTrackIsJoinedWith(i, "kVtxY") != 0
					|| fFitterConfig->GetTrackIsJoinedWith(i, "kVtxZ") != 0 )
			{
				sharedVertex = false;
			}
		}
	}

	if(nTracks == 1 || (nTracks > 1 && sharedVertex) )
	{
		// Are we allowed to move the vertex freely?
		bool canMove = (   !fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxX)
				&& !fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxY)
				&& !fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxZ)
		);

		// Fill vectors of each track's seed direction, and allowed energies
		std::vector<TVector3> directions;
		std::vector<double> minE, maxE, currE, stepE;
		std::vector<bool> fixE;
		for(unsigned int i = 0; i < nTracks; ++i)
		{
			TVector3 direction;
			direction.SetMagThetaPhi(1.0,
					fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirTh),
					fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirPhi)
			);
			directions.push_back(direction);

			minE.push_back(fFitterTrackParMap.GetMinValue(i, FitterParameterType::kEnergy));
			maxE.push_back(fFitterTrackParMap.GetMaxValue(i, FitterParameterType::kEnergy));
			currE.push_back(fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kEnergy));
			stepE.push_back(fFitterTrackParMap.GetStep(i, FitterParameterType::kEnergy));
			fixE.push_back(fFitterTrackParMap.GetIsFixed(i, FitterParameterType::kEnergy));
		}

		// Average the directions of all of the tracks.  We're not weighting by energy here, although perhaps we could
		TVector3 direction(0,0,0);
		for(unsigned int i = 0; i < nTracks; ++i)
		{
			direction += (directions.at(i));
		}
		direction *= 1.0 / directions.size();

		// Now find the allowed range for the vertex time
		double minT, maxT, currT, stepT;
		bool fixT;
		minT = fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxT);
		maxT = fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxT);
		currT = fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxT);
		stepT = 0.5; // ns
		fixT = (fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxT) || !WCSimAnalysisConfig::Instance()->GetUseTime());


		// Get the allowed ranges for the vertex positions
		TVector3 startVertex(fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxX),
				fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxY),
				fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxZ)
		);

		TVector3 mins(fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxX),
				fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxY),
				fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxZ));

		TVector3 maxes(fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxX),
				fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxY),
				fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxZ));

		// To work out the maximum and minimum distance we can move the vertex, we need to
		// figure out where the track would escape the detector
		std::vector<double> escapeDists;
		for(int i = 0; i < 3; ++i)
		{
			if(fabs(direction[i]) > 1e-6)
			{
				escapeDists.push_back( (maxes[i] - startVertex[i])/direction[i] );
				escapeDists.push_back( (mins[i] - startVertex[i])/direction[i] );
			}
		}
		assert(escapeDists.size() > 1);
		std::sort(escapeDists.begin(), escapeDists.end());
		double maxVal = static_cast<unsigned int>(-1);
		double minVal = -1.0 * maxVal;
		std::cout << "minVal = " << minVal << " and start maxVal = " << maxVal << std::endl;
		for(unsigned int i = 0; i < escapeDists.size(); ++i)
		{
			if(escapeDists[i] < 0 && escapeDists[i] > minVal){ minVal = escapeDists[i]; }
			if(escapeDists[i] > 0 && escapeDists[i] < maxVal){ maxVal = escapeDists[i]; }
		}
		std::cout << "MinVal = " << minVal << " and maxVal = " << maxVal << std::endl;


		// Calculate the step size by projecting the (x,y,z) step sizes along the average direction
		TVector3 stepSizes(fFitterTrackParMap.GetStep(0, FitterParameterType::kVtxX),
				fFitterTrackParMap.GetStep(0, FitterParameterType::kVtxY),
				fFitterTrackParMap.GetStep(0, FitterParameterType::kVtxZ)
		);
		double stepSize = sqrt(pow((stepSizes.X() * direction.X()), 2) + pow((stepSizes.Y() * direction.Y()), 2) + pow((stepSizes.Z() * direction.Z()), 2));
    stepSize = 50;
		std::cout << "Step size = " << stepSize << std::endl;
		//////////////////////////////////////////////////////


		TrackType::Type trackType = fFitterTrackParMap.GetTrackType(0);
		std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes(trackType);
		TrackAndType trackParX = std::make_pair(0, FitterParameterType::kVtxX);
		TrackAndType trackParY = std::make_pair(0, FitterParameterType::kVtxY);
		TrackAndType trackParZ = std::make_pair(0, FitterParameterType::kVtxZ);
		TrackAndType trackParT = std::make_pair(0, FitterParameterType::kVtxT);
		TrackAndType trackParTh = std::make_pair(0, FitterParameterType::kDirTh);
		TrackAndType trackParPhi = std::make_pair(0, FitterParameterType::kDirPhi);
		TrackAndType trackParE = std::make_pair(0, FitterParameterType::kEnergy);


		// Now do the actual Markov chain stepping
		std::vector<std::vector<double> > markovChain;

		std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
		std::vector<double> minVals = fFitterTrackParMap.GetMinValues();
		std::vector<double> maxVals = fFitterTrackParMap.GetMaxValues();
		std::vector<double> steps = fFitterTrackParMap.GetSteps();
		std::vector<bool> currentlyFixed = fFitterTrackParMap.GetCurrentlyFixed();

		// Build a vector of the initial parameters to be the first entry in the chain
		std::vector<double> startParams;
		startParams.push_back(0.0); // Vertex shift
		startParams.push_back(0.0); // Vertex time shift
		for(size_t iTrack = 0; iTrack < currE.size(); ++iTrack)
		{
			startParams.push_back(currE.at(iTrack));  // Energy of each track
		}
		markovChain.push_back(startParams);


		double * startX = &(startVals[0]);
		double old2LnL = WrapFunc(startX);


		int iTry = 0;
		TRandom3 ran;
		int accepted = 0;
		while(iTry < nTries)
		{
			std::cout << "Metropolis-Hastings along track iteration " << iTry << std::endl;
			std::vector<double> oldParams = markovChain.back();

			// Vary the distance
			double newDistance = ran.Gaus(oldParams.at(0), stepSize);
			TVector3 newVertex = startVertex + newDistance * direction;

			// Vary the time - compensate for the shift in vertex

			double adjustedTime = currT - newDistance/cInCmPerNs; // Purely accounting for the change in vertex

			double timeShift = ran.Gaus(oldParams.at(1), stepT); // Also shift the time by a random amount
			double newTime = adjustedTime + timeShift;

			// Vary the energy
			std::vector<double> newEnergies;
			for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
			{
				if(!fixE.at(iTrack))
				{
					newEnergies.push_back(ran.Gaus(oldParams.at(iTrack+2), stepE.at(iTrack)));
				}
			}

			// Make a single vector of all the track params to pass to WrapFunc
			std::vector<double> newParams;
			for(unsigned int i = 0; i < startVals.size(); ++i)
			{
				newParams.push_back(startVals.at(i));
			}
			newParams.at(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxX)) = newVertex.X();
			newParams.at(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxY)) = newVertex.Y();
			newParams.at(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxZ)) = newVertex.Z();
			newParams.at(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxT)) = newTime;

			// Only set the new energies if they're free
			for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
			{
				if(!fixE.at(iTrack))
				{
					newParams.at(fFitterTrackParMap.GetIndex(iTrack, FitterParameterType::kEnergy)) = newEnergies.at(iTrack);
				}
			}

			// Work out the likelihood
			double * x = &(newParams[0]);
			double new2LnL = WrapFunc(x);

			// Put the adjust parameters into a vector that we can append to the Markov chain
			std::vector<double> shiftedParams;
			shiftedParams.push_back(newDistance);
			shiftedParams.push_back(timeShift);
			for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
			{

				if(!fixE.at(iTrack))
				{

					shiftedParams.push_back(newParams.at(fFitterTrackParMap.GetIndex(iTrack, FitterParameterType::kEnergy)));
				}
			}
			std::cout << "Shifted params: " << " x = " << newDistance << "   time = " << timeShift << "   so new time = " << newTime << std::endl;

			// Keep the new parameters if they improve the log-likelihood
			if( new2LnL < old2LnL ) // they're -ve log likelihoods
			{
				markovChain.push_back(shiftedParams);
				accepted++;
			}
			else
			{
				double uniform = ran.Rndm();
				// Or keep them if they get lucky
				if(uniform < TMath::Exp(0.5*(old2LnL - new2LnL)))
				{
					markovChain.push_back(shiftedParams);
					accepted++;
				}
				// Otherwise retain the old ones
				else
				{
					markovChain.push_back(oldParams);
				}
			}
			++iTry;
		}


		// Give the chain some time to settle by getting rid of early entries
		// I'm only throwing away the first 1/4 because the seed should be right-ish to begin with
		std::vector<std::vector<double> > cleaned;
		for(unsigned int i = markovChain.size()/4 + 1; i < markovChain.size(); ++i)
		{
			cleaned.push_back(markovChain.at(i));
		}

		// Average the remaining entries
		std::vector<double> runningTotals;
		double numToAverage = 0.0;
		for(unsigned int iChain = 0; iChain < cleaned.size(); iChain += 5)
		{
			for(unsigned int iVariable = 0; iVariable < cleaned.at(iChain).size(); ++iVariable)
			{
				if(iVariable == runningTotals.size())
				{
					runningTotals.push_back(0.0);
				}
				runningTotals.at(iVariable) += markovChain.at(iChain).at(iVariable);
			}
			numToAverage++;

		}
		assert(numToAverage);

		double finalDistance = runningTotals.at(0)/numToAverage;
		double finalTimeShift = runningTotals.at(1)/numToAverage;

		TVector3 newVertex = startVertex + finalDistance * direction;
		double newTime = currT - finalDistance/cInCmPerNs + finalTimeShift;
		std::cout << "Shift track vertex by " << finalDistance << std::endl;
		fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxX), newVertex.X());
		fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxY), newVertex.Y());
		fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxZ), newVertex.Z());
		fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxT), newTime);

		for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
		{
			if(!fixE.at(iTrack))
			{
				fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(iTrack, FitterParameterType::kEnergy), runningTotals.at(2+iTrack));
			}
		}
		UpdateBestFits();
	}
	return;
}

void WCSimLikelihoodFitter::FitAlongTrack()
{
  const unsigned int nTracks = fFitterConfig->GetNumTracks();

  // Do all the tracks share the same vertex
  double sharedVertex = true;
  if(nTracks > 1)
  {
    sharedVertex = true;
    for(unsigned int i = 1; i < nTracks; ++i)
    {
      if(   fFitterConfig->GetTrackIsJoinedWith(i, "kVtxX") != 0
         || fFitterConfig->GetTrackIsJoinedWith(i, "kVtxY") != 0
         || fFitterConfig->GetTrackIsJoinedWith(i, "kVtxZ") != 0 )
      {
        sharedVertex = false;
      }
    }
  }

  if(nTracks == 1 || (nTracks > 1 && sharedVertex) )
  {

    bool canMove = (   !fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxX)
                    && !fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxY)
                    && !fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxZ)
                   );
      
    std::vector<TVector3> directions;
    std::vector<double> minE, maxE, currE, stepE;
    std::vector<bool> fixE;
    for(unsigned int i = 0; i < nTracks; ++i)
    {
      TVector3 direction;
      direction.SetMagThetaPhi(1.0, 
                               fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirTh), 
                               fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirPhi)
                              );
      directions.push_back(direction);

      minE.push_back(fFitterTrackParMap.GetMinValue(i, FitterParameterType::kEnergy));
      maxE.push_back(fFitterTrackParMap.GetMaxValue(i, FitterParameterType::kEnergy));
      currE.push_back(fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kEnergy));
      stepE.push_back(fFitterTrackParMap.GetStep(i, FitterParameterType::kEnergy));
      fixE.push_back(fFitterTrackParMap.GetIsFixed(i, FitterParameterType::kEnergy));
    }

    TVector3 direction(0,0,0);
    for(unsigned int i = 0; i < nTracks; ++i)
    {
      direction += (directions.at(i));
    }
    direction *= 1.0 / directions.size();
  
    double minT, maxT, currT, stepT;
    bool fixT;
    minT = fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxT);
    maxT = fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxT);
    currT = fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxT);
    stepT = fFitterTrackParMap.GetStep(0, FitterParameterType::kVtxT);
    fixT = (fFitterTrackParMap.GetIsFixed(0, FitterParameterType::kVtxT) || !WCSimAnalysisConfig::Instance()->GetUseTime());

    TVector3 startVertex(fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxX),
                         fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxY),
                         fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxZ)
                         );

    TVector3 stepSizes(fFitterTrackParMap.GetStep(0, FitterParameterType::kVtxX),
                       fFitterTrackParMap.GetStep(0, FitterParameterType::kVtxY),
                       fFitterTrackParMap.GetStep(0, FitterParameterType::kVtxZ)
                       );

    TVector3 mins(fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxX),
                  fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxY),
                  fFitterTrackParMap.GetMinValue(0, FitterParameterType::kVtxZ));
    
    TVector3 maxes(fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxX),
                   fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxY),
                   fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kVtxZ));

    std::vector<double> escapeDists;
    for(int i = 0; i < 3; ++i)
    {
      if(fabs(direction[i]) > 1e-6)
      {
        escapeDists.push_back( (maxes[i] - startVertex[i])/direction[i] );
        escapeDists.push_back( (mins[i] - startVertex[i])/direction[i] );
      }
    }
    assert(escapeDists.size() > 1);
    std::sort(escapeDists.begin(), escapeDists.end());
    double maxVal = static_cast<unsigned int>(-1);
    double minVal = -1.0 * maxVal;
    std::cout << "minVal = " << minVal << " and start maxVal = " << maxVal << std::endl;
    for(unsigned int i = 0; i < escapeDists.size(); ++i)
    {
      if(escapeDists[i] < 0 && escapeDists[i] > minVal){ minVal = escapeDists[i]; }
      if(escapeDists[i] > 0 && escapeDists[i] < maxVal){ maxVal = escapeDists[i]; }
    }
    std::cout << "MinVal = " << minVal << " and maxVal = " << maxVal << std::endl;

    double stepSize = sqrt(pow((stepSizes.X() * direction.X()), 2) + pow((stepSizes.Y() * direction.Y()), 2) + pow((stepSizes.Z() * direction.Z()), 2));
    stepSize = 250;
    std::cout << "Step size = " << stepSize << std::endl;
 

	  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

	  // Alternatively: use a different algorithm to check the minimizer works
	  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "GSLSimAn");


	  min->SetMaxFunctionCalls(2500);
	  min->SetMaxIterations(10);
	  min->SetPrintLevel(3);
	  //min->SetTolerance(0.);
    min->SetErrorDef(1.0);
	  min->SetStrategy(2);
	  std::cout << " Tolerance = " << min->Tolerance() << std::endl;

	  // Convert nTracks into number of parameters
	  unsigned int nParsNoConst = 1;
    // nParsNoConst += (int)(fixT);
    for(unsigned int i = 0; i < fixE.size(); ++i)
    {
      nParsNoConst =+ (int)(fixE.at(i));
    }
    const unsigned int nPars = nParsNoConst;


	  // Tell the minimizer the function to minimize
	  // We have to wrap it in this functor because it's a member function of this class
	  ROOT::Math::Functor func(this,&WCSimLikelihoodFitter::WrapFuncAlongTrack, nPars);


	  // Tell the minimizer the functor and variables to consider
	  min->SetFunction(func);
    if( canMove )
    {
	    min->SetLimitedVariable(0, "Distance from seed", minVal, stepSize, minVal, maxVal);
    }
    else
    {
      min->SetFixedVariable(0, "Distance from seed", 0);
    }
    for(unsigned int i = 0;  i < minE.size(); ++i)
    {
      TString name("Track %d energy", i);
      if(fixE.at(i) == false)
      {
        min->SetLimitedVariable(i+1, "Energy", currE.at(i), stepE.at(i), minE.at(i), maxE.at(i));
      }
      else
      {
        min->SetFixedVariable(i+1, "Energy", currE.at(i));
      }
    }
    
    if(fixT == false && WCSimAnalysisConfig::Instance()->GetUseTime())
    {
      min->SetLimitedVariable(minE.size()+1, "Time shift", 0, 1., -3, 3);
    }
    else
    {
      min->SetFixedVariable(minE.size()+1, "Time shift",0);
    }
    min->Minimize();

    double finalDistance = min->X()[0];
    TVector3 newVertex = startVertex + finalDistance * direction;
    std::cout << "Shift track vertex by " << finalDistance << std::endl;
    std::cout << "NewVertex = " << std::endl;
    newVertex.Print();
    fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxX), newVertex.X());
    fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxY), newVertex.Y());
    fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxZ), newVertex.Z());
    fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(0, FitterParameterType::kVtxT), fFitterTrackParMap.GetCurrentValue(0, FitterParameterType::kVtxT) + finalDistance/30.0 + min->X()[minE.size()+1]);
    for(unsigned int iTrack = 0; iTrack < minE.size(); ++iTrack)
    {
      fFitterTrackParMap.SetCurrentValue(fFitterTrackParMap.GetIndex(iTrack, FitterParameterType::kEnergy), min->X()[iTrack+1]);
    }
    UpdateBestFits();
  }
  else{ assert(0); }
	return;
}

void WCSimLikelihoodFitter::FitPiZeroEvent()
{
    fCalls = 0;
	  fFitterTrackParMap.Set();

    // Run seed
    SeedEvent();

    FixVertex();
    FreeEnergy();
    FreeDirection();
    FreeConversionLength();
    FitPiZero();

    FreeVertex();
    

   //  // Now fit the energy and conversion length
   //  FixVertex();
   //  FixDirection();
   //  FitPiZero();
   //  FreeVertex();
   //  FreeDirection();

    // Now freely fit the vertex and direction    
    FixEnergy();
    //FixConversionLength();
    FitPiZero();
    FreeEnergy();

    // And finally tidy up the energy
    FixVertex();
    FixDirection();
    FreeEnergy();
    FitPiZero();

    // Leave everything free at the end
    FreeVertex();
    FreeDirection();    
    FreeEnergy();
    FreeConversionLength();
    Fit();

}

void WCSimLikelihoodFitter::FixConversionLength()
{
  fFitterTrackParMap.FixConversionLength();
}

void WCSimLikelihoodFitter::FreeConversionLength()
{
  fFitterTrackParMap.FreeConversionLength();
}

Bool_t WCSimLikelihoodFitter::GetUsePiZeroMassConstraint() const
{
	return (fFitterConfig->GetNumTracks() == 2 && fFitterConfig->GetForcePiZeroMass());
}

// Get the energy of photon track 1 given the energy of photon track 0
// assuming they arise from the decay of a pi zero (and so have its invariant
// mass)
// \param x The array of current values that will be passed to/returned from WrapFuncPiZero
//          (which comes from the FitterTrackParMap)
Double_t WCSimLikelihoodFitter::GetPiZeroSecondTrackEnergy(const Double_t * x)
{

  int track0ThIndex = fFitterTrackParMap.GetIndex(std::make_pair(0, FitterParameterType::kDirTh));
  int track0PhiIndex = fFitterTrackParMap.GetIndex(std::make_pair(0, FitterParameterType::kDirPhi));
  int track1ThIndex = fFitterTrackParMap.GetIndex(std::make_pair(1, FitterParameterType::kDirTh));
  int track1PhiIndex = fFitterTrackParMap.GetIndex(std::make_pair(1, FitterParameterType::kDirPhi));

  TVector3 track0Dir, track1Dir;
  track0Dir.SetMagThetaPhi(1, x[track0ThIndex], x[track0PhiIndex]);
  track1Dir.SetMagThetaPhi(1, x[track1ThIndex], x[track1PhiIndex]);

  double track0Energy = x[fFitterTrackParMap.GetIndex(std::make_pair(0, FitterParameterType::kEnergy))];
  double massOfPiZeroMeV = 134.9766; // PDG 2014: http://pdg.lbl.gov/2014/listings/rpp2014-list-pi-zero.pdf

  double cosTheta01 = track0Dir.Dot(track1Dir);
                      // This is the dot product of unit vectors in the directions of track 0 and track 1
                      // i.e. the cosine of the angle between the two photons

  double track1Energy = 0.5 * massOfPiZeroMeV * massOfPiZeroMeV / (track0Energy * (1.0 - cosTheta01));
  if(track1Energy > 3000) { track1Energy = 3000; }

  std::cout << "Track 0 energy = " << track0Energy << "   Track 1 energy = " << track1Energy << "   Inv. mass = " << sqrt((2*track0Energy*track1Energy)*(1-cosTheta01)) << std::endl;
  return track1Energy;

}

void WCSimLikelihoodFitter::SetEvent(Int_t iEvent)
{
	fEvent = iEvent;
	WCSimInterface::Instance()->BuildEvent(iEvent);
	fRootEvent = WCSimInterface::Instance()->GetWCSimEvent(iEvent);
	fLikelihoodDigitArray = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(iEvent);
	std::cout << "There are " << fLikelihoodDigitArray->GetNDigits() << " digits" << std::endl;
	if(fTotalLikelihood != NULL)
	{
		delete fTotalLikelihood;
		fTotalLikelihood = NULL;
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}
	else
	{
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}
}

Bool_t WCSimLikelihoodFitter::CanFitEvent() const
{
  return (fLikelihoodDigitArray->GetNHits() > 0);
}
