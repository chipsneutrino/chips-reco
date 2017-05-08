/*
 *
 *  Created on: 21 Aug 2015
 *      Author: ajperch
 */
#include "WCSimAnalysisConfig.hh"
#include "WCSimChargePredictor.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimFitterInterface.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimOutputTree.hh"
#include "WCSimPiZeroSeed.hh"
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
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TMinuit.h"
#include "TPrincipal.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <string>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include <map>
#include <algorithm>

#include "WCSimPiZeroFitter.hh"


#ifndef REFLEX_DICTIONARY
ClassImp(WCSimPiZeroFitter)
#endif

/**
 * @todo Remove hardcoding of track type
 */
WCSimPiZeroFitter::WCSimPiZeroFitter(WCSimFitterConfig * config) : WCSimLikelihoodFitter(config), fPiZeroSeedGenerator(config)
{
  fFitterPlots = NULL;
  fTotalLikelihood = NULL;
  fRootEvent = NULL;
  fLikelihoodDigitArray = NULL;
  fTrueLikelihoodTracks = NULL;
  fCalls = 0;

  ResetEvent();
}

WCSimPiZeroFitter::~WCSimPiZeroFitter()
{
}

void WCSimPiZeroFitter::FitAfterFixingDirectionAndEnergy(WCSimPiZeroSeed * seed)
{
	SetStartingTracks(seed);

	// FreeVertex();
	// FreeEnergy();
	// FreeDirection();

	// FixVertex();
	FixDirection();
	FixEnergy();
	// FreeEnergy();
	FitPiZero("Simplex");

	FreeDirection();
	FreeVertex();
	FreeEnergy();
	Fit("Simplex");
}

void WCSimPiZeroFitter::FitAfterFixingEnergy(WCSimPiZeroSeed * seed)
{
  return;
	SetStartingTracks(seed);
	FreeVertex();
	FreeDirection();

	FixEnergy();
	FitPiZero("Simplex");
	FreeEnergy();

	FitPiZero("Simplex");
}


void WCSimPiZeroFitter::SetStartingTracks(WCSimPiZeroSeed * seed)
{

    // Set first track starting values
    if(CanSetParam(0, FitterParameterType::kVtxX)){
        fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxX, seed->GetTrack1()->GetX());
    }
    if(CanSetParam(0, FitterParameterType::kVtxY)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxY, seed->GetTrack1()->GetY());
    }
    if(CanSetParam(0, FitterParameterType::kVtxZ)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxZ, seed->GetTrack1()->GetZ());
    }
    if(CanSetParam(0, FitterParameterType::kVtxT)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxT, seed->GetTrack1()->GetT());
    }
    if(CanSetParam(0, FitterParameterType::kDirTh)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kDirTh, seed->GetTrack1()->GetTheta());
    }
    if(CanSetParam(0, FitterParameterType::kDirPhi)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kDirPhi, seed->GetTrack1()->GetPhi());
    }
    if(CanSetParam(0, FitterParameterType::kEnergy)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kEnergy, seed->GetTrack1()->GetE());
    }
    if(CanSetParam(0, FitterParameterType::kConversionDistance)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kConversionDistance, seed->GetTrack1()->GetConversionDistance());
    }

    // Set second track starting values
    if(CanSetParam(1, FitterParameterType::kVtxX)){
        fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxX, seed->GetTrack2()->GetX());
    }
    if(CanSetParam(1, FitterParameterType::kVtxY)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxY, seed->GetTrack2()->GetY());
    }
    if(CanSetParam(1, FitterParameterType::kVtxZ)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxZ, seed->GetTrack2()->GetZ());
    }
    if(CanSetParam(1, FitterParameterType::kVtxT)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxT, seed->GetTrack2()->GetT());
    }
    if(CanSetParam(1, FitterParameterType::kDirTh)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirTh, seed->GetTrack2()->GetTheta());
    }
    if(CanSetParam(1, FitterParameterType::kDirPhi)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirPhi, seed->GetTrack2()->GetPhi());
    }
    if(CanSetParam(1, FitterParameterType::kEnergy)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kEnergy, seed->GetTrack2()->GetE());
    }
    if(CanSetParam(1, FitterParameterType::kConversionDistance)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kConversionDistance, seed->GetTrack2()->GetConversionDistance());
    }

    return;

}

bool WCSimPiZeroFitter::CanSetParam(const unsigned int &iTrack, const FitterParameterType::Type &type)
{
    // Only see the parameters if they are not requested to be fixed:
    // Also check if any of the parameters are joined with any other tracks: we should
    // let the lower-numbered track win in this case, because it has the higher Hough peak
	return ((!fFitterTrackParMap.GetIsFixed(iTrack, type))
			 && (fFitterConfig->GetTrackIsJoinedWith(iTrack, type) >= iTrack));

}

UInt_t WCSimPiZeroFitter::GetNPars()
{
  // Do we know how to fit this number of tracks?
  unsigned int nPars = fFitterConfig->GetNumIndependentParameters();
  return nPars;
}

void WCSimPiZeroFitter::RunFits()
{
  std::cout << "WCSimLikelihoodFitter::RunFits() " << std::endl;

	UInt_t firstEvent = fFitterConfig->GetFirstEventToFit();
	UInt_t numEventsToFit = fFitterConfig->GetNumEventsToFit();


	for(UInt_t iEvent = firstEvent; iEvent < firstEvent + numEventsToFit ; ++iEvent)
	{
		FitEventNumber(iEvent);
		fOutputTree->SaveTree();
		fFitterPlots->SavePlots();
	}

}

void WCSimPiZeroFitter::FitEventNumber(Int_t iEvent) {
	std::cout << "WCSimPiZeroFitter::FitEventNumber " << iEvent << std::endl;
	fIsFirstCall = true;
	fEvent = iEvent;
	fCalls = 0;

	//////////////////////////////////////////////////////////
	// Build all the combinations of tracks to try as seeds:
	//////////////////////////////////////////////////////////
	std::vector<WCSimPiZeroSeed*> allSeedTracks = this->GetSeeds();
	std::vector<WCSimPiZeroSeed*>::iterator seedTrackItr = allSeedTracks.begin();

	// We try two minimiser algorithms:
	// 1. Fix the energy and direction of both tracks, run the fit to find a minimum
	//    Then free everything and re-run the fit, starting at this minimum
	// 2. As above, except it only fixes the energies and leaves the directions free
    SetEvent(iEvent);
	fFitterTrackParMap.Set();

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


    assert(allSeedTracks.size() > 0);
    assert(fLikelihoodDigitArray->GetNDigits() > 0);

    // Temp storage of seed variables for outputTree..
    std::vector<WCSimLikelihoodTrackBase*> fSeeds;

	while(seedTrackItr != allSeedTracks.end())
	{
		//Add the two photon tracks to the fSeeds vector!!!
		WCSimPiZeroSeed * seed = *seedTrackItr;
		fSeeds.push_back(WCSimLikelihoodTrackFactory::MakeTrack(TrackType::PhotonLike,
				  	  	  	  	  	  	  	  	  	  	  	  	  seed->GetTrack1()->GetX(), seed->GetTrack1()->GetY(), seed->GetTrack1()->GetZ(),
																  seed->GetTrack1()->GetT(),
																  seed->GetTrack1()->GetTheta(), seed->GetTrack1()->GetPhi(),
																  seed->GetTrack1()->GetE(), seed->GetTrack1()->GetConversionDistance()));

		fSeeds.push_back(WCSimLikelihoodTrackFactory::MakeTrack(TrackType::PhotonLike,
				  	  	  	  	  	  	  	  	  	  	  	  	  seed->GetTrack2()->GetX(), seed->GetTrack2()->GetY(), seed->GetTrack2()->GetZ(),
																  seed->GetTrack2()->GetT(),
																  seed->GetTrack2()->GetTheta(), seed->GetTrack2()->GetPhi(),
																  seed->GetTrack2()->GetE(), seed->GetTrack2()->GetConversionDistance()));

		// These will update fMinimum and the fBestFits if the minimiser improves on the previous best
		FitAfterFixingDirectionAndEnergy(*seedTrackItr);
		FitAfterFixingEnergy(*seedTrackItr);
		++seedTrackItr;
	}

	// Need to fill the SeedInfo in the outputTree...
  if(fOutputTree != NULL){
	  fOutputTree->SetSeed(fSeeds);
  }

	fTrueLikelihoodTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();

  if( fMinimum > 0 )
  {
	  FillPlots();
	  FillTree();
  }

  std::cout << "Fitted event number " << iEvent << std::endl;
  return;
}

std::vector<WCSimPiZeroSeed*> WCSimPiZeroFitter::GetSeeds()
{
	return fPiZeroSeedGenerator.GetSeeds(fEvent);
}

std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase*> > WCSimPiZeroFitter::SeedCheat()
{
  std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase*> > allSeeds;

  double vtxX = 40.79;
  double vtxY = 20.16;
  double vtxZ = -31.82;
  double vtxT = 907.00;

  std::map<FitterParameterType::Type, double> map50;
  std::map<FitterParameterType::Type, double> map250;
  map50[FitterParameterType::kConversionDistance] = 50.0;
  map250[FitterParameterType::kConversionDistance] = 250.0;

  WCSimLikelihoodTrackBase * track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 541.052, map50);
  WCSimLikelihoodTrackBase * track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9172, 0.5348, 221.180, map50);
  allSeeds.push_back(std::make_pair(track1, track2));

  track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 528.678, map250);
  track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9955, 0.5412, 352.048, map250);
  allSeeds.push_back(std::make_pair(track1, track2));

  return allSeeds;
  track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 541.052, map250);
  track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9172, 0.5348, 221.180, map250);

  allSeeds.push_back(std::make_pair(track1, track2));
  track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 528.678, map250);
  track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9955, 0.5412, 352.048, map50);
  allSeeds.push_back(std::make_pair(track1, track2));

  track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 541.052, map250);
  track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9172, 0.5348, 221.180, map50);
  allSeeds.push_back(std::make_pair(track1, track2));

  track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 528.678, map50);
  track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9955, 0.5412, 352.047, map250);
  allSeeds.push_back(std::make_pair(track1, track2));

  track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 541.052, map50);
  track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9172, 0.5348, 221.180, map250);
  allSeeds.push_back(std::make_pair(track1, track2));

  track1 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    2.3085, 0.5737, 528.678, map50);
  track2 = WCSimLikelihoodTrackFactory::MakeTrack(
    TrackType::PhotonLike,
    vtxX, vtxY, vtxZ, vtxT,
    1.9955, 0.5412, 352.047, map50);
  allSeeds.push_back(std::make_pair(track1, track2));

  return allSeeds;


}

void WCSimPiZeroFitter::SeedEvent()
{

	// First we fit a single electron track using the Hough transform.  This gets us a vertex
	// and a direction to use as the jumping off point to generate the seeds for the fitter
	bool cheat = false;
	std::vector<WCSimPiZeroSeed* > allSeedTracks;
	if(!cheat){
		allSeedTracks = fPiZeroSeedGenerator.GetSeeds(fEvent);
		/* |-> Fixes the conversion distance of each track to a pre-chosen value, and tries all combinations
		 	|-> Each combination iterates over the direction of the single track seed by perturbing the direction slightly
		     |--> For each perturbation of the single track, does a grid search to find the best
		 	 	  directions for the second track and saves one where the tracks are close in energy and
				  one where they're separate
		  Returns a vector of all the track pairs we want to try for the seed
		*/


		////////////////////////////////////////////////////////////////////////////////////
		// Now try all the seed tracks one at a time in the fitter, and pick the best result
		////////////////////////////////////////////////////////////////////////////////////
	}

	std::vector<WCSimPiZeroSeed*>::iterator seedTrackItr;
	seedTrackItr = allSeedTracks.begin();
	std::cout << "There are " << allSeedTracks.size() << " seeds" << std::endl;
	for(; seedTrackItr != allSeedTracks.end(); ++seedTrackItr)
	{
		std::cout << "Tracks used for seeding: " << std::endl;
		(*seedTrackItr)->Print();
	}
	return;
}
