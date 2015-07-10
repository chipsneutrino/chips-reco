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
#include "TMath.h"
#include "TMinuit.h"
#include "TStopwatch.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodFitter)
#endif

/**
 * @todo Remove hardcoding of track type
 */
WCSimLikelihoodFitter::WCSimLikelihoodFitter()
{
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
  unsigned int nPars = WCSimFitterConfig::Instance()->GetNumIndependentParameters();
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
  std::cout << "Function evaluated " << ++fCalls << " times for event " << fEvent << std::endl;
  UInt_t nTracks = WCSimFitterConfig::Instance()->GetNumTracks();

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

	  WCSimLikelihoodTrackBase * track = WCSimLikelihoodTrackFactory::MakeTrack(
			  	  	  	  	  	  fFitterTrackParMap.GetTrackType(iTrack),
			  	  	  	  	  	  x[fFitterTrackParMap.GetIndex(trackParX)],
		          						  x[fFitterTrackParMap.GetIndex(trackParY)],
		          						  x[fFitterTrackParMap.GetIndex(trackParZ)],
		          						  x[fFitterTrackParMap.GetIndex(trackParT)],
		          						  x[fFitterTrackParMap.GetIndex(trackParTh)],
		          						  x[fFitterTrackParMap.GetIndex(trackParPhi)],
		          						  x[fFitterTrackParMap.GetIndex(trackParE)]
		          						  );
	  tracksToFit.push_back(track);
    std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << "  " << TrackType::AsString(track->GetType()) << std::endl;
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

  std::cout << "deleting tracks" << std::endl;
  for(size_t iTrack = 0; iTrack < tracksToFit.size(); ++iTrack)
  {
    std::cout << iTrack << "/" << tracksToFit.size() << std::endl;
	  delete (tracksToFit.at(iTrack));
  }
  std::cout << "Clearing" << std::endl;
  tracksToFit.clear();
  std::cout << "Done" << std::endl;
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
	fEvent = iEvent;
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


	if(WCSimFitterInterface::Instance()->GetMakeFits() ) {
		WCSimInterface::Instance()->BuildEvent(iEvent);
    std::cout << "Fitting event " << iEvent << std::endl;
    
	  fFitterTrackParMap.Set();
    SeedEvent();
    
    FixVertex();
    FitEnergy();
    
    FixEnergy();
    FreeVertex();
    FitVertex();
    
    FixVertex();
    FreeEnergy();
    FitEnergy();

/*



		Minimize2LnL();
*/
		fTrueLikelihoodTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
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
	}

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

	std::vector<WCSimLikelihoodTrackBase*> *correctTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
  
	fTotalLikelihood->SetTracks(*correctTracks);
	fTotalLikelihood->Calc2LnL();
	std::vector<double> correctPredictedCharges = fTotalLikelihood->GetPredictedChargeVector();
	std::vector<double> correct2LnLs = fTotalLikelihood->GetTotal2LnLVector();
	fFitterTree->FillHitComparison(fEvent, fLikelihoodDigitArray, predictedCharges, correctPredictedCharges, measuredCharges, best2LnLs, correct2LnLs);


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

	UInt_t firstEvent = WCSimFitterConfig::Instance()->GetFirstEventToFit();
	UInt_t numEventsToFit = WCSimFitterConfig::Instance()->GetNumEventsToFit();


	for(UInt_t iEvent = firstEvent; iEvent < firstEvent + WCSimFitterConfig::Instance()->GetNumEventsToFit() ; ++iEvent)
	{
		FitEventNumber(iEvent);
		WCSimFitterInterface::Instance()->SaveResults();
	}
}

void WCSimLikelihoodFitter::RunSurfaces() {
	fRootEvent = WCSimInterface::Instance()->GetWCSimEvent(WCSimFitterConfig::Instance()->GetFirstEventToFit());
	fLikelihoodDigitArray = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(WCSimFitterConfig::Instance()->GetFirstEventToFit());
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
	WCSimFitterInterface::Instance()->SaveProfiles();
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
	  WCSimEmissionProfiles ep(track);
	  stoppingDistance = ep.GetStoppingDistance(track);
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
	WCSimReco * myReco = WCSimRecoFactory::Instance()->MakeReco(); // This calls new() - delete when done
	WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
	// myReco->RunFilter(recoEvent); // Might need to do this first, I'm not sure
	myReco->RunRecoRings(recoEvent);


	// Set the seed vertex coordinates using the Hough fit
	double seedX = recoEvent->GetVtxX();
	double seedY = recoEvent->GetVtxY();
	double seedZ = recoEvent->GetVtxZ();

	// Convert the direction into (theta, phi) coordinates
	Double_t recoDirX = recoEvent->GetDirX();
	Double_t recoDirY = recoEvent->GetDirY();
	Double_t recoDirZ = recoEvent->GetDirZ();
	Double_t seedTheta   = TMath::ACos(recoDirZ);
	Double_t seedPhi     = -999.9;
	if(recoDirY != 0.0){ seedPhi = TMath::ATan2(recoDirY,recoDirX); }// Ensure range is -pi to pi
	else{ seedPhi = (recoDirX < 0.0)? 0.5*TMath::Pi() : -0.5*TMath::Pi(); }

	// Now we need to loop over all the track parameters available to the fitter
	// and set them to the corresponding seed parameter
	for(int iTrack = 0; iTrack < WCSimFitterInterface::Instance()->GetNumTracks(); ++iTrack)
	{
    // Only see the parameters if they are not requested to be fixed:
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxX)){
		  fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxX, seedX);
    }
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxY)){
		  fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxY, seedY);
    }
    if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxZ)){
		  fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxZ, seedZ);
    }
    if(iTrack==0 && !fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kDirTh)){
		  fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kDirTh, seedTheta);
    }
    if(iTrack==0 && !fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kDirPhi)){
		  fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kDirPhi, seedPhi);
    }
	}
	delete myReco;
}


void WCSimLikelihoodFitter::FixVertex()
{
    fFitterTrackParMap.FixVertex();
}

void WCSimLikelihoodFitter::FreeVertex()
{
	fFitterTrackParMap.FreeVertex();
}

void WCSimLikelihoodFitter::FixEnergy()
{
	fFitterTrackParMap.FixEnergy();
}

void WCSimLikelihoodFitter::FreeEnergy()
{
	fFitterTrackParMap.FreeEnergy();
}

void WCSimLikelihoodFitter::FitEnergy()
{
  std::cout << "FITTING ENERGY" << std::endl;
	TrackAndType trackAndEnergy;
	if( WCSimFitterConfig::Instance()->GetNumTracks() > 2 )
	{
		std::cerr << "Error: Performing a grid search with more than two tracks will take an extremely long time" << std::endl;
		assert(WCSimFitterConfig::Instance()->GetNumTracks() <= 2);
	}

	// Get energies allowed for the first track
	// We have to make an emission profile using a temporary track set to this energy
	// Then get the energy bins for this track into a vector, and delete the track
	TrackAndType trackZeroEnergy(0, FitterParameterType::kEnergy);
	TrackType::Type trackZeroType = WCSimFitterConfig::Instance()->GetTrackType(0);

	std::cout << "Track zero type should be " << TrackType::AsString(trackZeroType) << std::endl;
	WCSimLikelihoodTrackBase * tempTrack = WCSimLikelihoodTrackFactory::MakeTrack(trackZeroType);
	tempTrack->SetE(fFitterTrackParMap.GetMinValue(0, FitterParameterType::kEnergy));

	std::cout << "Making empty emission profile" << std::endl;
	WCSimEmissionProfiles * tempProfileZero = new WCSimEmissionProfiles();
	std::cout << "Setting track 0" << std::endl;
	tempProfileZero->SetTrack(tempTrack);

	std::cout << "Getting profile energies" << std::endl;
	std::vector<Double_t> energyBinsZero = tempProfileZero->GetProfileEnergies();

	// Emission profiles can only be used for energies between (not equal to) the first and last bins
	assert(energyBinsZero.size() > 2);
	for(unsigned int i = 0; i < energyBinsZero.size(); ++i)
	{
		std::cout << "Energy bin " << i << " = " << energyBinsZero.at(i) << std::endl;
	}
	energyBinsZero.pop_back();
	energyBinsZero.erase(energyBinsZero.begin(), energyBinsZero.begin() + 1);

	std::cout << "after removing first and last..." << std::endl;
	for(unsigned int i = 0; i < energyBinsZero.size(); ++i)
	{
		std::cout << "Energy bin " << i << " = " << energyBinsZero.at(i) << std::endl;
	}
	delete tempTrack;
	delete tempProfileZero;

	// Now we have to do the same thing to get the energy bins for the second track
	TrackType::Type trackOneType = TrackType::Unknown;
	std::pair<UInt_t, FitterParameterType::Type> trackOneEnergy = std::make_pair(1,FitterParameterType::kEnergy);
	std::vector<Double_t> energyBinsOne;
	if( WCSimFitterConfig::Instance()->GetNumTracks() > 1)
	{
		trackOneType = WCSimFitterConfig::Instance()->GetTrackType(1);
		WCSimLikelihoodTrackBase * tempTrack = WCSimLikelihoodTrackFactory::MakeTrack(trackOneType);
		tempTrack->SetE(WCSimFitterConfig::Instance()->GetParStart(1, "kEnergy"));
		WCSimEmissionProfiles * tempProfileOne = new WCSimEmissionProfiles();
		tempProfileOne->SetTrack(tempTrack);
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
//	for(unsigned int iBin = 0; iBin < energyBinsZero.size(); iBin += 4) // LEIGH
	{
		std::cout << "Energy bin " << iBin << " is " << energyBinsZero.at(iBin) << std::endl;
    std::cout << "Setting initial energy to " << energyBinsZero.at(iBin) << std::endl;
		fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kEnergy, energyBinsZero.at(iBin));
		std::cout << "about to do the loop" << std::endl;

		unsigned int jBin = 0;
		do
		{
//		  std::cout << "jBin = " << jBin << std::endl;
		  std::cout << "Current energy bin " << iBin << "," << jBin << std::endl; // LEIGH
			if( WCSimFitterConfig::Instance()->GetNumTracks() > 1)
			{
				fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kEnergy, energyBinsOne.at(jBin));
			}
			std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
      for(int i = 0; i < startVals.size(); ++i)
      {
        std::cout << startVals.at(i) << std::endl;
      }
			Double_t * x = &(startVals[0]);
      std::cout << "GETTING TEMPMIN" << std::endl;
			double tempMin = WrapFunc(x);
      std::cout << "Energy = " << energyBinsZero.at(iBin) << "    -2LnL = " << tempMin << "   Best (till now) = " << best2LnL << std::endl;
			if( tempMin < best2LnL || isFirstLoop )
			{
				best2LnL = tempMin;
				isFirstLoop = false;
				bestEnergies.clear();
				std::cout << "Pushing back for iBin = " << iBin << std::endl;
				bestEnergies.push_back(energyBinsZero.at(iBin));
				if( WCSimFitterConfig::Instance()->GetNumTracks() == 2 )
				{
						std::cout << "And for jBin = " << std::endl;
						bestEnergies.push_back(energyBinsOne.at(jBin));
				}
			}
			jBin++;
//			jBin+=4; // LEIGH
		} while(jBin < energyBinsOne.size());


	}


	// Update everything with the best fit
	fMinimum = best2LnL;
	for(int iTrack = 0; iTrack < WCSimFitterInterface::Instance()->GetNumTracks(); ++iTrack)
	{
		fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kEnergy, bestEnergies.at(iTrack));
	}

	std::cout << "Grid search finished" << std::endl;
	std::cout << "Best 2 Ln(L) = " << fMinimum << std::endl;
  UpdateBestFits();

	return;
}

void WCSimLikelihoodFitter::FitVertex()
{
	// Set up the minimizer
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

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
			min->SetFixedVariable(i, names[i], startVals[i]);
		}
		else{
			isAnythingFree = true;
			min->SetLimitedVariable(i, names[i], startVals[i], steps[i], minVals[i], maxVals[i]);
		}
	}

	// Print the parameters we're using
	for(UInt_t j = 0; j < nPars; ++j)
	{
		std::cout << j << "   " << names[j] << "   " << startVals[j] << "   " << minVals[j] << "   " << maxVals[j] << std::endl;
	}




	  // Perform the minimization
	  fCalls = 0;
	  if( isAnythingFree )
	  {
  Double_t * x = &(startVals[0]);
  std::cout << "First" << std::endl;
  WrapFunc(x);
  std::cout << "Again" << std::endl;
  WrapFunc(x);
  std::cout << "Now min" << std::endl;
	    min->Minimize();
	    fMinimum = min->MinValue();
	    fStatus = min->Status();
	  }
	  // This is what we ought to do:
	  const Double_t * outPar = min->X();
	  // Future tracks need to start from this best fit
	  for( int i = 0; i < nPars; ++i)
	  {
	    fFitterTrackParMap.SetCurrentValue(i, outPar[i]);
	  }

	  // Minuit gives us a const double array, but (for now) we want to do the grid search and then modify it
	  std::cout << "Here's outPar..." << std::endl;
	  std::cout << outPar[0] << "  " << outPar[1] << std::endl;

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
	  for(UInt_t iTrack = 0; iTrack < WCSimFitterConfig::Instance()->GetNumTracks(); ++iTrack)
	  {
		  TrackAndType trackParX = std::make_pair(iTrack, FitterParameterType::kVtxX);
		  TrackAndType trackParY = std::make_pair(iTrack, FitterParameterType::kVtxY);
		  TrackAndType trackParZ = std::make_pair(iTrack, FitterParameterType::kVtxZ);
		  TrackAndType trackParT = std::make_pair(iTrack, FitterParameterType::kVtxT);
		  TrackAndType trackParTh = std::make_pair(iTrack, FitterParameterType::kDirTh);
		  TrackAndType trackParPhi = std::make_pair(iTrack, FitterParameterType::kDirPhi);
		  TrackAndType trackParE = std::make_pair(iTrack, FitterParameterType::kEnergy);

		  WCSimLikelihoodTrackBase * track = WCSimLikelihoodTrackFactory::MakeTrack(
		  	  	  	      fFitterTrackParMap.GetTrackType(iTrack),
	  	  	  	  	    fFitterTrackParMap.GetCurrentValue(trackParX),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParY),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParZ),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParT),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParTh),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParPhi),
  	  	  	  	  	  fFitterTrackParMap.GetCurrentValue(trackParE));
		  std::cout << "Best-fit track number " << iTrack << std::endl;
		  track->Print();
		  fBestFit.push_back(track);
	  }
}
