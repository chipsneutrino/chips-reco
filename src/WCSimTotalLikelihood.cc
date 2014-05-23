#include "WCSimChargeLikelihood.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"

#include <iostream>
#include <vector>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimTotalLikelihood)
#endif

///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimTotalLikelihood::WCSimTotalLikelihood( WCSimLikelihoodDigitArray * myLikelihoodDigitArray ) : 
  fLikelihoodDigitArray(myLikelihoodDigitArray), fChargeLikelihood(myLikelihoodDigitArray)
{  
  // Put this in the initialization list when the time part is ready 
  fTimeLikelihood( myLikelihoodDigitArray, fChargeLikelihood ); 
}

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimTotalLikelihood::~WCSimTotalLikelihood()
{
  //shouldn't we destroy the likelihood 
  //      arrray created in the initialisation list?
  //      (btw WCSimLikelihoodDigitArray doesn't have a copy constructor)
  //delete fLikelihoodDigitArray;
}

///////////////////////////////////////////////////////////////////////////
//  Add some tracks
///////////////////////////////////////////////////////////////////////////
void WCSimTotalLikelihood::SetTracks( std::vector<WCSimLikelihoodTrack> myTracks)
{
  this->ResetTracks();
  fTracks = myTracks;

  std::vector<WCSimLikelihoodTrack *> myTrackPtrs;
  std::vector<WCSimLikelihoodTrack>::iterator trackIter;
  for(trackIter = fTracks.begin(); trackIter != fTracks.end(); ++trackIter)
  {
    std::cout << "TrackIter -> " << std::endl;
    myTrackPtrs.push_back(&(*trackIter)); 
  }
  fChargeLikelihood.SetTracks(myTrackPtrs);
  myTimeLikelihood.SetTracks(myTrackPtrs);

  return;
}

///////////////////////////////////////////////////////////////////////////
//  Clear current tracks and remove from the charge and time likelihoods
///////////////////////////////////////////////////////////////////////////
void WCSimTotalLikelihood::ResetTracks()
{
  fTracks.clear();
  fChargeLikelihood.ClearTracks();
  fTimeLikelihood.ClearTracks();
}

//////////////////////////////////////////////////////////////////////////
//  Calculate -2 Ln(likelihood) for the two components and combine them
/////////////////////////////////////////////////////////////////////////
Double_t WCSimTotalLikelihood::Calc2LnL()
{
  // nb. at the moment we're just adding these together
  // may have to account for correlations somewhere
  return fChargeLikelihood.Calc2LnL() + fTimeLikelihood.Calc2LnL();
}
