#include "WCSimChargeLikelihood.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimAnalysisConfig.hh"

#include <iostream>
#include <vector>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimTotalLikelihood)
#endif

///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimTotalLikelihood::WCSimTotalLikelihood( WCSimLikelihoodDigitArray * myLikelihoodDigitArray ) : 
  fLikelihoodDigitArray(myLikelihoodDigitArray)
{  
  fChargeLikelihoodVector.push_back(
      WCSimChargeLikelihood(fLikelihoodDigitArray) );

  //FIXME new time likelihood will only ask for predicted charge later
  fTimeLikelihood = WCSimTimeLikelihood(
      fLikelihoodDigitArray, &(fChargeLikelihoodVector[0]) );
}

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimTotalLikelihood::~WCSimTotalLikelihood()
{
  fChargeLikelihoodVector.clear();
}

///////////////////////////////////////////////////////////////////////////
//  Add some tracks
///////////////////////////////////////////////////////////////////////////
void WCSimTotalLikelihood::SetTracks( std::vector<WCSimLikelihoodTrack> &myTracks)
{
  this->ResetTracks();
  fTracks = myTracks;

  std::vector<WCSimLikelihoodTrack *> myTrackPtrs;
  std::vector<WCSimLikelihoodTrack>::iterator trackIter;
  //FIXME This is a really ugly solution to the new idea
  //  of one Charge Likelihood object for one track...
  int i = 0;
  for(trackIter = fTracks.begin(); trackIter != fTracks.end(); ++trackIter)
  {
    std::cout << "TrackIter -> " << std::endl;
    myTrackPtrs.push_back(&(*trackIter)); 
    if (i > 0) {
      fChargeLikelihoodVector.push_back(
          WCSimChargeLikelihood(fLikelihoodDigitArray) );
    }
    fChargeLikelihoodVector[i].AddTrack(&(*trackIter));
    i++;
  }
  fTimeLikelihood.SetTracks(myTrackPtrs);

  return;
}

///////////////////////////////////////////////////////////////////////////
//  Clear current tracks and remove from the charge and time likelihoods
///////////////////////////////////////////////////////////////////////////
void WCSimTotalLikelihood::ResetTracks()
{
  for(unsigned int i = 0; i < fChargeLikelihoodVector.size(); i++) {
    fChargeLikelihoodVector[i].ClearTracks();
  }

  fTimeLikelihood.ClearTracks();
  fTracks.clear();
}

//////////////////////////////////////////////////////////////////////////
//  Calculate -2 Ln(likelihood) for the two components and combine them
/////////////////////////////////////////////////////////////////////////
Double_t WCSimTotalLikelihood::Calc2LnL()
{
  // nb. at the moment we're just adding these together
  // may have to account for correlations somewhere
  Double_t likelihood;

  //TODO: Loop over digits
  //for (iDigit in fLikelihoodDigitArray) {
    for (unsigned int iTrack; iTrack < fTracks.size(); iTrack++) {
      //TODO: calculate for one digit only
      likelihood += fChargeLikelihoodVector[iTrack].Calc2LnL();
    }
    if(WCSimAnalysisConfig::Instance()->GetUseChargeAndTime())
    {
      //TODO: calculate for one digit only
      //TODO: use charge predicted by charge likelihood!
      likelihood += fTimeLikelihood.Calc2LnL();
    }
  //} //for iDigit

  return likelihood;
}
