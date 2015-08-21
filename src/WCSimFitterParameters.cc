/*
 * WCSimFitterParameters.cc
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#include "WCSimFitterParameters.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrackParameterEnums.hh"
#include "TMath.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>


////////////////////////////////////////////////////////////////////
// FITTER PARAMETER
////////////////////////////////////////////////////////////////////

WCSimFitterParameter::WCSimFitterParameter() : fType(FitterParameterType::kUnknown), fIsFixed(0), fStart(0), fMin(0), fMax(0){
}

WCSimFitterParameter::WCSimFitterParameter(FitterParameterType::Type type,
		bool isFixed, double start, double min, double max) :
				fType(type), fIsFixed(isFixed), fStart(start),
				fMin(min), fMax(max){
	return;
}

WCSimFitterParameter::~WCSimFitterParameter() {
}

void WCSimFitterParameter::Print() {
  std::cout << "Fitter parameter information: " << std::endl;
  std::cout << "fType = " << FitterParameterType::AsString(fType) << std::endl;
  std::cout << "fMin = " << fMin << std::endl;
  std::cout << "fMax = " << fMax << std::endl;
  std::cout << "fStart = " << fStart << std::endl;
  std::cout << "fIsFixed = " << fIsFixed << std::endl << std::endl;
}
////////////////////////////////////////////////////////////////////
// SINGLE TRACK FITTER PARAMETERS
////////////////////////////////////////////////////////////////////

WCSimFitterSingleTrackParameters::WCSimFitterSingleTrackParameters() {
	SetDefaultParameters();
}

void WCSimFitterSingleTrackParameters::SetDefaultParameters(){
	WCSimFitterParameter parX(FitterParameterType::kVtxX, false, 0, 0, 2000);
	WCSimFitterParameter parY(FitterParameterType::kVtxY, false, 0, 0, 2000);
	WCSimFitterParameter parZ(FitterParameterType::kVtxZ, false, 0, 0, 2000);
	WCSimFitterParameter parT(FitterParameterType::kVtxZ, false, 0, 0, 2000);
	WCSimFitterParameter parTh(FitterParameterType::kDirTh, false, 0, 0, TMath::Pi());
	WCSimFitterParameter parPhi(FitterParameterType::kDirPhi, false, 0, 0, 2.0 * TMath::Pi());
	WCSimFitterParameter parE(FitterParameterType::kEnergy, false, 1500, 1000, 2000);
  WCSimFitterParameter parConv(FitterParameterType::kConversionDistance, true, 0, 0, 250);

  fParameters.clear();
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kVtxX,   parX               ) );
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kVtxY,   parY               ) );
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kVtxZ,   parZ               ) );
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kVtxT,   parT               ) );
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kDirTh,  parTh              ) );
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kDirPhi, parPhi             ) );
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kEnergy, parE               ) );
  fParameters.insert( std::map< FitterParameterType::Type, WCSimFitterParameter >::value_type ( FitterParameterType::kConversionDistance, parConv) );

}

WCSimFitterSingleTrackParameters::~WCSimFitterSingleTrackParameters() {
}

WCSimFitterParameter WCSimFitterSingleTrackParameters::GetParameter(
		FitterParameterType::Type type) {
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  return (*itr).second;
}
	
bool WCSimFitterSingleTrackParameters::GetParIsFixed(FitterParameterType type)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  return (*itr).second.GetIsFixed();
}

double WCSimFitterSingleTrackParameters::GetParMax(FitterParameterType type)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  return (*itr).second.GetMax();
}

double WCSimFitterSingleTrackParameters::GetParMin(FitterParameterType type)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  return (*itr).second.GetMin();
}

double WCSimFitterSingleTrackParameters::GetParStart(FitterParameterType type)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  return (*itr).second.GetStart();
}

void WCSimFitterSingleTrackParameters::SetParameter(
		FitterParameterType::Type type, bool isFixed, double start,
		double min, double max) {
	WCSimFitterParameter par(type, isFixed, start, min, max);
	fParameters[type] = par;
}

void WCSimFitterSingleTrackParameters::SetParMin(FitterParameterType::Type type, double min)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  (*itr).second.SetMin(min);
}

void WCSimFitterSingleTrackParameters::SetParMax(FitterParameterType::Type type, double max)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  (*itr).second.SetMax(max);
}

void WCSimFitterSingleTrackParameters::SetParStart(FitterParameterType::Type type, double start)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  (*itr).second.SetStart(start);
}

void WCSimFitterSingleTrackParameters::SetParRange(FitterParameterType::Type type, double min, double max)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  (*itr).second.SetMin(min);
  (*itr).second.SetMax(max);
}

void WCSimFitterSingleTrackParameters::SetParIsFixed(FitterParameterType::Type type, bool fixIt)
{
  std::map<FitterParameterType::Type, WCSimFitterParameter>::iterator itr;
  itr = fParameters.find(type);
  assert(itr != fParameters.end());
  (*itr).second.SetIsFixed(fixIt);
}

unsigned int WCSimFitterSingleTrackParameters::GetNumParameters() {
	return fParameters.size();
}


////////////////////////////////////////////////////////////////////
// FITTER PARAMETERS
////////////////////////////////////////////////////////////////////

WCSimFitterParameters::WCSimFitterParameters() : fNumTracks(0), fNumParameters(0){
	// TODO Auto-generated constructor stub

}

WCSimFitterParameters::~WCSimFitterParameters() {
	// TODO Auto-generated destructor stub
}

void WCSimFitterParameters::SetNumTracks(unsigned int nTracks)
{
  int tracksNeeded = fNumTracks;
  if(tracksNeeded < nTracks)
  {
    for(unsigned int toAdd = 0; toAdd < (nTracks - tracksNeeded) ; ++toAdd)
    {
      WCSimFitterSingleTrackParameters trackPars;
      AddTrack(trackPars);
      fTrackTypes.push_back(TrackType::Unknown);
   }
  }
  else if(tracksNeeded > nTracks)
  {
    fTrackPars.erase( fTrackPars.begin() + nTracks, fTrackPars.end() );
    fTrackTypes.erase( fTrackTypes.begin() + nTracks, fTrackTypes.end() );
  }


  fNumTracks = fTrackPars.size();
  assert(fNumTracks == nTracks);
  return;
}


void WCSimFitterParameters::AddTrack(WCSimFitterSingleTrackParameters trackPars) {
	fTrackPars.push_back(trackPars);
	fNumTracks = static_cast<int>(fTrackPars.size());
  fNumParameters += trackPars.GetNumParameters();
}

unsigned int WCSimFitterParameters::GetNumIndependentParameters() const{
	int toSubtract = 0;
	std::map<std::pair<unsigned int, unsigned int>, std::vector<FitterParameterType::Type> >::const_iterator mapItr;
	for(mapItr = fJoinedParams.begin(); mapItr != fJoinedParams.end() ; ++mapItr)
	{
		toSubtract += (*mapItr).second.size();
	}

	return fNumParameters - toSubtract;
}

unsigned int WCSimFitterParameters::GetNumParameters() const
{
	return fNumParameters;
}


WCSimFitterSingleTrackParameters * WCSimFitterParameters::GetTrackParameters(
		unsigned int trackNum) {
	assert(trackNum < fTrackPars.size());
	return &(fTrackPars.at(trackNum));
}

void WCSimFitterParameters::JoinParametersTogether(unsigned int track1, unsigned int track2,
		FitterParameterType::Type type) {
	assert(track1 < track2);
	std::pair<unsigned int, unsigned int> tracks = std::make_pair(track1, track2);
	std::map<std::pair<unsigned int, unsigned int>, std::vector<FitterParameterType::Type> >::iterator joinItr = fJoinedParams.find(tracks);
	if( joinItr != fJoinedParams.end() )
	{
		std::vector<FitterParameterType::Type> typeVec =(*joinItr).second;
		if( std::find( typeVec.begin(), typeVec.end(), type) == typeVec.end())
		{
			fJoinedParams[tracks].push_back(type);
		}
	}
	else
	{
		std::vector<FitterParameterType::Type> typeVec;
		typeVec.push_back(type);
		fJoinedParams[tracks] = typeVec;
	}
}

bool WCSimFitterParameters::GetJoinParametersTogether(unsigned int track1,
		unsigned int track2, FitterParameterType::Type type) {
	std::pair<unsigned int, unsigned int> tracks = std::make_pair(track1, track2);
	std::map<std::pair<unsigned int, unsigned int>, std::vector<FitterParameterType::Type> >::iterator joinItr = fJoinedParams.find(tracks);
	if( joinItr != fJoinedParams.end() )
	{
		std::vector<FitterParameterType::Type> typeVec =(*joinItr).second;
		if( std::find( typeVec.begin(), typeVec.end(), type) != typeVec.end())
		{
			return true;
		}
	}
	return false;
}

unsigned int WCSimFitterParameters::GetTrackIsJoinedWith(unsigned int track, FitterParameterType::Type type) {
	for(unsigned int iTrack = 0; iTrack < track; ++iTrack)
	{
		std::pair<unsigned int, unsigned int> tracks = std::make_pair(iTrack, track);
		std::map<std::pair<unsigned int, unsigned int>, std::vector<FitterParameterType::Type> >::iterator joinItr = fJoinedParams.find(tracks);
		if( joinItr != fJoinedParams.end() )
		{
			std::vector<FitterParameterType::Type> typeVec =(*joinItr).second;
			if( std::find( typeVec.begin(), typeVec.end(), type) != typeVec.end())
			{
				return iTrack;
			}
		}

	}
	// Return itself if not joined
	return track;
}

bool WCSimFitterParameters::GetIsParameterJoined(unsigned int track,
		FitterParameterType::Type type) {
	for(unsigned int iTrack = 0; iTrack < track; ++iTrack)
	{
		std::pair<unsigned int, unsigned int> tracks = std::make_pair(iTrack, track);
		std::map<std::pair<unsigned int, unsigned int>, std::vector<FitterParameterType::Type> >::iterator joinItr = fJoinedParams.find(tracks);
		if( joinItr != fJoinedParams.end() )
		{
			std::vector<FitterParameterType::Type> typeVec =(*joinItr).second;
			if( std::find( typeVec.begin(), typeVec.end(), type) != typeVec.end())
			{
				return true;
			}
		}

	}
	return false;
}

void WCSimFitterParameters::SetTrackType(unsigned int nTrack, TrackType::Type trackType)
{
  // std::cout << "WCSimFitterParameters::SetTrackType(" << nTrack << ", " << TrackType::AsString(trackType) << std::endl;
	fTrackTypes.at(nTrack) = trackType;
}

TrackType::Type WCSimFitterParameters::GetTrackType(const unsigned int &nTrack) const
{
	return fTrackTypes.at(nTrack);
}
