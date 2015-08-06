/*
 * WCSimFitterTrackParMap.cc
 *
 *  Created on: 27 May 2015
 *      Author: andy
 */

#include "WCSimFitterTrackParMap.hh"
#include "WCSimFitterInterface.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimTrackParameterEnums.hh"
#include <cassert>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimFitterTrackParMap)
#endif

WCSimFitterTrackParMap::WCSimFitterTrackParMap() {
	// TODO Auto-generated constructor stub

}

WCSimFitterTrackParMap::~WCSimFitterTrackParMap() {
	// TODO Auto-generated destructor stub
}

void WCSimFitterTrackParMap::Set() {
	fTrackAndTypeIndexMap.clear();
	unsigned int nTracks = WCSimFitterConfig::Instance()->GetNumTracks();
  //std::cout << "Number of tracks in ::Set is " << nTracks << std::endl;

	unsigned int arrayCounter = 0; // How many unique track parameters there are
	for (unsigned int iTrack = 0; iTrack < nTracks; ++iTrack) {

	  std::vector<FitterParameterType::Type> allTypes = 
        FitterParameterType::GetAllAllowedTypes(WCSimFitterInterface::Instance()->GetTrackType(iTrack));
		for (unsigned int iParam = 0; iParam < allTypes.size(); ++iParam) {
			// Is it joined?
			TrackAndType trackPar(iTrack, allTypes.at(iParam));
			unsigned int isJoinedWith =
					WCSimFitterConfig::Instance()->GetTrackIsJoinedWith(iTrack,
							FitterParameterType::AsString(trackPar.second).c_str());
			if (isJoinedWith == iTrack) // Itself, i.e. not joined with any other track
			{
				// Then set it to the value from the first track
        //std::cout << "Array counter = " << arrayCounter << "  " << "Independent params = " << WCSimFitterConfig::Instance()->GetNumIndependentParameters() << std::endl; 
				assert(arrayCounter	< WCSimFitterConfig::Instance()->GetNumIndependentParameters());
				fTrackAndTypeIndexMap[trackPar] = arrayCounter;
				arrayCounter++;
			}
			else
			{
				// Otherwise don't; set it to the next index value instead
				TrackAndType joinPar(isJoinedWith, allTypes.at(iParam));
				fTrackAndTypeIndexMap[trackPar] = fTrackAndTypeIndexMap[joinPar];
			}
		}
	}
	SetArrays();
}

void WCSimFitterTrackParMap::SetArrays()
{
  //std::cout << "SetArrays" << std::endl;
	ResizeVectors();

	Int_t iParam = 0;
	unsigned int numTracks = WCSimFitterConfig::Instance()->GetNumTracks();


	for(unsigned int jTrack = 0; jTrack < numTracks; ++jTrack)
	{
		fTypes.at(jTrack) = (WCSimFitterConfig::Instance()->GetTrackType(jTrack));
		std::vector<FitterParameterType::Type> paramTypes = FitterParameterType::GetAllAllowedTypes(fTypes.at(jTrack));
		std::vector<FitterParameterType::Type>::const_iterator typeItr = paramTypes.begin();

		for( ; typeItr != paramTypes.end(); ++typeItr)
		{
			const char * name  = FitterParameterType::AsString(*typeItr).c_str();

			if(WCSimFitterConfig::Instance()->GetIsParameterJoined(jTrack, name)){
				continue;
			}
			// Only set it if it's a unique parameter
			else{
				fCurrentValues[iParam] = WCSimFitterConfig::Instance()->GetParStart(jTrack,name);
				fMinValues[iParam]   = WCSimFitterConfig::Instance()->GetParMin(jTrack,name);
				fMaxValues[iParam]   = WCSimFitterConfig::Instance()->GetParMax(jTrack,name);
				fSteps[iParam] = (WCSimFitterConfig::Instance()->GetParMax(jTrack,name) - WCSimFitterConfig::Instance()->GetParMin(jTrack, name)) / 250.0;
				fAlwaysFixed[iParam] = WCSimFitterConfig::Instance()->GetIsFixedParameter(jTrack,name);
				fCurrentlyFixed[iParam] = fAlwaysFixed[iParam];
				TString parName = Form("%s_track%d", name, jTrack);
				fNames[iParam]  = std::string(parName.Data());
				fIsEnergy[iParam] = ((*typeItr) == FitterParameterType::kEnergy);
				iParam++;
			}
		}
	}
}
void WCSimFitterTrackParMap::ResizeVectors()
{
  //std::cout << "Resize vectors" << std::endl;
	unsigned int arraySize = fTrackAndTypeIndexMap.size();
	fCurrentValues.clear();
	fMinValues.clear();
	fMaxValues.clear();
	fCurrentlyFixed.clear();
	fAlwaysFixed.clear();
	fIsEnergy.clear();
	fSteps.clear();
	fNames.clear();
	fTypes.clear();

	fCurrentValues.resize(arraySize);
	fMinValues.resize(arraySize);
	fMaxValues.resize(arraySize);
	fCurrentlyFixed.resize(arraySize);
	fAlwaysFixed.resize(arraySize);
	fIsEnergy.resize(arraySize);
	fSteps.resize(arraySize);
	fNames.resize(arraySize);
	fTypes.resize(WCSimFitterInterface::Instance()->GetNumTracks());
  //std::cout << "WCSimFitterInterface::ResizeVectors() = " << WCSimFitterInterface::Instance()->GetNumTracks() << std::endl;


}

void WCSimFitterTrackParMap::FixDirection(int track) {
	FixOrFreeDirection(1, track);
}

void WCSimFitterTrackParMap::FreeDirection(int track) {
	FixOrFreeDirection(0, track);
}

void WCSimFitterTrackParMap::FixVertex(int track) {
	FixOrFreeVertex(1, track);
}

void WCSimFitterTrackParMap::FreeVertex(int track) {
	FixOrFreeVertex(0, track);
}

// Fixes the vertex (position, not direction) for a certain track (or all tracks
// if no argument given) if fixIt = true, otherwise frees it
void WCSimFitterTrackParMap::FixOrFreeDirection(bool fixIt, int track) {
	int firstTrack = (track != -1) ? track : 0;
	int lastTrack = (track != -1) ? track : WCSimFitterInterface::Instance()->GetNumTracks();

	for(int iTrack = firstTrack; iTrack < lastTrack; ++iTrack)
	{
		TrackAndType myDirTh(iTrack, FitterParameterType::kDirTh);
		TrackAndType myDirPhi(iTrack, FitterParameterType::kDirPhi);

    if(!fAlwaysFixed.at(fTrackAndTypeIndexMap[myDirTh])){
		  fCurrentlyFixed.at(fTrackAndTypeIndexMap[myDirTh]) = fixIt;
    }
    if(!fAlwaysFixed.at(fTrackAndTypeIndexMap[myDirPhi])){
		  fCurrentlyFixed.at(fTrackAndTypeIndexMap[myDirPhi]) = fixIt;
    }
	}
}


// Fixes the vertex (position, not direction) for a certain track (or all tracks
// if no argument given) if fixIt = true, otherwise frees it
void WCSimFitterTrackParMap::FixOrFreeVertex(bool fixIt, int track) {
	int firstTrack = (track != -1) ? track : 0;
	int lastTrack = (track != -1) ? track : WCSimFitterInterface::Instance()->GetNumTracks();

	for(int iTrack = firstTrack; iTrack < lastTrack; ++iTrack)
	{
		TrackAndType myVtxX(iTrack, FitterParameterType::kVtxX);
		TrackAndType myVtxY(iTrack, FitterParameterType::kVtxY);
		TrackAndType myVtxZ(iTrack, FitterParameterType::kVtxZ);
		TrackAndType myVtxT(iTrack, FitterParameterType::kVtxT);

    if(!fAlwaysFixed.at(fTrackAndTypeIndexMap[myVtxX])){
		  fCurrentlyFixed.at(fTrackAndTypeIndexMap[myVtxX]) = fixIt;
    }
    if(!fAlwaysFixed.at(fTrackAndTypeIndexMap[myVtxY])){
		  fCurrentlyFixed.at(fTrackAndTypeIndexMap[myVtxY]) = fixIt;
    }
    if(!fAlwaysFixed.at(fTrackAndTypeIndexMap[myVtxZ])){
		  fCurrentlyFixed.at(fTrackAndTypeIndexMap[myVtxZ]) = fixIt;
    }
    if(!fAlwaysFixed.at(fTrackAndTypeIndexMap[myVtxT])){
		  fCurrentlyFixed.at(fTrackAndTypeIndexMap[myVtxT]) = fixIt;
    }
	}
}

void WCSimFitterTrackParMap::FixEnergy(int track) {
	FixOrFreeEnergy(1, track);
}

void WCSimFitterTrackParMap::FreeEnergy(int track) {
	FixOrFreeEnergy(0, track);
}

void WCSimFitterTrackParMap::FixOrFreeEnergy(bool fixIt, int track) {
	int firstTrack = (track != -1) ? track : 0;
		int lastTrack = (track != -1) ? track : WCSimFitterInterface::Instance()->GetNumTracks();

		for(int iTrack = firstTrack; iTrack < lastTrack; ++iTrack)
		{
			TrackAndType myEnergy(iTrack, FitterParameterType::kEnergy);
			fCurrentlyFixed.at(fTrackAndTypeIndexMap[myEnergy]) = fixIt;
		}
}

std::vector<double> WCSimFitterTrackParMap::GetCurrentValues() {
	return fCurrentValues;
}

std::vector<double> WCSimFitterTrackParMap::GetMinValues() {
	return fMinValues;
}

std::vector<double> WCSimFitterTrackParMap::GetMaxValues() {
	return fMaxValues;
}

std::vector<bool> WCSimFitterTrackParMap::GetCurrentlyFixed() {
	return fCurrentlyFixed;
}

std::vector<bool> WCSimFitterTrackParMap::GetAlwaysFixed() {
	return fAlwaysFixed;
}

std::vector<bool> WCSimFitterTrackParMap::GetIsEnergy() {
	return fIsEnergy;
}

std::vector<double> WCSimFitterTrackParMap::GetSteps() {
	return fSteps;
}

std::vector<std::string> WCSimFitterTrackParMap::GetNames() {
	return fNames;
}

void WCSimFitterTrackParMap::SetCurrentValue(TrackAndType pairToSet,
		double value) {
  int index = GetIndex(pairToSet);
	fCurrentValues[index] = value;
}

void WCSimFitterTrackParMap::SetCurrentValue(int trackNum,
		FitterParameterType type, double value) {
	TrackAndType myPair(trackNum, type);
	SetCurrentValue(myPair, value);
}

void WCSimFitterTrackParMap::SetCurrentValue(int arrayIndex, double value) {
	assert(arrayIndex >= 0 && static_cast<unsigned int>(arrayIndex) < fCurrentValues.size());
	fCurrentValues.at(arrayIndex) = value;
}

int WCSimFitterTrackParMap::GetIndex(int track, FitterParameterType type) {
	TrackAndType myPair(track, type);
	return GetIndex(myPair);
}

int WCSimFitterTrackParMap::GetIndex(TrackAndType trackAndType) {
	std::map<TrackAndType, unsigned int>::iterator findIndex = fTrackAndTypeIndexMap.find(trackAndType);
	assert(findIndex != fTrackAndTypeIndexMap.end());
	return (*findIndex).second;
}

double WCSimFitterTrackParMap::GetMinValue(int track,
		FitterParameterType type) {
	TrackAndType myPair(track, type);
	return GetMinValue(myPair);
}

double WCSimFitterTrackParMap::GetMinValue(TrackAndType trackAndType) {
	return fMinValues[GetIndex(trackAndType)];
}

double WCSimFitterTrackParMap::GetMaxValue(int track,
		FitterParameterType type) {
	TrackAndType myPair(track, type);
	return GetMaxValue(myPair);
}

double WCSimFitterTrackParMap::GetMaxValue(TrackAndType trackAndType) {
	return fMaxValues[GetIndex(trackAndType)];
}

double WCSimFitterTrackParMap::GetCurrentValue(int track,
		FitterParameterType type) {
	TrackAndType myPair(track, type);
	return GetCurrentValue(myPair);
}

double WCSimFitterTrackParMap::GetCurrentValue(TrackAndType trackAndType) {
	int index = GetIndex(trackAndType);
	return fCurrentValues[index];
}

double WCSimFitterTrackParMap::GetStep(int track,
		FitterParameterType type) {
	TrackAndType myPair(track, type);
	return GetStep(myPair);
}

double WCSimFitterTrackParMap::GetStep(TrackAndType trackAndType) {
	int index = GetIndex(trackAndType);
	return fSteps[index];
}

TrackType::Type WCSimFitterTrackParMap::GetTrackType(
		int track) {
  //std::cout << "Track = " << track << "   " << "Size = " << fTypes.size() << std::endl;
	assert(track >= 0 && static_cast<unsigned int>(track) < fTypes.size());
	return fTypes.at(track);
}

bool WCSimFitterTrackParMap::GetIsFixed(int track, FitterParameterType type){
  TrackAndType myPair(track, type);
  return GetIsFixed(myPair);
}

bool WCSimFitterTrackParMap::GetIsFixed(TrackAndType trackAndType){
  int index = GetIndex(trackAndType);
  return fCurrentlyFixed[index];
}



