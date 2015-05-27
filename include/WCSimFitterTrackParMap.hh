/*
 * WCSimFitterTrackParMap.hh
 *
 *  Created on: 27 May 2015
 *      Author: andy
 */

#ifndef WCSIMFITTERTRACKPARMAP_HH_
#define WCSIMFITTERTRACKPARMAP_HH_
#include <map>
#include <string>
#include <vector>
#include <TObject.h>
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrackParameterEnums.hh"
typedef std::pair<UInt_t, FitterParameterType::Type> TrackAndType;

class WCSimFitterTrackParMap : public TObject {
public:
	WCSimFitterTrackParMap();
	virtual ~WCSimFitterTrackParMap();
	void Set();
	std::vector<double> GetCurrentValues();
	std::vector<double> GetMinValues();
	std::vector<double> GetMaxValues();
	std::vector<bool> GetCurrentlyFixed();
	std::vector<bool> GetAlwaysFixed();
	std::vector<bool> GetIsEnergy();
	std::vector<double> GetSteps();
	std::vector<std::string> GetNames();

	void SetCurrentValue(TrackAndType pairToSet, double value);
	void SetCurrentValue(int trackNum, FitterParameterType type, double value);
	void SetCurrentValue(int arrayIndex, double value);

	void FixVertex(int track = -1);
	void FreeVertex(int track = -1);
	void FixEnergy(int track = -1);
	void FreeEnergy(int track = -1);

	int GetIndex(int track, FitterParameterType type);
	int GetIndex(TrackAndType trackAndType);
	double GetMinValue(int track, FitterParameterType type);
	double GetMinValue(TrackAndType trackAndType);
	double GetMaxValue(int track, FitterParameterType type);
	double GetMaxValue(TrackAndType trackAndType);
	double GetCurrentValue(int track, FitterParameterType type);
	double GetCurrentValue(TrackAndType trackAndType);
	WCSimLikelihoodTrack::TrackType GetTrackType(int track);


private:
	void SetArrays();
	void ResizeVectors();
	void FixOrFreeVertex(bool fixIt, int track = -1);
	void FixOrFreeEnergy(bool fixIt, int track = -1);
	std::map<TrackAndType, UInt_t> fTrackAndTypeIndexMap;
	std::vector<double> fCurrentValues;
	std::vector<double> fMinValues;
	std::vector<double> fMaxValues;
	std::vector<bool> fCurrentlyFixed;
	std::vector<bool> fAlwaysFixed;
	std::vector<bool> fIsEnergy;
	std::vector<double> fSteps;
	std::vector<std::string> fNames;
	std::vector<WCSimLikelihoodTrack::TrackType> fTypes;
	int fNumParameters;

ClassDef(WCSimFitterTrackParMap,1)

};

#endif /* WCSIMFITTERTRACKPARMAP_HH_ */
