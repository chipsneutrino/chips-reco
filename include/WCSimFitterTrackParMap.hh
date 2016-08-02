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
class WCSimFitterConfig;
typedef std::pair<unsigned int, FitterParameterType::Type> TrackAndType;

class WCSimFitterTrackParMap : public TObject {
public:
	WCSimFitterTrackParMap(WCSimFitterConfig * config);
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
	void FixDirection(int track = -1);
	void FreeDirection(int track = -1);
	void FixEnergy(int track = -1);
	void FreeEnergy(int track = -1);
	void FixConversionLength(int track = -1);
	void FreeConversionLength(int track = -1);
	void FixTime(int track = -1);
	void FreeTime(int track = -1);

	unsigned int GetIndex(int track, FitterParameterType type);
	unsigned int GetIndex(TrackAndType trackAndType);
	double GetMinValue(int track, FitterParameterType type);
	double GetMinValue(TrackAndType trackAndType);
	double GetMaxValue(int track, FitterParameterType type);
	double GetMaxValue(TrackAndType trackAndType);
	double GetCurrentValue(int track, FitterParameterType type);
	double GetCurrentValue(TrackAndType trackAndType);
	double GetStep(int track, FitterParameterType type);
	double GetStep(TrackAndType trackAndType);
	TrackType::Type GetTrackType(int track);

  bool GetIsFixed(int track, FitterParameterType type);
  bool GetIsFixed(TrackAndType trackAndType);

  void Print();
private:
	void SetArrays();
	void ResizeVectors();
	void FixOrFreeVertex(bool fixIt, int track = -1);
	void FixOrFreeDirection(bool fixIt, int track = -1);
	void FixOrFreeEnergy(bool fixIt, int track = -1);
	void FixOrFreeConversionLength(bool fixIt, int track = -1);
	void FixOrFreeTime(bool fixIt, int track = -1);
	std::map<TrackAndType, unsigned int> fTrackAndTypeIndexMap;
	std::vector<double> fCurrentValues;
	std::vector<double> fMinValues;
	std::vector<double> fMaxValues;
	std::vector<bool> fCurrentlyFixed;
	std::vector<bool> fAlwaysFixed;
	std::vector<bool> fIsEnergy;
	std::vector<double> fSteps;
	std::vector<std::string> fNames;
	std::vector<TrackType::Type> fTypes;
	int fNumParameters;
	WCSimFitterConfig * fFitterConfig;

ClassDef(WCSimFitterTrackParMap,1)

};

#endif /* WCSIMFITTERTRACKPARMAP_HH_ */
