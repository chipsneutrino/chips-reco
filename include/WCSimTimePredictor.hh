/*
 * WCSimTimePredictor.hh
 *
 *  Created on: 7 Jul 2015
 *      Author: andy
 */
#ifndef INCLUDE_WCSIMTIMEPREDICTOR_HH_
#define INCLUDE_WCSIMTIMEPREDICTOR_HH_

#include <vector>
#include <TObject.h>
class WCSimLikelihoodDigitArray;
class WCSimLikelihoodTrackBase;
class WCSimLikelihoodDigit;
class WCSimEmissionProfileManager;
class WCSimLikelihoodTuner;

class WCSimTimePredictor : public TObject
{
public:
	WCSimTimePredictor();
	WCSimTimePredictor(WCSimLikelihoodDigitArray *myArray, WCSimEmissionProfileManager *myEmissionProfileManager);
	virtual ~WCSimTimePredictor();

	/**
		 * Add another track to calculate the likelihood for several particles at once
		 * @param myTrack Track object to add
		 */
	void AddTrack(WCSimLikelihoodTrackBase *myTrack);

	/**
		 * Set all the tracks that contribute to the likelihood at once
		 * @param myTrack Vector of all the track objects to consider
		 */
	void SetTracks(std::vector<WCSimLikelihoodTrackBase *> myTracks);

	/// Remove all the tracks currently loaded
	void ClearTracks();

	/**
		 * Replace the current event with hits from a different one
		 * @param myDigitArray New array of PMT responses
		 */
	void UpdateDigitArray(WCSimLikelihoodDigitArray *myDigitArray);

	double GetPredictedTime(unsigned int iDigit);
	double GetPredictedTime(WCSimLikelihoodDigit *myDigit);
	double GetHitTime(WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit, bool &isHit);
	double GetHitTime(WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit, bool &isHit,
					  const double &survivalDistance, const double &escapeDistance);

	std::vector<double> GetAllPredictedTimes();

private:
	double GetTravelDistance(WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit);
	double GetSurvivalDistance(WCSimLikelihoodTrackBase *myTrack);
	double GetEscapeDistance(WCSimLikelihoodTrackBase *myTrack);
	std::vector<WCSimLikelihoodTrackBase *> fTracks;
	WCSimLikelihoodDigitArray *fDigitArray;
	WCSimEmissionProfileManager *fEmissionProfileManager;
	WCSimLikelihoodTuner *fTuner;
	ClassDef(WCSimTimePredictor, 0)
};

#endif /* INCLUDE_WCSIMTIMEPREDICTOR_HH_ */
