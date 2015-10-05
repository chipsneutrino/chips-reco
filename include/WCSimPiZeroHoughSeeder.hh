/*
 * WCSimPiZeroHoughSeeder.hh
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#ifndef WCSIMPIZEROHOUGHSEEDER_HH_
#define WCSIMPIZEROHOUGHSEEDER_HH_

#include "WCSimPiZeroSeeder.hh"

#include "TVector3.h"

#include <vector>
#include <utility>

class WCSimLikelihoodTrackBase;
class WCSimFitterConfig;

class WCSimPiZeroHoughSeeder: public WCSimPiZeroSeeder
{
public:
	WCSimPiZeroHoughSeeder(WCSimFitterConfig * config);
	virtual ~WCSimPiZeroHoughSeeder();
	void MakeSeeds();

private:
	void RunHough();
	void FitUsingHoughResults();
	double FitTwoTracks();
	void SetStartingEnergies(WCSimRecoRing* ring1, WCSimRecoRing* ring2, WCSimLikelihoodTrackBase* track1, WCSimLikelihoodTrackBase* track2);

	TVector3 FindSharedVertex(std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> tracks,
							  double &convDist1, double &convDist2);

	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > fHoughResults;
	ClassDef(WCSimPiZeroHoughSeeder,0);

};

#endif /* WCSIMPIZEROHOUGHSEEDER_HH_ */
