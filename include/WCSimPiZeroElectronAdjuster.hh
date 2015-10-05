/*
 * WCSimPiZeroElectronAdjuster.hh
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#ifndef WCSIMPIZEROELECTRONADJUSTER_HH_
#define WCSIMPIZEROELECTRONADJUSTER_HH_

#include "WCSimPiZeroSeeder.hh"

class WCSimLikelihoodTrackBase;
class WCSimFitterConfig;
class WCSimPiZeroSeed;

class WCSimPiZeroElectronAdjuster: public WCSimPiZeroSeeder
{
public:
	WCSimPiZeroElectronAdjuster(WCSimFitterConfig * config, WCSimLikelihoodTrackBase * singleElectron, const double &minus2LnL);
	virtual ~WCSimPiZeroElectronAdjuster();
	void MakeSeeds();



private:

	std::vector<WCSimPiZeroSeed*> IterateOverConversionDistances();
	std::vector<WCSimPiZeroSeed*>  IterateOverFirstTrackPerturbations(const double &convDistTrack1, const double &convDistTrack2);
	std::vector<TVector3> GetFirstTrackDirectionsToTry();
	std::vector<WCSimPiZeroSeed*>  GridSearchOverSecondTrackDirection(const double &convDistTrack1, const double &convDistTrack2,
																	  const TVector3 &track1Dir,
																	  int numTimesAlready);

	double GetFirstTrackEnergyEstimator(const double &singleTrackFitE, const double &cosThetaSep); //> Take the energy from a single track electron fit and estimate the leading photon's energy

	// Calculate the energy of the second photon assuming we have pi0 -> gamma gamma and know the direction
	// of both photons and the energy of the first
	std::pair<double, double> GetPiZeroPhotonEnergies(const TVector3 &track1Dir, const double &track1Energy, const TVector3 &track2Dir );

	bool SimilarEnergy(WCSimPiZeroSeed* seed) const;
	bool SimilarEnergy(WCSimLikelihoodTrackBase * track1, WCSimLikelihoodTrackBase * track2) const;
	bool SimilarEnergy(const double &energy1, const double &energy2) const;

	WCSimLikelihoodTrackBase * fSingleElectronTrack;
	double fSingleElectron2LnL;

	ClassDef(WCSimPiZeroElectronAdjuster,0);
};

#endif /* WCSIMPIZEROELECTRONADJUSTER_HH_ */
