/*
 * WCSimPiZeroFitter.hh
 *
 *  Created on: 21 Aug 2015
 *      Author: ajperch
 */

#ifndef WCSIMPIZEROFITTER_HH_
#define WCSIMPIZEROFITTER_HH_

#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimReco.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRootEvent.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimFitterTrackParMap.hh"
#include "WCSimTrackParameterEnums.hh"

#include <vector>
#include <map>

struct SingleElectronSeed
{
	TVector3 fVertex;
	TVector3 fDirection;
	double fTime;
	double fEnergy;
};

class WCSimPiZeroFitter : public WCSimLikelihoodFitter
{
public:
	WCSimPiZeroFitter();
	virtual ~WCSimPiZeroFitter();
	void RunFits();

protected:
private:
	// For building seed track combinations
	void SeedWithSingleElectron();
	void FitSingleElectronSeed();  // Hough transform only gives  the vertex and direction, so fit the energy
	double WrapFuncSingleElectron(const double * x);


	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > IterateOverConversionDistances();
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > IterateOverFirstTrackPerturbations(const double &convDistTrack1, const double &convDistTrack2);
	std::vector<TVector3> GetFirstTrackDirectionsToTry();
	std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase *> > GridSearchOverSecondTrackDirection(const double &convDistTrack1, const double &convDistTrack2,
	        																							 const TVector3 &track1Dir,
	        																							 int numTimesAlready,
                                                         double &bestSimilar2LnL,
                                                         double &bestDifferent2LnL);

  double GetFirstTrackEnergyEstimator(const double &singleTrackFitE, const double &cosThetaSep); //> Take the energy from a single track electron fit and estimate the leading photon's energy

	// Calculate the energy of the second photon assuming we have pi0 -> gamma gamma and know the direction
	// of both photons and the energy of the first
  std::pair<double, double> GetPiZeroPhotonEnergies(const TVector3 &track1Dir, const double &track1Energy, const TVector3 &track2Dir );

	// For fitting using each seed track combination
	void FitAfterFixingDirectionAndEnergy(std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*>);
	void FitAfterFixingEnergy(std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*>);
	void SetStartingTracks(std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> tracks);

	bool CanSetParam(const unsigned int &iTrack, const FitterParameterType::Type &type);
	void FitEventNumber(Int_t iEvent);


	/**
	 * Work out how many free parameters Minuit needs to run over
	 * @return Number of parameters to use in Minuit
	 */
	UInt_t GetNPars();

	SingleElectronSeed fSingleElectronSeed;

  std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase*> > SeedCheat(); // Temporary function to save running the same grid search over and over
  
	int fCalls;
	ClassDef(WCSimPiZeroFitter,1);
};

#endif /* WCSIMPIZEROFITTER_HH_ */
