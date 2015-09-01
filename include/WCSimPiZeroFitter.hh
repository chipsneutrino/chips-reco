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
	void FitSingleElectronSeedEnergy();  // Hough transform only gives  the vertex and direction, so fit the energy
	double WrapFuncSeedTrackEnergy(const double * x);
	double GetFirstPhotonEnergyEstimate(const double &singleElectronEnergy);


	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > IterateOverConversionDistances();
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > IterateOverFirstTrackPerturbations(const double &convDistTrack1, const double &convDistTrack2);
	std::vector<TVector3> GetFirstTrackDirectionsToTry();
	std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase *> > GridSearchOverSecondTrackDirection(const double &convDistTrack1, const double &convDistTrack2,
	        																							 const TVector3 &track1Dir,
	        																							 int numTimesAlready);

	// Calculate the energy of the second photon assuming we have pi0 -> gamma gamma and know the direction
	// of both photons and the energy of the first
	Double_t GetPiZeroSecondTrackEnergy(const TVector3 &track1Dir, const double &track1Energy, const TVector3 &track2Dir );

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



	int fCalls;
	ClassDef(WCSimPiZeroFitter,1);
};

#endif /* WCSIMPIZEROFITTER_HH_ */
