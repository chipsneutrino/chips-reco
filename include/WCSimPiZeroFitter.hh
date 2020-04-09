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
#include "WCSimPiZeroSeedGenerator.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRootEvent.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimFitterTrackParMap.hh"
#include "WCSimTrackParameterEnums.hh"

#include <vector>
#include <map>

class WCSimFitterConfig;

class WCSimPiZeroFitter : public WCSimLikelihoodFitter
{
public:
	WCSimPiZeroFitter(WCSimFitterConfig *config);
	virtual ~WCSimPiZeroFitter();
	void RunFits();

protected:
private:
	std::vector<WCSimPiZeroSeed *> GetSeeds();

	// For building seed track combinations
	void SeedEvent();
	void SeedWithHoughTransform();
	void SeedWithSingleElectron();
	void FitSingleElectronSeed(); // Hough transform only gives  the vertex and direction, so fit the energy
	double WrapFuncSingleElectron(const double *x);

	// For fitting using each seed track combination
	void FitAfterFixingDirectionAndEnergy(WCSimPiZeroSeed *seed);
	void FitAfterFixingEnergy(WCSimPiZeroSeed *seed);
	void SetStartingTracks(WCSimPiZeroSeed *seed);

	bool CanSetParam(const unsigned int &iTrack, const FitterParameterType::Type &type);
	void FitEventNumber(Int_t iEvent);

	/**
		 * Work out how many free parameters Minuit needs to run over
		 * @return Number of parameters to use in Minuit
		 */
	UInt_t GetNPars();

	std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase *>> SeedCheat(); // Temporary function to save running the same grid search over and over

	int fCalls;
	WCSimPiZeroSeedGenerator fPiZeroSeedGenerator;
	ClassDef(WCSimPiZeroFitter, 0);
};

#endif /* WCSIMPIZEROFITTER_HH_ */
