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
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > IterateOverConversionDistances();
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > IterateOverFirstTrackPerturbations(const double &convDistTrack1, const double &convDistTrack2);
	std::vector<TVector3> GetFirstTrackDirectionsToTry();
	std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase *> GridSearchOverSecondTrackDirection(const double &convDistTrack1, const double &convDistTrack2,
	        																							 const TVector3 &track1Dir);

	// For fitting using each seed track combination
	void FitAfterFixingDirectionAndEnergy(std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*>);
	void FitAfterFixingEnergy(std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*>);

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
