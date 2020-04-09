/*
 * WCSimPiZeroSingleElectronSeeder.hh
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#pragma once

#include "WCSimPiZeroSeeder.hh"

class WCSimLikelihoodTrackBase;
class WCSimFitterConfig;

class WCSimPiZeroSingleElectronSeeder : public WCSimPiZeroSeeder
{
public:
	WCSimPiZeroSingleElectronSeeder(WCSimFitterConfig *config);
	virtual ~WCSimPiZeroSingleElectronSeeder();
	WCSimLikelihoodTrackBase *GetSeed();
	void MakeSeeds();
	double GetMinus2LnL();
	void SetEvent(const int &iEvent);

private:
	double FitSingleTrack();
	WCSimLikelihoodTrackBase *fSingleElectronSeed;
	double fMinus2LnL;
	ClassDef(WCSimPiZeroSingleElectronSeeder, 0);
};
