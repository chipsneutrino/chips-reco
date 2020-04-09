/*
 * WCSimPiZeroSeedGenerator.hh
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#pragma once

#include <vector>
#include "WCSimPiZeroHoughSeeder.hh"
#include "WCSimPiZeroSingleElectronSeeder.hh"
#include "TObject.h"

class WCSimFitterConfig;
class WCSimPiZeroSeed;

class WCSimPiZeroSeedGenerator : public TObject
{
public:
	WCSimPiZeroSeedGenerator(WCSimFitterConfig *config);
	virtual ~WCSimPiZeroSeedGenerator();

	std::vector<WCSimPiZeroSeed *> GetSeeds(const int &event);

	virtual unsigned int GetNumSeeds() const;

protected:
	void MakeSeeds();
	std::vector<WCSimPiZeroSeed *> fPiZeroSeeds;
	bool fMadeSeeds;
	void SetEvent(const int &event);

private:
	WCSimPiZeroHoughSeeder fPiZeroHoughSeeder;
	WCSimPiZeroSingleElectronSeeder fPiZeroSingleElectronSeeder;
	WCSimFitterConfig *fFitterConfig;
	int fEvent;
	ClassDef(WCSimPiZeroSeedGenerator, 0);
};
