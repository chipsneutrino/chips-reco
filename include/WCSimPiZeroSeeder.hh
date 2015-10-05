/*
 * WCSimPiZeroSeeder.hh
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#ifndef WCSIMPIZEROSEEDER_HH_
#define WCSIMPIZEROSEEDER_HH_

#include "WCSimLikelihoodFitter.hh"
#include <vector>
class WCSimPiZeroSeed;
class WCSimFitterConfig;

class WCSimPiZeroSeeder: public WCSimLikelihoodFitter
{
public:
	WCSimPiZeroSeeder(WCSimFitterConfig* config);
	virtual ~WCSimPiZeroSeeder();

	std::vector<WCSimPiZeroSeed*> GetSeeds();
	WCSimPiZeroSeed * GetSeed(unsigned int i);
	unsigned int GetNumSeeds();
	virtual void SetEvent(const int &event);

protected:
	virtual void MakeSeeds() = 0;
	std::vector<WCSimPiZeroSeed*> fPiZeroSeeds;
	bool fMadeSeeds;

	ClassDef(WCSimPiZeroSeeder,0);

};

#endif /* WCSIMPIZEROSEEDER_HH_ */
