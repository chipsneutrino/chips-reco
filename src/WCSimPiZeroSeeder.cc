/*
 * WCSimPiZeroSeeder.cc
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#include "WCSimFitterConfig.hh"
#include "WCSimPiZeroSeeder.hh"
#include "WCSimPiZeroSeed.hh"
#include <cassert>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimPiZeroSeeder)
#endif

WCSimPiZeroSeeder::WCSimPiZeroSeeder(WCSimFitterConfig * config) : WCSimLikelihoodFitter(config)
{
	// TODO Auto-generated constructor stub
	fMadeSeeds = false;
}

WCSimPiZeroSeeder::~WCSimPiZeroSeeder()
{
	// TODO Auto-generated destructor stub
}

std::vector<WCSimPiZeroSeed*> WCSimPiZeroSeeder::GetSeeds()
{
	if(!fMadeSeeds){
		MakeSeeds();
	}
	return fPiZeroSeeds;
}

WCSimPiZeroSeed* WCSimPiZeroSeeder::GetSeed(unsigned int i)
{
	if(!fMadeSeeds){
		MakeSeeds();
	}
	assert(i < fPiZeroSeeds.size());
	return fPiZeroSeeds.at(i);
}

unsigned int WCSimPiZeroSeeder::GetNumSeeds()
{
	if(fMadeSeeds){
		return fPiZeroSeeds.size();
	}
	return 0;
}

void WCSimPiZeroSeeder::SetEvent(const int& event)
{
	std::cout << "Setting event" << event << std::endl;
	WCSimLikelihoodFitter::SetEvent(event);
	fPiZeroSeeds.clear();
	fMadeSeeds = false;
}
