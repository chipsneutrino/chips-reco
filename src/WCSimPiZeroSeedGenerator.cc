/*
 * WCSimPiZeroSeedGenerator.cc
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#include "WCSimFitterConfig.hh"
#include "WCSimPiZeroSeed.hh"
#include "WCSimPiZeroSeedGenerator.hh"
#include "WCSimPiZeroElectronAdjuster.hh"
#include "WCSimPiZeroHoughSeeder.hh"
#include "WCSimPiZeroSingleElectronSeeder.hh"
#include <vector>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimPiZeroSeedGenerator)
#endif

WCSimPiZeroSeedGenerator::WCSimPiZeroSeedGenerator(WCSimFitterConfig * config) : fPiZeroHoughSeeder(config), fPiZeroSingleElectronSeeder(config)
{
	// TODO Auto-generated constructor stub
	fFitterConfig = config;
	fMadeSeeds = false;
	fEvent = -999;
}

WCSimPiZeroSeedGenerator::~WCSimPiZeroSeedGenerator()
{
	// TODO Auto-generated destructor stub
}

std::vector<WCSimPiZeroSeed*> WCSimPiZeroSeedGenerator::GetSeeds(const int &event)
{
	std::cout << "About to call SetEvent" << std::endl;
	SetEvent(event);
	if(!fMadeSeeds){  MakeSeeds(); }
	std::cout << "Returning fPiZeroSeeds" << std::endl;
	return fPiZeroSeeds;
}

unsigned int WCSimPiZeroSeedGenerator::GetNumSeeds() const
{
	if(!fMadeSeeds) { return 0; }
	return fPiZeroSeeds.size();
}

void WCSimPiZeroSeedGenerator::MakeSeeds()
{
	if(fMadeSeeds){ return; }
	fPiZeroSeeds.clear();

	fPiZeroHoughSeeder.MakeSeeds();
	fPiZeroSingleElectronSeeder.MakeSeeds();

	bool singleElectronIsBest = true;
	WCSimLikelihoodTrackBase * singleElectron = fPiZeroSingleElectronSeeder.GetSeed();
	double single2LnL = fPiZeroSingleElectronSeeder.GetMinus2LnL();

	std::cout << "WCSimPiZeroSeedGenerator::MakeSeeds: single electron 2lnL = " << single2LnL << std::endl;
	for(unsigned int iSeed = 0; iSeed < fPiZeroHoughSeeder.GetNumSeeds(); ++iSeed)
	{
		fPiZeroSeeds.push_back(fPiZeroHoughSeeder.GetSeed(iSeed));
		std::cout << "WCSimPiZeroSeedGenerator::MakeSeeds: hough seed 2lnL = " << fPiZeroSeeds.back()->GetMinus2LnL() << std::endl;
		if(fPiZeroSeeds.back()->GetMinus2LnL() < single2LnL)
		{
			singleElectronIsBest = false;
		}
	}


	std::cout << "Was single electron best? " << singleElectronIsBest << std::endl;
	if(singleElectronIsBest)
	{
		std::cout << "Adjusting single electron seeds" << std::endl;
		WCSimPiZeroElectronAdjuster adj(fFitterConfig, singleElectron, single2LnL);
		adj.SetEvent(fEvent);
		adj.MakeSeeds();
		std::cout << "There were " << adj.GetNumSeeds() << " adjusted seeds" << std::endl;
		for(unsigned int iSeed = 0; iSeed < adj.GetNumSeeds(); ++iSeed)
		{
			fPiZeroSeeds.push_back(adj.GetSeed(iSeed));
		}
	}
	return;
}

void WCSimPiZeroSeedGenerator::SetEvent(const int& event)
{
	if(fEvent == event)
	{
		return;
	}
	else
	{
		fPiZeroHoughSeeder.SetEvent(event);
		fPiZeroSingleElectronSeeder.SetEvent(event);
		fEvent = event;
		fMadeSeeds = false;
		fPiZeroSeeds.clear();
	}
}
