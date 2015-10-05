/*
 * WCSimPiZeroSeed.cc
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#include "WCSimPiZeroSeed.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrackFactory.hh"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimPiZeroSeed)
#endif

WCSimPiZeroSeed::WCSimPiZeroSeed(WCSimLikelihoodTrackBase* track1,
		WCSimLikelihoodTrackBase* track2, double minus2LnL)
{
	// These are always new'ed so this class can take responsibility for deleting them
	fTrack1 = WCSimLikelihoodTrackFactory::MakeTrack(track1->GetType(),
													 track1->GetX(), track1->GetY(), track1->GetZ(),
													 track1->GetT(),
													 track1->GetTheta(), track1->GetPhi(),
													 track1->GetE(),
													 track1->GetConversionDistance());
	fTrack2 = WCSimLikelihoodTrackFactory::MakeTrack(track2->GetType(),
													 track2->GetX(), track2->GetY(), track2->GetZ(),
													 track2->GetT(),
													 track2->GetTheta(), track2->GetPhi(),
													 track2->GetE(),
													 track2->GetConversionDistance());
	fMinus2LnL = minus2LnL;
}

// Copy constructor
WCSimPiZeroSeed::WCSimPiZeroSeed(const WCSimPiZeroSeed& other)
{
	WCSimLikelihoodTrackBase * track1 = other.GetTrack1();
	WCSimLikelihoodTrackBase * track2 = other.GetTrack2();
	fTrack1 = WCSimLikelihoodTrackFactory::MakeTrack(track1->GetType(),
													 track1->GetX(), track1->GetY(), track1->GetZ(),
													 track1->GetT(),
													 track1->GetTheta(), track1->GetPhi(),
													 track1->GetE(),
													 track1->GetConversionDistance());
	fTrack2 = WCSimLikelihoodTrackFactory::MakeTrack(track2->GetType(),
													 track2->GetX(), track2->GetY(), track2->GetZ(),
													 track2->GetT(),
													 track2->GetTheta(), track2->GetPhi(),
													 track2->GetE(),
													 track2->GetConversionDistance());
	fMinus2LnL = other.GetMinus2LnL();
}

// Copy assignment
WCSimPiZeroSeed& WCSimPiZeroSeed::operator =(const WCSimPiZeroSeed& rhs)
{
	delete fTrack1;
	delete fTrack2;
	WCSimLikelihoodTrackBase * track1 = rhs.GetTrack1();
	WCSimLikelihoodTrackBase * track2 = rhs.GetTrack2();
	fTrack1 = WCSimLikelihoodTrackFactory::MakeTrack(track1->GetType(),
													 track1->GetX(), track1->GetY(), track1->GetZ(),
													 track1->GetT(),
													 track1->GetTheta(), track1->GetPhi(),
													 track1->GetE(),
													 track1->GetConversionDistance());
	fTrack2 = WCSimLikelihoodTrackFactory::MakeTrack(track2->GetType(),
													 track2->GetX(), track2->GetY(), track2->GetZ(),
													 track2->GetT(),
													 track2->GetTheta(), track2->GetPhi(),
													 track2->GetE(),
													 track2->GetConversionDistance());
	fMinus2LnL = rhs.GetMinus2LnL();
	return *this;
}

bool operator <(const WCSimPiZeroSeed& lhs, const WCSimPiZeroSeed& rhs)
{
	return lhs.GetMinus2LnL() < rhs.GetMinus2LnL();
}

bool operator >(const WCSimPiZeroSeed& lhs, const WCSimPiZeroSeed& rhs)
{
	return lhs.GetMinus2LnL() > rhs.GetMinus2LnL();
}

bool operator <=(const WCSimPiZeroSeed& lhs, const WCSimPiZeroSeed& rhs)
{
	return !(lhs > rhs);
}

bool operator >=(const WCSimPiZeroSeed& lhs, const WCSimPiZeroSeed& rhs)
{
	return !(lhs < rhs);
}

WCSimPiZeroSeed::~WCSimPiZeroSeed()
{
	delete fTrack1;
	delete fTrack2;
}

WCSimLikelihoodTrackBase* WCSimPiZeroSeed::GetTrack1() const
{
	return fTrack1;
}

WCSimLikelihoodTrackBase* WCSimPiZeroSeed::GetTrack2() const
{
	return fTrack2;
}

double WCSimPiZeroSeed::GetMinus2LnL() const
{
	return fMinus2LnL;
}

void WCSimPiZeroSeed::Print()
{
	std::cout << " *** WCSimPiZeroSeed::Print() *** " << std::endl;
	std::cout << "Track 1: " << std::endl;
	fTrack1->Print();
	std::cout << "Track 2: " << std::endl;
	fTrack2->Print();
	std::cout << "These have: -2 * ln(Likelihood) = " << fMinus2LnL << std::endl;
}
