#include "WCSimHitComparison.hh"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimHitPrediction)
ClassImp(WCSimSingleHitComparison)
ClassImp(WCSimHitComparison)
#endif

WCSimHitPrediction::WCSimHitPrediction()
{
	fPredictedCharge = 0.0;
	fPredictedTime = 0.0;
	fTotal2LnL = 0.0;
	fCharge2LnL = 0.0;
	fTime2LnL = 0.0;
	fHit2LnL = 0.0;
	fCutoff2LnL = 0.0;
}

WCSimHitPrediction::WCSimHitPrediction(const WCSimHitPrediction &other)
{
	fPredictedCharge = other.fPredictedCharge;
	fPredictedTime = other.fPredictedTime;
	fTotal2LnL = other.fTotal2LnL;
	fCharge2LnL = other.fCharge2LnL;
	fTime2LnL = other.fTime2LnL;
	fHit2LnL = other.fHit2LnL;
}

WCSimHitPrediction &WCSimHitPrediction::operator=(const WCSimHitPrediction &rhs)
{
	if (&rhs != this)
	{
		fPredictedCharge = rhs.fPredictedCharge;
		fPredictedTime = rhs.fPredictedTime;
		fTotal2LnL = rhs.fTotal2LnL;
		fCharge2LnL = rhs.fCharge2LnL;
		fTime2LnL = rhs.fTime2LnL;
		fHit2LnL = rhs.fHit2LnL;
		fCutoff2LnL = rhs.fCutoff2LnL;
	}
	return *this;
}

WCSimHitPrediction::~WCSimHitPrediction()
{
	// Empty, nothing to delete here
}

WCSimSingleHitComparison::WCSimSingleHitComparison(const WCSimSingleHitComparison &other) : fDigit(other.fDigit), fBestFit(other.fBestFit), fCorrectFit(other.fCorrectFit)
{
	// Empty
}

WCSimSingleHitComparison &WCSimSingleHitComparison::operator=(const WCSimSingleHitComparison &rhs)
{
	if (&rhs != this)
	{
		fDigit = rhs.fDigit;
		fBestFit = rhs.fBestFit;
		fCorrectFit = rhs.fCorrectFit;
	}
	return *this;
}

WCSimSingleHitComparison::~WCSimSingleHitComparison()
{
	// Empty, nothing to delete here
}

WCSimHitComparison::WCSimHitComparison()
{
	// Empty, nothing to delete here
}

WCSimHitComparison::~WCSimHitComparison()
{
	// Empty, nothing to delete here
}

WCSimHitComparison::WCSimHitComparison(const WCSimHitComparison &other) : fSingleHitComparisons(other.fSingleHitComparisons)
{
	// Empty
}

WCSimHitComparison &WCSimHitComparison::operator=(const WCSimHitComparison &rhs)
{
	if (&rhs != this)
	{
		fSingleHitComparisons = rhs.fSingleHitComparisons;
	}
	return *this;
}

void WCSimHitComparison::push_back(WCSimSingleHitComparison hitComparison)
{
	fSingleHitComparisons.push_back(hitComparison);
}

WCSimSingleHitComparison WCSimHitComparison::WCSimHitComparison::at(int i)
{
	return fSingleHitComparisons[i];
}

u_int WCSimHitComparison::size()
{
	return fSingleHitComparisons.size();
}