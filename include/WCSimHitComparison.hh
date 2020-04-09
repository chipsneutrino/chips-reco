#ifndef WCSIMHITCOMPARISON_HH
#define WCSIMHITCOMPARISON_HH

#include <vector>
#include "TObject.h"
#include "WCSimLikelihoodDigit.hh"

/**
 * @brief Container class to hold predicted charge, time and likelihood info
 */
class WCSimHitPrediction : public TObject
{
public:
	WCSimHitPrediction();
	WCSimHitPrediction(double predQ, double predT, double total2LnL, double charge2LnL, double time2LnL,
					   double hit2LnL) : fPredictedCharge(predQ), fPredictedTime(predT), fTotal2LnL(total2LnL), fCharge2LnL(charge2LnL), fTime2LnL(time2LnL), fHit2LnL(hit2LnL), fCutoff2LnL(total2LnL - charge2LnL - time2LnL)
	{
	}

	WCSimHitPrediction(const WCSimHitPrediction &other);
	WCSimHitPrediction &operator=(const WCSimHitPrediction &rhs);
	~WCSimHitPrediction();

	double GetPredictedCharge() const
	{
		return fPredictedCharge;
	}
	double GetPredictedTime() const
	{
		return fPredictedTime;
	}
	double GetTotal2LnL() const
	{
		return fTotal2LnL;
	}
	double GetCharge2LnL() const
	{
		return fCharge2LnL;
	}
	double GetTime2LnL() const
	{
		return fTime2LnL;
	}
	double GetHit2LnL() const
	{
		return fHit2LnL;
	}
	double GetCutoff2LnL() const
	{
		return fCutoff2LnL;
	}

private:
	double fPredictedCharge;
	double fPredictedTime;
	double fTotal2LnL;
	double fCharge2LnL;
	double fTime2LnL;
	double fHit2LnL;
	double fCutoff2LnL;

	ClassDef(WCSimHitPrediction, 1);
};

class WCSimSingleHitComparison
{

public:
	WCSimSingleHitComparison(WCSimLikelihoodDigit digit, WCSimHitPrediction bestFit, WCSimHitPrediction correctFit) : fDigit(digit), fBestFit(bestFit), fCorrectFit(correctFit)
	{
	}

	WCSimSingleHitComparison() : fDigit(
									 WCSimLikelihoodDigit(-9999.9, -9999.9, -9999.9, -9999.9, 0, -999, -9999.9, -9999.9, -9999.9, "",
														  0x0, 0, 0)),
								 fBestFit(WCSimHitPrediction(0, 0, 0, 0, 0, 0)), fCorrectFit(
																					 WCSimHitPrediction(0, 0, 0, 0, 0, 0))
	{
	}

	WCSimSingleHitComparison(const WCSimSingleHitComparison &other);
	WCSimSingleHitComparison &operator=(const WCSimSingleHitComparison &rhs);
	~WCSimSingleHitComparison();

	double GetFitPredictedCharge() const
	{
		return fBestFit.GetPredictedCharge();
	}

	double GetFitPredictedTime() const
	{
		return fBestFit.GetPredictedTime();
	}

	double GetFitTotal2LnL() const
	{
		return fBestFit.GetTotal2LnL();
	}
	double GetFitCharge2LnL() const
	{
		return fBestFit.GetCharge2LnL();
	}
	double GetFitTime2LnL() const
	{
		return fBestFit.GetTime2LnL();
	}
	double GetFitHit2LnL() const
	{
		return fBestFit.GetHit2LnL();
	}
	double GetFitCutoff2LnL() const
	{
		return fBestFit.GetCutoff2LnL();
	}

	double GetCorrectPredictedCharge() const
	{
		return fCorrectFit.GetPredictedCharge();
	}

	double GetCorrectPredictedTime() const
	{
		return fCorrectFit.GetPredictedTime();
	}
	double GetCorrectTotal2LnL() const
	{
		return fCorrectFit.GetTotal2LnL();
	}
	double GetCorrectCharge2LnL() const
	{
		return fCorrectFit.GetCharge2LnL();
	}
	double GetCorrectTime2LnL() const
	{
		return fCorrectFit.GetTime2LnL();
	}
	double GetCorrectHit2LnL() const
	{
		return fCorrectFit.GetHit2LnL();
	}
	double GetCorrectCutoff2LnL() const
	{
		return fCorrectFit.GetCutoff2LnL();
	}

	double GetHitQ() const
	{
		return fDigit.GetQ();
	}
	double GetHitT() const
	{
		return fDigit.GetT();
	}
	double GetPMTX() const
	{
		return fDigit.GetX();
	}
	double GetPMTY() const
	{
		return fDigit.GetY();
	}
	double GetPMTZ() const
	{
		return fDigit.GetZ();
	}
	double GetTubeID() const
	{
		return fDigit.GetTubeId();
	}
	TVector3 GetPMTPos() const
	{
		return fDigit.GetPos();
	}

private:
	WCSimLikelihoodDigit fDigit; // For PMT (x,y,z,q,t)
	WCSimHitPrediction fBestFit;
	WCSimHitPrediction fCorrectFit;

	ClassDef(WCSimSingleHitComparison, 1);
};

typedef std::vector<WCSimSingleHitComparison> WCSimHitComparison;
#endif
