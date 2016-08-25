#include "WCSimHitComparison.hh"

#ifndef REFLEX_DICTIONARY
    ClassImp(WCSimHitPrediction)
    ClassImp(WCSimSingleHitComparison)
#endif

    WCSimHitPrediction::WCSimHitPrediction(const WCSimHitPrediction& other)
    {
            fPredictedCharge = other.fPredictedCharge;
            fPredictedTime = other.fPredictedTime;
            fTotal2LnL = other.fTotal2LnL;
            fCharge2LnL = other.fCharge2LnL;
            fTime2LnL = other.fTime2LnL;
            fHit2LnL  = other.fHit2LnL;
    }

    WCSimHitPrediction& WCSimHitPrediction::operator=(
                const WCSimHitPrediction &rhs
                )
    {
        if(&rhs != this)
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

    WCSimSingleHitComparison::WCSimSingleHitComparison( 
            const WCSimSingleHitComparison& other
            ) : 
        fDigit(other.fDigit), fBestFit(other.fBestFit), fCorrectFit(other.fCorrectFit)
    {
        // Empty
    }

    WCSimSingleHitComparison& WCSimSingleHitComparison::operator=(
            const WCSimSingleHitComparison& rhs
            )
    {
        if(&rhs != this)
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
