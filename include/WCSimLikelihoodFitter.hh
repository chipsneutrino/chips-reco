#ifndef WCSIMLIKELIHOODFITTER_H
#define WCSIMLIKELIHOODFITTER_H

#include "WCSimRootEvent.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimChargeLikelihood.hh"
#include <vector>

//class WCSimChargeLikelihood;

class WCSimLikelihoodFitter
{
    public:
        WCSimLikelihoodFitter( WCSimRootEvent*);
        virtual ~WCSimLikelihoodFitter();
        void Minimize2LnL();
        double Calc2LnL();
        double Calc2LnL(WCSimLikelihoodTrack * myTrack);
        double Charge2LnL();
        double Charge2LnL(WCSimLikelihoodTrack * myTrack);
        double CalcMuDirect();
        double CalcMuIndirect();
        double Time2LnL();
    protected:
    private:
        WCSimChargeLikelihood *myLikelihood;
        WCSimRootEvent * fRootEvent;
        WCSimLikelihoodTrack * fTrack;
        std::vector<WCSimRecoDigit* > fRecoDigitList;
        WCSimLikelihoodDigitArray * fLikelihoodDigitArray;
};

#endif // WCSIMLIKELIHOODFITTER_H
