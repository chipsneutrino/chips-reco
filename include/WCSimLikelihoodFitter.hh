#ifndef WCSIMLIKELIHOODFITTER_H
#define WCSIMLIKELIHOODFITTER_H

#include "WCSimRootEvent.hh"
#include "WCSimChargeLikelihood.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodTrack.hh"
#include <vector>

class WCSimLikelihoodFitter
{
    public:
        WCSimLikelihoodFitter( WCSimRootEvent*);
        virtual ~WCSimLikelihoodFitter();
        void Minimize2LnL();
        double Calc2LnL();
        double Charge2LnL();
        double CalcMuDirect();
        double CalcMuIndirect();
        double Time2LnL();
    protected:
    private:
        WCSimRootEvent * fRootEvent;
        WCSimLikelihoodTrack * fTrack;
        std::vector<WCSimRecoDigit* > fRecoDigitList;
        WCSimLikelihoodDigitArray * fLikelihoodDigitArray;
};

#endif // WCSIMLIKELIHOODFITTER_H
