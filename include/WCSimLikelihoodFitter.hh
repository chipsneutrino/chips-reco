#ifndef WCSIMLIKELIHOODFITTER_H
#define WCSIMLIKELIHOODFITTER_H

#include "WCSimRootEvent.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimReco.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTotalLikelihood.hh"
#include <vector>
#include <map>

//class WCSimChargeLikelihood;

class WCSimLikelihoodFitter
{
    public:
        WCSimLikelihoodFitter( WCSimRootEvent*);
        virtual ~WCSimLikelihoodFitter();
        Int_t Minimize2LnL(Int_t nTracks);
        Int_t GetNPars(Int_t nTracks);
        Double_t WrapFunc(const Double_t * x);
        Double_t GetMinimum();
        void SeedParams( WCSimReco * myReco );
        std::vector<WCSimLikelihoodTrack> GetBestFit();
        WCSimLikelihoodTrack GetSeedParams();

    protected:
    private:
        WCSimTotalLikelihood * fTotalLikelihood;
        WCSimRootEvent * fRootEvent;
        std::vector<WCSimRecoDigit* > fRecoDigitList;
        WCSimLikelihoodDigitArray * fLikelihoodDigitArray;
        WCSimLikelihoodTrack::TrackType fType;
        std::map<Int_t, Int_t> fParMap;
        std::vector<WCSimLikelihoodTrack> fBestFit;
        WCSimLikelihoodTrack RescaleParams(Double_t x, Double_t y, Double_t z, Double_t t,
                                           Double_t th, Double_t phi, Double_t E,
                                           WCSimLikelihoodTrack::TrackType type);
        

        Double_t fMinimum;
        Double_t fSeedVtxX;
        Double_t fSeedVtxY;
        Double_t fSeedVtxZ;
        Double_t fSeedTheta;
        Double_t fSeedPhi;
        Bool_t fIsFirstCall;
	
    ClassDef(WCSimLikelihoodFitter,1);
};

#endif // WCSIMLIKELIHOODFITTER_H
