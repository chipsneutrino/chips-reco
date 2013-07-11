#ifndef WCSIMLIKELIHOODTUNER_H
#define WCSIMLIKELIHOODTUNER_H

#include "TFile.h"
#include "TString.h"

#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimRootGeom.hh"

class WCSimLikelihoodTuner
{
    public:
        WCSimLikelihoodTuner();
        virtual ~WCSimLikelihoodTuner();
   
        void LoadEmissionProfiles(WCSimLikelihoodTrack::TrackType myType);
        void LoadEmissionProfiles(WCSimLikelihoodTrack * myTracK); 
        Double_t TransmissionFunction(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        Double_t Efficiency(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        Double_t SolidAngle(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        Double_t ScatteringTable(Double_t s);

        Double_t CalculateIndIntegrals(Int_t i, WCSimLikelihoodTrack * myTrack);
        Double_t CalculateChIntegrals(Int_t i, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

        void CalculateCoefficients(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        std::vector<Double_t> CalculateCoefficientsVector(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
   
        void TabulateIndirectIntegrals( WCSimLikelihoodTrack::TrackType myType, TString filename);
        void TabulateDirectIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename);
        void TabulateIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename);

        
    
    protected:
    private:
        Double_t fDirCoeffs[3];
        Double_t fIndCoeffs[3];
        TString * fProfileLocation;
        TFile * fProfiles;
        WCSimLikelihoodTrack::TrackType fIsOpen;
        TObjArray * fHistArray;
        TObjArray * fAngHistArray;
        TObjArray * fFluxArray;
};

#endif // WCSIMLIKELIHOODTUNER_H
