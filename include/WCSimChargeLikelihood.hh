#ifndef WCSIMCHARGELIKELIHOOD_H
#define WCSIMCHARGELIKELIHOOD_H

#include <vector>
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimRecoEvent.hh"
#include "TMath.h"

class WCSimChargeLikelihood
{
    public:
        WCSimChargeLikelihood( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack);
        virtual ~WCSimChargeLikelihood();
        void Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack );

        Double_t Calc2LnL();
        Double_t ChargeExpectation();
        Double_t GetLightFlux();
        Double_t GetMuIndirect();
        Double_t GetMuDirect();
        Double_t LookupChIntegrals(int i); 
        Double_t LookupIndIntegrals(int i); 
        Double_t ScatteringTable();
        void GetTrackParameters();
        Double_t CalculateExactLikelihood();
        Double_t GetR0Interval();
        Double_t GetCosTheta0Interval();

    protected:
    private:
        std::vector<Double_t> fCoeffsCh;
        std::vector<Double_t> fCoeffsInd;
        // Tuning parameters from a parabolic fit to solid angle * transmission * PMT acceptance
        // all as a function of distance along the track that light was emitted from


        /////////////////////////////////////////////
        // VARIABLE NAMES FOLLOW arXiv.0902.2222v2 //
        /////////////////////////////////////////////

        // Use these for integral lookup tables:
        Double_t fEnergy; // particle energy
        Double_t fR0;     // vertex to PMT distance
        Double_t fCosTheta0; // angle to PMT as viewed from vertex
        Double_t fEta;    // angle of incidence of emitted light at the PMT

        Double_t fEnergyInterval; // steps to increment variables by
        Double_t fR0Interval;             // from one bin of the lookup table
        Double_t fCosTheta0Interval;      // to the next
        Double_t fChIntegral[100][1000][100];
        Double_t fIndIntegral[100];

      // Use these for the scattering table, ie. to relate the strength of scattered to direct light
        Double_t fRadius;   // Distance from centre of tank to the source
        Double_t fAngle;    // Angle between source, centre of tank, and PMT
        Double_t fTheta;    // Angle between source direction and ray to PMT (same as theta0)
        Double_t fPhi;       // Angle between plane containing tank centre, PMT and source, and plane
                            // containing the track and the tank centre

      // The track and event parameters for which we calculate the likelihood
        WCSimLikelihoodTrack * fTrack;
        WCSimLikelihoodDigitArray * fDigitArray;
        WCSimLikelihoodDigit * fDigit;
        WCSimLikelihoodTuner * fTuner;

      // The fitted functions are defined using various variables that relate the track to the 
      // PMT hit in question.  I calculate these in GetTrackParameters() and set a flag when
      // this has been done
      Bool_t fGotTrackParameters;

      // Do we want to calculate or look up the integrals for the likelihood?
      // Calculating is fine for drawing/testing, but look them up when we do the fit
      Bool_t fCalculateIntegrals;

};

#endif // WCSIMCHARGELIKELIHOOD_H
