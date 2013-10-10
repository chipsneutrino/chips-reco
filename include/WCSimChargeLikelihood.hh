#ifndef WCSIMCHARGELIKELIHOOD_H
#define WCSIMCHARGELIKELIHOOD_H

#include <vector>
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimRecoEvent.hh"

#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"

class WCSimChargeLikelihood
{
    public:
        WCSimChargeLikelihood( WCSimLikelihoodDigitArray * myDigitArray, Bool_t calculateIntegrals = false);
        virtual ~WCSimChargeLikelihood();
        void Initialize( WCSimLikelihoodDigitArray * myDigitArray, Bool_t calculateIntegrals = false );
        void AddTrack( WCSimLikelihoodTrack * myTrack);
        void SetTracks( std::vector<WCSimLikelihoodTrack*> myTrack );
        void UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray);

        Double_t Calc2LnL();
        Double_t CalculateExactLikelihood(); // Calc2LnL can (optionally) use tabulated integrals, this always calculates them by hand
        Double_t ChargeExpectation(WCSimLikelihoodTrack * myTrack);
        Double_t GetLightFlux(WCSimLikelihoodTrack * myTrack);
        Double_t GetMuIndirect(WCSimLikelihoodTrack * myTrack);
        Double_t GetMuDirect(WCSimLikelihoodTrack * myTrack);
        Double_t ScatteringTable();
        void GetTrackParameters(WCSimLikelihoodTrack * myTrack);

        Double_t DigitizerMinus2LnL( const Double_t &myUndigiHit, const Double_t &myDigiHit );
        Double_t DigitizerLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit );
        Double_t DigitizePickerLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit );
        Double_t DigitizeGausExpoLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit );
        Double_t DigitizeGausLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit );
        Double_t DigitizerExpectation( const Double_t & myUndigiHit );

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

      // Use these for the scattering table, ie. to relate the strength of scattered to direct light
        Double_t fRadius;   // Distance from centre of tank to the source
        Double_t fAngle;    // Angle between source, centre of tank, and PMT
        Double_t fTheta;    // Angle between source direction and ray to PMT (same as theta0)
        Double_t fPhi;       // Angle between plane containing tank centre, PMT and source, and plane
                            // containing the track and the tank centre

      // The track and event parameters for which we calculate the likelihood
        std::vector<WCSimLikelihoodTrack *> fTracks;
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

	    // The digitizer samples the 1pe distribution repeatedly, then applies a threshold function,
	    // then multiplies by an efficiency term
	    Double_t fEfficiency;

      // For < 10 predicted photons, we look up the digitizer probability from a 2D histogram;
      TH2D * fDigiPDF;


};

#endif // WCSIMCHARGELIKELIHOOD_H
