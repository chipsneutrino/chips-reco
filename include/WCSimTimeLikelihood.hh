/**
 * \class WCSimTimeLikelihood
 * This class is used to calculate the contribution to the
 * total log-ikelihood due to the timing information of the PMT
 * hits.
 */

#ifndef WCSIMTIMELIKELIHOOD_H
#define WCSIMTIMELIKELIHOOD_H

#include <vector>
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimChargeLikelihood.hh"
#include "WCSimRecoEvent.hh"

#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"

class WCSimTimeLikelihood : public TObject
{
  public:
    /**
     * Constructor
     * @param myDigitArray The PMT responses for a single event
     * @param myChargeLikelihood Charge likelihood object to get charge prediction from
     */
    //TODO: do we need charge likelihood in the constructor?
    WCSimTimeLikelihood( WCSimLikelihoodDigitArray * myDigitArray, WCSimChargeLikelihood *myChargeLikelihood);

    virtual ~WCSimTimeLikelihood();

    /**
     * Add another track to calculate the likelihood for several particles at once
     * @param myTrack Track object to add
     */
    void AddTrack( WCSimLikelihoodTrack * myTrack);

    /**
     * Set all the tracks that contribute to the likelihood at once
     * @param myTrack Vector of all the track objects to consider
     */
    void SetTracks( std::vector<WCSimLikelihoodTrack*> myTrack );

    /// Remove all the tracks currently loaded
    void ClearTracks();

    /**
     * Replace the current event with hits from a different one
     * @param myDigitArray New array of PMT responses
     */
    void UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray);

    /**
     * Calculate -2 log(likelihood) for the current PMT response,
     * given the current set of hypothesised tracks
     * @return -2 log(likelihood)
     */
    Double_t Calc2LnL();

    Double_t CorrectedTime( WCSimLikelihoodTrack * myTrack, Double_t primaryTime );
    Double_t GetPredictedCharge( WCSimLikelihoodTrack * myTrack );

    //TODO: time likelihood is not additive! work out how it behaves for >1 tracks
    Double_t TimeLikelihood( WCSimLikelihoodTrack * myTrack, Double_t correctedTime );
    Double_t TimeLikelihood( WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit* myDigit, Double_t correctedTime );

    void GetExternalVariables( const char *fName );
    //void GetLikelihoodParameters(); //???

    //TODO: update comment
    /**
     * Something something (anything?)
     * @param myTrack The charged particle track
     */
    void GetTrackParameters(WCSimLikelihoodTrack * myTrack);

  protected:
  private:
    /**
     * Called by constructor.  Initializes vectors to hold the track integrals, and
     * the WCSimLikelihoodTuner and WCSimDigitizerLikelihood member variables
     * @param myDigitArray PMT responses for this event
     */
    void Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimChargeLikelihood *myChargeLikelihood);

    //XXX: debug stuff
    TFile *fDebugFile;
    TH1F *fDebugHist9;

    //Energy of the track lepton
    Double_t fEnergy;

    //Energy dependent variables extracted from file
    Double_t fTrackMidpoint;  //track midpoint (middle)
    //coefficients to get the above value
    Double_t fTrackMidpointCoeffs[2];
    //table of tertiary parameters (for energy values -> charge bins)
    Double_t fTertiaryParameters[7][5][4];
    //table of cuts defining the charge bins
    Int_t fSizeChargeCuts;  //size of cuts array
    Int_t fNumChargeCuts;   //for convenience: size-1
    Double_t fChargeCuts[1000]; //probably not more

    //table of secondary parameters (for charge bins)
    Double_t fSecondaryParameters[7][5];

    //The time likelihood function
    static double fFunctionForm(double *x, double *par);
    TF1 *fLikelihoodFunction;
    //norm factors of the likelihood function array (one for each charge bin)
    Double_t *fNormFuncParams;

    //The charge dependent parameter function
    TF1 *fChargeParameterFunction;
    //The energy dependent parameter function
    TF1 *fEnergyParameterFunction;

    //charge likelihood object to get the predicted charge
    WCSimChargeLikelihood *fChargeLikelihood;

    // The track and event parameters for which we calculate the likelihood
    std::vector<WCSimLikelihoodTrack *>   fTracks;     ///< Vector of simultaneous tracks contributing to the likelihood
    WCSimLikelihoodDigitArray           * fDigitArray; ///< Response for all the detector's PMTs for this event
    WCSimLikelihoodDigit                * fDigit;      ///< The PMT being considered

    Bool_t fGotTrackParameters;


    //TODO: needed?
    ClassDef(WCSimTimeLikelihood,0)
};

#endif // WCSIMTIMELIKELIHOOD_H
