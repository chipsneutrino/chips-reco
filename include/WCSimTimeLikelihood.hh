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

class WCSimTimeLikelihood
{
  public:
    WCSimTimeLikelihood( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack, WCSimChargeLikelihood *myChargeLikelihood);
    virtual ~WCSimTimeLikelihood();
    void Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack, WCSimChargeLikelihood *myChargeLikelihood);

    Double_t Calc2LnL();
    Double_t CorrectedTime( Double_t primaryTime );
    Double_t TimeLikelihood( Double_t correctedTime );
    Double_t TimeLikelihood( WCSimLikelihoodDigit* myDigit, Double_t correctedTime );
    Double_t GetPredictedCharge();
    void GetExternalVariables( const char *fName );
    void GetLikelihoodParameters();
    void GetTrackParameters();

  protected:
  private:

    //XXX: debug 
    TFile *fDebugFile;
    TH1F *fDebugHist9;

    //PMT-vertex variables
    //Double_t fPmtX, fPmtY, fPmtZ; //pmt position
    //Double_t fTrkX, fTrkY, fTrkZ; //track vertex position
    //Double_t fMidX, fMidY, fMidZ; //track midpoint position
    //Double_t fPhi, fTheta;        //track direction
    //Double_t fTrkT;               //track time

    //Energy of the track lepton
    Double_t fEnergy;

    //Energy dependent variables extracted from file
    Double_t fTrackMidpoint;  //track midpoint (middle)
    //table of secondary parameters (for charge bins)
    Double_t fSecondaryParameters[7][5];
    //table of cuts defining the charge bins
    Int_t fSizeChargeCuts;  //size of cuts array
    Int_t fNumChargeCuts;   //for convenience: size-1
    Double_t fChargeCuts[1000]; //probably not more

    //The time likelihood function
    static double fFunctionForm(double *x, double *par);
    TF1 *fLikelihoodFunction;

    //The charge dependent parameter function
    TF1 *fChargeParameterFunction;

    //charge likelihood object to get the predicted charge
    WCSimChargeLikelihood *fChargeLikelihood;

    // The track and event parameters for which we calculate the likelihood
    WCSimLikelihoodTrack * fTrack;
    WCSimLikelihoodDigitArray * fDigitArray;
    WCSimLikelihoodDigit * fDigit;

    //TODO: do I really need this and the GetTrackParameters function?
    Bool_t fGotTrackParameters;

};

#endif // WCSIMTIMELIKELIHOOD_H
