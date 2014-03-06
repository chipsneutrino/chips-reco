#ifndef WCSIMLIKELIHOODTUNER_H
#define WCSIMLIKELIHOODTUNER_H

#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"

#include "WCSimLikelihoodTrack.hh"
#include "WCSimRootGeom.hh"



class WCSimLikelihoodTuner
{
    public:
        WCSimLikelihoodTuner(Bool_t calculateIntegrals);
        WCSimLikelihoodTuner(Double_t xMax, Double_t yMax, Double_t zMax, Bool_t calculateIntegrals);
        void Initialize();
        virtual ~WCSimLikelihoodTuner();
        void UpdateDigitArray(WCSimLikelihoodDigitArray * myDigitArray);
   
        void LoadEmissionProfiles(WCSimLikelihoodTrack::TrackType myType);
        void LoadEmissionProfiles(WCSimLikelihoodTrack * myTracK); 
        Double_t TransmissionFunction(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        Double_t Efficiency(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
//		Double_t QuantumEfficiency(WCSimLikelihoodTrack * myTrack);
        Double_t SolidAngle(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        Double_t ScatteringTable(Double_t s);
        std::vector<Double_t> CalculateJ( Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        void CalculateCutoff(WCSimLikelihoodTrack * myTrack);

		    void LoadTabulatedIntegrals(WCSimLikelihoodTrack * myTrack);
        Int_t GetEBin(Double_t energy);
        Int_t GetSBin(Double_t sMax);
        Int_t GetESBin(Double_t energy, Double_t sMax);
        Int_t GetIntegralBin(Double_t energy, Double_t s);
        Int_t GetIntegralBin(Double_t energy, Double_t R0, Double_t cosTheta0, Double_t s);
        Double_t LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, int sPower); 
        std::vector<Double_t> LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit); 
        Double_t LookupIndIntegrals(WCSimLikelihoodTrack * myTrack, Int_t sPower); 
        std::vector<Double_t> LookupIndIntegrals(WCSimLikelihoodTrack * myTrack); 

        Double_t CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack, int i);
        std::vector<Double_t> CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack);
        Double_t CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, int i);
        std::vector<Double_t> CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

        void CalculateCoefficients(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        std::vector<Double_t> CalculateCoefficientsVector(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
   
        void TabulateIndirectIntegrals( WCSimLikelihoodTrack::TrackType myType, TString filename);
        void TabulateDirectIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename);
        void TabulateIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename);

        
    
    protected:
    private:
        Double_t fDirCoeffs[3];
        Double_t fIndCoeffs[3];
        Double_t fExtent[3];
		    Double_t fAverageQE;

        Bool_t fConstrainExtent;

        // Work out where to cut off integrals, and cache the last one
        WCSimLikelihoodTrack * fLastCutoff;
        Double_t fCutoffIntegral;

        // Read in from the emission profile files
        TString * fProfileLocation;
        TFile * fProfiles;
        WCSimLikelihoodTrack::TrackType fIsOpen;
        TObjArray * fHistArray;
        TObjArray * fAngHistArray;
        TObjArray * fFluxArray;
        TH1D * fWhichHisto;

	  
	  	// These are used to speed up reading in the pre-computed integrals
	  	// by ensuring that if the appropriate file is already open, and the appropriate
	  	// tree has already been loaded, we don't waste time getting them again
	  	Bool_t fCalculateIntegrals;
	  	TFile * fRhoIntegralFile;
	  	TFile * fRhoGIntegralFile;
	  	TTree * fRhoIntegralTree;
	  	TTree * fRhoGIntegralTree;
      Int_t fIntegralEnergyBin;
      Double_t fIntegralSMax;
      
      
      // The binning scheme for the table of indirect integrals      
      UInt_t fNBinsRho;
      Int_t fNSBinsRho;
      Int_t fNEBinsRho;
      Double_t fSMinRho;
      Double_t fSMaxRho;
      Double_t fEMinRho;
      Double_t fEMaxRho;


      // The binning scheme for the table of direct integrals
      UInt_t fNBinsRhoG;
      Int_t fNR0Bins;
		  Int_t fNCosTheta0Bins;
		  Int_t fNSBins;
		  Int_t fNEBins;
      Double_t fR0Min;
		  Double_t fR0Max;
		  Double_t fEMin;
		  Double_t fEMax;
      Double_t fSMin;
		  Double_t fSMax;
      Double_t fCosTheta0Min;
      Double_t fCosTheta0Max;
      
      // The table of integrals
		  Double_t * fRhoIntegrals;
		  Double_t * fRhoGIntegrals;

	  	WCSimLikelihoodTrack::TrackType fIntegralParticleType;
	  
};

#endif // WCSIMLIKELIHOODTUNER_H
