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
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimRecoEvent.hh"

#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"

class WCSimTimeLikelihood: public TObject {
	public:
		/**
		 * Constructor
		 * @param myDigitArray The PMT responses for a single event
		 * @param myChargeLikelihood Charge likelihood object to get charge prediction from
		 */
		WCSimTimeLikelihood(WCSimLikelihoodDigitArray * myDigitArray);

		virtual ~WCSimTimeLikelihood();

		/**
		 * Add another track to calculate the likelihood for several particles at once
		 * @param myTrack Track object to add
		 */
		void AddTrack(WCSimLikelihoodTrackBase * myTrack);

		/**
		 * Set all the tracks that contribute to the likelihood at once
		 * @param myTracks Vector of all the track objects to consider
		 */
		void SetTracks(std::vector<WCSimLikelihoodTrackBase*> myTracks);

		/// Remove all the tracks currently loaded
		void ClearTracks();

		/**
		 * Replace the current event with hits from a different one
		 * @param myDigitArray New array of PMT responses
		 */
		void UpdateDigitArray(WCSimLikelihoodDigitArray * myDigitArray);

		/**
		 * Calculate -2 log(likelihood) for the given PMT response,
		 * given the current set of hypothesised tracks
		 * @param myDigit Current PMT to give result for
		 * @param predictedCharges Vector of charge predictions for fTracks
		 * @return -2 log(likelihood)
		 */
		Double_t Calc2LnL(WCSimLikelihoodDigit* myDigit, std::vector<Double_t> predictedCharges);

		Double_t CorrectedTime(Int_t trackIndex, Double_t primaryTime);

		//TODO: time likelihood is not additive! work out how it behaves for >1 tracks
		Double_t TimeLikelihood(Int_t trackIndex, Double_t correctedTime, std::vector<Double_t> predictedCharges);
		Double_t TimeLikelihood(Int_t trackIndex, WCSimLikelihoodDigit* myDigit, Double_t correctedTime,
				std::vector<Double_t> predictedCharges);

		void GetExternalVariables(const char *fName);
		//void GetLikelihoodParameters(); //???

		//TODO: update comment
		/**
		 * Something something (anything?)
		 * @param trackIndex Index of the charged particle track in fTracks vector
		 */
		void GetTrackParameters(Int_t trackIndex);

	protected:
	private:
		/**
		 * Called by constructor. Constructs parameter functions.
		 * @param myDigitArray PMT responses for this event
		 * @param myChargeLikelihood charge likelihood object to get the charge prediction
		 */
		void Initialize(WCSimLikelihoodDigitArray * myDigitArray);

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

		// The track and event parameters for which we calculate the likelihood
		std::vector<WCSimLikelihoodTrackBase *> fTracks; ///< Vector of simultaneous tracks contributing to the likelihood
		WCSimLikelihoodDigitArray * fDigitArray; ///< Response for all the detector's PMTs for this event
		WCSimLikelihoodDigit * fDigit;      ///< The PMT being considered

		/// Flag to check if GetTrackParameters() has been called to calculate
		/// the integral coefficients. Shows index of the track in fTracks
		/// vector or is negative if track parameters have not been calculated
		Int_t fGotTrackParameters;

		//TODO: needed?
		ClassDef(WCSimTimeLikelihood,0)
};

#endif // WCSIMTIMELIKELIHOOD_H
