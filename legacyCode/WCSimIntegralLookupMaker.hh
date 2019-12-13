/*
 * WCSimIntegralLookupMaker.hh
 *
 *  Created on: 9 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUPMAKER_HH_
#define WCSIMINTEGRALLOOKUPMAKER_HH_

#include "WCSimEmissionProfileManager.hh"
#include "WCSimIntegralLookup.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimLikelihoodTrack.hh"
#include "THnSparse.h"
#include "TObject.h"
#include "TString.h"

class TH1F;

class WCSimIntegralLookupMaker: public TObject {
	public:
		WCSimIntegralLookupMaker(TrackType::Type particle, int nR0Bins, double R0Min, double R0Max, int nCosTh0Bins,
				double cosTh0Min, double cosTh0Max);
		WCSimIntegralLookupMaker();

		virtual ~WCSimIntegralLookupMaker();
		void Run(TString fileName);

	protected:
		void SetBins(const int &nEBins, const double &eMin, const double &eMax, const int &nSBins, const double &sMin,
				const double &sMax, const int &nR0Bins, const double &R0Min, const double &R0Max,
				const int &nCosTh0Bins, const double &cosTh0Min, const double &cosTh0Max);
		void CheckBins(const int nBins, const double min, const double max);

		void MakeLookupTables();
		void MakeCutoffS();
		void MakeRhoTables();
		void MakeRhoGTables();

		void SaveLookupTables(TString fileName);

		// Binning information
		int fNEBins;
		double fEMin;
		double fEMax;

		int fNR0Bins;
		double fR0Min;
		double fR0Max;

		int fNCosTh0Bins;
		double fCosTh0Min;
		double fCosTh0Max;

		TrackType::Type fType;
		WCSimEmissionProfileManager fEmissionProfileManager;

		// The integral histograms:
		WCSimIntegralLookupHists fIntegrals;
		THnSparseF * fRhoInt;
		THnSparseF * fRhoSInt;
		THnSparseF * fRhoSSInt;

		THnSparseF * fRhoGInt;
		THnSparseF * fRhoGSInt;
		THnSparseF * fRhoGSSInt;

	private:
		int fNSBins;
		double fSMin;
		double fSMax;

		TH1F * fCutoffS;

		ClassDef(WCSimIntegralLookupMaker,1)

};

#endif /* WCSIMINTEGRALLOOKUPMAKER_HH_ */