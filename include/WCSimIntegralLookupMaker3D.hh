/*
 * WCSimIntegralLookupMaker3D.hh
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUPMAKER3D_HH_
#define WCSIMINTEGRALLOOKUPMAKER3D_HH_
#include "WCSimEmissionProfileManager.hh"
#include "WCSimIntegralLookup3D.hh"
#include "WCSimTrackParameterEnums.hh"
#include "TObject.h"
class TH1F;
class TH3F;

class WCSimIntegralLookupMaker3D: public TObject {
	public:
		WCSimIntegralLookupMaker3D(TrackType::Type particle, int nR0Bins, double R0Min, double R0Max, int nCosTh0Bins,
				double cosTh0Min, double cosTh0Max);

		virtual ~WCSimIntegralLookupMaker3D();
		void Run(TString fileName);

	protected:
		void SetBins(const int &nEBins, const double &eMin, const double &eMax, const int &nSBins, const double &sMin,
				const double &sMax, const int &nR0Bins, const double &R0Min, const double &R0Max,
				const int &nCosTh0Bins, const double &cosTh0Min, const double &cosTh0Max);
		void CheckBins(const int nBins, const double min, const double max);

		void MakeLookupTables();
		void MakeRhoTables();
		void MakeRhoGTables();

		void SmoothLookupTables();
		void SmoothRhoTables();
		void SmoothRhoGTables();

		/**
		 * @brief Function to take a TGraph and return a smoothed version of it.
		 *        The method is to construct a TH1D using the graph points, call
		 *        TH1::Smooth and then fill a new graph with the smoothed histogram
		 *        contents
		 *
		 * @param graph The graph to be smoothed
		 *
		 * @return A smoothed graph
		 */
		TGraph SmoothGraph(TGraph * graph, int nTimes = 1);

		/**
		 * @brief Replace all the (x,y) points of a TGraph with those of a different one
		 *        Will clear the first graph and then set all the points to those
		 *        from the second graph, leaving all the titles, draw options etc. intact
		 *
		 * @param original The graph whose points are to be replaced
		 * @param replacement The graph holding the replacements points
		 */
		void ReplaceGraphPoints(TGraph * original, const TGraph& replacement);

		void MakeSplines();
		void MakeRhoSplines();
		void MakeRhoGSplines();

		void SaveLookupTables(TString fileName);

	private:
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
		WCSimIntegralLookupHistArray fIntegrals;
		TH1F * fRhoInt;
		TH1F * fRhoSInt;
		TH1F * fRhoSSInt;

		TH3F * fRhoGInt;
		TH3F * fRhoGSInt;
		TH3F * fRhoGSSInt;

		ClassDef(WCSimIntegralLookupMaker3D,1)
};

#endif /* WCSIMINTEGRALLOOKUPMAKER3D_HH_ */
