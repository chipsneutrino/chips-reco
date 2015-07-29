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
#include "WCSimIntegralLookup.hh"
#include "TObject.h"
class TH1F;
class TH3F;

class WCSimIntegralLookupMaker3D : public TObject {
public:
	WCSimIntegralLookupMaker3D( TrackType::Type particle,
							  int nR0Bins, double R0Min, double R0Max,
							  int nCosTh0Bins, double cosTh0Min, double cosTh0Max );

	virtual ~WCSimIntegralLookupMaker3D();
	void Run(TString fileName);

protected:
	void SetBins(const int &nEBins, const double &eMin, const double &eMax,
				 const int &nSBins, const double &sMin, const double &sMax,
				 const int &nR0Bins, const double &R0Min, const double &R0Max,
				 const int &nCosTh0Bins, const double &cosTh0Min, const double &cosTh0Max);

	void MakeLookupTables();
	void MakeRhoTables();
	void MakeRhoGTables();

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
	WCSimIntegralLookupHists3D fIntegrals;
	TH1F * fRhoInt;
	TH1F * fRhoSInt;
	TH1F * fRhoSSInt;

	TH3F * fRhoGInt;
	TH3F * fRhoGSInt;
	TH3F * fRhoGSSInt;

  ClassDef(WCSimIntegralLookupMaker3D,1)
};

#endif /* WCSIMINTEGRALLOOKUPMAKER3D_HH_ */
