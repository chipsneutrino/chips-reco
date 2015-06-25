/*
 * WCSimIntegralLookupMaker3D.hh
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUPMAKER3D_HH_
#define WCSIMINTEGRALLOOKUPMAKER3D_HH_
#include "WCSimIntegralLookupMaker.hh"
#include "WCSimTrackParameterEnums.hh"

class WCSimIntegralLookupMaker3D : public WCSimIntegralLookupMaker {
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
	
  
  ClassDef(WCSimIntegralLookupMaker,1)
};

#endif /* WCSIMINTEGRALLOOKUPMAKER3D_HH_ */
