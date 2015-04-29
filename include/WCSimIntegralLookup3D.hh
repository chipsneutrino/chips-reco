/*
 * WCSimIntegralLookup3D.hh
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUP3D_HH_
#define WCSIMINTEGRALLOOKUP3D_HH_
#include "WCSimIntegralLookup.hh"

class WCSimIntegralLookup3D : public WCSimIntegralLookup {
public:

	WCSimIntegralLookup3D( TString fileName );
	virtual ~WCSimIntegralLookup3D();

	double GetRhoIntegral( const double &E, const double &s );
	double GetRhoSIntegral( const double &E, const double &s );
	double GetRhoSSIntegral( const double &E, const double &s );

	double GetRhoGIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSSIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);

  THnSparseF * GetRhoIntegralHist();

private:

	double GetIntegral1D( const double &E, THnSparseF * hist);
	double GetIntegral3D( const double &E, const double &R0, const double &cosTh0, THnSparseF * hist);
	double GetCutoffS( const double &E ) const;

	WCSimIntegralLookupHists fIntegrals;
	TFile * fHistFile;

};


#endif /* WCSIMINTEGRALLOOKUP3D_HH_ */
