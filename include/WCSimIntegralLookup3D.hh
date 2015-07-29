/*
 * WCSimIntegralLookup3D.hh
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUP3D_HH_
#define WCSIMINTEGRALLOOKUP3D_HH_
#include "WCSimIntegralLookup.hh"
class TH1F;
class TH3F;

class WCSimIntegralLookupHists3D{
public:
	WCSimIntegralLookupHists3D();
	WCSimIntegralLookupHists3D( TH1F * rhoInt,  TH1F * rhoSInt,  TH1F * rhoSSInt,
			   	   	   	   	  TH3F * rhoGInt, TH3F * rhoGSInt, TH3F * rhoGSSInt,
			   	   	   	   	  TH1F * cutoffS );

	TH1F * GetRhoInt() const;
	TH1F * GetRhoSInt() const;
	TH1F * GetRhoSSInt() const;
	TH3F * GetRhoGInt() const;
	TH3F * GetRhoGSInt() const;
	TH3F * GetRhoGSSInt() const;
	TH1F * GetCutoffS() const;

	void SetHists( TH1F * rhoInt,  TH1F * rhoSInt,  TH1F * rhoSSInt,
				   TH3F * rhoGInt, TH3F * rhoGSInt, TH3F * rhoGSSInt,
				   TH1F * cutoffS );

private:
	TH1F * fRhoInt;
	TH1F * fRhoSInt;
	TH1F * fRhoSSInt;

	TH3F * fRhoGInt;
	TH3F * fRhoGSInt;
	TH3F * fRhoGSSInt;

	TH1F * fCutoffS;
};

class WCSimIntegralLookup3D {
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
