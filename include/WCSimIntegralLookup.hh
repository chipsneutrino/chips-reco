/*
 * WCSimIntegralLookup.hh
 *
 *  Created on: 9 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUP_HH_
#define WCSIMINTEGRALLOOKUP_HH_

#include "WCSimEmissionProfiles.hh"
#include "THnSparse.h"
#include <exception>

class TH1F;
class TFile;


class WCSimIntegralLookupHists{
public:
	WCSimIntegralLookupHists();
	WCSimIntegralLookupHists( THnSparseF * rhoInt,  THnSparseF * rhoSInt,  THnSparseF * rhoSSInt,
			   	   	   	   	  THnSparseF * rhoGInt, THnSparseF * rhoGSInt, THnSparseF * rhoGSSInt,
			   	   	   	   	  TH1F * cutoffS );

	THnSparseF * GetRhoInt() const;
	THnSparseF * GetRhoSInt() const;
	THnSparseF * GetRhoSSInt() const;
	THnSparseF * GetRhoGInt() const;
	THnSparseF * GetRhoGSInt() const;
	THnSparseF * GetRhoGSSInt() const;
	TH1F * GetCutoffS() const;

	void SetHists( THnSparseF * rhoInt,  THnSparseF * rhoSInt,  THnSparseF * rhoSSInt,
				   THnSparseF * rhoGInt, THnSparseF * rhoGSInt, THnSparseF * rhoGSSInt,
				   TH1F * cutoffS );

private:
	THnSparseF * fRhoInt;
	THnSparseF * fRhoSInt;
	THnSparseF * fRhoSSInt;

	THnSparseF * fRhoGInt;
	THnSparseF * fRhoGSInt;
	THnSparseF * fRhoGSSInt;

	TH1F * fCutoffS;
};

class WCSimIntegralLookup {
public:

	WCSimIntegralLookup( TString fileName );
	WCSimIntegralLookup();
	virtual ~WCSimIntegralLookup();

	virtual double GetRhoIntegral( const double &E, const double &s);
	virtual double GetRhoSIntegral( const double &E, const double &s);
	virtual double GetRhoSSIntegral( const double &E, const double &s);

	virtual double GetRhoGIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);
	virtual double GetRhoGSIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);
	virtual double GetRhoGSSIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);

  virtual THnSparseF * GetRhoIntegralHist();

private:

	double GetIntegral2D( const double &E, const double &s, THnSparseF * hist);
	double GetIntegral4D( const double &E, const double &s, const double &R0, const double &cosTh0, THnSparseF * hist);
	double GetCutoffS( const double &E ) const;

	WCSimIntegralLookupHists fIntegrals;
	TFile * fHistFile;

};

#endif /* WCSIMINTEGRALLOOKUP_HH_ */
