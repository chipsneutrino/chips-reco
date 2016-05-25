/*
 * WCSimIntegralLookup3D.hh
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUP3D_HH_
#define WCSIMINTEGRALLOOKUP3D_HH_
#include "WCSimIntegralLookup.hh"
#include "TH1F.h"
#include <vector>
#include "TSpline.h"
#include "TObject.h"
#include <iostream>
class TH3F;

class WCSimIntegralLookupHistArray : public TObject{
public:
	WCSimIntegralLookupHistArray();
    WCSimIntegralLookupHistArray(const WCSimIntegralLookupHistArray& other);

	/**
	 * Constructor for when you don't already have an array of integrals
	 * graphs and splines and are going to build them yourself
	 * @param R0Max Maximum value of R0 (in cm) ie. the starting radius of the track
	 * @param R0Min Minimum value of R0 (in cm)
	 * @param R0Bins Number of bins in R0
	 * @param cosTh0Max Maximum value (usually 1) of cosTheta0, the initial angle from track direction to PMT
	 * @param cosTh0Min Minimum value of cosTheta0 (usually -1)
	 * @param cosTh0Bins Number of bins in cosTheta0
	 */
	WCSimIntegralLookupHistArray( double R0Max, double R0Min, int R0Bins,
								  double cosTh0Max, double cosTh0Min, int cosTh0Bins);
    WCSimIntegralLookupHistArray& operator=(const WCSimIntegralLookupHistArray& rhs);
    virtual ~WCSimIntegralLookupHistArray();

	void SetBins( double R0Max, double R0Min, int R0Bins,
				   double cosTh0Max, double cosTh0Min, int cosTh0Bins);

	TGraph * GetRhoInt() const{ return fRhoInt; }
	TGraph * GetRhoSInt() const{ return fRhoSInt; }
	TGraph * GetRhoSSInt() const{ return fRhoSSInt; }

	TGraph * GetRhoGInt(double R0, double cosTh0) const;
	TGraph * GetRhoGSInt(double R0, double cosTh0) const;
	TGraph * GetRhoGSSInt(double R0, double cosTh0) const;
	
    TSpline3 * GetRhoGSpline(double R0, double cosTh0) const;
	TSpline3 * GetRhoGSSpline(double R0, double cosTh0) const;
	TSpline3 * GetRhoGSSSpline(double R0, double cosTh0) const;

	void SetRhoG(double R0, double cosTh0, TGraph gr);
	void SetRhoGS(double R0, double cosTh0, TGraph gr);
	void SetRhoGSS(double R0, double cosTh0, TGraph gr);
	void SetRhoSSpline(TSpline3 * spline);
	void SetRhoSSSpline(TSpline3 * spline);

	void SetRhoGSpline(double R0, double cosTh0, TSpline3 * spline);
	void SetRhoGSSpline(double R0, double cosTh0, TSpline3 * spline);
	void SetRhoGSSSpline(double R0, double cosTh0, TSpline3 * spline);

	double GetRhoIntegral(const double E, const int sPower);
	double GetRhoGIntegral(const double E, const double R0, const double cosTh0, const int sPower);
    void Verify();
    void ClearGraphs();

private:
	void ResetArrays();
	void SetGraph(int bin, std::vector<TGraph*>& graphVec, TGraph gr);
	void SetSpline(int bin, std::vector<TSpline3*>& splineVec, TSpline3 * spline);

	unsigned long int GetArrayIndex(double R0, double cosTheta0) const;
	std::vector<TGraph*>   fRhoGArr;
	std::vector<TSpline3*> fRhoGSplineArr;
	std::vector<TGraph*>   fRhoGSArr;
	std::vector<TSpline3*> fRhoGSSplineArr;
	std::vector<TGraph*>   fRhoGSSArr;
	std::vector<TSpline3*> fRhoGSSSplineArr;
	double fR0Min;
	double fR0Max;
	int    fR0Bins;
	double fCosTheta0Min;
	double fCosTheta0Max;
	int    fCosTheta0Bins;

	TGraph * fRhoInt;
	TGraph * fRhoSInt;
	TGraph * fRhoSSInt;
	TSpline3 * fRhoSSpline;
	TSpline3 * fRhoSSSpline;

	ClassDef(WCSimIntegralLookupHistArray,1);
};

class WCSimIntegralLookup3D : public TObject{
public:

	WCSimIntegralLookup3D( TString fileName );
	virtual ~WCSimIntegralLookup3D();

	double GetRhoIntegral( const double &E, const double &s );
	double GetRhoSIntegral( const double &E, const double &s );
	double GetRhoSSIntegral( const double &E, const double &s );

	double GetRhoGIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSSIntegral( const double &E, const double &s, const double &R0, const double &cosTh0);

    void SaveIntegrals(const double E, const double s, const double R0, const double cosTh0);

private:

	double GetCutoffS( const double &E ) const;

	WCSimIntegralLookupHistArray * fIntegrals;
	TFile * fHistFile;

    ClassDef(WCSimIntegralLookup3D, 1);
};


#endif /* WCSIMINTEGRALLOOKUP3D_HH_ */
