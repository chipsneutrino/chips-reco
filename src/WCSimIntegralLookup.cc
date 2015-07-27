/*
 * WCSimIntegralLookup.cc
 *
 *  Created on: 9 Mar 2015
 *      Author: ajperch
 */

#include "WCSimIntegralLookup.hh"
#include "THnSparse.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TString.h"
#include <stdexcept>

///////////////////////////////////////////////////////////////////////////////////////////////
// Class to hold all the histograms so that the reader and writer code don't get out of sync
///////////////////////////////////////////////////////////////////////////////////////////////
WCSimIntegralLookupHists::WCSimIntegralLookupHists() {
	fRhoInt = 0x0;
	fRhoSInt = 0x0;
	fRhoSSInt = 0x0;

	fRhoGInt = 0x0;
	fRhoGSInt = 0x0;
	fRhoGSSInt = 0x0;

	fCutoffS = 0x0;
}

WCSimIntegralLookupHists::WCSimIntegralLookupHists(THnSparseF* rhoInt, THnSparseF* rhoSInt, THnSparseF* rhoSSInt,
		THnSparseF* rhoGInt, THnSparseF* rhoGSInt, THnSparseF* rhoGSSInt, TH1F* cutoffS) {
	SetHists(rhoInt, rhoSInt, rhoSSInt, rhoGInt, rhoGSInt, rhoGSSInt, cutoffS);

}

THnSparseF* WCSimIntegralLookupHists::GetRhoInt() const {
	return fRhoInt;
}

THnSparseF* WCSimIntegralLookupHists::GetRhoSInt() const {
	return fRhoSInt;
}

THnSparseF* WCSimIntegralLookupHists::GetRhoSSInt() const {
	return fRhoSSInt;
}

THnSparseF* WCSimIntegralLookupHists::GetRhoGInt() const {
	return fRhoGInt;
}

THnSparseF* WCSimIntegralLookupHists::GetRhoGSInt() const {
  // std::cout << fRhoGSInt << std::endl;
  // std::cout << fRhoGSInt->GetEntries() << std::endl;
  // std::cout << "Returning fRhoGSInt" << std::endl;
	return fRhoGSInt;
}

THnSparseF* WCSimIntegralLookupHists::GetRhoGSSInt() const {
	return fRhoGSSInt;
}

TH1F* WCSimIntegralLookupHists::GetCutoffS() const {
	return fCutoffS;
}

void WCSimIntegralLookupHists::SetHists(THnSparseF* rhoInt, THnSparseF* rhoSInt, THnSparseF* rhoSSInt,
		THnSparseF* rhoGInt, THnSparseF* rhoGSInt, THnSparseF* rhoGSSInt, TH1F* cutoffS) {

	fRhoInt = rhoInt;
	fRhoSInt = rhoSInt;
	fRhoSSInt = rhoSSInt;

	fRhoGInt = rhoGInt;
	fRhoGSInt = rhoGSInt;
	fRhoGSSInt = rhoGSSInt;

	fCutoffS = cutoffS;
  // std::cout << "Finished setting hists" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Now the code used to read the hists
//////////////////////////////////////////////////////////////////////////////////////////////

WCSimIntegralLookup::~WCSimIntegralLookup() {
	// TODO Auto-generated destructor stub
}

WCSimIntegralLookup::WCSimIntegralLookup()
{
}

WCSimIntegralLookup::WCSimIntegralLookup(TString fileName) {

	fHistFile = 0x0;
	fHistFile = new TFile(fileName.Data(), "READ");
    if (!fHistFile)
    {

        throw std::runtime_error("Could not open file");
    }
    if( fHistFile->IsZombie())
    {
    	throw std::runtime_error("Integral lookup TFile is a zombie");
    }

	THnSparseF * rhoInt = 0x0;
	THnSparseF * rhoSInt = 0x0;
	THnSparseF * rhoSSInt = 0x0;
	THnSparseF * rhoGInt = 0x0;
	THnSparseF * rhoGSInt = 0x0;
	THnSparseF * rhoGSSInt = 0x0;
	TH1F * cutoffS = 0x0;

	fHistFile->GetObject("fRhoInt",rhoInt);
	fHistFile->GetObject("fRhoSInt",rhoSInt);
	fHistFile->GetObject("fRhoSSInt",rhoSSInt);
	fHistFile->GetObject("fRhoGInt",rhoGInt);
	fHistFile->GetObject("fRhoGSInt",rhoGSInt);
	fHistFile->GetObject("fRhoGSSInt",rhoGSSInt);
	fHistFile->GetObject("fCutoffS",cutoffS);

	if( rhoInt == 0x0 ) { throw std::runtime_error("Could not get rhoInt from file"); }
	if( rhoSInt == 0x0 ) { throw std::runtime_error("Could not get rhoSInt from file"); }
	if( rhoSSInt == 0x0 ) { throw std::runtime_error("Could not get rhoSSInt from file"); }
	if( rhoGInt == 0x0 ) { throw std::runtime_error("Could not get rhoGInt from file"); }
	if( rhoGSInt == 0x0 ) { throw std::runtime_error("Could not get rhoGSInt from file"); }
	if( rhoGSSInt == 0x0 ) { throw std::runtime_error("Could not get rhoGSSInt from file"); }
	if( cutoffS == 0x0 ) { throw std::runtime_error("Could not get cutoffS from file"); }

	fIntegrals.SetHists(rhoInt, rhoSInt, rhoSSInt, rhoGInt, rhoGSInt, rhoGSSInt, cutoffS);

}

double WCSimIntegralLookup::GetRhoIntegral(const double &E, const double &s) {
	return GetIntegral2D(E, s, fIntegrals.GetRhoInt());
}

double WCSimIntegralLookup::GetRhoSIntegral(const double &E, const double &s) {
	return GetIntegral2D(E, s, fIntegrals.GetRhoSInt());
}

double WCSimIntegralLookup::GetRhoSSIntegral(const double &E, const double &s) {
	return GetIntegral2D(E, s, fIntegrals.GetRhoSSInt());
}

double WCSimIntegralLookup::GetRhoGIntegral(const double &E, const double &s, const double &R0, const double &cosTh0) {
	return GetIntegral4D(E, s, R0, cosTh0, fIntegrals.GetRhoGInt());
}

double WCSimIntegralLookup::GetRhoGSIntegral(const double &E, const double &s, const double &R0, const double &cosTh0) {
	return GetIntegral4D(E, s, R0, cosTh0, fIntegrals.GetRhoGSInt());
}

double WCSimIntegralLookup::GetRhoGSSIntegral(const double &E, const double &s, const double &R0, const double &cosTh0) {
	return GetIntegral4D(E, s, R0, cosTh0, fIntegrals.GetRhoGSSInt());
}

THnSparseF * WCSimIntegralLookup::GetRhoIntegralHist()
{
  return fIntegrals.GetRhoInt();
}

double WCSimIntegralLookup::GetIntegral2D(const double &E, const double &s, THnSparseF * hist) {
	double cutoff = GetCutoffS(E);
	double x[2] = {E, s};
	if( s > cutoff) { x[1] = cutoff; }

	int bin = hist->GetBin(x, false); // Don't allocate it if it's empty

	if( bin != -1){ return ( hist->GetBinContent(bin) ); }
	else { return 0; }
}

double WCSimIntegralLookup::GetIntegral4D(const double &E, const double &s, const double &R0, const double &cosTh0, THnSparseF * hist) {
	double cutoff = GetCutoffS(E);
	double x[4] = {E, s, R0, cosTh0};
	if( s > cutoff ) { x[1] = cutoff; }

	int bin = hist->GetBin(x, false); // false -> don't allocate bin if it's empty
	if( bin != -1) { return (hist->GetBinContent(bin) ); }
	else { return 0; }
}

double WCSimIntegralLookup::GetCutoffS(const double& E) const {
	return ( fIntegrals.GetCutoffS()->GetBinContent(fIntegrals.GetCutoffS()->GetXaxis()->FindBin(E)) );
}


