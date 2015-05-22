/*
 * WCSimIntegralLookup3D.cc
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#include <stdexcept>
#include "WCSimIntegralLookup3D.hh"
#include <TFile.h>
#include <TH1F.h>
#include <THnSparse.h>

WCSimIntegralLookup3D::WCSimIntegralLookup3D(TString fileName) {

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

	fHistFile->GetObject("fRhoInt",rhoInt);
	fHistFile->GetObject("fRhoSInt",rhoSInt);
	fHistFile->GetObject("fRhoSSInt",rhoSSInt);
	fHistFile->GetObject("fRhoGInt",rhoGInt);
	fHistFile->GetObject("fRhoGSInt",rhoGSInt);
	fHistFile->GetObject("fRhoGSSInt",rhoGSSInt);

	if( rhoInt == 0x0 ) { throw std::runtime_error("Could not get rhoInt from file"); }
	if( rhoSInt == 0x0 ) { throw std::runtime_error("Could not get rhoSInt from file"); }
	if( rhoSSInt == 0x0 ) { throw std::runtime_error("Could not get rhoSSInt from file"); }
	if( rhoGInt == 0x0 ) { throw std::runtime_error("Could not get rhoGInt from file"); }
	if( rhoGSInt == 0x0 ) { throw std::runtime_error("Could not get rhoGSInt from file"); }
	if( rhoGSSInt == 0x0 ) { throw std::runtime_error("Could not get rhoGSSInt from file"); }

	fIntegrals.SetHists(rhoInt, rhoSInt, rhoSSInt, rhoGInt, rhoGSInt, rhoGSSInt, 0);
}

WCSimIntegralLookup3D::~WCSimIntegralLookup3D() {
}

double WCSimIntegralLookup3D::GetRhoIntegral(const double& E, const double &s) {
	return GetIntegral1D(E, fIntegrals.GetRhoInt());
}

double WCSimIntegralLookup3D::GetRhoSIntegral(const double& E, const double &s) {
	return GetIntegral1D(E, fIntegrals.GetRhoSInt());
}

double WCSimIntegralLookup3D::GetRhoSSIntegral(const double& E, const double &s) {
	return GetIntegral1D(E, fIntegrals.GetRhoSSInt());
}

double WCSimIntegralLookup3D::GetRhoGIntegral(const double& E, const double &s, const double& R0, const double& cosTh0) {
	return GetIntegral3D(E, R0, cosTh0, fIntegrals.GetRhoGInt());
}

double WCSimIntegralLookup3D::GetRhoGSIntegral(const double& E, const double &s, const double& R0, const double& cosTh0) {
	return GetIntegral3D(E, R0, cosTh0, fIntegrals.GetRhoGSInt());
}

double WCSimIntegralLookup3D::GetRhoGSSIntegral(const double& E, const double &s, const double& R0, const double& cosTh0) {
	return GetIntegral3D(E, R0, cosTh0, fIntegrals.GetRhoGSSInt());
}

THnSparseF* WCSimIntegralLookup3D::GetRhoIntegralHist() {
	return fIntegrals.GetRhoInt();
}

double WCSimIntegralLookup3D::GetIntegral1D(const double& E, THnSparseF* hist) {
	double x[1] = {E};
	int bin = hist->GetBin(x, false); // Don't allocate it if it's empty
	if( bin != -1){ return ( hist->GetBinContent(bin) ); }
	else { return 0; }
}

double WCSimIntegralLookup3D::GetIntegral3D(const double& E, const double& R0, const double& cosTh0, THnSparseF* hist) {
	double x[3] = {E, R0, cosTh0};
	int bin = hist->GetBin(x, false); // false -> don't allocate bin if it's empty

//	if(0.7 < cosTh0 && cosTh0 < 0.8){
//		std::cout << "Getting integral for:  E = " << E << "  R0 = " << R0 << "   cosTh0 = " << cosTh0 << "  " << std::endl;
//		std::cout << " bin = " << bin;
//		if(bin != -1) {std::cout << "  int = " << hist->GetBinContent(bin) << std::endl;}
//	}
	if( bin != -1) { return (hist->GetBinContent(bin) ); }
	else { return 0; }
}
