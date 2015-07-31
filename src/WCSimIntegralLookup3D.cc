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
#include <TH3F.h>

//////////////////////////////////////////////////////////////////////////////////////
// First the class used to hold the histograms (compare to WCSimIntegralLookupHists3D //
// in WCSimIntegralLookup.hh/cc														//
//////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// Class to hold all the histograms so that the reader and writer code don't get out of sync
///////////////////////////////////////////////////////////////////////////////////////////////
WCSimIntegralLookupHists3D::WCSimIntegralLookupHists3D() {
	fRhoInt = 0x0;
	fRhoSInt = 0x0;
	fRhoSSInt = 0x0;

	fRhoGInt = 0x0;
	fRhoGSInt = 0x0;
	fRhoGSSInt = 0x0;

	fCutoffS = 0x0;
}

WCSimIntegralLookupHists3D::WCSimIntegralLookupHists3D(TH1F* rhoInt, TH1F* rhoSInt, TH1F* rhoSSInt,
		TH3F* rhoGInt, TH3F* rhoGSInt, TH3F* rhoGSSInt, TH1F* cutoffS) {
	SetHists(rhoInt, rhoSInt, rhoSSInt, rhoGInt, rhoGSInt, rhoGSSInt, cutoffS);

}

TH1F* WCSimIntegralLookupHists3D::GetRhoInt() const {
	return fRhoInt;
}

TH1F* WCSimIntegralLookupHists3D::GetRhoSInt() const {
	return fRhoSInt;
}

TH1F* WCSimIntegralLookupHists3D::GetRhoSSInt() const {
	return fRhoSSInt;
}

TH3F* WCSimIntegralLookupHists3D::GetRhoGInt() const {
	return fRhoGInt;
}

TH3F* WCSimIntegralLookupHists3D::GetRhoGSInt() const {
  // std::cout << fRhoGSInt << std::endl;
  // std::cout << fRhoGSInt->GetEntries() << std::endl;
  // std::cout << "Returning fRhoGSInt" << std::endl;
	return fRhoGSInt;
}

TH3F* WCSimIntegralLookupHists3D::GetRhoGSSInt() const {
	return fRhoGSSInt;
}

TH1F* WCSimIntegralLookupHists3D::GetCutoffS() const {
	return fCutoffS;
}

void WCSimIntegralLookupHists3D::SetHists(TH1F* rhoInt, TH1F* rhoSInt, TH1F* rhoSSInt,
		TH3F* rhoGInt, TH3F* rhoGSInt, TH3F* rhoGSSInt, TH1F* cutoffS) {

	fRhoInt = rhoInt;
	fRhoSInt = rhoSInt;
	fRhoSSInt = rhoSSInt;

	fRhoGInt = rhoGInt;
	fRhoGSInt = rhoGSInt;
	fRhoGSSInt = rhoGSSInt;

	fCutoffS = cutoffS;
  // std::cout << "Finished setting hists" << std::endl;
}




//////////////////////////////////////////////////////////////////////////////////////
// Now the main class to do the looking-up                                          //
//////////////////////////////////////////////////////////////////////////////////////


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

	TH1F * rhoInt = 0x0;
	TH1F * rhoSInt = 0x0;
	TH1F * rhoSSInt = 0x0;
	TH3F * rhoGInt = 0x0;
	TH3F * rhoGSInt = 0x0;
	TH3F * rhoGSSInt = 0x0;

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

  std::cout << "Setting hists" << std::endl;
	fIntegrals.SetHists(rhoInt, rhoSInt, rhoSSInt, rhoGInt, rhoGSInt, rhoGSSInt, 0);
  std::cout << "Set hists" << std::endl;
}

WCSimIntegralLookup3D::~WCSimIntegralLookup3D() {
  fHistFile->Close();
  delete fHistFile;
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

TH1F* WCSimIntegralLookup3D::GetRhoIntegralHist() {
	return fIntegrals.GetRhoInt();
}

double WCSimIntegralLookup3D::GetIntegral1D(const double& E, TH1F* hist) {
  // std::cout << "GetIntegral1D" << std::endl;
  return hist->Interpolate(E);
}

double WCSimIntegralLookup3D::GetIntegral3D(const double& E, const double& R0, const double& cosTh0, TH3F* hist) {
  // std::cout << "GetIntegral3D" << std::endl;
  // std::cout << "hist = " << hist << "  " << hist->GetEntries() << std::endl;
	// std::cout << "Getting integral for:  E = " << E << "  R0 = " << R0 << "   cosTh0 = " << cosTh0 << "  " << std::endl;
  // std::cout << "It has " << hist->GetNdimensions() << " dimensions" << std::endl;
  // std::cout << "Axis 0 goes from " << hist->GetXaxis()->GetXmin() << " to " << hist->GetXaxis()->GetXmax() << std::endl;
  // std::cout << "Axis 1 goes from " << hist->GetYaxis()->GetXmin() << " to " << hist->GetYaxis()->GetXmax() << std::endl;
  // std::cout << "Axis 2 goes from " << hist->GetZaxis()->GetXmin() << " to " << hist->GetZaxis()->GetXmax() << std::endl;
  // std::cout << "Hist bin = " << hist->GetXaxis()->FindBin(E) << ", " << hist->GetYaxis()->FindBin(R0) << ", " << hist->GetZaxis()->FindBin(cosTh0) << std::endl;

  // Sometimes cosTheta is too close to +/- 1 that it's before the centre of the
  // first bin or beyond the centre of the last bin, so TH3::Interpolate doesn't work
  int binE = hist->GetXaxis()->FindBin(E);
  float lowE = hist->GetXaxis()->GetBinLowEdge(binE);
  float hiE = hist->GetXaxis()->GetBinUpEdge(binE);

  int lowBin = hist->FindBin(lowE, R0, cosTh0);
  int hiBin = hist->FindBin(hiE, R0, cosTh0);

  float interp = (hist->GetBinContent(lowBin) * (hiE-E) + hist->GetBinContent(hiBin) * (E-lowE)) / (hiE - lowE);
  return interp;

//	if(0.7 < cosTh0 && cosTh0 < 0.8){
	//  std::cout << "Getting integral for:  E = " << E << "  R0 = " << R0 << "   cosTh0 = " << cosTh0 << "  " << std::endl;
	//  std::cout << " bin = " << bin << std::endl;
	//  if(bin != -1) {std::cout << "  int = " << hist->GetBinContent(bin) << std::endl;}
  //  std::cout << "Done" << std::endl;
	//if( bin != -1) { return (hist->GetBinContent(bin) ); }
	//else { return 0; }
}

