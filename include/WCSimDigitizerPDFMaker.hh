/*
 * WCSimDigitizerPDFMaker.h
 *
 *  Created on: 5 Jan 2015
 *      Author: andy
 */

#pragma once

#include "WCSimCHIPSPMT.hh"
#include "WCSimSK1pePMT.hh"
#include "WCSimTOTPMT.hh"

class TH1D;
class TH2D;
class TRandom3;
class WCSimSK1pePMT;
class WCSimCHIPSPMT;
class WCSimTOTPMT;

class WCSimDigitizerPDFMaker
{
public:
	WCSimDigitizerPDFMaker();
	virtual ~WCSimDigitizerPDFMaker();

	void SetBinning(int bins = 1000, int min = 0, int max = 10);
	int GetBins() const;
	int GetMin() const;
	int GetMax() const;

	void SetNumThrows(int num);
	int GetNumThrows() const;

	void SetType(int type);
	int GetType() const;

	void Run();

private:
	void MakeHisto();
	void LoopDigitize();
	void Digitize();
	void FillEmptyBins();
	void NormHists();
	void CreateDigiNormHists();
	void CreateLnHists();
	void SaveHistogram();

	int ThrowPoisson();

	TRandom3 *fRandom;

	int fNumZero;
	double fMu;
	int fNumThrows;
	int fType; // 0=sk1pe, 1=pmtSim, 2=tot
	int fChargeMin;
	int fChargeMax;
	int fNumChargeBins;

	TH2D *fProbRawHisto;
	TH2D *fProbPoissonHisto;
	TH1D *fPoisson;

	TH2D *fProbRawHisto_digiNorm;
	TH2D *fProbPoissonHisto_digiNorm;
	TH2D *fProbRawHisto_ln;
	TH2D *fProbPoissonHisto_ln;
	TH2D *fProbRawHisto_digiNorm_ln;
	TH2D *fProbPoissonHisto_digiNorm_ln;
	
	WCSimSK1pePMT *fSK1pePMT;
	WCSimCHIPSPMT *fCHIPSPMT;
	WCSimTOTPMT *fTOTPMT;
};
