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

	void SetBinning(int bins = 1000, double min = 0, double max = 10);
	int GetBins() const;
	double GetMin() const;
	double GetMax() const;

	void SetMu(double mu);
	double GetMu() const;

	void SetNumThrows(int num);
	int GetNumThrows() const;

	void SetType(int type);
	int GetType() const;

	void SetTotPoisson(bool set);
	bool GetTotPoisson() const;

	void Run();

private:
	void MakeHisto();
	void LoopDigitize();
	void Digitize();
	void FillEmptyBins();
	void NormHistogram();
	void lnHist();
	void SaveHistogram();

	int ThrowPoisson();

	TRandom3 *fRandom;

	int fNumZero;
	double fMu;
	int fNumThrows;
	int fType; // 0=sk1pe, 1=pmtSim, 2=tot

	TH2D *fProbHisto;
	double fChargeMin;
	double fChargeMax;
	int fNumChargeBins;
	bool totPoisson;

	TH1D *fDebug;

	WCSimSK1pePMT *fSK1pePMT;
	WCSimCHIPSPMT *fCHIPSPMT;
	WCSimTOTPMT *fTOTPMT;
};
