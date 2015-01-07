/*
 * WCSimDigitizerPDFMaker.h
 *
 *  Created on: 5 Jan 2015
 *      Author: andy
 */

#ifndef WCSIMDIGITIZERPDFMAKER_H_
#define WCSIMDIGITIZERPDFMAKER_H_

class TH1D;
class TH2D;
class TRandom3;

class WCSimDigitizerPDFMaker {
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

	void Run();


private:
	void MakeHisto();
	void LoopDigitize();
	void Digitize();
	void NormHistogram();
	void SaveHistogram();

	void BuildQPE0();
	int ThrowPoisson();
	double Rn1pe();
	void Threshold(double& pe,int& iflag);

	double fEfficiency;
	double fMu;
	int fNumThrows;

	TH2D * fProbHisto;
	double fChargeMin;
	double fChargeMax;
	int fNumChargeBins;

	TRandom3 * fRandom;
	double fQPE0[501];

  TH1D * fDebug;


};

#endif /* WCSIMDIGITIZERPDFMAKER_H_ */
