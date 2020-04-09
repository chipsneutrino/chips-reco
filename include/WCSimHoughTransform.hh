#pragma once

#include "TObject.h"

class TH2D;

class WCSimHoughTransform : public TObject
{

public:
	WCSimHoughTransform(Int_t phiBins, Int_t costhetaBins);
	~WCSimHoughTransform();

	void Fill(Double_t nx, Double_t ny, Double_t nz, Double_t weight);
	void Fill(Double_t phi, Double_t costheta, Double_t weight);
	void Reset();

	Int_t GetEntries()
	{
		return fEntries;
	}
	Double_t GetIntegral()
	{
		return fIntegral;
	}

	Double_t GetEntry(Int_t binX, Int_t binY)
	{
		return fHoughArray[binX * fYbins + binY];
	}
	Double_t GetX(Int_t binx)
	{
		return -180.0 + 360.0 * (binx + 0.5) / fXbins;
	}
	Double_t GetY(Int_t biny)
	{
		return -1.0 + 2.0 * (biny + 0.5) / fYbins;
	}

	Int_t GetPhiBins()
	{
		return fXbins;
	}
	Int_t GetThetaBins()
	{
		return fYbins;
	}

	Int_t GetBinX(Double_t x)
	{
		return (Int_t)(fXbins * (x + 180.0) / 360.0);
	}
	Int_t GetBinY(Double_t y)
	{
		return (Int_t)(fYbins * (y + 1.0) / 2.0);
	}

	void FindPeak(Double_t &height);
	void FindPeak(Double_t &phi, Double_t &costheta, Double_t &height);
	void FindPeak(Double_t &hx, Double_t &hy, Double_t &hz, Double_t &height);

	TH2D *GetTH2D(const char *histname);

	TH2D *GetRotatedTH2D(const char *histname, Double_t deltaPhi = 0.0);

private:
	Int_t fEntries;
	Double_t fIntegral;

	Int_t fXbins;
	Int_t fYbins;

	Double_t *fHoughArray;

	ClassDef(WCSimHoughTransform, 0)
};
