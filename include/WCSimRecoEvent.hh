#pragma once

#include "TObject.h"

#include <vector>

class WCSimRecoDigit;
class WCSimRecoVertex;
class WCSimRecoRing;
class WCSimHoughTransformArray;

class WCSimRecoEvent : public TObject
{
public:
	WCSimRecoEvent();
	~WCSimRecoEvent();

	void Reset();

	void SetHeader(Int_t run, Int_t event, Int_t trigger);

	Int_t GetRun()
	{
		return fRunNum;
	}
	Int_t GetEvent()
	{
		return fEventNum;
	}
	Int_t GetTrigger()
	{
		return fTriggerNum;
	}

	void AddDigit(WCSimRecoDigit *digit);
	void AddVetoDigit(WCSimRecoDigit *digit);
	void AddFilterDigit(WCSimRecoDigit *digit);
	void AddRing(WCSimRecoRing *track);

	WCSimRecoDigit *GetDigit(Int_t n);
	Int_t GetNDigits();

	WCSimRecoDigit *GetFilterDigit(Int_t n);
	Int_t GetNFilterDigits();

	WCSimRecoRing *GetRing(Int_t n);
	Int_t GetNRings();

	WCSimRecoDigit *GetVetoDigit(Int_t n);
	Int_t GetNVetoDigits();

	WCSimRecoRing *GetPrimaryRing();

	void SetHoughArray(WCSimHoughTransformArray *array);
	WCSimHoughTransformArray *GetHoughArray();

	void SetVertex(Double_t x, Double_t y, Double_t z, Double_t t);
	void SetDirection(Double_t px, Double_t py, Double_t pz);
	void SetConeAngle(Double_t angle);
	void SetTrackLength(Double_t length);
	void SetVtxFOM(Double_t fom, Int_t nsteps, Bool_t pass);
	void SetVtxStatus(Int_t status);

	WCSimRecoVertex *GetVertex();

	Double_t GetVtxX();
	Double_t GetVtxY();
	Double_t GetVtxZ();
	Double_t GetVtxTime();

	Double_t GetDirX();
	Double_t GetDirY();
	Double_t GetDirZ();

	Double_t GetConeAngle();
	Double_t GetTrackLength();

	Double_t GetVtxFOM();
	Int_t GetVtxIterations();
	Bool_t GetVtxPass();
	Int_t GetVtxStatus();

	Bool_t FoundVertex();
	Bool_t FoundDirection();
	Bool_t FoundRings();

	void SetFilterDone()
	{
		fIsFilterDone = 1;
	}
	Bool_t IsFilterDone()
	{
		return fIsFilterDone;
	}

	void SetVertexFinderDone()
	{
		fIsVertexFinderDone = 1;
	}
	Bool_t IsVertexFinderDone()
	{
		return fIsVertexFinderDone;
	}

	void SetRingFinderDone()
	{
		fIsRingFinderDone = 1;
	}
	Bool_t IsRingFinderDone()
	{
		return fIsRingFinderDone;
	}

	std::vector<WCSimRecoDigit *> *GetDigitList()
	{
		return fDigitList;
	}
	std::vector<WCSimRecoDigit *> *GetVetoDigitList()
	{
		return fVetoDigitList;
	}
	std::vector<WCSimRecoDigit *> *GetFilterDigitList()
	{
		return fFilterDigitList;
	}
	std::vector<WCSimRecoRing *> *GetRingList()
	{
		return fRingList;
	}

	void PrintDigitList(const char *filename);
	void PrintFilterDigitList(const char *filename);

	void PrintEvent();

private:
	void ClearDigits();
	void ClearVetoDigits();
	void ClearFilterDigits();
	void ClearRings();

	Int_t fRunNum;
	Int_t fEventNum;
	Int_t fTriggerNum;

	std::vector<WCSimRecoDigit *> *fDigitList;
	std::vector<WCSimRecoDigit *> *fVetoDigitList;
	std::vector<WCSimRecoDigit *> *fFilterDigitList;
	std::vector<WCSimRecoRing *> *fRingList;
	WCSimHoughTransformArray *fHoughArray;

	WCSimRecoVertex *fVertex;

	Bool_t fIsFilterDone;
	Bool_t fIsVertexFinderDone;
	Bool_t fIsRingFinderDone;

	ClassDef(WCSimRecoEvent, 0)
};
