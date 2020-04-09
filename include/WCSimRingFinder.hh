#ifndef WCSIMRINGFINDER_HH
#define WCSIMRINGFINDER_HH

#include "TObject.h"

#include <vector>

class WCSimRecoDigit;
class WCSimRecoEvent;
class WCSimRecoVertex;
class WCSimRecoRing;
class WCSimHoughTransform;
class WCSimHoughTransformArray;

class WCSimRingFinder : public TObject
{
public:
	static WCSimRingFinder *Instance();

	static void UseRecoVertex();

	static void HoughX(Int_t x);
	static void HoughY(Int_t y);
	static void HoughPoints(Int_t n);
	static void ConeAngleBins(Int_t bins);
	static void ConeAngleMin(Double_t min);
	static void ConeAngleMax(Double_t max);

	static void PrintParameters();
	void RunPrintParameters();

	void SetHoughX(Int_t x)
	{
		fHoughX = x;
	}
	void SetHoughY(Int_t y)
	{
		fHoughY = y;
	}
	void SetHoughPoints(Int_t n)
	{
		fHoughPoints = n;
	}
	void SetConeAngleBins(Int_t bins)
	{
		fConeAngleBins = bins;
	}
	void SetConeAngleMin(Double_t min)
	{
		fConeAngleMin = min;
	}
	void SetConeAngleMax(Double_t max)
	{
		fConeAngleMax = max;
	}

	void SetUsingRecoVertex()
	{
		fUseVertex = 1;
	}
	Bool_t UsingRecoVertex()
	{
		return fUseVertex;
	}
	void SetUsingTSpectrum2()
	{
		fUseTSpectrum2 = 1;
	}
	Bool_t UsingTSpectrum2()
	{
		return fUseTSpectrum2;
	}

	std::vector<WCSimRecoRing *> *Run(WCSimRecoVertex *vtx);
	std::vector<WCSimRecoRing *> *Run(WCSimRecoEvent *evt, WCSimRecoVertex *vtx);
	std::vector<WCSimRecoRing *> *Run(std::vector<WCSimRecoDigit *> *digitlist, WCSimRecoVertex *vtx);

	WCSimHoughTransform *HoughTransform(WCSimRecoEvent *evt, WCSimRecoVertex *vtx, Double_t angle);
	WCSimHoughTransform *HoughTransform(std::vector<WCSimRecoDigit *> *digitlist, WCSimRecoVertex *vtx,
										Double_t angle);

	WCSimHoughTransformArray *HoughTransformArray(WCSimRecoEvent *evt, WCSimRecoVertex *vtx);
	WCSimHoughTransformArray *HoughTransformArray(std::vector<WCSimRecoDigit *> *digitlist, WCSimRecoVertex *vtx);

	WCSimHoughTransformArray *GetHoughTransformArray()
	{
		return fHoughTransformArray;
	}

private:
	WCSimRingFinder();
	~WCSimRingFinder();

	void Reset();

	// use reco vertex
	Bool_t fUseVertex;
	Bool_t fUseTSpectrum2;

	// hough transform parameters
	Int_t fHoughX;
	Int_t fHoughY;
	Int_t fHoughPoints;
	Int_t fConeAngleBins;
	Double_t fConeAngleMin;
	Double_t fConeAngleMax;

	// hough transform object
	WCSimHoughTransform *fHoughTransform;

	// hough transform array
	WCSimHoughTransformArray *fHoughTransformArray;

	ClassDef(WCSimRingFinder, 0)
};

#endif
