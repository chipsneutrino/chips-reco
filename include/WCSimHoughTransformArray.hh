#pragma once

#include "TObject.h"

#include <vector>

class WCSimHoughTransform;

class WCSimHoughTransformArray : public TObject
{

public:
	WCSimHoughTransformArray(Int_t bins, Double_t min, Double_t max, Int_t xbins, Int_t ybins);
	~WCSimHoughTransformArray();

	void BuildArray(Int_t bins, Double_t min, Double_t max, Int_t xbins, Int_t ybins);
	void PrintArray();

	Int_t GetBins()
	{
		return fConeAngleBins;
	}
	Int_t GetConeAngleBins()
	{
		return fConeAngleBins;
	}
	Double_t GetConeAngleMin()
	{
		return fConeAngleMin;
	}
	Double_t GetConeAngleMax()
	{
		return fConeAngleMax;
	}
	Double_t GetConeAngle(Int_t bin);
	Double_t GetXBins()
	{
		return fHoughX;
	}
	Double_t GetYBins()
	{
		return fHoughY;
	}

	Double_t GetAngle(Int_t bin);
	Int_t FindBin(Double_t angle);

	WCSimHoughTransform *GetHoughTransform(Int_t nAngle);

	void FindPeak(Int_t &bin, Double_t &height);
	void FindPeak(Double_t &phi, Double_t &costheta, Double_t &angle, Double_t &height);
	void FindPeak(Double_t &hx, Double_t &hy, Double_t &hz, Double_t &angle, Double_t &height);
	void FindPeak(std::vector<Double_t> &hx, std::vector<Double_t> &hy, std::vector<Double_t> &hz,
				  std::vector<Double_t> &angle, std::vector<Double_t> &height);
	void AngleMaxima(std::vector<std::vector<Float_t>> &fHoughArrayMaxAngle);
	void FitTSpectrum2(std::vector<Double_t> &houghDirX, std::vector<Double_t> &houghDirY,
					   std::vector<Double_t> &houghDirZ, std::vector<Double_t> &houghAngle,
					   std::vector<Double_t> &houghHeight);
	void FitTSpectrum3(std::vector<Double_t> &houghDirX, std::vector<Double_t> &houghDirY,
					   std::vector<Double_t> &houghDirZ, std::vector<Double_t> &houghAngle,
					   std::vector<Double_t> &houghHeight);

	// Leigh
	void FitMultiPeaksSmooth(std::vector<Double_t> &houghDirX, std::vector<Double_t> &houghDirY,
							 std::vector<Double_t> &houghDirZ, Double_t seedDirX, Double_t seedDirY, Double_t seedDirZ,
							 std::vector<Double_t> &houghAngle, std::vector<Double_t> &houghHeight);

	void Reset();

private:
	void DeleteArray();

	Int_t fHoughX;			// # bins in Hough X axis
	Int_t fHoughY;			// # bins in Hough y axis
	Int_t fConeAngleBins;	//	# bins in Hough cone angle
	Double_t fConeAngleMin; // min cone angle
	Double_t fConeAngleMax; // max cone angle

	std::vector<WCSimHoughTransform *> fHoughArray; // vector of Hough transforms, stepping through cone angle

	ClassDef(WCSimHoughTransformArray, 0)
};
