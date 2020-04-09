#pragma once

#include "TObject.h"

class WCSimRecoRing : public TObject
{
public:
	WCSimRecoRing(Double_t vtxx, Double_t vtxy, Double_t vtxz, Double_t dirx, Double_t diry, Double_t dirz,
				  Double_t angle, Double_t height);
	~WCSimRecoRing();

	Double_t GetVtxX()
	{
		return fVtxX;
	}
	Double_t GetVtxY()
	{
		return fVtxY;
	}
	Double_t GetVtxZ()
	{
		return fVtxZ;
	}

	Double_t GetDirX()
	{
		return fDirX;
	}
	Double_t GetDirY()
	{
		return fDirY;
	}
	Double_t GetDirZ()
	{
		return fDirZ;
	}

	Double_t GetAngle()
	{
		return fAngle;
	}
	Double_t GetHeight()
	{
		return fHeight;
	}

private:
	Double_t fVtxX;
	Double_t fVtxY;
	Double_t fVtxZ;

	Double_t fDirX;
	Double_t fDirY;
	Double_t fDirZ;

	Double_t fAngle;
	Double_t fHeight;

	ClassDef(WCSimRecoRing, 0)
};
