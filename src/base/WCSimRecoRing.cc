#include "WCSimRecoRing.hh"

#include "WCSimRecoObjectTable.hh"

ClassImp (WCSimRecoRing)

WCSimRecoRing::WCSimRecoRing(Double_t vtxx, Double_t vtxy, Double_t vtxz, Double_t dirx, Double_t diry, Double_t dirz,
		Double_t angle, Double_t height) {
	fVtxX = vtxx;
	fVtxY = vtxy;
	fVtxZ = vtxz;
	fDirX = dirx;
	fDirY = diry;
	fDirZ = dirz;

	fAngle = angle;
	fHeight = height;

	WCSimRecoObjectTable::Instance()->NewRing();
}

WCSimRecoRing::~WCSimRecoRing() {
	WCSimRecoObjectTable::Instance()->DeleteRing();
}

