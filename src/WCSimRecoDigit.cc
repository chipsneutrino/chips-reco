#include "WCSimRecoDigit.hh"

#include "WCSimRecoObjectTable.hh"

ClassImp (WCSimRecoDigit)

WCSimRecoDigit::WCSimRecoDigit(Int_t region, Int_t tubeID, Double_t x, Double_t y, Double_t z, Double_t rawT,
		Double_t rawQ, Double_t calT, Double_t calQ) {
	fRegion = region;
	fTubeID = tubeID;

	fX = x;
	fY = y;
	fZ = z;

	fRawTime = rawT;
	fRawQPEs = rawQ;

	fCalTime = calT;
	fCalQPEs = calQ;

	fIsFiltered = 0;

	WCSimRecoObjectTable::Instance()->NewDigit();
}

WCSimRecoDigit::~WCSimRecoDigit() {
	WCSimRecoObjectTable::Instance()->DeleteDigit();
}
