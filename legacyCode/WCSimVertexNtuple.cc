#include "WCSimVertexNtuple.hh"
#include "WCSimVertexGeometry.hh"
#include "WCSimVertexFinder.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <iostream>
#include <cmath>

ClassImp (WCSimVertexNtuple)

WCSimVertexNtuple::WCSimVertexNtuple() :
		WCSimNtuple(), fRecoFile(0), fRecoTree(0) {

}

WCSimVertexNtuple::~WCSimVertexNtuple() {

}

void WCSimVertexNtuple::Reset() {
	std::cout << " *** WCSimVertexNtuple::Reset() *** " << std::endl;

	fRecoFile = 0;
	fRecoTree = 0;

	return;
}

void WCSimVertexNtuple::Fill(WCSimTrueEvent* trueEvent, WCSimRecoEvent* recoEvent) {
	this->ResetVariables();
	this->WriteVariables(trueEvent, recoEvent);
}

void WCSimVertexNtuple::ResetVariables() {
	// reset variables
	// ===============
	fRunNum = -1;
	fEventNum = -1;
	fTriggerNum = -1;

	fVtxX = 0.0;
	fVtxY = 0.0;
	fVtxZ = 0.0;
	fVtxTime = 0.0;
	fDirX = 0.0;
	fDirY = 0.0;
	fDirZ = 0.0;

	fTrueVtxX = 0.0;
	fTrueVtxY = 0.0;
	fTrueVtxZ = 0.0;
	fTrueDirX = 0.0;
	fTrueDirY = 0.0;
	fTrueDirZ = 0.0;
	fTrueLength = 0.0;

	index = -1;
	last = -1;

	fIsFiltered = 0;

	fDigitX = 0.0;
	fDigitY = 0.0;
	fDigitZ = 0.0;
	fDigitT = 0.0;
	fDigitQ = 0.0;

	fConeAngle = 0.0;
	fAngle = 0.0;
	fZenith = 0.0;
	fAzimuth = 0.0;
	fSolidAngle = 0.0;

	fDistPoint = 0.0;
	fDistTrack = 0.0;
	fDistPhoton = 0.0;
	fDistScatter = 0.0;

	fDeltaTime = 0.0;
	fDeltaSigma = 0.0;

	fDeltaAngle = 0.0;
	fDeltaPoint = 0.0;
	fDeltaTrack = 0.0;
	fDeltaPhoton = 0.0;
	fDeltaScatter = 0.0;

	fPointPath = 0.0;
	fExtendedPath = 0.0;
	fPointResidual = 0.0;
	fExtendedResidual = 0.0;

	// --- debug ---
	fExtendedResidualCorrected = 0.0;

	return;
}

void WCSimVertexNtuple::WriteVariables(WCSimTrueEvent* fTrueEvent, WCSimRecoEvent* fRecoEvent) {
	std::cout << " *** WCSimVertexNtuple::WriteVariables() *** " << std::endl;

	// run vertex finder
	// =================
	// WCSimVertexFinder::Instance()->Reset();
	// WCSimRecoVertex* fRecoVertex = WCSimVertexFinder::Instance()->Run(fRecoEvent,
	//                                                                  fTrueEvent->GetVtxX(),
	//                                                                  fTrueEvent->GetVtxY(),
	//	                                                              fTrueEvent->GetVtxZ(),
	//	                                                              fTrueEvent->GetDirX(),
	//	                                                              fTrueEvent->GetDirY(),
	//	                                                              fTrueEvent->GetDirZ());
	//

	// get true information
	// ====================

	fTrueVtxX = fTrueEvent->GetVtxX();
	fTrueVtxY = fTrueEvent->GetVtxY();
	fTrueVtxZ = fTrueEvent->GetVtxZ();
	fTrueDirX = fTrueEvent->GetDirX();
	fTrueDirY = fTrueEvent->GetDirY();
	fTrueDirZ = fTrueEvent->GetDirZ();
	fTrueLength = fTrueEvent->GetLength();

	// get reconstructed vertex
	// ========================
	WCSimRecoVertex* fRecoVertex = (WCSimRecoVertex*) (fRecoEvent->GetVertex());

	// check vertex
	// ============
	if (fRecoVertex == 0)
		return;

	// calculate vertex residuals
	// ==========================
	WCSimVertexGeometry::Instance()->CalcResiduals(fRecoEvent, fRecoVertex);

	// header information
	// ==================
	fRunNum = fRecoEvent->GetRun();
	fEventNum = fRecoEvent->GetEvent();
	fTriggerNum = fRecoEvent->GetTrigger();

	// get vertex position
	// ===================
	fVtxX = fRecoVertex->GetX();
	fVtxY = fRecoVertex->GetY();
	fVtxZ = fRecoVertex->GetZ();
	fVtxTime = fRecoVertex->GetTime();

	fDirX = fRecoVertex->GetDirX();
	fDirY = fRecoVertex->GetDirY();
	fDirZ = fRecoVertex->GetDirZ();

	// get digits used in fit
	// ======================
	last = WCSimVertexGeometry::Instance()->GetNDigits();

	for (Int_t idigit = 0; idigit < last; idigit++) {
		index = idigit;

		fIsFiltered = WCSimVertexGeometry::Instance()->IsFiltered(idigit);

		fDigitX = WCSimVertexGeometry::Instance()->GetDigitX(idigit);
		fDigitY = WCSimVertexGeometry::Instance()->GetDigitY(idigit);
		fDigitZ = WCSimVertexGeometry::Instance()->GetDigitZ(idigit);
		fDigitT = WCSimVertexGeometry::Instance()->GetDigitT(idigit);
		fDigitQ = WCSimVertexGeometry::Instance()->GetDigitQ(idigit);

		fConeAngle = WCSimVertexGeometry::Instance()->GetConeAngle(idigit);
		fAngle = WCSimVertexGeometry::Instance()->GetAngle(idigit);
		fZenith = WCSimVertexGeometry::Instance()->GetZenith(idigit);
		fAzimuth = WCSimVertexGeometry::Instance()->GetAzimuth(idigit);
		fSolidAngle = WCSimVertexGeometry::Instance()->GetSolidAngle(idigit);

		fDistPoint = WCSimVertexGeometry::Instance()->GetDistPoint(idigit);
		fDistTrack = WCSimVertexGeometry::Instance()->GetDistTrack(idigit);
		fDistPhoton = WCSimVertexGeometry::Instance()->GetDistPhoton(idigit);
		fDistScatter = WCSimVertexGeometry::Instance()->GetDistScatter(idigit);

		fDeltaTime = WCSimVertexGeometry::Instance()->GetDeltaTime(idigit);
		fDeltaSigma = WCSimVertexGeometry::Instance()->GetDeltaSigma(idigit);

		fDeltaAngle = WCSimVertexGeometry::Instance()->GetDeltaAngle(idigit);
		fDeltaPoint = WCSimVertexGeometry::Instance()->GetDeltaPoint(idigit);
		fDeltaTrack = WCSimVertexGeometry::Instance()->GetDeltaTrack(idigit);
		fDeltaPhoton = WCSimVertexGeometry::Instance()->GetDeltaPhoton(idigit);
		fDeltaScatter = WCSimVertexGeometry::Instance()->GetDeltaScatter(idigit);

		fPointPath = WCSimVertexGeometry::Instance()->GetPointPath(idigit);
		fExtendedPath = WCSimVertexGeometry::Instance()->GetExtendedPath(idigit);
		fPointResidual = WCSimVertexGeometry::Instance()->GetPointResidual(idigit);
		fExtendedResidual = WCSimVertexGeometry::Instance()->GetExtendedResidual(idigit);

		// --- debug ---
		fExtendedResidualCorrected = fExtendedResidual
				+ WCSimVertexGeometry::Instance()->GetDeltaCorrection(idigit, fTrueLength);

		this->WriteToFile();
	}

	return;
}

void WCSimVertexNtuple::WriteToFile() {
	TDirectory* tmpd = 0;

	if (fRecoFile == 0) {
		tmpd = gDirectory;
		std::cout << " *** WCSimVertexNtuple::OpenFile(...) *** " << std::endl;
		std::cout << "  opening file: " << fNtpFileName.Data() << std::endl;
		fRecoFile = new TFile(fNtpFileName.Data(), "recreate");
		fRecoTree = new TTree("ntuple", "my analysis ntuple");
		fRecoTree->Branch("run", &fRunNum, "run/I");
		fRecoTree->Branch("event", &fEventNum, "event/I");
		fRecoTree->Branch("trigger", &fTriggerNum, "trigger/I");
		fRecoTree->Branch("vtxX", &fVtxX, "vtxX/D");
		fRecoTree->Branch("vtxY", &fVtxY, "vtxY/D");
		fRecoTree->Branch("vtxZ", &fVtxZ, "vtxZ/D");
		fRecoTree->Branch("vtxTime", &fVtxTime, "vtxTime/D");
		fRecoTree->Branch("dirX", &fDirX, "dirX/D");
		fRecoTree->Branch("dirY", &fDirY, "dirY/D");
		fRecoTree->Branch("dirZ", &fDirZ, "dirZ/D");
		fRecoTree->Branch("trueVtxX", &fTrueVtxX, "trueVtxX/D");
		fRecoTree->Branch("trueVtxY", &fTrueVtxY, "trueVtxY/D");
		fRecoTree->Branch("trueVtxZ", &fTrueVtxZ, "trueVtxZ/D");
		fRecoTree->Branch("trueDirX", &fTrueDirX, "trueDirX/D");
		fRecoTree->Branch("trueDirY", &fTrueDirY, "trueDirY/D");
		fRecoTree->Branch("trueDirZ", &fTrueDirZ, "trueDirZ/D");
		fRecoTree->Branch("trueLength", &fTrueLength, "trueLength/D");
		fRecoTree->Branch("index", &index, "index/I");
		fRecoTree->Branch("last", &last, "last/I");
		fRecoTree->Branch("filter", &fIsFiltered, "filter/I");
		fRecoTree->Branch("digitX", &fDigitX, "digitX/D");
		fRecoTree->Branch("digitY", &fDigitY, "digitY/D");
		fRecoTree->Branch("digitZ", &fDigitZ, "digitZ/D");
		fRecoTree->Branch("digitT", &fDigitT, "digitT/D");
		fRecoTree->Branch("digitQ", &fDigitQ, "digitQ/D");
		fRecoTree->Branch("coneAngle", &fConeAngle, "coneAngle/D");
		fRecoTree->Branch("angle", &fAngle, "angle/D");
		fRecoTree->Branch("zenith", &fZenith, "zenith/D");
		fRecoTree->Branch("azimuth", &fAzimuth, "azimuth/D");
		fRecoTree->Branch("solidAngle", &fSolidAngle, "solidAngle/D");
		fRecoTree->Branch("distPoint", &fDistPoint, "distPoint/D");
		fRecoTree->Branch("distTrack", &fDistTrack, "distTrack/D");
		fRecoTree->Branch("distPhoton", &fDistPhoton, "distPhoton/D");
		fRecoTree->Branch("distScatter", &fDistScatter, "distScatter/D");
		fRecoTree->Branch("deltaTime", &fDeltaTime, "deltaTime/D");
		fRecoTree->Branch("deltaSigma", &fDeltaSigma, "deltaSigma/D");
		fRecoTree->Branch("deltaAngle", &fDeltaAngle, "deltaAngle/D");
		fRecoTree->Branch("deltaPoint", &fDeltaPoint, "deltaPoint/D");
		fRecoTree->Branch("deltaTrack", &fDeltaTrack, "deltaTrack/D");
		fRecoTree->Branch("deltaPhoton", &fDeltaPhoton, "deltaPhoton/D");
		fRecoTree->Branch("deltaScatter", &fDeltaScatter, "deltaScatter/D");
		fRecoTree->Branch("pointPath", &fPointPath, "pointPath/D");
		fRecoTree->Branch("extendedPath", &fExtendedPath, "extendedPath/D");
		fRecoTree->Branch("pointResidual", &fPointResidual, "pointResidual/D");
		fRecoTree->Branch("extendedResidual", &fExtendedResidual, "extendedResidual/D");

		// --- debug ---
		fRecoTree->Branch("extendedResidualCorrected", &fExtendedResidualCorrected, "extendedResidualCorrected/D");

		gDirectory = tmpd;
	}

	if (fRecoFile) {
		tmpd = gDirectory;
		fRecoFile->cd();
		fRecoTree->Fill();
		gDirectory = tmpd;
	}

	return;
}

void WCSimVertexNtuple::CloseFile() {
	TDirectory* tmpd = 0;

	if (fRecoFile) {
		tmpd = gDirectory;
		std::cout << " *** WCSimVertexNtuple::CloseFile(...) *** " << std::endl;
		std::cout << "  closing file: " << fNtpFileName.Data() << std::endl;
		fRecoFile->cd();
		fRecoTree->Write();
		fRecoFile->Close();
		gDirectory = tmpd;
	}

	return;
}
