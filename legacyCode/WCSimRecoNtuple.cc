#include "WCSimRecoNtuple.hh"

#include "WCSimVertexFinder.hh"

#include "WCSimGeometry.hh"
#include "WCSimRecoVertex.hh"
#include "WCSimRecoRing.hh"
#include "WCSimTrueTrack.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <iostream>
#include <cmath>

ClassImp (WCSimRecoNtuple)

WCSimRecoNtuple::WCSimRecoNtuple() :
		WCSimNtuple(), fRecoFile(0), fRecoTree(0) {

}

WCSimRecoNtuple::~WCSimRecoNtuple() {

}

void WCSimRecoNtuple::Reset() {
	std::cout << " *** WCSimRecoNtuple::Reset() *** " << std::endl;

	fRecoFile = 0;
	fRecoTree = 0;

	return;
}

void WCSimRecoNtuple::Fill(WCSimTrueEvent* trueEvent, WCSimRecoEvent* recoEvent) {
	this->ResetVariables();
	this->WriteVariables(trueEvent, recoEvent);
}

void WCSimRecoNtuple::ResetVariables() {
	// reset variables
	// ===============
	fRunNum = -1;
	fEventNum = -1;
	fTriggerNum = -1;

	fTruePDG = 0;
	fTrueP = 0.0;
	fTrueE = 0.0;
	fTrueG4VtxX = 0.0;
	fTrueG4VtxY = 0.0;
	fTrueG4VtxZ = 0.0;
	fTrueVtxX = 0.0;
	fTrueVtxY = 0.0;
	fTrueVtxZ = 0.0;
	fTrueDirX = 0.0;
	fTrueDirY = 0.0;
	fTrueDirZ = 0.0;
	fTrueVtxTraceFwd = -99999.9;
	fTrueVtxTraceBwd = -99999.9;
	fTrueVtxDistToEdge = -99999.9;

	for (int i = 0; i < 1000; i++) {
		fTrueTrackPDG[i] = 99999;
		fTrueTrackParentPDG[i] = 99999;
		fTrueTrackVtxX[i] = 0.0;
		fTrueTrackVtxY[i] = 0.0;
		fTrueTrackVtxZ[i] = 0.0;
		fTrueTrackDirX[i] = 0.0;
		fTrueTrackDirY[i] = 0.0;
		fTrueTrackDirZ[i] = 0.0;
	}

	fTrueG4EndX = 0.0;
	fTrueG4EndY = 0.0;
	fTrueG4EndZ = 0.0;
	fTrueEndX = 0.0;
	fTrueEndY = 0.0;
	fTrueEndZ = 0.0;
	fTrueEndDistToEdge = -99999.9;
	fTrueLength = 0.0;

	fNDigits = 0;
	fNFilterDigits = 0;
	fNRings = 0;

	fFoundVertex = 0;
	fFoundDirection = 0;
	fFoundRings = 0;

	fRecoVtxX = 0.0;
	fRecoVtxY = 0.0;
	fRecoVtxZ = 0.0;
	fRecoVtxTime = 0.0;
	fRecoDirX = 0.0;
	fRecoDirY = 0.0;
	fRecoDirZ = 0.0;
	fRecoConeAngle = 0.0;
	fRecoTrackLength = 0.0;
	fRecoVtxFOM = 0.0;
	fRecoVtxSteps = 0;
	fRecoVtxPass = 0;
	fRecoVtxStatus = 0;

	fRecoVtxTraceFwd = -99999.9;
	fRecoVtxTraceBwd = -99999.9;
	fRecoVtxDistToEdge = -99999.9;

	fRecoDelta = -99999.9;
	fRecoDeltaTrans = -99999.9;
	fRecoDeltaLong = -99999.9;
	fRecoDeltaAngle = -99999.9;

	fRecoRingVtxX = 0.0;
	fRecoRingVtxY = 0.0;
	fRecoRingVtxZ = 0.0;
	fRecoRingDirX = 0.0;
	fRecoRingDirY = 0.0;
	fRecoRingDirZ = 0.0;
	fRecoRingAngle = 0.0;
	fRecoRingHeight = 0.0;
	fRecoRingHeightFrac = 0.0;

	fRecoRingDeltaAngle = -99999.9;

	// vertex reconstruction
	// =====================
	fSimplePositionX = 0.0;
	fSimplePositionY = 0.0;
	fSimplePositionZ = 0.0;
	fSimplePositionDelta = -99999.9;
	fSimplePositionDeltaTrans = -99999.9;
	fSimplePositionDeltaLong = -99999.9;

	fPointPositionX = 0.0;
	fPointPositionY = 0.0;
	fPointPositionZ = 0.0;
	fPointPositionFOM = 0.0;
	fPointPositionSteps = 0;
	fPointPositionPass = 0;
	fPointPositionDelta = -99999.9;
	fPointPositionDeltaTrans = -99999.9;
	fPointPositionDeltaLong = -99999.9;

	fSimpleDirectionX = 0.0;
	fSimpleDirectionY = 0.0;
	fSimpleDirectionZ = 0.0;
	fSimpleDirectionDeltaAngle = -99999.9;

	fPointDirectionX = 0.0;
	fPointDirectionY = 0.0;
	fPointDirectionZ = 0.0;
	fPointDirectionFOM = 0.0;
	fPointDirectionSteps = 0;
	fPointDirectionPass = 0;
	fPointDirectionDeltaAngle = -99999.9;

	fPointVtxX = 0.0;
	fPointVtxY = 0.0;
	fPointVtxZ = 0.0;
	fPointDirX = 0.0;
	fPointDirY = 0.0;
	fPointDirZ = 0.0;
	fPointConeAngle = 0.0;
	fPointTrackLength = 0.0;
	fPointVtxFOM = 0.0;
	fPointVtxSteps = 0;
	fPointVtxPass = 0;
	fPointVtxDelta = -99999.9;
	fPointVtxDeltaTrans = -99999.9;
	fPointVtxDeltaLong = -99999.9;
	fPointVtxDeltaAngle = -99999.9;

	fExtendedVtxX = 0.0;
	fExtendedVtxY = 0.0;
	fExtendedVtxZ = 0.0;
	fExtendedDirX = 0.0;
	fExtendedDirY = 0.0;
	fExtendedDirZ = 0.0;
	fExtendedConeAngle = 0.0;
	fExtendedTrackLength = 0.0;
	fExtendedVtxFOM = 0.0;
	fExtendedVtxSteps = 0;
	fExtendedVtxPass = 0;
	fExtendedVtxDelta = -99999.9;
	fExtendedVtxDeltaTrans = -99999.9;
	fExtendedVtxDeltaLong = -99999.9;
	fExtendedVtxDeltaAngle = -99999.9;

	// --- debug ---
	fExtendedVtxTimeParam0 = 0.0;
	fExtendedVtxConeParam0 = 0.0;
	fExtendedVtxConeParam1 = 0.0;
	fExtendedVtxConeParam2 = 0.0;

	return;
}

void WCSimRecoNtuple::WriteVariables(WCSimTrueEvent* fTrueEvent, WCSimRecoEvent* fRecoEvent) {
	std::cout << " *** WCSimRecoNtuple::WriteVariables() *** " << std::endl;

	// header information
	// ==================
	fRunNum = fRecoEvent->GetRun();
	fEventNum = fRecoEvent->GetEvent();
	fTriggerNum = fRecoEvent->GetTrigger();

	// truth information
	// =================
	fTruePDG = fTrueEvent->GetPDG();
	fTrueP = fTrueEvent->GetMomentum();
	fTrueE = fTrueEvent->GetEnergy();

	fTrueG4VtxX = fTrueEvent->GetG4VtxX();
	fTrueG4VtxY = fTrueEvent->GetG4VtxY();
	fTrueG4VtxZ = fTrueEvent->GetG4VtxZ();
	fTrueG4EndX = fTrueEvent->GetG4EndX();
	fTrueG4EndY = fTrueEvent->GetG4EndY();
	fTrueG4EndZ = fTrueEvent->GetG4EndZ();

	fTrueVtxX = fTrueEvent->GetVtxX();
	fTrueVtxY = fTrueEvent->GetVtxY();
	fTrueVtxZ = fTrueEvent->GetVtxZ();
	fTrueEndX = fTrueEvent->GetEndX();
	fTrueEndY = fTrueEvent->GetEndY();
	fTrueEndZ = fTrueEvent->GetEndZ();
	fTrueDirX = fTrueEvent->GetDirX();
	fTrueDirY = fTrueEvent->GetDirY();
	fTrueDirZ = fTrueEvent->GetDirZ();

	fTrueLength = fTrueEvent->GetLength();

	if (fTruePDG != 0) {
		fTrueVtxTraceFwd = WCSimGeometry::Instance()->ForwardProjectionToEdge(fTrueVtxX, fTrueVtxY, fTrueVtxZ,
				fTrueDirX, fTrueDirY, fTrueDirZ);
		fTrueVtxTraceBwd = WCSimGeometry::Instance()->BackwardProjectionToEdge(fTrueVtxX, fTrueVtxY, fTrueVtxZ,
				fTrueDirX, fTrueDirY, fTrueDirZ);
		fTrueVtxDistToEdge = WCSimGeometry::Instance()->DistanceToEdge(fTrueVtxX, fTrueVtxY, fTrueVtxZ);

		fTrueEndDistToEdge = WCSimGeometry::Instance()->DistanceToEdge(fTrueEndX, fTrueEndY, fTrueEndZ);
	}

	// Andy P addition: get track truth information too
	int nTracks = fTrueEvent->GetNTracks();
	for (int i = 0; i < nTracks; i++) {
		WCSimTrueTrack * myTrack = fTrueEvent->GetTrack(i);
		fTrueTrackVtxX[i] = myTrack->GetVtxX();
		fTrueTrackVtxY[i] = myTrack->GetVtxY();
		fTrueTrackVtxZ[i] = myTrack->GetVtxZ();
		fTrueTrackDirX[i] = myTrack->GetDirX();
		fTrueTrackDirY[i] = myTrack->GetDirY();
		fTrueTrackDirZ[i] = myTrack->GetDirZ();
		fTrueTrackPDG[i] = myTrack->GetTrackPDG();
		fTrueTrackParentPDG[i] = myTrack->GetParentPDG();
	}

	// reconstructed variables
	// =======================
	fNDigits = fRecoEvent->GetNDigits();
	fNFilterDigits = fRecoEvent->GetNFilterDigits();
	fNRings = fRecoEvent->GetNRings();

	fFoundVertex = fRecoEvent->FoundVertex();
	fRecoVtxX = fRecoEvent->GetVtxX();
	fRecoVtxY = fRecoEvent->GetVtxY();
	fRecoVtxZ = fRecoEvent->GetVtxZ();
	fRecoVtxTime = fRecoEvent->GetVtxTime();

	fFoundDirection = fRecoEvent->FoundDirection();
	fRecoDirX = fRecoEvent->GetDirX();
	fRecoDirY = fRecoEvent->GetDirY();
	fRecoDirZ = fRecoEvent->GetDirZ();

	fRecoConeAngle = fRecoEvent->GetConeAngle();
	fRecoTrackLength = fRecoEvent->GetTrackLength();

	fRecoVtxFOM = fRecoEvent->GetVtxFOM();
	fRecoVtxSteps = fRecoEvent->GetVtxIterations();
	fRecoVtxPass = fRecoEvent->GetVtxPass();
	fRecoVtxStatus = fRecoEvent->GetVtxStatus();

	if (fFoundVertex && fFoundDirection) {
		fRecoVtxTraceFwd = WCSimGeometry::Instance()->ForwardProjectionToEdge(fRecoVtxX, fRecoVtxY, fRecoVtxZ,
				fRecoDirX, fRecoDirY, fRecoDirZ);
		fRecoVtxTraceBwd = WCSimGeometry::Instance()->BackwardProjectionToEdge(fRecoVtxX, fRecoVtxY, fRecoVtxZ,
				fRecoDirX, fRecoDirY, fRecoDirZ);
		fRecoVtxDistToEdge = WCSimGeometry::Instance()->DistanceToEdge(fRecoVtxX, fRecoVtxY, fRecoVtxZ);
	}

	fFoundRings = fRecoEvent->FoundRings();

	if (fFoundRings) {
		WCSimRecoRing* myPrimaryRing = (WCSimRecoRing*) (fRecoEvent->GetPrimaryRing());
		fRecoRingVtxX = myPrimaryRing->GetVtxX();
		fRecoRingVtxY = myPrimaryRing->GetVtxY();
		fRecoRingVtxZ = myPrimaryRing->GetVtxZ();
		fRecoRingDirX = myPrimaryRing->GetDirX();
		fRecoRingDirY = myPrimaryRing->GetDirY();
		fRecoRingDirZ = myPrimaryRing->GetDirZ();
		fRecoRingAngle = myPrimaryRing->GetAngle();
		fRecoRingHeight = myPrimaryRing->GetHeight();
		std::cout << "fRecoRingHeight  = " << fRecoRingHeight << std::endl;
		fRecoRingHeightFrac = 0.0;

		if (fRecoRingHeight > 0.0) {
			fRecoRingHeightFrac = fRecoRingHeight / (double) fNFilterDigits;
			std::cout << " fRecoRingHeightFrac = " << fRecoRingHeightFrac << std::endl;
		}

	}

	// compare with truth
	// ==================
	if (fTruePDG != 0 && fFoundVertex && fFoundDirection) {
		Double_t dx = fRecoVtxX - fTrueVtxX;
		Double_t dy = fRecoVtxY - fTrueVtxY;
		Double_t dz = fRecoVtxZ - fTrueVtxZ;
		Double_t ds = sqrt(dx * dx + dy * dy + dz * dz);
		Double_t px = fTrueDirX;
		Double_t py = fTrueDirY;
		Double_t pz = fTrueDirZ;

		Double_t epsilon = 1.0e-7;

		if (ds > epsilon) {
			px = dx / ds;
			py = dy / ds;
			pz = dz / ds;
		}

		Double_t costh = px * fTrueDirX + py * fTrueDirY + pz * fTrueDirZ;
		Double_t sinth = 0.0;
		if (costh < 1.0 - epsilon)
			sinth = sqrt(1.0 - costh * costh);
		else
			costh = 1.0;

		Double_t costh2 = fTrueDirX * fRecoDirX + fTrueDirY * fRecoDirY + fTrueDirZ * fRecoDirZ;
		Double_t angle = 0.0;
		if (costh2 < 1.0 - epsilon)
			angle = (180.0 / TMath::Pi()) * acos(costh2);
		else
			angle = 0.0;

		fRecoDelta = ds;
		fRecoDeltaTrans = ds * sinth;
		fRecoDeltaLong = ds * costh;
		fRecoDeltaAngle = angle;
	}

	if (fTruePDG != 0 && fFoundRings) {
		Double_t epsilon = 1.0e-7;
		Double_t costh2 = fTrueDirX * fRecoRingDirX + fTrueDirY * fRecoRingDirY + fTrueDirZ * fRecoRingDirZ;
		Double_t angle = 0.0;
		if (costh2 < 1.0 - epsilon)
			angle = (180.0 / TMath::Pi()) * acos(costh2);
		else
			angle = 0.0;

		fRecoRingDeltaAngle = angle;
	}

	// reconstruction chain
	// =====================
	if (fFoundVertex && fFoundDirection) {
		Double_t dx = 0.0;
		Double_t dy = 0.0;
		Double_t dz = 0.0;
		Double_t ds = 0.0;
		Double_t px = 0.0;
		Double_t py = 0.0;
		Double_t pz = 0.0;
		Double_t qx = 0.0;
		Double_t qy = 0.0;
		Double_t qz = 0.0;

		Double_t epsilon = 1.0e-7;

		Double_t costh = 0.0;
		Double_t sinth = 0.0;
		Double_t costh2 = 0.0;
		Double_t angle = 0.0;

		// simple position
		WCSimRecoVertex* simplePos = WCSimVertexFinder::Instance()->GetSimplePosition();
		fSimplePositionX = simplePos->GetX();
		fSimplePositionY = simplePos->GetY();
		fSimplePositionZ = simplePos->GetZ();

		dx = fSimplePositionX - fTrueVtxX;
		dy = fSimplePositionY - fTrueVtxY;
		dz = fSimplePositionZ - fTrueVtxZ;
		ds = sqrt(dx * dx + dy * dy + dz * dz);
		px = fTrueDirX;
		py = fTrueDirY;
		pz = fTrueDirZ;

		if (ds > epsilon) {
			px = dx / ds;
			py = dy / ds;
			pz = dz / ds;
		}

		costh = px * fTrueDirX + py * fTrueDirY + pz * fTrueDirZ;
		sinth = 0.0;
		if (costh < 1.0 - epsilon)
			sinth = sqrt(1.0 - costh * costh);
		else
			costh = 1.0;

		fSimplePositionDelta = ds;
		fSimplePositionDeltaTrans = ds * sinth;
		fSimplePositionDeltaLong = ds * costh;

		// point position
		WCSimRecoVertex* pointPos = WCSimVertexFinder::Instance()->GetPointPosition();
		fPointPositionX = pointPos->GetX();
		fPointPositionY = pointPos->GetY();
		fPointPositionZ = pointPos->GetZ();
		fPointPositionFOM = pointPos->GetFOM();
		fPointPositionSteps = pointPos->GetIterations();
		fPointPositionPass = pointPos->GetPass();

		dx = fPointPositionX - fTrueVtxX;
		dy = fPointPositionY - fTrueVtxY;
		dz = fPointPositionZ - fTrueVtxZ;
		ds = sqrt(dx * dx + dy * dy + dz * dz);
		px = fTrueDirX;
		py = fTrueDirY;
		pz = fTrueDirZ;

		if (ds > epsilon) {
			px = dx / ds;
			py = dy / ds;
			pz = dz / ds;
		}

		costh = px * fTrueDirX + py * fTrueDirY + pz * fTrueDirZ;
		sinth = 0.0;
		if (costh < 1.0 - epsilon)
			sinth = sqrt(1.0 - costh * costh);
		else
			costh = 1.0;

		fPointPositionDelta = ds;
		fPointPositionDeltaTrans = ds * sinth;
		fPointPositionDeltaLong = ds * costh;

		// simple direction
		WCSimRecoVertex* simpleDir = WCSimVertexFinder::Instance()->GetSimpleDirection();
		fSimpleDirectionX = simpleDir->GetDirX();
		fSimpleDirectionY = simpleDir->GetDirY();
		fSimpleDirectionZ = simpleDir->GetDirZ();

		qx = fSimpleDirectionX;
		qy = fSimpleDirectionY;
		qz = fSimpleDirectionZ;

		costh2 = qx * fTrueDirX + qy * fTrueDirY + qz * fTrueDirZ;
		if (costh2 < 1.0 - epsilon)
			angle = (180.0 / TMath::Pi()) * acos(costh2);
		else
			angle = 0.0;

		fSimpleDirectionDeltaAngle = angle;

		// point direction
		WCSimRecoVertex* pointDir = WCSimVertexFinder::Instance()->GetPointDirection();
		fPointDirectionX = pointDir->GetDirX();
		fPointDirectionY = pointDir->GetDirY();
		fPointDirectionZ = pointDir->GetDirZ();
		fPointDirectionFOM = pointDir->GetFOM();
		fPointDirectionSteps = pointDir->GetIterations();
		fPointDirectionPass = pointDir->GetPass();

		qx = fPointDirectionX;
		qy = fPointDirectionY;
		qz = fPointDirectionZ;

		costh2 = qx * fTrueDirX + qy * fTrueDirY + qz * fTrueDirZ;
		if (costh2 < 1.0 - epsilon)
			angle = (180.0 / TMath::Pi()) * acos(costh2);
		else
			angle = 0.0;

		fPointDirectionDeltaAngle = angle;

		// point vertex
		WCSimRecoVertex* pointVtx = WCSimVertexFinder::Instance()->GetPointVertex();
		fPointVtxX = pointVtx->GetX();
		fPointVtxY = pointVtx->GetY();
		fPointVtxZ = pointVtx->GetZ();
		fPointDirX = pointVtx->GetDirX();
		fPointDirY = pointVtx->GetDirY();
		fPointDirZ = pointVtx->GetDirZ();
		fPointConeAngle = pointVtx->GetConeAngle();
		fPointTrackLength = pointVtx->GetTrackLength();
		fPointVtxFOM = pointVtx->GetFOM();
		fPointVtxSteps = pointVtx->GetIterations();
		fPointVtxPass = pointVtx->GetPass();

		dx = fPointVtxX - fTrueVtxX;
		dy = fPointVtxY - fTrueVtxY;
		dz = fPointVtxZ - fTrueVtxZ;
		ds = sqrt(dx * dx + dy * dy + dz * dz);
		qx = fPointDirX;
		qy = fPointDirY;
		qz = fPointDirZ;
		px = fTrueDirX;
		py = fTrueDirY;
		pz = fTrueDirZ;

		if (ds > epsilon) {
			px = dx / ds;
			py = dy / ds;
			pz = dz / ds;
		}

		costh = px * fTrueDirX + py * fTrueDirY + pz * fTrueDirZ;
		sinth = 0.0;
		if (costh < 1.0 - epsilon)
			sinth = sqrt(1.0 - costh * costh);
		else
			costh = 1.0;

		costh2 = qx * fTrueDirX + qy * fTrueDirY + qz * fTrueDirZ;
		if (costh2 < 1.0 - epsilon)
			angle = (180.0 / TMath::Pi()) * acos(costh2);
		else
			angle = 0.0;

		fPointVtxDelta = ds;
		fPointVtxDeltaTrans = ds * sinth;
		fPointVtxDeltaLong = ds * costh;
		fPointVtxDeltaAngle = angle;

		// extended
		WCSimRecoVertex* extendedVtx = WCSimVertexFinder::Instance()->GetExtendedVertex();
		fExtendedVtxX = extendedVtx->GetX();
		fExtendedVtxY = extendedVtx->GetY();
		fExtendedVtxZ = extendedVtx->GetZ();
		fExtendedDirX = extendedVtx->GetDirX();
		fExtendedDirY = extendedVtx->GetDirY();
		fExtendedDirZ = extendedVtx->GetDirZ();
		fExtendedConeAngle = extendedVtx->GetConeAngle();
		fExtendedTrackLength = extendedVtx->GetTrackLength();
		fExtendedVtxFOM = extendedVtx->GetFOM();
		fExtendedVtxSteps = extendedVtx->GetIterations();
		fExtendedVtxPass = extendedVtx->GetPass();

		dx = fExtendedVtxX - fTrueVtxX;
		dy = fExtendedVtxY - fTrueVtxY;
		dz = fExtendedVtxZ - fTrueVtxZ;
		ds = sqrt(dx * dx + dy * dy + dz * dz);
		qx = fExtendedDirX;
		qy = fExtendedDirY;
		qz = fExtendedDirZ;
		px = fTrueDirX;
		py = fTrueDirY;
		pz = fTrueDirZ;

		if (ds > epsilon) {
			px = dx / ds;
			py = dy / ds;
			pz = dz / ds;
		}

		costh = px * fTrueDirX + py * fTrueDirY + pz * fTrueDirZ;
		sinth = 0.0;
		if (costh < 1.0 - epsilon)
			sinth = sqrt(1.0 - costh * costh);
		else
			costh = 1.0;

		costh2 = qx * fTrueDirX + qy * fTrueDirY + qz * fTrueDirZ;
		if (costh2 < 1.0 - epsilon)
			angle = (180.0 / TMath::Pi()) * acos(costh2);
		else
			angle = 0.0;

		fExtendedVtxDelta = ds;
		fExtendedVtxDeltaTrans = ds * sinth;
		fExtendedVtxDeltaLong = ds * costh;
		fExtendedVtxDeltaAngle = angle;

		// --- debug ---
		fExtendedVtxTimeParam0 = WCSimVertexFinder::Instance()->fTimeParam0;
		fExtendedVtxConeParam0 = WCSimVertexFinder::Instance()->fConeParam0;
		fExtendedVtxConeParam1 = WCSimVertexFinder::Instance()->fConeParam1;
		fExtendedVtxConeParam2 = WCSimVertexFinder::Instance()->fConeParam2;
	}

	// write variables
	// ===============
	this->WriteToFile();

	return;
}

void WCSimRecoNtuple::WriteToFile() {
	TDirectory* tmpd = 0;

	if (fRecoFile == 0) {
		tmpd = gDirectory;
		std::cout << " *** WCSimRecoNtuple::OpenFile() *** " << std::endl;
		std::cout << "  opening file: " << fNtpFileName.Data() << std::endl;
		fRecoFile = new TFile(fNtpFileName.Data(), "recreate");
		fRecoTree = new TTree("ntuple", "my analysis ntuple");

//    fRecoTree->Branch("run",&fRunNum,"run/I");
		fRecoTree->Branch("eventID", &fEventNum, "eventID/I");
//    fRecoTree->Branch("trigger",&fTriggerNum,"trigger/I");

//    fRecoTree->Branch("truePDG",&fTruePDG,"truePDG/I");
//    fRecoTree->Branch("trueP",&fTrueP,"trueP/D");
//    fRecoTree->Branch("trueE",&fTrueE,"trueE/D"); 

//		fRecoTree->Branch("trueTrackPDG",&fTrueTrackPDG,"trueTrackPDG[1000]/I");
//		fRecoTree->Branch("trueTrackParentPDG",&fTrueTrackParentPDG,"trueTrackParentPDG[1000]/I");
//		fRecoTree->Branch("trueTrackVtxX",&fTrueTrackVtxX,"trueTrackVtxX[1000]/D");
//		fRecoTree->Branch("trueTrackVtxY",&fTrueTrackVtxY,"trueTrackVtxY[1000]/D");
//		fRecoTree->Branch("trueTrackVtxZ",&fTrueTrackVtxZ,"trueTrackVtxZ[1000]/D");
//		fRecoTree->Branch("trueTrackDirX",&fTrueTrackDirX,"trueTrackDirX[1000]/D");
//		fRecoTree->Branch("trueTrackDirY",&fTrueTrackDirY,"trueTrackDirY[1000]/D");
//		fRecoTree->Branch("trueTrackDirZ",&fTrueTrackDirZ,"trueTrackDirZ[1000]/D");

//    fRecoTree->Branch("trueG4VtxX",&fTrueG4VtxX,"trueG4VtxX/D");
//    fRecoTree->Branch("trueG4VtxY",&fTrueG4VtxY,"trueG4VtxY/D");
//    fRecoTree->Branch("trueG4VtxZ",&fTrueG4VtxZ,"trueG4VtxZ/D");  
//    fRecoTree->Branch("trueG4EndX",&fTrueG4EndX,"trueG4EndX/D");
//    fRecoTree->Branch("trueG4EndY",&fTrueG4EndY,"trueG4EndY/D");
//    fRecoTree->Branch("trueG4EndZ",&fTrueG4EndZ,"trueG4EndZ/D");  
//
//    fRecoTree->Branch("trueVtxX",&fTrueVtxX,"trueVtxX/D");
//    fRecoTree->Branch("trueVtxY",&fTrueVtxY,"trueVtxY/D");
//    fRecoTree->Branch("trueVtxZ",&fTrueVtxZ,"trueVtxZ/D");   
//    fRecoTree->Branch("trueEndX",&fTrueEndX,"trueEndX/D");
//    fRecoTree->Branch("trueEndY",&fTrueEndY,"trueEndY/D");
//    fRecoTree->Branch("trueEndZ",&fTrueEndZ,"trueEndZ/D"); 
//    fRecoTree->Branch("trueDirX",&fTrueDirX,"trueDirX/D");
//    fRecoTree->Branch("trueDirY",&fTrueDirY,"trueDirY/D");
//    fRecoTree->Branch("trueDirZ",&fTrueDirZ,"trueDirZ/D");

//    fRecoTree->Branch("trueLength",&fTrueLength,"trueLength/D");

//    fRecoTree->Branch("trueVtxTraceFwd",&fTrueVtxTraceFwd,"trueVtxTraceFwd/D");
//    fRecoTree->Branch("trueVtxTraceBwd",&fTrueVtxTraceBwd,"trueVtxTraceBwd/D");
//    fRecoTree->Branch("trueVtxDistToEdge",&fTrueVtxDistToEdge,"trueVtxDistToEdge/D");
//    fRecoTree->Branch("trueEndDistToEdge",&fTrueEndDistToEdge,"trueEndDistToEdge/D");
//
//    fRecoTree->Branch("nDigits",&fNDigits,"nDigits/I");
//    fRecoTree->Branch("nFilterDigits",&fNFilterDigits,"nFilterDigits/I");
//    fRecoTree->Branch("nRings",&fNRings,"nRings/I");
//    fRecoTree->Branch("foundVertex",&fFoundVertex,"foundVertex/I");
//    fRecoTree->Branch("foundDirection",&fFoundDirection,"foundDirection/I");
//    fRecoTree->Branch("foundRings",&fFoundRings,"foundRings/I");
//
		fRecoTree->Branch("recoVtxX", &fRecoVtxX, "recoVtxX/D");
		fRecoTree->Branch("recoVtxY", &fRecoVtxY, "recoVtxY/D");
		fRecoTree->Branch("recoVtxZ", &fRecoVtxZ, "recoVtxZ/D");
//    fRecoTree->Branch("recoVtxTime",&fRecoVtxTime,"recoVtxTime/D");
		fRecoTree->Branch("recoDirX", &fRecoDirX, "recoDirX/D");
		fRecoTree->Branch("recoDirY", &fRecoDirY, "recoDirY/D");
		fRecoTree->Branch("recoDirZ", &fRecoDirZ, "recoDirZ/D");
//    fRecoTree->Branch("recoConeAngle",&fRecoConeAngle,"recoConeAngle/D");
//    fRecoTree->Branch("recoTrackLength",&fRecoTrackLength,"recoTrackLength/D");
//    fRecoTree->Branch("recoVtxFOM",&fRecoVtxFOM,"recoVtxFOM/D");
//    fRecoTree->Branch("recoVtxSteps",&fRecoVtxSteps,"recoVtxSteps/I");  
//    fRecoTree->Branch("recoVtxPass",&fRecoVtxPass,"recoVtxPass/I");
//    fRecoTree->Branch("recoVtxStatus",&fRecoVtxStatus,"recoVtxStatus/I");
//    fRecoTree->Branch("recoVtxTraceFwd",&fRecoVtxTraceFwd,"recoVtxTraceFwd/D");
//    fRecoTree->Branch("recoVtxTraceBwd",&fRecoVtxTraceBwd,"recoVtxTraceBwd/D");
//    fRecoTree->Branch("recoVtxDistToEdge",&fRecoVtxDistToEdge,"recoVtxDistToEdge/D");
//
		fRecoTree->Branch("recoRingVtxX", &fRecoRingVtxX, "recoRingVtxX/D");
		fRecoTree->Branch("recoRingVtxY", &fRecoRingVtxY, "recoRingVtxY/D");
		fRecoTree->Branch("recoRingVtxZ", &fRecoRingVtxZ, "recoRingVtxZ/D");
		fRecoTree->Branch("recoRingDirX", &fRecoRingDirX, "recoRingDirX/D");
		fRecoTree->Branch("recoRingDirY", &fRecoRingDirY, "recoRingDirY/D");
		fRecoTree->Branch("recoRingDirZ", &fRecoRingDirZ, "recoRingDirZ/D");
		fRecoTree->Branch("recoRingAngle", &fRecoRingAngle, "recoRingAngle/D");
		fRecoTree->Branch("recoRingHeight", &fRecoRingHeight, "recoRingHeight/D");
		fRecoTree->Branch("recoRingHeightFrac", &fRecoRingHeightFrac, "recoRingHeightFrac/D");

//    fRecoTree->Branch("recoDelta",&fRecoDelta,"recoDelta/D");
//    fRecoTree->Branch("recoDeltaTrans",&fRecoDeltaTrans,"recoDeltaTrans/D");
//    fRecoTree->Branch("recoDeltaLong",&fRecoDeltaLong,"recoDeltaLong/D");
//    fRecoTree->Branch("recoDeltaAngle",&fRecoDeltaAngle,"recoDeltaAngle/D");
//    fRecoTree->Branch("recoRingDeltaAngle",&fRecoRingDeltaAngle,"recoRingDeltaAngle/D");
//
//    fRecoTree->Branch("simplePositionX",&fSimplePositionX,"simplePositionX/D");
//    fRecoTree->Branch("simplePositionY",&fSimplePositionY,"simplePositionY/D");
//    fRecoTree->Branch("simplePositionZ",&fSimplePositionZ,"simplePositionZ/D");
//    fRecoTree->Branch("simplePositionDelta",&fSimplePositionDelta,"simplePositionDelta/D");
//    fRecoTree->Branch("simplePositionDeltaTrans",&fSimplePositionDeltaTrans,"simplePositionDeltaTrans/D");
//    fRecoTree->Branch("simplePositionDeltaLong",&fSimplePositionDeltaLong,"simplePositionDeltaLong/D");
//
//    fRecoTree->Branch("pointPositionX",&fPointPositionX,"pointPositionX/D");
//    fRecoTree->Branch("pointPositionY",&fPointPositionY,"pointPositionY/D");
//    fRecoTree->Branch("pointPositionZ",&fPointPositionZ,"pointPositionZ/D");
//    fRecoTree->Branch("pointPositionFOM",&fPointPositionFOM,"pointPositionFOM/D");  
//    fRecoTree->Branch("pointPositionSteps",&fPointPositionSteps,"pointPositionSteps/I");
//    fRecoTree->Branch("pointPositionPass",&fPointPositionPass,"pointPositionPass/I");
//    fRecoTree->Branch("pointPositionDelta",&fPointPositionDelta,"pointPositionDelta/D");
//    fRecoTree->Branch("pointPositionDeltaTrans",&fPointPositionDeltaTrans,"pointPositionDeltaTrans/D");
//    fRecoTree->Branch("pointPositionDeltaLong",&fPointPositionDeltaLong,"pointPositionDeltaLong/D");
//
//    fRecoTree->Branch("simpleDirectionX",&fSimpleDirectionX,"simpleDirectionX/D");
//    fRecoTree->Branch("simpleDirectionY",&fSimpleDirectionY,"simpleDirectionY/D");
//    fRecoTree->Branch("simpleDirectionZ",&fSimpleDirectionZ,"simpleDirectionZ/D");
//    fRecoTree->Branch("simpleDirectionDeltaAngle",&fSimpleDirectionDeltaAngle,"simpleDirectionDeltaAngle/D");
//
//    fRecoTree->Branch("pointDirectionX",&fPointDirectionX,"pointDirectionX/D");
//    fRecoTree->Branch("pointDirectionY",&fPointDirectionY,"pointDirectionY/D");
//    fRecoTree->Branch("pointDirectionZ",&fPointDirectionZ,"pointDirectionZ/D");
//    fRecoTree->Branch("pointDirectionFOM",&fPointDirectionFOM,"pointDirectionFOM/D");
//    fRecoTree->Branch("pointDirectionSteps",&fPointDirectionSteps,"pointDirSteps/I");
//    fRecoTree->Branch("pointDirectionPass",&fPointDirectionPass,"pointDirPass/I");
//    fRecoTree->Branch("pointDirectionDeltaAngle",&fPointDirectionDeltaAngle,"pointDirectionDeltaAngle/D");
//
//    fRecoTree->Branch("pointVtxX",&fPointVtxX,"pointVtxX/D");
//    fRecoTree->Branch("pointVtxY",&fPointVtxY,"pointVtxY/D");
//    fRecoTree->Branch("pointVtxZ",&fPointVtxZ,"pointVtxZ/D");
//    fRecoTree->Branch("pointDirX",&fPointDirX,"pointDirX/D");
//    fRecoTree->Branch("pointDirY",&fPointDirY,"pointDirY/D");
//    fRecoTree->Branch("pointDirZ",&fPointDirZ,"pointDirZ/D");
//    fRecoTree->Branch("pointConeAngle",&fPointConeAngle,"pointConeAngle/D");
//    fRecoTree->Branch("pointTrackLength",&fPointTrackLength,"pointTrackLength/D");
//    fRecoTree->Branch("pointVtxFOM",&fPointVtxFOM,"pointVtxFOM/D");
//    fRecoTree->Branch("pointVtxSteps",&fPointVtxSteps,"pointVtxSteps/I");
//    fRecoTree->Branch("pointVtxPass",&fPointVtxPass,"pointVtxPass/I");
//    fRecoTree->Branch("pointVtxDelta",&fPointVtxDelta,"pointVtxDelta/D");
//    fRecoTree->Branch("pointVtxDeltaTrans",&fPointVtxDeltaTrans,"pointVtxDeltaTrans/D");
//    fRecoTree->Branch("pointVtxDeltaLong",&fPointVtxDeltaLong,"pointVtxDeltaLong/D");
//    fRecoTree->Branch("pointVtxDeltaAngle",&fPointVtxDeltaAngle,"pointVtxDeltaAngle/D");
//
//    fRecoTree->Branch("extendedVtxX",&fExtendedVtxX,"extendedVtxX/D");
//    fRecoTree->Branch("extendedVtxY",&fExtendedVtxY,"extendedVtxY/D");
//    fRecoTree->Branch("extendedVtxZ",&fExtendedVtxZ,"extendedVtxZ/D");
//    fRecoTree->Branch("extendedDirX",&fExtendedDirX,"extendedDirX/D");
//    fRecoTree->Branch("extendedDirY",&fExtendedDirY,"extendedDirY/D");
//    fRecoTree->Branch("extendedDirZ",&fExtendedDirZ,"extendedDirZ/D"); 
//    fRecoTree->Branch("extendedConeAngle",&fExtendedConeAngle,"extendedConeAngle/D");
//    fRecoTree->Branch("extendedTrackLength",&fExtendedTrackLength,"extendedTrackLength/D");
//    fRecoTree->Branch("extendedVtxFOM",&fExtendedVtxFOM,"extendedVtxFOM/D");
//    fRecoTree->Branch("extendedVtxSteps",&fExtendedVtxSteps,"extendedVtxSteps/I");
//    fRecoTree->Branch("extendedVtxPass",&fExtendedVtxPass,"extendedVtxPass/I");
//    fRecoTree->Branch("extendedVtxDelta",&fExtendedVtxDelta,"extendedVtxDelta/D");
//    fRecoTree->Branch("extendedVtxDeltaTrans",&fExtendedVtxDeltaTrans,"extendedVtxDeltaTrans/D");
//    fRecoTree->Branch("extendedVtxDeltaLong",&fExtendedVtxDeltaLong,"extendedVtxDeltaLong/D");
//    fRecoTree->Branch("extendedVtxDeltaAngle",&fExtendedVtxDeltaAngle,"extendedVtxDeltaAngle/D");
//
//    // --- debug ---
//    fRecoTree->Branch("extendedVtxTimeParam0",&fExtendedVtxTimeParam0,"extendedVtxTimeParam0/D");
//    fRecoTree->Branch("extendedVtxConeParam0",&fExtendedVtxConeParam0,"extendedVtxConeParam0/D");
//    fRecoTree->Branch("extendedVtxConeParam1",&fExtendedVtxConeParam1,"extendedVtxConeParam1/D");
		fRecoTree->Branch("extendedVtxConeParam2", &fExtendedVtxConeParam2, "extendedVtxConeParam2/D"); // muon-ness parameter
//
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

void WCSimRecoNtuple::CloseFile() {
	TDirectory* tmpd = 0;

	if (fRecoFile) {
		tmpd = gDirectory;
		std::cout << " *** WCSimRecoNtuple::CloseFile() *** " << std::endl;
		std::cout << "  closing file: " << fNtpFileName.Data() << std::endl;
		fRecoFile->cd();
		fRecoTree->Write();
		fRecoFile->Close();
		gDirectory = tmpd;
	}

	return;
}
