#include "WCSimTrueEvent.hh"
#include "WCSimTrueTrack.hh"

#include <iostream>
#include <cmath>

ClassImp (WCSimTrueEvent)

WCSimTrueEvent::WCSimTrueEvent() {
	fTrackList = new std::vector<WCSimTrueTrack*>;

	this->Reset();
}

WCSimTrueEvent::~WCSimTrueEvent() {
	delete fTrackList;
}

void WCSimTrueEvent::SetHeader(Int_t ipdg, Double_t g4vx, Double_t g4vy, Double_t g4vz, Double_t g4ex, Double_t g4ey,
		Double_t g4ez, Double_t vx, Double_t vy, Double_t vz, Double_t ex, Double_t ey, Double_t ez, Double_t px,
		Double_t py, Double_t pz, Double_t trkE, Double_t trkP) {
	fIpdg = ipdg;

	fG4VtxX = g4vx;
	fG4VtxY = g4vy;
	fG4VtxZ = g4vz;

	fG4EndX = g4ex;
	fG4EndY = g4ey;
	fG4EndZ = g4ez;

	fVtxX = vx;
	fVtxY = vy;
	fVtxZ = vz;

	fEndX = ex;
	fEndY = ey;
	fEndZ = ez;

	fDirX = px;
	fDirY = py;
	fDirZ = pz;

	fTrkE = trkE;
	fTrkP = trkP;

	fLength = sqrt(
			(fEndX - fVtxX) * (fEndX - fVtxX) + (fEndY - fVtxY) * (fEndY - fVtxY) + (fEndZ - fVtxZ) * (fEndZ - fVtxZ));
}

void WCSimTrueEvent::Reset() {
	fIpdg = 0;

	fG4VtxX = -99999.9;
	fG4VtxY = -99999.9;
	fG4VtxZ = -99999.9;

	fG4EndX = -99999.9;
	fG4EndY = -99999.9;
	fG4EndZ = -99999.9;

	fVtxX = -99999.9;
	fVtxY = -99999.9;
	fVtxZ = -99999.9;

	fEndX = -99999.9;
	fEndY = -99999.9;
	fEndZ = -99999.9;

	fDirX = 0.0;
	fDirY = 0.0;
	fDirZ = 0.0;

	fTrkE = 0.0;
	fTrkP = 0.0;

	fLength = 0.0;

	this->ClearTracks();

	return;
}

Int_t WCSimTrueEvent::GetNTracks() {
	return fTrackList->size();
}

WCSimTrueTrack* WCSimTrueEvent::GetTrack(Int_t itrack) {
	return (WCSimTrueTrack*) (fTrackList->at(itrack));
}

void WCSimTrueEvent::AddTrack(WCSimTrueTrack* myTrack) {
	fTrackList->push_back(myTrack);
}

void WCSimTrueEvent::ClearTracks() {
	fTrackList->clear();
}

void WCSimTrueEvent::PrintEvent() {
	std::cout << " *** WCSimTrueEvent::PrintEvent() *** " << std::endl;

//	// Andy P add: print track info too if it's an e-, e+, mu- or mu+
//	for(int i = 0; i < this->GetNTracks(); i++)
//	{
//		WCSimTrueTrack * track = this->GetTrack(i);
//		if(abs(track->GetTrackPDG()) == 11 || abs(track->GetTrackPDG()) == 13)
//		{
//			track->PrintTrack();
//		}
//	}

	std::cout << " * Ipdg = " << fIpdg << std::endl << " * G4VtxX = " << fG4VtxX << std::endl << " * G4VtxY = "
			<< fG4VtxY << std::endl << " * G4VtxZ = " << fG4VtxZ << std::endl << " * G4EndX = " << fG4EndX << std::endl
			<< " * G4EndY = " << fG4EndY << std::endl << " * G4EndZ = " << fG4EndZ << std::endl << " * VtxX = " << fVtxX
			<< std::endl << " * VtxY = " << fVtxY << std::endl << " * VtxZ = " << fVtxZ << std::endl << " * EndX = "
			<< fEndX << std::endl << " * EndY = " << fEndY << std::endl << " * EndZ = " << fEndZ << std::endl
			<< " * DirX = " << fDirX << std::endl << " * DirY = " << fDirY << std::endl << " * DirZ = " << fDirZ
			<< std::endl << " * TrkE = " << fTrkP << std::endl << " * TrkP = " << fTrkE << std::endl << " * Length = "
			<< fLength << std::endl << " *********************************** " << std::endl;

	return;
}

