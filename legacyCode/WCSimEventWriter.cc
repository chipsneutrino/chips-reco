#include "WCSimEventWriter.hh"
#include "WCSimGeometry.hh"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TClonesArray.h"

#include <iostream>
#include <cmath>

ClassImp (WCSimEventWriter)

WCSimEventWriter::WCSimEventWriter() :
		fScale(1.0), fFileDigits(0), fTreeDigits(0), fFileTracks(0), fTreeTracks(0), fWriteOutGeometry(0), fWriteOutTracks(
				0), fWriteOutDigits(0), fChain(0), fChainGeom(0), fEvent(0), fGeometry(0) {
	this->Initialize();
}

WCSimEventWriter::~WCSimEventWriter() {

}

void WCSimEventWriter::Initialize() {
	// create event chain
	if (!fChain) {
		fChain = new TChain("wcsimT", "chain");
		fChain->SetBranchAddress("wcsimrootevent", &fEvent);
	}

	// create geometry chain
	if (!fChainGeom) {
		fChainGeom = new TChain("wcsimGeoT", "chainGeom");
		fChainGeom->SetBranchAddress("wcsimrootgeom", &fGeometry);
	}
}

void WCSimEventWriter::AddFile(const char* file) {
	std::cout << " adding: " << file << std::endl;

	fChain->Add(file);
	fChainGeom->Add(file);

	std::cout << "   ... total entries=" << fChain->GetEntries() << std::endl;

	// bug fix
	// =======
	// if( fChain->GetEntries()>0 ){
	//   std::cout << " applying fix: GetBranch(\"wcsimrootevent\")->SetAutoDelete(true) " << std::endl;
	//   fChain->GetBranch("wcsimrootevent")->SetAutoDelete(true);
	// }

	return;
}

void WCSimEventWriter::Run(Int_t entries) {
	std::cout << " *** WCSimEventWriter::Run() *** " << std::endl;

	// check for events
	if (fChain->GetEntries() <= 0) {
		return;
	}

	// check for geometry
	if (fChainGeom->GetEntries() <= 0) {
		return;
	}

	// Load Geometry
	// =============
	// get first entry of geometry tree
	fChainGeom->GetEntry(0);

	// build geometry
	WCSimGeometry::BuildGeometry (fGeometry);

	if (fWriteOutGeometry) {
		WCSimGeometry::WriteGeometry();
	}

	// convert [cm] to [m]
	fScale = 0.01;

	// Loop over entries
	// =================
	Int_t nLast = entries;
	if (nLast > fChain->GetEntries())
		nLast = fChain->GetEntries();
	if (nLast < 0)
		nLast = fChain->GetEntries();

	std::cout << "  looping over " << nLast << " entries " << std::endl;

	// loop over entries
	for (Int_t nEvent = 0; nEvent < nLast; nEvent++) {
		std::cout << " " << nEvent << "/" << fChain->GetEntries() << " " << std::endl;

		if (fEvent)
			delete fEvent;
		fEvent = 0;
		fChain->GetEntry(nEvent);

		for (Int_t nTrigger = 0; nTrigger < fEvent->GetNumberOfEvents(); nTrigger++) {
			WCSimRootTrigger* myTrigger = (WCSimRootTrigger*) (fEvent->GetTrigger(nTrigger));

			// Event Header
			// ============
			WCSimRootEventHeader* myHeader = (WCSimRootEventHeader*) (myTrigger->GetHeader());
			fRun = myHeader->GetRun();
			fEventNumber = myHeader->GetEvtNum();
			fSubEventNumber = myHeader->GetSubEvtNumber();

			// Get Array of True Tracks
			// ========================
			TClonesArray* fTrackArray = (TClonesArray*) (myTrigger->GetTracks());

			for (Int_t nTrack = 0; nTrack < 1 + fTrackArray->GetLast(); nTrack++) {
				WCSimRootTrack* myTrack = (WCSimRootTrack*) (fTrackArray->At(nTrack));

				fTrackNumber = nTrack;
				fPdg = myTrack->GetIpnu();
				fFlag = myTrack->GetFlag();
				fM = myTrack->GetM();
				fP = myTrack->GetP();
				fE = myTrack->GetE();
				fTrkDirX = myTrack->GetDir(0);
				fTrkDirY = myTrack->GetDir(1);
				fTrkDirZ = myTrack->GetDir(2);
				fVtxX = fScale * myTrack->GetStart(0);
				fVtxY = fScale * myTrack->GetStart(1);
				fVtxZ = fScale * myTrack->GetStart(2);
				fVtxTime = myTrack->GetTime();
				fEndX = fScale * myTrack->GetStop(0);
				fEndY = fScale * myTrack->GetStop(1);
				fEndZ = fScale * myTrack->GetStop(2);
				fPdgParent = myTrack->GetParenttype();

				if (fWriteOutTracks)
					this->WriteThisTrack();
			}

			// Get Array of Digitized Hits
			// ===========================
			TClonesArray* fDigiHitArray = (TClonesArray*) (myTrigger->GetCherenkovDigiHits());

			for (Int_t nDigit = 0; nDigit < 1 + fDigiHitArray->GetLast(); nDigit++) {
				WCSimRootCherenkovDigiHit* myDigit = (WCSimRootCherenkovDigiHit*) (fDigiHitArray->At(nDigit));

				fTube = myDigit->GetTubeId();
				fQPEs = myDigit->GetQ();
				fTime = myDigit->GetT();

				fRegion = WCSimGeometry::Instance()->GetRegion(fTube);
				fXpos = fScale * WCSimGeometry::Instance()->GetX(fTube);
				fYpos = fScale * WCSimGeometry::Instance()->GetY(fTube);
				fZpos = fScale * WCSimGeometry::Instance()->GetZ(fTube);

				if (fWriteOutDigits)
					this->WriteThisDigit();
			}

		} // nTrigger
	} // nEvent

	  // Write out data
	  // ==============
	this->CloseFiles();

	return;
}

void WCSimEventWriter::WriteThisDigit() {
	if (!fWriteOutDigits)
		return;

	TDirectory* tmpd = 0;

	if (!fFileDigits) {
		tmpd = gDirectory;
		fFileDigits = new TFile("wcsim.digits.root", "recreate");
		fTreeDigits = new TTree("ntuple", "list of wc digits");
		fTreeDigits->Branch("run", &fRun, "run/I");
		fTreeDigits->Branch("event", &fEventNumber, "event/I");
		fTreeDigits->Branch("subevent", &fSubEventNumber, "subevent/I");
		fTreeDigits->Branch("tube", &fTube, "tube/I");
		fTreeDigits->Branch("x", &fXpos, "x/D");
		fTreeDigits->Branch("y", &fYpos, "y/D");
		fTreeDigits->Branch("z", &fZpos, "z/D");
		fTreeDigits->Branch("q", &fQPEs, "q/D");
		fTreeDigits->Branch("t", &fTime, "t/D");
		fTreeDigits->SetAutoSave(100);
		gDirectory = tmpd;
	}

	tmpd = gDirectory;
	fFileDigits->cd();
	fTreeDigits->Fill();
	gDirectory = tmpd;
}

void WCSimEventWriter::WriteThisTrack() {
	if (!fWriteOutTracks)
		return;

	TDirectory* tmpd = 0;

	if (!fFileTracks) {
		tmpd = gDirectory;
		fFileTracks = new TFile("wcsim.tracks.root", "recreate");
		fTreeTracks = new TTree("ntuple", "list of tracks");
		fTreeTracks->Branch("run", &fRun, "run/I");
		fTreeTracks->Branch("event", &fEventNumber, "event/I");
		fTreeTracks->Branch("subevent", &fSubEventNumber, "subevent/I");
		fTreeTracks->Branch("track", &fTrackNumber, "track/I");
		fTreeTracks->Branch("pdg", &fPdg, "pdg/I");
		fTreeTracks->Branch("flag", &fFlag, "flag/I");
		fTreeTracks->Branch("M", &fM, "M/D");
		fTreeTracks->Branch("P", &fP, "P/D");
		fTreeTracks->Branch("E", &fE, "E/D");
		fTreeTracks->Branch("trkdirX", &fTrkDirX, "trkdirX/D");
		fTreeTracks->Branch("trkdirY", &fTrkDirY, "trkdirY/D");
		fTreeTracks->Branch("trkdirZ", &fTrkDirZ, "trkdirZ/D");
		fTreeTracks->Branch("vtxX", &fVtxX, "vtxX/D");
		fTreeTracks->Branch("vtxY", &fVtxY, "vtxY/D");
		fTreeTracks->Branch("vtxZ", &fVtxZ, "vtxZ/D");
		fTreeTracks->Branch("vtxTime", &fVtxTime, "vtxTime/D");
		fTreeTracks->Branch("endX", &fEndX, "endX/D");
		fTreeTracks->Branch("endY", &fEndY, "endY/D");
		fTreeTracks->Branch("endZ", &fEndZ, "endZ/D");
		fTreeTracks->Branch("pdgparent", &fPdgParent, "pdgparent/I");

		fTreeTracks->SetAutoSave(100);
		gDirectory = tmpd;
	}

	tmpd = gDirectory;
	fFileTracks->cd();
	fTreeTracks->Fill();
	gDirectory = tmpd;
}

void WCSimEventWriter::CloseFiles() {
	TDirectory* tmpd = 0;

	// Write out Tracks
	// ================
	if (fWriteOutTracks) {
		tmpd = gDirectory;
		fFileTracks->cd();
		fTreeTracks->Write();
		fFileTracks->Close();
		gDirectory = tmpd;
	}

	// Write out Digits
	// ================
	if (fWriteOutDigits) {
		tmpd = gDirectory;
		fFileDigits->cd();
		fTreeDigits->Write();
		fFileDigits->Close();
		gDirectory = tmpd;
	}
}

