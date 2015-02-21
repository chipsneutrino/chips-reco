/*
 * WCSimFitterTree.cc
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#include "WCSimLikelihoodTrack.hh"
#include "WCSimFitterTree.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimRecoSummary.hh"
#include "WCSimTrueEvent.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>

WCSimFitterTree::WCSimFitterTree() {
	// TODO Auto-generated constructor stub
	fFitTree = 0x0;
	fTrueTree = 0x0;
	fRecoSummaryTree = 0x0;
	fWCSimTree = 0x0;
	fGeoTree = 0x0;
//	fWCSimLikelihoodRecoEvent = 0x0;

	TTimeStamp ts;
	unsigned int year, month, day, hour, minute, second;
	ts.GetDate(true, 0, &year, &month, &day);
	ts.GetTime(true, 0, &hour, &minute, &second);
	fSaveFileName.Form("fitterTree_%d%d%d_%d%d.root", year, month, day, hour, minute);

	fGeometry = 0x0;
	fWCSimRootEvent = 0x0;
	fRecoSummary = new WCSimRecoSummary();

	fEvent = 0;

	fFitTrackNum = -999;
	f2LnL = -9999.99;
	fFitVtxX = -9999.99;
	fFitVtxY = -9999.99;
	fFitVtxZ = -9999.99;
	fFitVtxT = -9999.99;
	fFitDirTheta = -9999.99;
	fFitDirPhi = -9999.99;
	fFitEnergy = -9999.99;
	fFitPDG = -999;

	// The truth variables:
	fTrueTrackNum = -999;
	fTrueVtxX = -9999.99;
	fTrueVtxY = -9999.99;
	fTrueVtxZ = -9999.99;
	fTrueVtxT = -9999.99;
	fTrueDirTheta = -9999.99;
	fTrueDirPhi = -9999.99;
	fTrueEnergy = -9999.99;
	fTruePDG = -999;

}

WCSimFitterTree::~WCSimFitterTree() {
	// TODO Auto-generated destructor stub
	delete fRecoSummary;
}

void WCSimFitterTree::MakeTree() {

	fSaveFile = new TFile(fSaveFileName.Data(), "RECREATE");

	if(fTrueTree != 0x0)
	{
		delete fTrueTree;
	}
	if(fFitTree != 0x0)
	{
		delete fFitTree;
	}
	if(fGeoTree != 0x0)
	{
		delete fGeoTree;
	}
	if(fGeometry != 0x0)
	{
		delete fGeometry;
	}
	if(fRecoSummaryTree != 0x0)
	{
		delete fRecoSummaryTree;
	}
	if(fWCSimTree != 0x0)
	{
		delete fWCSimTree;
	}

	fTrueTree = new TTree("fTruthTree","Truth");
	fFitTree = new TTree("fFitTree","Fit");
	fGeoTree = new TTree("wcsimGeoT","Geometry");
	fRecoSummaryTree = new TTree("fRecoSummaryTree","wcsimReco");
	fWCSimTree = new TTree("wcsimT","wcsimT");


	fTrueTree->Branch("event", &fEvent);
	fTrueTree->Branch("track", &fTrueTrackNum);
	fTrueTree->Branch("vtxX", &fTrueVtxX);
	fTrueTree->Branch("vtxY", &fTrueVtxY);
	fTrueTree->Branch("vtxZ", &fTrueVtxZ);
	fTrueTree->Branch("vtxT", &fTrueVtxT);
	fTrueTree->Branch("theta", &fTrueDirTheta);
	fTrueTree->Branch("phi", &fTrueDirPhi);
	fTrueTree->Branch("energy", &fTrueEnergy);
	fTrueTree->Branch("pdg", &fTruePDG);

	fFitTree->Branch("event", &fEvent);
	fTrueTree->Branch("track", &fFitTrackNum);
	fFitTree->Branch("2LnL", &f2LnL);
	fFitTree->Branch("vtxX", &fFitVtxX);
	fFitTree->Branch("vtxY", &fFitVtxY);
	fFitTree->Branch("vtxZ", &fFitVtxZ);
	fFitTree->Branch("vtxT", &fFitVtxT);
	fFitTree->Branch("theta", &fFitDirTheta);
	fFitTree->Branch("phi", &fFitDirPhi);
	fFitTree->Branch("energy", &fFitEnergy);
	fFitTree->Branch("pdg", &fFitPDG);


    fGeometry = new WCSimRootGeom();
	fGeoTree->Branch("wcsimrootgeom", "WCSimRootGeom", &fGeometry, 64000,0);
	fWCSimTree->Branch("wcsimrootevent", "WCSimRootEvent", &fWCSimRootEvent, 64000,0);

	fRecoSummaryTree->Branch("WCSimRecoSummary","WCSimRecoSummary", &fRecoSummary, 64000, 0);

	fEvent = 0;

}

TTree* WCSimFitterTree::GetFitTree() {
	return fFitTree;
}

TTree* WCSimFitterTree::GetTruthTree() {
	return fTrueTree;
}

TTree* WCSimFitterTree::GetRecoTree() {
	return fRecoSummaryTree;
}

TTree* WCSimFitterTree::GetGeoTree() {
	return fGeoTree;
}

void WCSimFitterTree::Fill(Int_t iEvent, std::vector<WCSimLikelihoodTrack> bestFitTracks,
		std::vector<WCSimLikelihoodTrack*> trueTracks, Double_t minimum) {

	fEvent = iEvent;

	// Fill the truth and fit tree entries
	if( fTrueTree == 0x0 && fFitTree == 0x0 && fGeoTree == 0x0)
	{
		MakeTree();
	}

	// Now fill the truth and best-fit trees themselves
	std::vector<WCSimLikelihoodTrack>::iterator bfIter;
	f2LnL = minimum;
	for(bfIter = bestFitTracks.begin(); bfIter != bestFitTracks.end(); ++bfIter)
	{
		fFitTrackNum = std::distance(bestFitTracks.begin(), bfIter);
		WCSimLikelihoodTrack track = (*bfIter);
		FillFitTrack(track, minimum);
	}

	std::vector<WCSimLikelihoodTrack*>::iterator truIter;
	for( truIter = trueTracks.begin(); truIter != trueTracks.end(); ++truIter)
	{
		fTrueTrackNum = std::distance(trueTracks.begin(), truIter);
		WCSimLikelihoodTrack track = *(*truIter);
		FillTrueTrack(track);
	}

	// Fill the geometry tree
	if(fGeoTree->GetEntries() == 0)
	{
      
      *fGeometry = *(WCSimGeometry::Instance()->GetWCSimGeometry());
      WCSimGeometry::Instance()->PrintGeometry();

		fGeoTree->Fill();
	}

	// Make the reco summary object for the event display, and fill its tree
	if(fRecoSummaryTree != 0x0)
	{
		MakeRecoSummary( bestFitTracks );
		fRecoSummaryTree->Fill();
	}

	// Copy the root event into the file for the event display to use
	fWCSimRootEvent = WCSimInterface::Instance()->GetWCSimEvent(iEvent);
	fWCSimTree->Fill();

}

void WCSimFitterTree::FillFitTrack(WCSimLikelihoodTrack track, Double_t twoLnL) {

	track.Print();
	fFitVtxX     = track.GetX();
	fFitVtxY     = track.GetY();
	fFitVtxZ     = track.GetZ();
	fFitVtxT     = track.GetT();
	fFitDirTheta = track.GetTheta();
	fFitDirPhi   = track.GetPhi();
	fFitEnergy   = track.GetE();
	fFitPDG      = track.GetPDG();

	if(fFitTree == 0x0)
	{
		MakeTree();
	}
	fFitTree->Fill();
}

void WCSimFitterTree::SaveTree() {

	TDirectory* tmpd = 0;
	tmpd = gDirectory;
	std::cout << " *** WCSimFitter::SaveTree() *** " << std::endl;
	if(fSaveFile == NULL)
	{
		fSaveFile = new TFile(fSaveFileName.Data(), "UPDATE");
	}

	fTrueTree->SetDirectory(fSaveFile);
	fFitTree->SetDirectory(fSaveFile);
	fWCSimTree->SetDirectory(fSaveFile);
	fGeoTree->SetDirectory(fSaveFile);
	fRecoSummaryTree->SetDirectory(fSaveFile);

	// Was it open?
	fSaveFile->cd();
	fTrueTree->Write();
	fFitTree->Write();

	fGeoTree->Write();
	fWCSimTree->Write();
	fRecoSummaryTree->Write();

	fSaveFile->Close();
	fSaveFile = NULL;
	tmpd->cd();
	return;
}


void WCSimFitterTree::FillTrueTrack(WCSimLikelihoodTrack track) {
	fTrueVtxX     = track.GetX();
	fTrueVtxY     = track.GetY();
	fTrueVtxZ     = track.GetZ();
	fTrueVtxT     = track.GetT();
	fTrueDirTheta = track.GetTheta();
	fTrueDirPhi   = track.GetPhi();
	fTrueEnergy   = track.GetE();
	fTruePDG      = track.GetPDG();

	if(fTrueTree == 0x0)
	{
		MakeTree();
	}
	fTrueTree->Fill();
}

void WCSimFitterTree::SetSaveFileName(TString saveName) {
	fSaveFileName = saveName;
}

TString WCSimFitterTree::GetSaveFileName() const
{
	return fSaveFileName;
}

void WCSimFitterTree::MakeRecoSummary(
		std::vector<WCSimLikelihoodTrack> bestFitTracks) {
	fRecoSummary->ResetValues();
	std::sort(bestFitTracks.begin(), bestFitTracks.end(), WCSimLikelihoodTrack::EnergyGreaterThanOrEqual);

	for(unsigned int iTrack = 0; iTrack < bestFitTracks.size() ; ++iTrack )
	{
		WCSimLikelihoodTrack * track = &(bestFitTracks.at(iTrack));
		if(iTrack == 0)
		{
			fRecoSummary->SetVertex( 10*track->GetX(), 10*track->GetY(), 10*track->GetZ() ); // Convert cm to mm
		}
		fRecoSummary->AddPrimary(track->GetPDG(), track->GetE(), track->GetDir());
	}
	return;
}
