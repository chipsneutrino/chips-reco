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
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimRecoSummary.hh"
#include "WCSimTrueEvent.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>

WCSimFitterTree::WCSimFitterTree(const TString &saveFileName) : 
                fSaveFile(0x0),
                fTrueTree(0x0), fFitTree(0x0), fRecoSummaryTree(0x0),
                fWCSimTree(0x0), fGeoTree(0x0), fHitComparisonTree(0x0),
                fRecoFailureTree(0x0), fGeometry(0x0), fWCSimRootEvent(0x0),
                fFailedEvent(0), fEvent(0), fFitTrackNum(-999), f2LnL(-999.9),
                fFitVtxX(-999.9), fFitVtxY(-999.9), fFitVtxZ(-999.9),
                fFitVtxT(-999.9), fFitDirTheta(-999.9), fFitDirPhi(-999.9),
                fFitEnergy(-999.9), fFitConversionDistance(0), 
                fFitEscapes(false), fFitPDG(-999),
                fTrueTrackNum(-999),
                fTrueVtxX(-999.9), fTrueVtxY(-999.9), fTrueVtxZ(-999.9),
                fTrueVtxT(-999.9), fTrueDirTheta(-999.9), fTrueDirPhi(-999.9),
                fTrueEnergy(-999.9), fTrueEscapes(false), 
                fTruePDG(-999), fHitComparison(0x0){

	// TODO Auto-generated constructor stub
	fSaveFileName.Form("%s_tree.root", saveFileName.Data());
	fRecoSummary = new WCSimRecoSummary();
    fSaveFile = 0x0;
}

WCSimFitterTree::~WCSimFitterTree() {
	// TODO Auto-generated destructor stub
	delete fRecoSummary;
}

void WCSimFitterTree::MakeTree() {
    std::cout << "fSaveFileName = " << fSaveFileName << " and file is " << fSaveFile << std::endl;
	fSaveFile->cd();
	if(fTrueTree != 0x0)
	{
		delete fTrueTree;
		fTrueTree = 0x0;
	}
	if(fFitTree != 0x0)
	{
		delete fFitTree;
		fFitTree = 0x0;
	}
	if(fGeoTree != 0x0)
	{
		delete fGeoTree;
		fGeoTree = 0x0;
	}
	if(fGeometry != 0x0)
	{
		delete fGeometry;
		fGeometry = 0x0;
	}
	if(fRecoSummaryTree != 0x0)
	{
		delete fRecoSummaryTree;
		fRecoSummaryTree = 0x0;
	}
	if(fWCSimTree != 0x0)
	{
		delete fWCSimTree;
		fWCSimTree = 0x0;
	}
	if(fHitComparisonTree != 0x0)
	{
		delete fHitComparisonTree;
		fHitComparisonTree = 0x0;
	}

	fTrueTree = new TTree("fTruthTree","Truth");
	fFitTree = new TTree("fFitTree","Fit");
	fGeoTree = new TTree("wcsimGeoT","Geometry");
	fRecoSummaryTree = new TTree("fRecoSummaryTree","wcsimReco");
	fWCSimTree = new TTree("wcsimT","wcsimT");
  fHitComparisonTree = new TTree("fHitComparisonTree","Hit comparison");
  fRecoFailureTree = new TTree("fRecoFailureTree","Reco failures");
  
  fHitComparison = new HitComparison();

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
	fTrueTree->Branch("escapes", &fTrueEscapes);
  

	fFitTree->Branch("event", &fEvent);
	fFitTree->Branch("track", &fFitTrackNum);
	fFitTree->Branch("2LnL", &f2LnL);
	fFitTree->Branch("vtxX", &fFitVtxX);
	fFitTree->Branch("vtxY", &fFitVtxY);
	fFitTree->Branch("vtxZ", &fFitVtxZ);
	fFitTree->Branch("vtxT", &fFitVtxT);
	fFitTree->Branch("theta", &fFitDirTheta);
	fFitTree->Branch("phi", &fFitDirPhi);
	fFitTree->Branch("energy", &fFitEnergy);
	fFitTree->Branch("pdg", &fFitPDG);
	fFitTree->Branch("escapes", &fFitEscapes);
  fFitTree->Branch("convDist",&fFitConversionDistance);


  fGeometry = new WCSimRootGeom();
	fGeoTree->Branch("wcsimrootgeom", "WCSimRootGeom", &fGeometry, 64000,0);
	fWCSimTree->Branch("wcsimrootevent", "WCSimRootEvent", &fWCSimRootEvent, 64000,0);

	fRecoSummaryTree->Branch("WCSimRecoSummary","WCSimRecoSummary", &fRecoSummary, 64000, 0);

  fHitComparisonTree->Branch("eventID",&fEvent);
  fHitComparisonTree->Branch("pmtID",&(fHitComparison->pmtID));
  fHitComparisonTree->Branch("pmtX",&(fHitComparison->pmtX));
  fHitComparisonTree->Branch("pmtY",&(fHitComparison->pmtY));
  fHitComparisonTree->Branch("pmtZ",&(fHitComparison->pmtZ));
  fHitComparisonTree->Branch("trueT",&(fHitComparison->trueT));
  fHitComparisonTree->Branch("predT",&(fHitComparison->predT));
  fHitComparisonTree->Branch("correctPredT",&(fHitComparison->correctPredT));
  fHitComparisonTree->Branch("predQ",&(fHitComparison->predQ));
  fHitComparisonTree->Branch("trueQ",&(fHitComparison->trueQ));
  fHitComparisonTree->Branch("correctPredQ",&(fHitComparison->correctPredQ));
  fHitComparisonTree->Branch("minus2LnL",&(fHitComparison->minus2LnL));
  fHitComparisonTree->Branch("correctMinus2LnL",&(fHitComparison->correctMinus2LnL));
  fHitComparisonTree->Branch("charge2LnL",&(fHitComparison->charge2LnL));
  fHitComparisonTree->Branch("correctCharge2LnL",&(fHitComparison->correctCharge2LnL));

  fHitComparisonTree->Branch("time2LnL",&(fHitComparison->time2LnL));
  fHitComparisonTree->Branch("correctTime2LnL",&(fHitComparison->correctTime2LnL));

  fRecoFailureTree->Branch("failedEvent",&fFailedEvent);
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

TTree* WCSimFitterTree::GetHitComparisonTree() {
  return fHitComparisonTree;
}

void WCSimFitterTree::Fill(Int_t iEvent,
						   std::vector<WCSimLikelihoodTrackBase*> bestFitTracks,
						   std::vector<Bool_t> bestFitEscapes,
						   std::vector<WCSimLikelihoodTrackBase*> trueTracks,
						   std::vector<Bool_t> trueTrackEscapes,
						   Double_t minimum) {

	fEvent = iEvent;

	// Fill the truth and fit tree entries
	if( fTrueTree == 0x0 && fFitTree == 0x0 && fGeoTree == 0x0)
	{
		MakeTree();
	}

	// Now fill the truth and best-fit trees themselves
	std::vector<WCSimLikelihoodTrackBase*>::iterator bfIter;
	f2LnL = minimum;
	for(bfIter = bestFitTracks.begin(); bfIter != bestFitTracks.end(); ++bfIter)
	{
		fFitTrackNum = std::distance(bestFitTracks.begin(), bfIter);
		Bool_t escapes = bestFitEscapes.at(fFitTrackNum);
		FillFitTrack((*bfIter), escapes, minimum);
  }

	std::vector<WCSimLikelihoodTrackBase*>::iterator truIter;
	for( truIter = trueTracks.begin(); truIter != trueTracks.end(); ++truIter)
	{
		fTrueTrackNum = std::distance(trueTracks.begin(), truIter);
		Bool_t escapes = trueTrackEscapes.at(fTrueTrackNum);
		WCSimLikelihoodTrackBase* track = (*truIter);
    std::cout << "Fill true track...";
		FillTrueTrack(track, escapes);
    std::cout << " ... done" << std::endl;
	}

	// Fill the geometry tree
	if(fGeoTree->GetEntries() == 0)
	{
      std::cout << "Fill geo tree ..." ;
      *fGeometry = *(WCSimGeometry::Instance()->GetWCSimGeometry());
      WCSimGeometry::Instance()->PrintGeometry();
      
		fGeoTree->Fill();
    std::cout << "... done" << std::endl;
	}

	// Make the reco summary object for the event display, and fill its tree
	if(fRecoSummaryTree != 0x0)
	{
    std::cout << "Make reco summary ...";
		MakeRecoSummary( bestFitTracks );
    std::cout << "... made it - fill tree ..." ;
		fRecoSummaryTree->Fill();
    std::cout << "... done" << std::endl;
	}

	// Copy the root event into the file for the event display to use
  std::cout << "Make root event ... " ;
	fWCSimRootEvent = WCSimInterface::Instance()->GetWCSimEvent(iEvent);
  std::cout << " ... made - filling tree ... " ;
	fWCSimTree->Fill();
  std::cout << " ... done" << std::endl;
  SaveTree();
  std::cout << "Saved tree" << std::endl;

}

void WCSimFitterTree::FillFitTrack(WCSimLikelihoodTrackBase * track, Bool_t escapes, Double_t twoLnL) {
  std::cout << "Filling fit track" << std::endl;
	track->Print();
	fFitVtxX     = track->GetX();
	fFitVtxY     = track->GetY();
	fFitVtxZ     = track->GetZ();
	fFitVtxT     = track->GetT();
	fFitDirTheta = track->GetTheta();
	fFitDirPhi   = track->GetPhi();
	fFitEnergy   = track->GetE();
	fFitPDG      = track->GetPDG();
  fFitConversionDistance = track->GetConversionDistance();
	fFitEscapes = escapes;

	if(fFitTree == 0x0)
	{
		MakeTree();
	}
	fFitTree->Fill();
}

void WCSimFitterTree::SaveTree() {

	TDirectory* tmpd = 0;
	tmpd = gDirectory;
	fSaveFile->cd();
	std::cout << " *** WCSimFitter::SaveTree() *** " << std::endl;
	std::cout << "Save file name = " << fSaveFileName << std::endl;
	//std::cout << "fTrueTree = " << fTrueTree << std::endl;
	//TrueTree->Print();

	assert(fSaveFile->IsOpen());
	//std::cout << "fTrueTree = " << fTrueTree << std::endl;
	//fTrueTree->Print();
	fTrueTree->Write("",TObject::kOverwrite);
	fFitTree->Write("",TObject::kOverwrite);

	fGeoTree->Write("",TObject::kOverwrite);
	fWCSimTree->Write("",TObject::kOverwrite);
	fRecoSummaryTree->Write("",TObject::kOverwrite);

	fHitComparisonTree->Write("",TObject::kOverwrite);
	fRecoFailureTree->Write("",TObject::kOverwrite);
	tmpd->cd();
	return;
}


void WCSimFitterTree::FillTrueTrack(WCSimLikelihoodTrackBase * track, Bool_t escapes) {
	std::cout << "Filling true tree" << std::endl;
	track->Print();
	fTrueVtxX     = track->GetX();
	fTrueVtxY     = track->GetY();
	fTrueVtxZ     = track->GetZ();
	fTrueVtxT     = track->GetT();
	fTrueDirTheta = track->GetTheta();
	fTrueDirPhi   = track->GetPhi();
	fTrueEnergy   = track->GetE();
	fTruePDG      = track->GetPDG();

	fTrueEscapes = escapes;

	if(fTrueTree == 0x0)
	{
		MakeTree();
	}
	fTrueTree->Fill();
	std::cout << "Showing fTrueTree(0) " << std::endl;
	fTrueTree->Show(0);
}

void WCSimFitterTree::SetSaveFileName(TString saveName) {
	fSaveFileName = saveName;
}

TString WCSimFitterTree::GetSaveFileName() const
{
	return fSaveFileName;
}

void WCSimFitterTree::MakeRecoSummary(
		std::vector<WCSimLikelihoodTrackBase*> bestFitTracks) {
	fRecoSummary->ResetValues();

  std::cout << "There are " << bestFitTracks.size() << " best-fit tracks" << std::endl;
	std::sort(bestFitTracks.begin(), bestFitTracks.end(), WCSimLikelihoodTrackBase::EnergyGreaterThanOrEqualPtrs);
	for(unsigned int iTrack = 0; iTrack < bestFitTracks.size() ; ++iTrack )
	{
    std::cout << "Make track " << iTrack << std::endl;
		WCSimLikelihoodTrackBase * track = (bestFitTracks.at(iTrack));
    std::cout << "track:" << track->GetX() << ", " << track->GetY() << ", " << track->GetZ() << "    ";
    std::cout << track->GetType() << std::endl;
    track->Print();
		if(iTrack == 0)
		{
			fRecoSummary->SetVertex( 10*track->GetX(), 10*track->GetY(), 10*track->GetZ() ); // Convert cm to mm
      fRecoSummary->SetVertexT(track->GetT());
		}
    std::cout << "Making reco summary" << std::endl;
		fRecoSummary->AddPrimary(track->GetPDG(), track->GetE(), track->GetDir());
	}
	return;
}

void WCSimFitterTree::FillHitComparison(
    const int &event,
    WCSimLikelihoodDigitArray* digitArray,
		const std::vector<double>& predictedCharges,
		const std::vector<double>& correctPredictedCharges,
		const std::vector<double>& measuredCharges,
    const std::vector<double>& predictedTimes,
    const std::vector<double>& correctPredictedTimes,
		const std::vector<double>& total2LnLs,
		const std::vector<double>& correct2LnLs,
		const std::vector<double>& q2LnLs,
		const std::vector<double>& correctQ2LnLs,
		const std::vector<double>& t2LnLs,
    const std::vector<double>& correctT2LnLs) {
  assert(predictedCharges.size() == measuredCharges.size() && predictedCharges.size() == total2LnLs.size() && predictedTimes.size() == predictedCharges.size());
  fEvent = event;
  for(unsigned int iPMT = 0; iPMT < predictedCharges.size(); ++iPMT)
  {
    WCSimLikelihoodDigit * digit = digitArray->GetDigit(iPMT);

    fHitComparison->Set(iPMT, digit->GetX(), digit->GetY(), digit->GetZ(), 
                        digit->GetT(),
                        predictedTimes.at(iPMT),
                        correctPredictedTimes.at(iPMT),
                        measuredCharges.at(iPMT),
                        predictedCharges.at(iPMT),
                        correctPredictedCharges.at(iPMT),
                        total2LnLs.at(iPMT),
                        correct2LnLs.at(iPMT),
                        q2LnLs.at(iPMT),
                        correctQ2LnLs.at(iPMT),
                        t2LnLs.at(iPMT),
                        correctT2LnLs.at(iPMT));
    fHitComparisonTree->Fill();
  }
}

void WCSimFitterTree::FillRecoFailures(const int &event)
{
  fFailedEvent = event;
  fRecoFailureTree->Fill();

}

void WCSimFitterTree::MakeSaveFile()
{
    std::cout << " *** WCSimFitterTree::MakeSaveFile *** " << std::endl;
    if(fSaveFile != NULL)
    {
        std::cerr << "Save file already exists, so you can't make it again" << std::endl;
        std::cout << "Name is " << fSaveFileName.Data() << std::endl;
        return;
    }
	TDirectory * tmpd = gDirectory;

	fSaveFile = new TFile(fSaveFileName.Data(), "CREATE");
    std::cout << "  Plots to be saved in file: " << fSaveFileName.Data() << std::endl;

    tmpd->cd();
}
