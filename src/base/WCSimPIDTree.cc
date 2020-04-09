/*
 * WCSimPIDTree.cc
 *
 *  Created on: 17 Jul 2017
 *      Author: jtingey
 */

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TString.h"
#include "TStopwatch.h"

#include "WCSimOutputTree.hh"
#include "WCSimPIDTree.hh"

#include <iostream>
#include <cstdlib>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimPIDTree);
#endif

///////////////////////////////////////////////////////////////
//  METHODS FOR PIDTreeEntry CLASS
//////////////////////////////////////////////////////////////

WCSimPIDTree::WCSimPIDTree(TTree *resultsTree)
{
	std::cout << " *** WCSimPIDTree::WCSimPIDTree() *** " << std::endl;

	//Set up the TruthInfo
	fTruthInfo = new TruthInfo();
	fBti = resultsTree->GetBranch("TruthInfo");
	fBti->SetAddress(&fTruthInfo);

	//Set up the PidInfo_el
	fPidInfo_el = new PidInfo();
	fBpi_el = resultsTree->GetBranch("PidInfo_ElectronLike");
	fBpi_el->SetAddress(&fPidInfo_el);

	//Set up the PidInfo_mu
	fPidInfo_mu = new PidInfo();
	fBpi_mu = resultsTree->GetBranch("PidInfo_MuonLike");
	fBpi_mu->SetAddress(&fPidInfo_mu);

	std::cout << "fBti Entries -> " << fBti->GetEntries() << std::endl;
	std::cout << "fBpi_el Entries -> " << fBpi_el->GetEntries() << std::endl;
	std::cout << "fBpi_mu Entries -> " << fBpi_mu->GetEntries() << std::endl;

	// Make the TTree
	MakeOutputTree();

	annNueCCQEvsNumuCCQE = -999;
	annNueCCQEvsNC = -999;
}

WCSimPIDTree::~WCSimPIDTree()
{
	// Empty
}

bool WCSimPIDTree::GetEntry(int entry)
{
	// Load the entry into the trees
	fBti->GetEntry(entry);
	fBpi_el->GetEntry(entry);
	fBpi_mu->GetEntry(entry);

	// Check the events are the same...
	if (!(fPidInfo_el->GetTotalQ() == fPidInfo_mu->GetTotalQ()))
	{
		//std::cout << "Error: Charges don't match for entry -> " << entry << std::endl;
		return false;
	}

	//Useful variables to have
	nHits = fPidInfo_el->GetNHits();
	totalQ = fPidInfo_el->GetTotalQ();
	veto = fPidInfo_el->GetVeto();

	if (totalQ <= 0.0 || nHits <= 0)
	{
		//std::cout << "Error: No Hits!" << std::endl;
		return false;
	}

	//Truth variables
	trueType = fTruthInfo->GetType();
	trueBeamPDG = fTruthInfo->GetBeamPDG();
	trueCCEvent = fTruthInfo->IsCC();
	trueNCEvent = fTruthInfo->IsNC();
	trueQEEvent = fTruthInfo->IsQE();
	trueResEvent = fTruthInfo->IsRes();
	trueDISEvent = fTruthInfo->IsDIS();
	trueCohEvent = fTruthInfo->IsCoherent();
	trueNueElectronElasticEvent = fTruthInfo->IsNueElectronElastic();
	trueInverseMuonDecayEvent = fTruthInfo->IsInverseMuonDecay();
	trueBeamE = fTruthInfo->GetBeamE();

	// PID Combined Variables
	deltaCharge2LnL = fPidInfo_el->GetCharge2LnL() - fPidInfo_mu->GetCharge2LnL();
	deltaTime2LnL = fPidInfo_el->GetTime2LnL() - fPidInfo_mu->GetTime2LnL();
	deltaCharge2LnLOverNHits = deltaCharge2LnL / (nHits);

	// PID Hit variables

	fracHitsUpstream = fPidInfo_el->GetFracHitsUpstream();
	fracHitsDownstream = fPidInfo_el->GetFracHitsDownstream();
	fracHitsInBottom = fPidInfo_el->GetFracHitsInBottom();
	fracHitsInTop = fPidInfo_el->GetFracHitsInTop();
	fracHitsAboveMid = fPidInfo_el->GetFracHitsAboveMid();
	fracHitsBelowMid = fPidInfo_el->GetFracHitsBelowMid();
	fracQUpstream = fPidInfo_el->GetFracQUpstream();
	fracQDownstream = fPidInfo_el->GetFracQDownstream();
	fracQInBottom = fPidInfo_el->GetFracQInBottom();
	fracQInTop = fPidInfo_el->GetFracQInTop();
	fracQAboveMid = fPidInfo_el->GetFracQAboveMid();
	fracQBelowMid = fPidInfo_el->GetFracQBelowMid();

	// Muon fit variables
	charge2LnL_mu = fPidInfo_mu->GetCharge2LnL();
	time2LnL_mu = fPidInfo_mu->GetTime2LnL();
	recoE_mu = fPidInfo_mu->GetEnergy();
	recoEOverQ_mu = fPidInfo_mu->GetEnergy() / (totalQ);

	totalQOutsideRing_mu = fPidInfo_mu->GetTotalQOutsideRing();
	totalQInRing_mu = fPidInfo_mu->GetTotalQInRing();
	totalQInRingHole_mu = fPidInfo_mu->GetTotalQInRingHole();
	fracQOutsideRing_mu = fPidInfo_mu->GetFracTotalQOutsideRing();
	fracQInRing_mu = fPidInfo_mu->GetFracTotalQInRing();
	fracQInRingHole_mu = fPidInfo_mu->GetFracTotalQInRingHole();

	nHitsOutsideRing_mu = fPidInfo_mu->GetNHitsOutsideRing();
	nHitsInRing_mu = fPidInfo_mu->GetNHitsInRing();
	nHitsInRingHole_mu = fPidInfo_mu->GetNHitsInRingHole();
	fracHitsOutsideRing_mu = fPidInfo_mu->GetFracNHitsOutsideRing();
	fracHitsInRing_mu = fPidInfo_mu->GetFracNHitsInRing();
	fracHitsInRingHole_mu = fPidInfo_mu->GetFracNHitsInRingHole();
	totalPredQ_mu = fPidInfo_mu->GetPredQInRing() + fPidInfo_mu->GetPredQOutsideRing() + fPidInfo_mu->GetPredQInRingHole();
	totalPredQOutsideRing_mu = fPidInfo_mu->GetPredQOutsideRing();
	totalPredQInRing_mu = fPidInfo_mu->GetPredQInRing();
	totalPredQInRingHole_mu = fPidInfo_mu->GetPredQInRingHole();
	fracPredQOutsideRing_mu = fPidInfo_mu->GetFracPredQOutsideRing();
	fracPredQInRing_mu = fPidInfo_mu->GetFracPredQInRing();
	fracPredQInRingHole_mu = fPidInfo_mu->GetFracPredQInRingHole();
	predictedChargeOverTotalCharge_mu = totalPredQ_mu / (totalQ);

	vtxRho_mu = fPidInfo_mu->GetVtxRho();
	endRho_mu = fPidInfo_mu->GetEndRho();
	vtxX_mu = fPidInfo_mu->GetVtxX();
	vtxY_mu = fPidInfo_mu->GetVtxY();
	vtxZ_mu = fPidInfo_mu->GetVtxZ();
	endX_mu = fPidInfo_mu->GetEndX();
	endY_mu = fPidInfo_mu->GetEndY();
	endZ_mu = fPidInfo_mu->GetEndZ();
	dirX_mu = fPidInfo_mu->GetDirX();
	dirY_mu = fPidInfo_mu->GetDirY();
	dirZ_mu = fPidInfo_mu->GetDirZ();
	escapes_mu = fPidInfo_mu->Escapes();

	// Electron fit variables
	charge2LnL_el = fPidInfo_el->GetCharge2LnL();
	time2LnL_el = fPidInfo_el->GetTime2LnL();
	recoE_el = fPidInfo_el->GetEnergy();
	recoEOverQ_el = fPidInfo_el->GetEnergy() / (totalQ);

	totalQOutsideRing_el = fPidInfo_el->GetTotalQOutsideRing();
	totalQInRing_el = fPidInfo_el->GetTotalQInRing();
	totalQInRingHole_el = fPidInfo_el->GetTotalQInRingHole();
	fracQOutsideRing_el = fPidInfo_el->GetFracTotalQOutsideRing();
	fracQInRing_el = fPidInfo_el->GetFracTotalQInRing();
	fracQInRingHole_el = fPidInfo_el->GetFracTotalQInRingHole();

	nHitsOutsideRing_el = fPidInfo_el->GetNHitsOutsideRing();
	nHitsInRing_el = fPidInfo_el->GetNHitsInRing();
	nHitsInRingHole_el = fPidInfo_el->GetNHitsInRingHole();
	fracHitsOutsideRing_el = fPidInfo_el->GetFracNHitsOutsideRing();
	fracHitsInRing_el = fPidInfo_el->GetFracNHitsInRing();
	fracHitsInRingHole_el = fPidInfo_el->GetFracNHitsInRingHole();
	totalPredQ_el = fPidInfo_el->GetPredQInRing() + fPidInfo_el->GetPredQOutsideRing() + fPidInfo_el->GetPredQInRingHole();
	totalPredQOutsideRing_el = fPidInfo_el->GetPredQOutsideRing();
	totalPredQInRing_el = fPidInfo_el->GetPredQInRing();
	totalPredQInRingHole_el = fPidInfo_el->GetPredQInRingHole();
	fracPredQOutsideRing_el = fPidInfo_el->GetFracPredQOutsideRing();
	fracPredQInRing_el = fPidInfo_el->GetFracPredQInRing();
	fracPredQInRingHole_el = fPidInfo_el->GetFracPredQInRingHole();
	predictedChargeOverTotalCharge_el = totalPredQ_el / (totalQ);

	vtxRho_el = fPidInfo_el->GetVtxRho();
	endRho_el = fPidInfo_el->GetEndRho();
	vtxX_el = fPidInfo_el->GetVtxX();
	vtxY_el = fPidInfo_el->GetVtxY();
	vtxZ_el = fPidInfo_el->GetVtxZ();
	endX_el = fPidInfo_el->GetEndX();
	endY_el = fPidInfo_el->GetEndY();
	endZ_el = fPidInfo_el->GetEndZ();
	dirX_el = fPidInfo_el->GetDirX();
	dirY_el = fPidInfo_el->GetDirY();
	dirZ_el = fPidInfo_el->GetDirZ();
	escapes_el = fPidInfo_el->Escapes();

	// Apply the preselection cuts, keep them loose...
	preselected = !(nHits < 50 || recoE_el < 550 || recoE_mu < 550 || recoE_el > 4950 || recoE_mu > 4950 || vtxRho_el > 1100 || vtxZ_el > 550 || vtxZ_el < -550);

	return true;
}

void WCSimPIDTree::SetNueVsNumu(float NueCCQEvsNumuCCQE)
{
	annNueCCQEvsNumuCCQE = NueCCQEvsNumuCCQE;
}

void WCSimPIDTree::SetNueVsNC(float NueCCQEvsNC)
{
	annNueCCQEvsNC = NueCCQEvsNC;
}

void WCSimPIDTree::FillTree()
{
	PIDOutputTree->Fill();
}

void WCSimPIDTree::Clear()
{
	fTruthInfo->Clear();
	fPidInfo_el->Clear();
	fPidInfo_mu->Clear();
}

void WCSimPIDTree::MakeOutputTree()
{
	PIDOutputTree = new TTree("PIDTree_ann", "PIDTree_ann");

	// PID Variables
	PIDOutputTree->Branch("annNueCCQEvsNumuCCQE", &annNueCCQEvsNumuCCQE);
	PIDOutputTree->Branch("annNueCCQEvsNC", &annNueCCQEvsNC);

	//Preselection
	PIDOutputTree->Branch("preselected", &preselected);

	//Truth Variables
	PIDOutputTree->Branch("trueType", &trueType);
	PIDOutputTree->Branch("trueBeamPDG", &trueBeamPDG);
	PIDOutputTree->Branch("trueCCEvent", &trueCCEvent);
	PIDOutputTree->Branch("trueNCEvent", &trueNCEvent);
	PIDOutputTree->Branch("trueQEEvent", &trueQEEvent);
	PIDOutputTree->Branch("trueResEvent", &trueResEvent);
	PIDOutputTree->Branch("trueDISEvent", &trueDISEvent);
	PIDOutputTree->Branch("trueCohEvent", &trueCohEvent);
	PIDOutputTree->Branch("trueNueElectronElasticEvent", &trueNueElectronElasticEvent);
	PIDOutputTree->Branch("trueInverseMuonDecayEvent", &trueInverseMuonDecayEvent);
	PIDOutputTree->Branch("trueBeamE", &trueBeamE);

	// Combined variables
	PIDOutputTree->Branch("veto", &veto);
	PIDOutputTree->Branch("nHits", &nHits);
	PIDOutputTree->Branch("totalQ", &totalQ);

	PIDOutputTree->Branch("deltaCharge2LnL", &deltaCharge2LnL);
	PIDOutputTree->Branch("deltaTime2LnL", &deltaTime2LnL);
	PIDOutputTree->Branch("deltaCharge2LnLOverNHits", &deltaCharge2LnLOverNHits);

	PIDOutputTree->Branch("fracHitsUpstream", &fracHitsUpstream);
	PIDOutputTree->Branch("fracHitsDownstream", &fracHitsDownstream);
	PIDOutputTree->Branch("fracHitsInBottom", &fracHitsInBottom);
	PIDOutputTree->Branch("fracHitsInTop", &fracHitsInTop);
	PIDOutputTree->Branch("fracHitsAboveMid", &fracHitsAboveMid);
	PIDOutputTree->Branch("fracHitsBelowMid", &fracHitsBelowMid);
	PIDOutputTree->Branch("fracQUpstream", &fracQUpstream);
	PIDOutputTree->Branch("fracQDownstream", &fracQDownstream);
	PIDOutputTree->Branch("fracQInBottom", &fracQInBottom);
	PIDOutputTree->Branch("fracQInTop", &fracQInTop);
	PIDOutputTree->Branch("fracQAboveMid", &fracQAboveMid);
	PIDOutputTree->Branch("fracQBelowMid", &fracQBelowMid);

	// Electron fit variables
	PIDOutputTree->Branch("charge2LnL_el", &charge2LnL_el);
	PIDOutputTree->Branch("time2LnL_el", &time2LnL_el);
	PIDOutputTree->Branch("recoE_el", &recoE_el);
	PIDOutputTree->Branch("recoEOverQ_el", &recoEOverQ_el);

	PIDOutputTree->Branch("fracHitsOutsideRing_el", &fracHitsOutsideRing_el);
	PIDOutputTree->Branch("fracHitsInRing_el", &fracHitsInRing_el);
	PIDOutputTree->Branch("fracHitsInRingHole_el", &fracHitsInRingHole_el);
	PIDOutputTree->Branch("fracQOutsideRing_el", &fracQOutsideRing_el);
	PIDOutputTree->Branch("fracQInRing_el", &fracQInRing_el);
	PIDOutputTree->Branch("fracQInRingHole_el", &fracQInRingHole_el);
	PIDOutputTree->Branch("fracPredQOutsideRing_el", &fracPredQOutsideRing_el);
	PIDOutputTree->Branch("fracPredQInRing_el", &fracPredQInRing_el);
	PIDOutputTree->Branch("fracPredQInRingHole_el", &fracPredQInRingHole_el);
	PIDOutputTree->Branch("predictedChargeOverTotalCharge_el", &predictedChargeOverTotalCharge_el);

	PIDOutputTree->Branch("vtxRho_el", &vtxRho_el);
	PIDOutputTree->Branch("endRho_el", &endRho_el);
	PIDOutputTree->Branch("vtxX_el", &vtxX_el);
	PIDOutputTree->Branch("vtxY_el", &vtxY_el);
	PIDOutputTree->Branch("vtxZ_el", &vtxZ_el);
	PIDOutputTree->Branch("endX_el", &endX_el);
	PIDOutputTree->Branch("endY_el", &endY_el);
	PIDOutputTree->Branch("endZ_el", &endZ_el);
	PIDOutputTree->Branch("dirX_el", &dirX_el);
	PIDOutputTree->Branch("dirY_el", &dirY_el);
	PIDOutputTree->Branch("dirZ_el", &dirZ_el);
	PIDOutputTree->Branch("escapes_el", &escapes_el);

	// Muon fit variables
	PIDOutputTree->Branch("charge2LnL_mu", &charge2LnL_mu);
	PIDOutputTree->Branch("time2LnL_mu", &time2LnL_mu);
	PIDOutputTree->Branch("recoE_mu", &recoE_mu);
	PIDOutputTree->Branch("recoEOverQ_mu", &recoEOverQ_mu);

	PIDOutputTree->Branch("fracHitsOutsideRing_mu", &fracHitsOutsideRing_mu);
	PIDOutputTree->Branch("fracHitsInRing_mu", &fracHitsInRing_mu);
	PIDOutputTree->Branch("fracHitsInRingHole_mu", &fracHitsInRingHole_mu);
	PIDOutputTree->Branch("fracQOutsideRing_mu", &fracQOutsideRing_mu);
	PIDOutputTree->Branch("fracQInRing_mu", &fracQInRing_mu);
	PIDOutputTree->Branch("fracQInRingHole_mu", &fracQInRingHole_mu);
	PIDOutputTree->Branch("fracPredQOutsideRing_mu", &fracPredQOutsideRing_mu);
	PIDOutputTree->Branch("fracPredQInRing_mu", &fracPredQInRing_mu);
	PIDOutputTree->Branch("fracPredQInRingHole_mu", &fracPredQInRingHole_mu);
	PIDOutputTree->Branch("predictedChargeOverTotalCharge_mu", &predictedChargeOverTotalCharge_mu);

	PIDOutputTree->Branch("vtxRho_mu", &vtxRho_mu);
	PIDOutputTree->Branch("endRho_mu", &endRho_mu);
	PIDOutputTree->Branch("vtxX_mu", &vtxX_mu);
	PIDOutputTree->Branch("vtxY_mu", &vtxY_mu);
	PIDOutputTree->Branch("vtxZ_mu", &vtxZ_mu);
	PIDOutputTree->Branch("endX_mu", &endX_mu);
	PIDOutputTree->Branch("endY_mu", &endY_mu);
	PIDOutputTree->Branch("endZ_mu", &endZ_mu);
	PIDOutputTree->Branch("dirX_mu", &dirX_mu);
	PIDOutputTree->Branch("dirY_mu", &dirY_mu);
	PIDOutputTree->Branch("dirZ_mu", &dirZ_mu);
	PIDOutputTree->Branch("escapes_mu", &escapes_mu);
}
