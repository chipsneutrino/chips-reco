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
ClassImp (WCSimPIDTree);
#endif

///////////////////////////////////////////////////////////////
//  METHODS FOR PIDTreeEntry CLASS
//////////////////////////////////////////////////////////////

WCSimPIDTree::WCSimPIDTree(TTree * tree_mu, TTree * tree_el) {

	assert(tree_mu != NULL && tree_el != NULL);
	assert(tree_mu->GetEntriesFast() == tree_el->GetEntriesFast());
	SetTreeBranches(tree_mu, tree_el);

	pidInfo_el = new PidInfo();
	truthInfo_el = new TruthInfo();
	seedInfo_el = new SeedInfo();

	pidInfo_mu = new PidInfo();
	truthInfo_mu = new TruthInfo();
	seedInfo_mu = new SeedInfo();

	MakeOutputTree();

	annNueCCQEvsNumuCCQE = -999;
	annNueCCQEvsNC = -999;

}

WCSimPIDTree::~WCSimPIDTree() {
	// Empty
}

void WCSimPIDTree::GetEntry(int entry) {
	//Load the entry into the trees
	bpi_mu->GetEntry(entry);
	bti_mu->GetEntry(entry);
	bsi_mu->GetEntry(entry);
	bpi_el->GetEntry(entry);
	bti_el->GetEntry(entry);
	bsi_el->GetEntry(entry);

	//Useful variables to have
	totalQ = (pidInfo_mu->GetTotalQ() == pidInfo_el->GetTotalQ()) ? pidInfo_mu->GetTotalQ() : -999;
	veto = (pidInfo_mu->GetVeto() == pidInfo_el->GetVeto()) ? pidInfo_mu->GetVeto() : 1;
	nHits = (pidInfo_mu->GetNHits() == pidInfo_el->GetNHits()) ? pidInfo_mu->GetNHits() : -999;

	//Truth variables
	trueType = (truthInfo_mu->GetType() == truthInfo_el->GetType()) ? truthInfo_mu->GetType() : -999;
	trueBeamPDG = (truthInfo_mu->GetBeamPDG() == truthInfo_el->GetBeamPDG()) ? truthInfo_mu->GetBeamPDG() : -999;
	trueCCEvent = (truthInfo_mu->IsCC() == truthInfo_el->IsCC()) ? truthInfo_mu->IsCC() : -999;
	trueNCEvent = (truthInfo_mu->IsNC() == truthInfo_el->IsNC()) ? truthInfo_mu->IsNC() : -999;
	trueQEEvent = (truthInfo_mu->IsQE() == truthInfo_el->IsQE()) ? truthInfo_mu->IsQE() : -999;
	trueResEvent = (truthInfo_mu->IsRes() == truthInfo_el->IsRes()) ? truthInfo_mu->IsRes() : -999;
	trueDISEvent = (truthInfo_mu->IsDIS() == truthInfo_el->IsDIS()) ? truthInfo_mu->IsDIS() : -999;
	trueCohEvent = (truthInfo_mu->IsCoherent() == truthInfo_el->IsCoherent()) ? truthInfo_mu->IsCoherent() : -999;
	trueNueElectronElasticEvent =
			(truthInfo_mu->IsNueElectronElastic() == truthInfo_el->IsNueElectronElastic()) ?
					truthInfo_mu->IsNueElectronElastic() : -999;
	trueInverseMuonDecayEvent =
			(truthInfo_mu->IsInverseMuonDecay() == truthInfo_el->IsInverseMuonDecay()) ?
					truthInfo_mu->IsInverseMuonDecay() : -999;

	//These are the variables that are actually used in the PID.
	//if(nHits > 0 && totalQ > 0 && !veto)
	if (nHits > 0) {
		// Combined variables
		deltaCharge2LnL = pidInfo_el->GetCharge2LnL() - pidInfo_mu->GetCharge2LnL();
		deltaTime2LnL = pidInfo_el->GetTime2LnL() - pidInfo_mu->GetTime2LnL();
		deltaCharge2LnLOverNHits = deltaCharge2LnL / (nHits);

		nHitsUpstream =
				(pidInfo_mu->GetNHitsUpstream() == pidInfo_el->GetNHitsUpstream()) ?
						pidInfo_mu->GetNHitsUpstream() : -999;
		nHitsDownstream =
				(pidInfo_mu->GetNHitsDownstream() == pidInfo_el->GetNHitsDownstream()) ?
						pidInfo_mu->GetNHitsDownstream() : -999;
		fracHitsUpstream =
				(pidInfo_mu->GetFracHitsUpstream() == pidInfo_el->GetFracHitsUpstream()) ?
						pidInfo_mu->GetFracHitsUpstream() : -999;
		fracHitsDownstream =
				(pidInfo_mu->GetFracHitsDownstream() == pidInfo_el->GetFracHitsDownstream()) ?
						pidInfo_mu->GetFracHitsDownstream() : -999;
		totalQUpstream =
				(pidInfo_mu->GetTotalQUpstream() == pidInfo_el->GetTotalQUpstream()) ?
						pidInfo_mu->GetTotalQUpstream() : -999;
		totalQDownstream =
				(pidInfo_mu->GetTotalQDownstream() == pidInfo_el->GetTotalQDownstream()) ?
						pidInfo_mu->GetTotalQDownstream() : -999;
		fracQUpstream =
				(pidInfo_mu->GetFracQUpstream() == pidInfo_el->GetFracQUpstream()) ?
						pidInfo_mu->GetFracQUpstream() : -999;
		fracQDownstream =
				(pidInfo_mu->GetFracQDownstream() == pidInfo_el->GetFracQDownstream()) ?
						pidInfo_mu->GetFracQDownstream() : -999;

		nRings = (seedInfo_mu->GetNumRings() == seedInfo_el->GetNumRings()) ? seedInfo_mu->GetNumRings() : -999;
		firstRingHeight =
				(seedInfo_mu->GetRingHeight(0) == seedInfo_el->GetRingHeight(0)) ? seedInfo_mu->GetRingHeight(0) : -999;
		lastRingHeight =
				(seedInfo_mu->GetRingHeight(nRings - 1) == seedInfo_el->GetRingHeight(nRings - 1)) ?
						seedInfo_mu->GetRingHeight(nRings - 1) : -999;

		// Muon fit variables
		charge2LnL_mu = pidInfo_mu->GetCharge2LnL();
		time2LnL_mu = pidInfo_mu->GetTime2LnL();
		recoE_mu = pidInfo_mu->GetEnergy();
		recoEOverQ_mu = pidInfo_mu->GetEnergy() / (totalQ);

		totalQOutsideRing_mu = pidInfo_mu->GetTotalQOutsideRing();
		totalQInRing_mu = pidInfo_mu->GetTotalQInRing();
		totalQInRingHole_mu = pidInfo_mu->GetTotalQInRingHole();
		fracQOutsideRing_mu = pidInfo_mu->GetFracTotalQOutsideRing();
		fracQInRing_mu = pidInfo_mu->GetFracTotalQInRing();
		fracQInRingHole_mu = pidInfo_mu->GetFracTotalQInRingHole();

		nHitsOutsideRing_mu = pidInfo_mu->GetNHitsOutsideRing();
		nHitsInRing_mu = pidInfo_mu->GetNHitsInRing();
		nHitsInRingHole_mu = pidInfo_mu->GetNHitsInRingHole();
		fracHitsOutsideRing_mu = pidInfo_mu->GetFracNHitsOutsideRing();
		fracHitsInRing_mu = pidInfo_mu->GetFracNHitsInRing();
		fracHitsInRingHole_mu = pidInfo_mu->GetFracNHitsInRingHole();
		totalPredQ_mu = pidInfo_mu->GetPredQInRing() + pidInfo_mu->GetPredQOutsideRing()
				+ pidInfo_mu->GetPredQInRingHole();
		totalPredQOutsideRing_mu = pidInfo_mu->GetPredQOutsideRing();
		totalPredQInRing_mu = pidInfo_mu->GetPredQInRing();
		totalPredQInRingHole_mu = pidInfo_mu->GetPredQInRingHole();
		fracPredQOutsideRing_mu = pidInfo_mu->GetFracPredQOutsideRing();
		fracPredQInRing_mu = pidInfo_mu->GetFracPredQInRing();
		fracPredQInRingHole_mu = pidInfo_mu->GetFracPredQInRingHole();
		predictedChargeOverTotalCharge_mu = totalPredQ_mu / (totalQ);

		vtxRho_mu = pidInfo_mu->GetVtxRho();
		endRho_mu = pidInfo_mu->GetEndRho();
		vtxX_mu = pidInfo_mu->GetVtxX();
		vtxY_mu = pidInfo_mu->GetVtxY();
		vtxZ_mu = pidInfo_mu->GetVtxZ();
		endX_mu = pidInfo_mu->GetEndX();
		endY_mu = pidInfo_mu->GetEndY();
		endZ_mu = pidInfo_mu->GetEndZ();
		dirX_mu = pidInfo_mu->GetDirX();
		dirY_mu = pidInfo_mu->GetDirY();
		dirZ_mu = pidInfo_mu->GetDirZ();
		escapes_mu = pidInfo_mu->Escapes();

		// Electron fit variables
		charge2LnL_el = pidInfo_el->GetCharge2LnL();
		time2LnL_el = pidInfo_el->GetTime2LnL();
		recoE_el = pidInfo_el->GetEnergy();
		recoEOverQ_el = pidInfo_el->GetEnergy() / (totalQ);

		totalQOutsideRing_el = pidInfo_el->GetTotalQOutsideRing();
		totalQInRing_el = pidInfo_el->GetTotalQInRing();
		totalQInRingHole_el = pidInfo_el->GetTotalQInRingHole();
		fracQOutsideRing_el = pidInfo_el->GetFracTotalQOutsideRing();
		fracQInRing_el = pidInfo_el->GetFracTotalQInRing();
		fracQInRingHole_el = pidInfo_el->GetFracTotalQInRingHole();

		nHitsOutsideRing_el = pidInfo_el->GetNHitsOutsideRing();
		nHitsInRing_el = pidInfo_el->GetNHitsInRing();
		nHitsInRingHole_el = pidInfo_el->GetNHitsInRingHole();
		fracHitsOutsideRing_el = pidInfo_el->GetFracNHitsOutsideRing();
		fracHitsInRing_el = pidInfo_el->GetFracNHitsInRing();
		fracHitsInRingHole_el = pidInfo_el->GetFracNHitsInRingHole();
		totalPredQ_el = pidInfo_el->GetPredQInRing() + pidInfo_el->GetPredQOutsideRing()
				+ pidInfo_el->GetPredQInRingHole();
		totalPredQOutsideRing_el = pidInfo_el->GetPredQOutsideRing();
		totalPredQInRing_el = pidInfo_el->GetPredQInRing();
		totalPredQInRingHole_el = pidInfo_el->GetPredQInRingHole();
		fracPredQOutsideRing_el = pidInfo_el->GetFracPredQOutsideRing();
		fracPredQInRing_el = pidInfo_el->GetFracPredQInRing();
		fracPredQInRingHole_el = pidInfo_el->GetFracPredQInRingHole();
		predictedChargeOverTotalCharge_el = totalPredQ_el / (totalQ);

		vtxRho_el = pidInfo_el->GetVtxRho();
		endRho_el = pidInfo_el->GetEndRho();
		vtxX_el = pidInfo_el->GetVtxX();
		vtxY_el = pidInfo_el->GetVtxY();
		vtxZ_el = pidInfo_el->GetVtxZ();
		endX_el = pidInfo_el->GetEndX();
		endY_el = pidInfo_el->GetEndY();
		endZ_el = pidInfo_el->GetEndZ();
		dirX_el = pidInfo_el->GetDirX();
		dirY_el = pidInfo_el->GetDirY();
		dirZ_el = pidInfo_el->GetDirZ();
		escapes_el = pidInfo_el->Escapes();

		preselected = !(nHits < 50 || recoE_el < 550 || recoE_mu < 550 || recoE_el > 4950 || recoE_mu > 4950
				|| vtxRho_el > 1100 || vtxZ_el > 900 || vtxZ_el < -900 || veto || fracHitsDownstream < 0 || nRings < 0
				|| firstRingHeight < 0 || lastRingHeight < 0);
	} else {
		// Combined variables
		deltaCharge2LnL = -999;
		deltaTime2LnL = -999;
		deltaCharge2LnLOverNHits = -999;

		nHitsUpstream = -999;
		nHitsDownstream = -999;
		fracHitsUpstream = -999;
		fracHitsDownstream = -999;
		totalQUpstream = -999;
		totalQDownstream = -999;
		fracQUpstream = -999;
		fracQDownstream = -999;

		nRings = -999;
		firstRingHeight = -999;
		lastRingHeight = -999;

		// Muon fit variables
		charge2LnL_mu = -999;
		time2LnL_mu = -999;
		recoE_mu = -999;
		recoEOverQ_mu = -999;

		totalQOutsideRing_mu = -999;
		totalQInRing_mu = -999;
		totalQInRingHole_mu = -999;
		fracQOutsideRing_mu = -999;
		fracQInRing_mu = -999;
		fracQInRingHole_mu = -999;

		nHitsOutsideRing_mu = -999;
		nHitsInRing_mu = -999;
		nHitsInRingHole_mu = -999;
		fracHitsOutsideRing_mu = -999;
		fracHitsInRing_mu = -999;
		fracHitsInRingHole_mu = -999;
		totalPredQ_mu = -999;
		totalPredQOutsideRing_mu = -999;
		totalPredQInRing_mu = -999;
		totalPredQInRingHole_mu = -999;
		fracPredQOutsideRing_mu = -999;
		fracPredQInRing_mu = -999;
		fracPredQInRingHole_mu = -999;
		predictedChargeOverTotalCharge_mu = -999;

		vtxRho_mu = -999;
		endRho_mu = -999;
		vtxX_mu = -999;
		vtxY_mu = -999;
		vtxZ_mu = -999;
		endX_mu = -999;
		endY_mu = -999;
		endZ_mu = -999;
		dirX_mu = -999;
		dirY_mu = -999;
		dirZ_mu = -999;
		escapes_mu = 1;

		// Electron fit variables
		charge2LnL_el = -999;
		time2LnL_el = -999;
		recoE_el = -999;
		recoEOverQ_el = -999;

		totalQOutsideRing_el = -999;
		totalQInRing_el = -999;
		totalQInRingHole_el = -999;
		fracQOutsideRing_el = -999;
		fracQInRing_el = -999;
		fracQInRingHole_el = -999;

		nHitsOutsideRing_el = -999;
		nHitsInRing_el = -999;
		nHitsInRingHole_el = -999;
		fracHitsOutsideRing_el = -999;
		fracHitsInRing_el = -999;
		fracHitsInRingHole_el = -999;
		totalPredQ_el = -999;
		totalPredQOutsideRing_el = -999;
		totalPredQInRing_el = -999;
		totalPredQInRingHole_el = -999;
		fracPredQOutsideRing_el = -999;
		fracPredQInRing_el = -999;
		fracPredQInRingHole_el = -999;
		predictedChargeOverTotalCharge_el = -999;

		vtxRho_el = -999;
		endRho_el = -999;
		vtxX_el = -999;
		vtxY_el = -999;
		vtxZ_el = -999;
		endX_el = -999;
		endY_el = -999;
		endZ_el = -999;
		dirX_el = -999;
		dirY_el = -999;
		dirZ_el = -999;
		escapes_el = 1;

		preselected = false;
	}
}

void WCSimPIDTree::SetNueVsNumu(float NueCCQEvsNumuCCQE) {
	annNueCCQEvsNumuCCQE = NueCCQEvsNumuCCQE;
}

void WCSimPIDTree::SetNueVsNC(float NueCCQEvsNC) {
	annNueCCQEvsNC = NueCCQEvsNC;
}

void WCSimPIDTree::FillTree() {
	//Fill the PIDOutputTree
	PIDOutputTree->Fill();
}

void WCSimPIDTree::SetTreeBranches(TTree * tree_mu, TTree * tree_el) {
	///// Muon Like  ...
	bpi_mu = tree_mu->GetBranch("PidInfo");
	bpi_mu->SetAddress(&pidInfo_mu);
	bti_mu = tree_mu->GetBranch("TruthInfo");
	bti_mu->SetAddress(&truthInfo_mu);
	bsi_mu = tree_mu->GetBranch("SeedInfo");
	bsi_mu->SetAddress(&seedInfo_mu);

	///// Electron Like ...
	bpi_el = tree_el->GetBranch("PidInfo");
	bpi_el->SetAddress(&pidInfo_el);
	bti_el = tree_el->GetBranch("TruthInfo");
	bti_el->SetAddress(&truthInfo_el);
	bsi_el = tree_el->GetBranch("SeedInfo");
	bsi_el->SetAddress(&seedInfo_el);
}

void WCSimPIDTree::Clear() {
	pidInfo_mu->Clear();
	truthInfo_mu->Clear();
	seedInfo_mu->Clear();

	pidInfo_el->Clear();
	truthInfo_el->Clear();
	seedInfo_el->Clear();
}

void WCSimPIDTree::MakeOutputTree() {
	PIDOutputTree = new TTree("PIDTree_ann", "PIDTree_ann");

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

	// Combined variables
	PIDOutputTree->Branch("deltaCharge2LnL", &deltaCharge2LnL);
	PIDOutputTree->Branch("deltaTime2LnL", &deltaTime2LnL);
	PIDOutputTree->Branch("deltaCharge2LnLOverNHits", &deltaCharge2LnLOverNHits);

	PIDOutputTree->Branch("veto", &veto);
	PIDOutputTree->Branch("nHits", &nHits);
	PIDOutputTree->Branch("nHitsUpstream", &nHitsUpstream);
	PIDOutputTree->Branch("nHitsDownstream", &nHitsDownstream);
	PIDOutputTree->Branch("fracHitsUpstream", &fracHitsUpstream);
	PIDOutputTree->Branch("fracHitsDownstream", &fracHitsDownstream);
	PIDOutputTree->Branch("totalQ", &totalQ);
	PIDOutputTree->Branch("totalQUpstream", &totalQUpstream);
	PIDOutputTree->Branch("totalQDownstream", &totalQDownstream);
	PIDOutputTree->Branch("fracQUpstream", &fracQUpstream);
	PIDOutputTree->Branch("fracQDownstream", &fracQDownstream);

	PIDOutputTree->Branch("nRings", &nRings);
	PIDOutputTree->Branch("firstRingHeight", &firstRingHeight);
	PIDOutputTree->Branch("lastRingHeight", &lastRingHeight);

	// Muon fit variables
	PIDOutputTree->Branch("charge2LnL_mu", &charge2LnL_mu);
	PIDOutputTree->Branch("time2LnL_mu", &time2LnL_mu);
	PIDOutputTree->Branch("recoE_mu", &recoE_mu);
	PIDOutputTree->Branch("recoEOverQ_mu", &recoEOverQ_mu);

	PIDOutputTree->Branch("totalQOutsideRing_mu", &totalQOutsideRing_mu);
	PIDOutputTree->Branch("totalQInRing_mu", &totalQInRing_mu);
	PIDOutputTree->Branch("totalQInRingHole_mu", &totalQInRingHole_mu);
	PIDOutputTree->Branch("fracQOutsideRing_mu", &fracQOutsideRing_mu);
	PIDOutputTree->Branch("fracQInRing_mu", &fracQInRing_mu);
	PIDOutputTree->Branch("fracQInRingHole_mu", &fracQInRingHole_mu);

	PIDOutputTree->Branch("nHitsOutsideRing_mu", &nHitsOutsideRing_mu);
	PIDOutputTree->Branch("nHitsInRing_mu", &nHitsInRing_mu);
	PIDOutputTree->Branch("nHitsInRingHole_mu", &nHitsInRingHole_mu);
	PIDOutputTree->Branch("fracHitsOutsideRing_mu", &fracHitsOutsideRing_mu);
	PIDOutputTree->Branch("fracHitsInRing_mu", &fracHitsInRing_mu);
	PIDOutputTree->Branch("fracHitsInRingHole_mu", &fracHitsInRingHole_mu);
	PIDOutputTree->Branch("totalPredQ_mu", &totalPredQ_mu);
	PIDOutputTree->Branch("totalPredQOutsideRing_mu", &totalPredQOutsideRing_mu);
	PIDOutputTree->Branch("totalPredQInRing_mu", &totalPredQInRing_mu);
	PIDOutputTree->Branch("totalPredQInRingHole_mu", &totalPredQInRingHole_mu);
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

	// Electron fit variables
	PIDOutputTree->Branch("charge2LnL_el", &charge2LnL_el);
	PIDOutputTree->Branch("time2LnL_el", &time2LnL_el);
	PIDOutputTree->Branch("recoE_el", &recoE_el);
	PIDOutputTree->Branch("recoEOverQ_el", &recoEOverQ_el);

	PIDOutputTree->Branch("totalQOutsideRing_el", &totalQOutsideRing_el);
	PIDOutputTree->Branch("totalQInRing_el", &totalQInRing_el);
	PIDOutputTree->Branch("totalQInRingHole_el", &totalQInRingHole_el);
	PIDOutputTree->Branch("fracQOutsideRing_el", &fracQOutsideRing_el);
	PIDOutputTree->Branch("fracQInRing_el", &fracQInRing_el);
	PIDOutputTree->Branch("fracQInRingHole_el", &fracQInRingHole_el);

	PIDOutputTree->Branch("nHitsOutsideRing_el", &nHitsOutsideRing_el);
	PIDOutputTree->Branch("nHitsInRing_el", &nHitsInRing_el);
	PIDOutputTree->Branch("nHitsInRingHole_el", &nHitsInRingHole_el);
	PIDOutputTree->Branch("fracHitsOutsideRing_el", &fracHitsOutsideRing_el);
	PIDOutputTree->Branch("fracHitsInRing_el", &fracHitsInRing_el);
	PIDOutputTree->Branch("fracHitsInRingHole_el", &fracHitsInRingHole_el);
	PIDOutputTree->Branch("totalPredQ_el", &totalPredQ_el);
	PIDOutputTree->Branch("totalPredQOutsideRing_el", &totalPredQOutsideRing_el);
	PIDOutputTree->Branch("totalPredQInRing_el", &totalPredQInRing_el);
	PIDOutputTree->Branch("totalPredQInRingHole_el", &totalPredQInRingHole_el);
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

	// PID Variables
	PIDOutputTree->Branch("annNueCCQEvsNumuCCQE", &annNueCCQEvsNumuCCQE);
	PIDOutputTree->Branch("annNueCCQEvsNC", &annNueCCQEvsNC);

	//Preselection
	PIDOutputTree->Branch("preselected", &preselected);
}
