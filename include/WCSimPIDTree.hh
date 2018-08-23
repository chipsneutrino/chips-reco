#ifndef WCSIMPIDTREE_HH
#define WCSIMPIDTREE_HH

/*
 * WCSimPIDTree.hh
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

#include "TObject.h"

#include "WCSimOutputTree.hh"

#include <iostream>
#include <cstdlib>

class TTree;
class TruthInfo;
class PidInfo;
class SeedInfo;

class WCSimPIDTree: public TObject {
	public:
		WCSimPIDTree(TTree * tree_mu, TTree * tree_el);
		~WCSimPIDTree();

		// NEED TO ADD ALL THE STUFF HERE!!!SS
		bool GetVeto() {
			return veto;
		}
		int GetNHits() {
			return nHits;
		}

		int GetNHitsUpstream() {
			return nHitsUpstream;
		}
		int GetNHitsDownstream() {
			return nHitsDownstream;
		}
		int GetNHitsInBottom() {
			return nHitsInBottom;
		}
		int GetNHitsInTop() {
			return nHitsInTop;
		}

		float GetTotalQ() {
			return totalQ;
		}
		float GetFracHitsDownstream() {
			return fracHitsDownstream;
		}

		float GetDeltaCharge2LnL() {
			return deltaCharge2LnL;
		}
		float GetDeltaTime2LnL() {
			return deltaTime2LnL;
		}
		float GetDeltaCharge2LnLOverNHits() {
			return deltaCharge2LnLOverNHits;
		}

		float GetFracHitsOutsideRing_mu() {
			return fracHitsOutsideRing_mu;
		}
		float GetFracHitsInRing_mu() {
			return fracHitsInRing_mu;
		}
		float GetFracHitsInRingHole_mu() {
			return fracHitsInRingHole_mu;
		}
		float GetFracHitsOutsideRing_el() {
			return fracHitsOutsideRing_el;
		}
		float GetFracHitsInRing_el() {
			return fracHitsInRing_el;
		}
		float GetFracHitsInRingHole_el() {
			return fracHitsInRingHole_el;
		}

		float GetFracQOutsideRing_mu() {
			return fracQOutsideRing_mu;
		}
		float GetFracQInRing_mu() {
			return fracQInRing_mu;
		}
		float GetFracQInRingHole_mu() {
			return fracQInRingHole_mu;
		}
		float GetFracQOutsideRing_el() {
			return fracQOutsideRing_el;
		}
		float GetFracQInRing_el() {
			return fracQInRing_el;
		}
		float GetFracQInRingHole_el() {
			return fracQInRingHole_el;
		}

		float GetFracPredQOutsideRing_mu() {
			return fracPredQOutsideRing_mu;
		}
		float GetFracPredQOutsideRing_el() {
			return fracPredQOutsideRing_el;
		}

		float GetPredictedChargeOverTotalCharge_mu() {
			return predictedChargeOverTotalCharge_mu;
		}
		float GetPredictedChargeOverTotalCharge_el() {
			return predictedChargeOverTotalCharge_el;
		}

		float GetRecoEOverQ_mu() {
			return recoEOverQ_mu;
		}
		float GetRecoEOverQ_el() {
			return recoEOverQ_el;
		}

		float GetDirX_el() {
			return dirX_el;
		}

		float GetNRings() {
			return nRings;
		}
		float GetFirstRingHeight() {
			return firstRingHeight;
		}
		float GetLastRingHeight() {
			return lastRingHeight;
		}

		void SetNueVsNumu(float NueCCQEvsNumuCCQE);
		void SetNueVsNC(float NueCCQEvsNC);

		bool IsPreselected() {
			return preselected;
		}

		int GetTrueType() {
			return trueType;
		}
		int GetTrueBeamPDG() {
			return trueBeamPDG;
		}
		int GetTrueCCEvent() {
			return trueCCEvent;
		}
		int GetTrueNCEvent() {
			return trueNCEvent;
		}
		int GetTrueQEEvent() {
			return trueQEEvent;
		}
		int GetTrueResEvent() {
			return trueResEvent;
		}
		int GetTrueDISEvent() {
			return trueDISEvent;
		}
		int GetTrueCohEvent() {
			return trueCohEvent;
		}
		int GetTrueNueElectronElasticEvent() {
			return trueNueElectronElasticEvent;
		}
		int GetTrueInverseMuonDecayEvent() {
			return trueInverseMuonDecayEvent;
		}

		void GetEntry(int entry);
		void FillTree();

		void SetTreeBranches(TTree * tree_mu, TTree * tree_el);
		void Clear();
		void MakeOutputTree();

		TTree * GetOutputTree() {
			return PIDOutputTree;
		}

	private:
		//Truth Variables
		Int_t trueType;
		Int_t trueBeamPDG;
		Bool_t trueCCEvent;
		Bool_t trueNCEvent;
		Bool_t trueQEEvent;
		Bool_t trueResEvent;
		Bool_t trueDISEvent;
		Bool_t trueCohEvent;
		Bool_t trueNueElectronElasticEvent;
		Bool_t trueInverseMuonDecayEvent;

		//Combined Variables
		Float_t deltaCharge2LnL;
		Float_t deltaTime2LnL;
		Float_t deltaCharge2LnLOverNHits;

		Bool_t veto;
		Int_t nHits;
		Int_t nHitsUpstream;
		Int_t nHitsDownstream;
		Int_t nHitsInBottom;
		Int_t nHitsInTop;
		Float_t fracHitsUpstream;
		Float_t fracHitsDownstream;
		Float_t fracHitsInBottom;
		Float_t fracHitsInTop;
		Float_t totalQ;
		Float_t totalQUpstream;
		Float_t totalQDownstream;
		Float_t totalQInBottom;
		Float_t totalQInTop;
		Float_t fracQUpstream;
		Float_t fracQDownstream;
		Float_t fracQInBottom;
		Float_t fracQInTop;

		Int_t nRings;
		Float_t firstRingHeight;
		Float_t lastRingHeight;

		//Muon Fit Variables
		Float_t charge2LnL_mu;
		Float_t time2LnL_mu;
		Float_t recoE_mu;
		Float_t recoEOverQ_mu;

		Float_t totalQOutsideRing_mu;
		Float_t totalQInRing_mu;
		Float_t totalQInRingHole_mu;
		Float_t fracQOutsideRing_mu;
		Float_t fracQInRing_mu;
		Float_t fracQInRingHole_mu;

		Int_t nHitsOutsideRing_mu;
		Int_t nHitsInRing_mu;
		Int_t nHitsInRingHole_mu;
		Float_t fracHitsOutsideRing_mu;
		Float_t fracHitsInRing_mu;
		Float_t fracHitsInRingHole_mu;
		Float_t totalPredQ_mu;
		Float_t totalPredQOutsideRing_mu;
		Float_t totalPredQInRing_mu;
		Float_t totalPredQInRingHole_mu;
		Float_t fracPredQOutsideRing_mu;
		Float_t fracPredQInRing_mu;
		Float_t fracPredQInRingHole_mu;
		Float_t predictedChargeOverTotalCharge_mu;

		Float_t vtxRho_mu;
		Float_t endRho_mu;
		Float_t vtxX_mu;
		Float_t vtxY_mu;
		Float_t vtxZ_mu;
		Float_t endX_mu;
		Float_t endY_mu;
		Float_t endZ_mu;
		Float_t dirX_mu;
		Float_t dirY_mu;
		Float_t dirZ_mu;
		Bool_t escapes_mu;

		//Electron Fit Variables
		Float_t charge2LnL_el;
		Float_t time2LnL_el;
		Float_t recoE_el;
		Float_t recoEOverQ_el;

		Float_t totalQOutsideRing_el;
		Float_t totalQInRing_el;
		Float_t totalQInRingHole_el;
		Float_t fracQOutsideRing_el;
		Float_t fracQInRing_el;
		Float_t fracQInRingHole_el;

		Int_t nHitsOutsideRing_el;
		Int_t nHitsInRing_el;
		Int_t nHitsInRingHole_el;
		Float_t fracHitsOutsideRing_el;
		Float_t fracHitsInRing_el;
		Float_t fracHitsInRingHole_el;
		Float_t totalPredQ_el;
		Float_t totalPredQOutsideRing_el;
		Float_t totalPredQInRing_el;
		Float_t totalPredQInRingHole_el;
		Float_t fracPredQOutsideRing_el;
		Float_t fracPredQInRing_el;
		Float_t fracPredQInRingHole_el;
		Float_t predictedChargeOverTotalCharge_el;

		Float_t vtxRho_el;
		Float_t endRho_el;
		Float_t vtxX_el;
		Float_t vtxY_el;
		Float_t vtxZ_el;
		Float_t endX_el;
		Float_t endY_el;
		Float_t endZ_el;
		Float_t dirX_el;
		Float_t dirY_el;
		Float_t dirZ_el;
		Bool_t escapes_el;

		// PID variables
		Float_t annNueCCQEvsNumuCCQE;
		Float_t annNueCCQEvsNC;

		// Preselection variable
		Bool_t preselected;

		//Declaration of leaves types
		PidInfo *pidInfo_el;
		TruthInfo *truthInfo_el;
		SeedInfo *seedInfo_el;

		PidInfo *pidInfo_mu;
		TruthInfo *truthInfo_mu;
		SeedInfo *seedInfo_mu;

		// Branches
		TBranch *bpi_el;
		TBranch *bti_el;
		TBranch *bsi_el;
		TBranch *bpi_mu;
		TBranch *bti_mu;
		TBranch *bsi_mu;

		//The output tree
		TTree * PIDOutputTree;

		ClassDef(WCSimPIDTree,2)

};

#endif /* WCSIMPIDTREE_HH */
