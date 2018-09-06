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

class WCSimPIDTree: public TObject {
	public:
		WCSimPIDTree(TTree * resultsTree);
		~WCSimPIDTree();

		void MakeOutputTree();

		bool GetEntry(int entry);
		void FillTree();

		void Clear();

		TTree * GetOutputTree() {
			return PIDOutputTree;
		}
		
		void SetNueVsNumu(float NueCCQEvsNumuCCQE);
		void SetNueVsNC(float NueCCQEvsNC);

		bool IsPreselected() {
			return preselected;
		}

		// Get methods
		//Combined Variables
		float GetDeltaCharge2LnL() {
			return deltaCharge2LnL;
		}
		float GetDeltaTime2LnL() {
			return deltaTime2LnL;
		}
		float GetDeltaCharge2LnLOverNHits() {
			return deltaCharge2LnLOverNHits;
		}

		bool GetVeto() {
			return veto;
		}
		int GetNHits() {
			return nHits;
		}
		float GetTotalQ() {
			return totalQ;
		}

		float GetFracHitsUpstream() {
			return fracHitsUpstream;
		}
		float GetFracHitsDownstream() {
			return fracHitsDownstream;
		}
		float GetFracHitsInBottom() {
			return fracHitsInBottom;
		}
		float GetFracHitsInTop() {
			return fracHitsInTop;
		}
		float GetFracHitsAboveMid() {
			return fracHitsAboveMid;
		}
		float GetFracHitsBelowMid() {
			return fracHitsBelowMid;
		}

		float GetFracQUpstream() {
			return fracQUpstream;
		}
		float GetFracQDownstream() {
			return fracQDownstream;
		}
		float GetFracQInBottom() {
			return fracQInBottom;
		}
		float GetFracQInTop() {
			return fracQInTop;
		}
		float GetFracQAboveMid() {
			return fracQAboveMid;
		}
		float GetFracQBelowMid() {
			return fracQBelowMid;
		}

		//Electron Fit Variables
		float GetCharge2LnL_el() {
			return charge2LnL_el;
		}
		float GetTime2LnL_el() {
			return time2LnL_el;
		}
		float GetRecoE_el() {
			return recoE_el;
		}
		float GetRecoEOverQ_el() {
			return recoEOverQ_el;
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

		float GetFracQOutsideRing_el() {
			return fracQOutsideRing_el;
		}
		float GetFracQInRing_el() {
			return fracQInRing_el;
		}
		float GetFracQInRingHole_el() {
			return fracQInRingHole_el;
		}

		float GetFracPredQOutsideRing_el() {
			return fracPredQOutsideRing_el;
		}
		float GetFracPredQInRing_el() {
			return fracPredQInRing_el;
		}
		float GetFracPredQInRingHole_el() {
			return fracPredQInRingHole_el;
		}
		float GetPredictedChargeOverTotalCharge_el() {
			return predictedChargeOverTotalCharge_el;
		}

		float GetVtxRho_el() {
			return vtxRho_el;
		}
		float GetEndRho_el() {
			return endRho_el;
		}
		float GetVtxX_el() {
			return vtxX_el;
		}
		float GetVtxY_el() {
			return vtxY_el;
		}
		float GetVtxZ_el() {
			return vtxZ_el;
		}
		float GetEndX_el() {
			return endX_el;
		}
		float GetEndY_el() {
			return endY_el;
		}
		float GetEndZ_el() {
			return endZ_el;
		}
		float GetDirX_el() {
			return dirX_el;
		}
		float GetDirY_el() {
			return dirY_el;
		}
		float GetDirZ_el() {
			return dirZ_el;
		}

		bool GetEscapes_el() {
			return escapes_el;
		}

		//Muon Fit Variables
		float GetCharge2LnL_mu() {
			return charge2LnL_mu;
		}
		float GetTime2LnL_mu() {
			return time2LnL_mu;
		}
		float GetRecoE_mu() {
			return recoE_mu;
		}
		float GetRecoEOverQ_mu() {
			return recoEOverQ_mu;
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

		float GetFracQOutsideRing_mu() {
			return fracQOutsideRing_mu;
		}
		float GetFracQInRing_mu() {
			return fracQInRing_mu;
		}
		float GetFracQInRingHole_mu() {
			return fracQInRingHole_mu;
		}

		float GetFracPredQOutsideRing_mu() {
			return fracPredQOutsideRing_mu;
		}
		float GetFracPredQInRing_mu() {
			return fracPredQInRing_mu;
		}
		float GetFracPredQInRingHole_mu() {
			return fracPredQInRingHole_mu;
		}
		float GetPredictedChargeOverTotalCharge_mu() {
			return predictedChargeOverTotalCharge_mu;
		}

		float GetVtxRho_mu() {
			return vtxRho_mu;
		}
		float GetEndRho_mu() {
			return endRho_mu;
		}
		float GetVtxX_mu() {
			return vtxX_mu;
		}
		float GetVtxY_mu() {
			return vtxY_mu;
		}
		float GetVtxZ_mu() {
			return vtxZ_mu;
		}
		float GetEndX_mu() {
			return endX_mu;
		}
		float GetEndY_mu() {
			return endY_mu;
		}
		float GetEndZ_mu() {
			return endZ_mu;
		}
		float GetDirX_mu() {
			return dirX_mu;
		}
		float GetDirY_mu() {
			return dirY_mu;
		}
		float GetDirZ_mu() {
			return dirZ_mu;
		}

		bool GetEscapes_mu() {
			return escapes_mu;
		}

		// What type of event is it...
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
		Int_t nHitsAboveMid;
		Int_t nHitsBelowMid;
		Float_t fracHitsUpstream;
		Float_t fracHitsDownstream;
		Float_t fracHitsInBottom;
		Float_t fracHitsInTop;
		Float_t fracHitsAboveMid;
		Float_t fracHitsBelowMid;		
		Float_t totalQ;
		Float_t totalQUpstream;
		Float_t totalQDownstream;
		Float_t totalQInBottom;
		Float_t totalQInTop;
		Float_t totalQAboveMid;
		Float_t totalQBelowMid;
		Float_t fracQUpstream;
		Float_t fracQDownstream;
		Float_t fracQInBottom;
		Float_t fracQInTop;
		Float_t fracQAboveMid;
		Float_t fracQBelowMid;

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
		TruthInfo *fTruthInfo;
		PidInfo *fPidInfo_el;
		PidInfo *fPidInfo_mu;

		// Branches
		TBranch *fBti;
		TBranch *fBpi_el;
		TBranch *fBpi_mu;
		
		//The output tree
		TTree * PIDOutputTree;

		ClassDef(WCSimPIDTree,2)

};

#endif /* WCSIMPIDTREE_HH */
