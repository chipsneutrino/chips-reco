#ifndef WCSIMOUTPUTTREE_HH
#define WCSIMOUTPUTTREE_HH

/*
 * WCSimOutputTree.hh
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimHitComparison.hh"

#include "TObject.h"

#include <string>
#include <vector>


class TFile;
class TString;
class TTree;
class WCSimRootGeom;
class WCSimLikelihoodRecoEvent;
class WCSimLikelihoodDigitArray;
class WCSimTrueEvent;
class WCSimRootEvent;
class WCSimRecoSummary;


class HitInfo : public TObject{

    public:
        HitInfo();
        HitInfo(
                bool veto, 
                int NHits, int NHitsUpstream, 
                float totalQ, float totalQUpstream
                );
        HitInfo(const HitInfo& other);
        HitInfo& operator=(const HitInfo& rhs);
        ~HitInfo();

        bool GetVeto() const{ return fVeto; }

        int GetNHits() const{ return fNHits; }
        int GetNHitsUpstream() const{ return fNHitsUpstream; }
        int GetNHitsDownstream() const{ return fNHitsDownstream; }
        float GetFracHitsUpstream() const{ return fFracHitsUpstream; }
        float GetFracHitsDownstream() const{ return fFracHitsDownstream; }

        float GetTotalQ() const{ return fTotalQ; }
        float GetTotalQUpstream() const{ return fTotalQUpstream; }
        float GetTotalQDownstream() const{ return fTotalQDownstream; }
        float GetFracQUpstream() const{ return fFracQUpstream; }
        float GetFracQDownstream() const{ return fFracQDownstream; }

    private:

        bool fVeto;

        int fNHits;
        int fNHitsUpstream;
        int fNHitsDownstream;
        float fFracHitsUpstream;
        float fFracHitsDownstream;

        float fTotalQ;
        float fTotalQUpstream;
        float fTotalQDownstream;
        float fFracQUpstream;
        float fFracQDownstream;

        ClassDef(HitInfo,1)

};


class TruthInfo : public TObject{

    public:
        TruthInfo(); 
        TruthInfo(int type, int beamPDG, float beamEnergy, float leadPDG, float leadEnergy);
        TruthInfo(const TruthInfo& other);
        TruthInfo& operator=(const TruthInfo& rhs);
        ~TruthInfo();

        int GetType(){ return fType; }
        float GetBeamE(){ return fBeamEnergy; }
        float GetBeamPDG(){ return fBeamPDG; }
        int GetLeadPDG(){ return fLeadPDG; }
        int GetLeadEnergy(){ return fLeadEnergy; }
        bool IsCC(){ return fIsCC; }
        bool IsNC(){ return fIsNC; }
        bool IsQE(){ return fIsQE; }
        bool IsRes(){ return fIsRes; }
        bool IsDIS(){ return fIsDIS; }
        bool IsCoherent(){ return fIsCoherent; }
        bool IsNueElectronElastic(){ return fIsNueElectronElastic; }
        bool IsInverseMuonDecay(){ return fIsInverseMuonDecay; }
        bool IsOther(){ return fIsOther; }


    private:

        int fType;  // Interaction type code from WCSimTruthSummary
        float fBeamPDG;  // PDG code of the beam particle (usually a neutrino)
        float fBeamEnergy;  // Energy of the incident beam particle
        int fLeadPDG;   // PDG code of the highest energy final state particle
        int fLeadEnergy; // Energy of the highest-energy final state particle
        bool fIsCC;
        bool fIsNC;
        bool fIsQE;
        bool fIsRes;
        bool fIsDIS;
        bool fIsCoherent;
        bool fIsNueElectronElastic;
        bool fIsInverseMuonDecay;
        bool fIsOther; // A CC event that doesn't fall into any of the above
                       // categories - sometimes there isn't a nuance code
                       // for the type of event Genie has made
        
        ClassDef(TruthInfo,1)
};


class RecoInfo : public TObject{

    public:

        RecoInfo();
        RecoInfo(const RecoInfo& other);
        RecoInfo& operator=(const RecoInfo& other);
        ~RecoInfo();

        void SetLikelihoods(float total, float hit, float charge, float time);
        
        void SetTotalQ(float in, float out, float hole);
        void SetNHits(int in, int out, int hole);
        void SetPredQ(float in, float out, float hole);

        void SetVtx(float x, float y, float z);
        void SetEnd(float x, float y, float z);
        void SetDir(float x, float y, float z);

        void SetEscapes(bool escapes);


        float GetTotal2LnL() const{ return fTotal2LnL; }
        float GetHit2LnL() const{ return fHit2LnL; }
        float GetCharge2LnL() const{ return fCharge2LnL; }
        float GetTime2LnL() const{ return fTime2LnL; }
        float GetCutoff2LnL() const{ return fCutoff2LnL; }

        float GetTotalQInRing() const{ return fTotalQInRing; }
        float GetTotalQOutsideRing() const{ return fTotalQOutsideRing; }
        float GetTotalQInRingHole() const{ return fTotalQInRingHole; }

        int GetNHitsInRing(){ return fNHitsInRing; }
        int GetNHitsOutsideRing(){ return fNHitsOutsideRing; }
        int GetNHitsInRingHole(){ return fNHitsInRingHole; }

        float GetPredQInRing() const{ return fPredQInRing; }
        float GetPredQOutsideRing() const{ return fPredQOutsideRing; }
        float GetPredQInRingHole() const{ return fPredQInRingHole; }

        float GetFracTotalQInRing() const{ return fFracTotalQInRing; }
        float GetFracTotalQOutsideRing() const{ return fFracTotalQOutsideRing; }
        float GetFracTotalQInRingHole() const{ return fFracTotalQInRingHole; }

        int GetFracNHitsInRing(){ return fFracNHitsInRing; }
        int GetFracNHitsOutsideRing(){ return fFracNHitsOutsideRing; }
        int GetFracNHitsInRingHole(){ return fFracNHitsInRingHole; }

        float GetFracPredQInRing() const{ return fFracPredQInRing; }
        float GetFracPredQOutsideRing() const{ return fFracPredQOutsideRing; }
        float GetFracPredQInRingHole() const{ return fFracPredQInRingHole; }

        float GetPredictedOverTotalCharge() const{ return fPredictedOverTotalCharge; }

        float GetVtxX() const{ return fVtxX; }
        float GetVtxY() const{ return fVtxY; }
        float GetVtxZ() const{ return fVtxZ; }
        float GetVtxRho() const{ return fVtxRho; }

        float GetEndX() const{ return fEndX; }
        float GetEndY() const{ return fEndY; }
        float GetEndZ() const{ return fEndZ; }
        float GetEndRho() const{ return fEndRho; }

        float GetDirX() const{ return fDirX; }
        float GetDirY() const{ return fDirY; }
        float GetDirZ() const{ return fDirZ; }

        bool Escapes() const{ return fEscapes; }
    private:
        
        float fTotal2LnL;
        float fHit2LnL;
        float fCharge2LnL;
        float fTime2LnL;
        float fCutoff2LnL;

        float fTotalQInRing;
        float fTotalQOutsideRing;
        float fTotalQInRingHole;

        int fNHitsInRing;
        int fNHitsOutsideRing;
        int fNHitsInRingHole;

        float fPredQInRing;
        float fPredQOutsideRing;
        float fPredQInRingHole;

        float fFracTotalQInRing;
        float fFracTotalQOutsideRing;
        float fFracTotalQInRingHole;

        int fFracNHitsInRing;
        int fFracNHitsOutsideRing;
        int fFracNHitsInRingHole;

        float fFracPredQInRing;
        float fFracPredQOutsideRing;
        float fFracPredQInRingHole;

        float fPredictedOverTotalCharge;

        float fVtxX;
        float fVtxY;
        float fVtxZ;
        float fVtxRho;

        float fEndX;
        float fEndY;
        float fEndZ;
        float fEndRho;

        float fDirX;
        float fDirY;
        float fDirZ;

        bool fEscapes;

        ClassDef(RecoInfo, 1)

};



class WCSimOutputTree : public TObject{
public:
	WCSimOutputTree(const TString &saveFileName);
	virtual ~WCSimOutputTree();

	void MakeTree();
	void SaveTree();
	void SetSaveFileName(TString saveName);
	TString GetSaveFileName() const;

    void Fill(
              bool failed,
              Int_t iEvent,
              std::string recoType,
              WCSimRecoSummary& summ,
              WCSimHitComparison& hitComp,
              HitInfo& hitInfo,
              RecoInfo& recoInfo,
              TruthInfo& truthInfo
              );

    void MakeSaveFile();
    void SetInputFile(const TString &location);


private:
    void DeletePointersIfExist();
    void BuildUID();

	TFile * fSaveFile;
	TString fSaveFileName;
    TTree * fResultsTree;
    TTree * fGeoTree;

    WCSimRootGeom * fGeometry;

    std::string fUID;
    std::string fInputFile;
    int fInputEvent;
	WCSimRecoSummary * fRecoSummary;
    WCSimHitComparison * fHitComparison;
    HitInfo * fHitInfo;
    RecoInfo * fRecoInfo;
    TruthInfo * fTruthInfo;
    std::string fRecoType;
    int fEvent;
    bool fFailed;

    ClassDef(WCSimOutputTree, 1)
};

#endif /* WCSIMOUTPUTTREE_HH */
