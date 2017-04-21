#ifndef WCSIMOUTPUTTREE_HH
#define WCSIMOUTPUTTREE_HH

/*
 * WCSimOutputTree.hh
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimHitComparison.hh"
#include "WCSimLikelihoodTrackFactory.hh"

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


class SeedInfo : public TObject{

	//For now can just cope with 1 seed (can have multiple tracks) per event.
	//This should be expanded so that you can store multiple seeds per event.

	public:
		SeedInfo();
		SeedInfo(std::vector<WCSimLikelihoodTrackBase*> tracks);
		SeedInfo(WCSimLikelihoodTrackBase* track);
		SeedInfo(const SeedInfo& other);
		SeedInfo& operator=(const SeedInfo& rhs);
    	~SeedInfo();


    	int GetNumTracks(){ return fTracks.size(); }
    	std::vector<WCSimLikelihoodTrackBase*> GetTracks(){ return fTracks; }
    	WCSimLikelihoodTrackBase* GetTrack(int p){
            if(p < fTracks.size()){
                return fTracks[p]; // Index starts at 0
            }
            else{
                std::cerr << "GetTrack(index) out of range [0..." << fTracks.size()-1 << "]" << std::endl;
                return WCSimLikelihoodTrackFactory::MakeTrack(TrackType::Unknown, -999, -999, -999, -999,
                											  -999, -999, -999, -999);
            }
        }

	private:
    	std::vector<WCSimLikelihoodTrackBase*> fTracks;
    	ClassDef(SeedInfo,1)
};


class TruthInfo : public TObject{

    public:
        TruthInfo(); 
        TruthInfo(int type, int beamPDG, float beamEnergy, int nPrimaries, std::vector<int> primaryPDGs, 
                  std::vector<double> primaryEnergies, std::vector<TVector3> primaryDirs);
        TruthInfo(const TruthInfo& other);
        TruthInfo& operator=(const TruthInfo& rhs);
        ~TruthInfo();

        void SetVtxTime(float t);
        void SetVtx(float x, float y, float z);
        void SetBeamDir(float x, float y, float z);

        int GetType(){ return fType; }
        float GetBeamE(){ return fBeamEnergy; }
        int GetBeamPDG(){ return fBeamPDG; }

        int GetPrimaryPDG(int p){ 
            if(p < fNPrimaries){
                return fPrimaryPDGs[p]; // Index starts at 0   
            }
            else{ 
                std::cerr << "GetPrimaryPDG(index) out of range [0..." << fNPrimaries-1 << "]" << std::endl;
                return -999;
            }
        } 
        double GetPrimaryEnergy(int p){
            if(p < fNPrimaries){
                return fPrimaryEnergies[p]; // Index starts at 0   
            }
            else{ 
                std::cerr << "GetPrimaryEnergies(index) out of range [0..." << fNPrimaries-1 << "]" << std::endl;
                return -999;
            }
        } // Index starts at 0
        TVector3 GetPrimaryDir(int p){ 
            if(p < fNPrimaries){
                return fPrimaryDirs[p]; // Index starts at 0   
            }
            else{ 
                std::cerr << "GetPrimaryDir(index) out of range [0..." << fNPrimaries-1 << "]" << std::endl;
                return TVector3(-999, -999, -999);
            }
        } // Index starts at 0 

        int GetNPrimaries(){ return fNPrimaries; }
        std::vector<int> GetPrimaryPDGs(){ return fPrimaryPDGs; }
        std::vector<double> GetPrimaryEnergies(){ return fPrimaryEnergies; }  
        std::vector<TVector3> GetPrimaryDirs(){ return fPrimaryDirs; }    

        bool IsCC(){ return fIsCC; }
        bool IsNC(){ return fIsNC; }
        bool IsQE(){ return fIsQE; }
        bool IsRes(){ return fIsRes; }
        bool IsDIS(){ return fIsDIS; }
        bool IsCoherent(){ return fIsCoherent; }
        bool IsNueElectronElastic(){ return fIsNueElectronElastic; }
        bool IsInverseMuonDecay(){ return fIsInverseMuonDecay; }
        bool IsOther(){ return fIsOther; }

        float GetVtxTime() const{ return fVtxTime; }
        float GetVtxX() const{ return fVtxX; }
        float GetVtxY() const{ return fVtxY; }
        float GetVtxZ() const{ return fVtxZ; }

        float GetBeamDirX() const{ return fBeamDirX; }
        float GetBeamDirY() const{ return fBeamDirY; }
        float GetBeamDirZ() const{ return fBeamDirZ; }

    private:

        int fType;  // Interaction type code from WCSimTruthSummary
        int fBeamPDG;  // PDG code of the beam particle (usually a neutrino)
        float fBeamEnergy;  // Energy of the incident beam particle

        int fNPrimaries;
        std::vector<int> fPrimaryPDGs; //Vector of primary PDG's
        std::vector<double> fPrimaryEnergies; //Vector of primary energies
        std::vector<TVector3> fPrimaryDirs; //Vector of primary directions <TVector3>

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

        float fVtxTime;
        float fVtxX;
        float fVtxY;
        float fVtxZ;
        float fBeamDirX;
        float fBeamDirY;
        float fBeamDirZ;
        
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

        void SetEnergy(float E);
        void SetVtxTime(float t);
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

        float  GetFracNHitsInRing(){ return fFracNHitsInRing; }
        float  GetFracNHitsOutsideRing(){ return fFracNHitsOutsideRing; }
        float  GetFracNHitsInRingHole(){ return fFracNHitsInRingHole; }

        float GetFracPredQInRing() const{ return fFracPredQInRing; }
        float GetFracPredQOutsideRing() const{ return fFracPredQOutsideRing; }
        float GetFracPredQInRingHole() const{ return fFracPredQInRingHole; }

        float GetPredictedOverTotalCharge() const{ return fPredictedOverTotalCharge; }

        float GetEnergy() const{ return fEnergy; }
        float GetVtxTime() const{ return fVtxTime; }
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

        float fFracNHitsInRing;
        float fFracNHitsOutsideRing;
        float fFracNHitsInRingHole;

        float fFracPredQInRing;
        float fFracPredQOutsideRing;
        float fFracPredQInRingHole;

        float fPredictedOverTotalCharge;

        float fEnergy;
        float fVtxTime;
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

	void SetSeed(std::vector<WCSimLikelihoodTrackBase*> tracks);
	void SetSeed(WCSimLikelihoodTrackBase* track);

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
    SeedInfo * fSeedInfo;
    std::string fRecoType;
    int fEvent;
    bool fFailed;

    ClassDef(WCSimOutputTree, 2)
};

#endif /* WCSIMOUTPUTTREE_HH */
