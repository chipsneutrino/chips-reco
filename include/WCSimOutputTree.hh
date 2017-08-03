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
#include "WCSimRecoRing.hh"

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

	//As likelihoods can not be compared between different types of seed, ntrack = nseeds.

	public:
		SeedInfo();
		SeedInfo(std::string type, std::vector<WCSimLikelihoodTrackBase*> tracks, int slices, std::vector<WCSimRecoRing*> ringVec, std::vector<double> ringTime);
		SeedInfo(const SeedInfo& other);
		SeedInfo& operator=(const SeedInfo& rhs);
    	~SeedInfo();

    	std::string GetSeedType(){ return fType; }

    	int GetNumTracks(){ return fTracks.size(); }
    	std::vector<WCSimLikelihoodTrackBase*> GetTracks(){ return fTracks; }
    	WCSimLikelihoodTrackBase* GetTrack(int p){
            if(p < fTracks.size()){
                return fTracks[p]; // Index starts at 0
            }
            else{
                std::cerr << "GetTrack(index) out of range [0..." << fTracks.size()-1 << "]" << std::endl;
                return NULL;
                //return WCSimLikelihoodTrackFactory::MakeTrack(TrackType::Unknown, -999, -999, -999, -999, -999, -999, -999, -999);
            }
        }

    	int GetNumSlices(){ return fNSlices; }
    	int GetNumRings(){
    		if( fRingHeight.size() == fRingTime.size()){
    			return fRingHeight.size();
    		}
    		else{
    			std::cout << "fRingVec and fRingTime are not the same size!!!" << std::endl;
    			return -999;
    		}
    	}

    	std::vector<double> GetRingHeights(){ return fRingHeight; }
    	double GetRingHeight(int p){
            if(p < fRingHeight.size()){ return fRingHeight[p]; }
            else{
                std::cerr << "GetRing(index) out of range [0..." << fRingHeight.size()-1 << "]" << std::endl;
                return -999;
            }
        }

    	std::vector<double> GetRingTimes(){ return fRingTime; }
    	double GetRingTime(int p){
            if(p < fRingTime.size()){ return fRingTime[p]; }
            else{
                std::cerr << "GetRingTime(index) out of range [0..." << fRingTime.size()-1 << "]" << std::endl;
                return -999;
            }
        }

    	std::vector<double> GetRingAngles(){ return fRingAngle; }
    	double GetRingAngle(int p){
            if(p < fRingAngle.size()){ return fRingAngle[p]; }
            else{
                std::cerr << "GetRingAngle(index) out of range [0..." << fRingAngle.size()-1 << "]" << std::endl;
                return -999;
            }
        }

    	std::vector<TVector3> GetRingVtxs(){ return fRingVtx; }
    	TVector3 GetRingVtx(int p){
            if(p < fRingVtx.size()){ return fRingVtx[p]; }
            else{
                std::cerr << "GetRingVtx(index) out of range [0..." << fRingVtx.size()-1 << "]" << std::endl;
                return TVector3(-999, -999, -999);
            }
        }

    	std::vector<TVector3> GetRingDirs(){ return fRingDir; }
    	TVector3 GetRingDir(int p){
            if(p < fRingDir.size()){ return fRingDir[p]; }
            else{
                std::cerr << "GetRingDir(index) out of range [0..." << fRingDir.size()-1 << "]" << std::endl;
                return TVector3(-999, -999, -999);
            }
        }

	private:
    	std::string fType; // Records the type of seed used to produce the passed on tracks.
    	std::vector<WCSimLikelihoodTrackBase*> fTracks; // Stores a vector of the tracks actually passed to the fitter.

    	int fNSlices; // Number of slices produced by the seeder.
    	std::vector<double> fRingHeight; // Stores the ring heights produces by the seeder.
    	std::vector<double> fRingTime; // Stores the vertex times of the rings. These are the same for all rings belonging to the same slice.
    	std::vector<double> fRingAngle; // Stores the ring angles
    	std::vector<TVector3> fRingVtx; // Stores TVector3 objects containing the ring vtx positions
    	std::vector<TVector3> fRingDir; // Stores Tvector3 objects containing the ring dirs.

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

	void SetSeed(std::string type, std::vector<WCSimLikelihoodTrackBase*> tracks, int slices, std::vector<WCSimRecoRing*> ringVec, std::vector<double> ringTime);

    void Fill(
              bool failed,
              Int_t iEvent,
              std::string recoType,
              WCSimRecoSummary& summ,
              //WCSimHitComparison& hitComp,
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
    //WCSimHitComparison * fHitComparison;
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
