/*
 * WCSimOutputTree.cc
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#include "WCSimAnalysisConfig.hh"
#include "WCSimGeometry.hh"
#include "WCSimHitComparison.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimOutputTree.hh"
#include "WCSimRecoSummary.hh"
#include "WCSimTrueEvent.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TMD5.h"
#include "TString.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>

#ifndef REFLEX_DICTIONARY
ClassImp(EventHeader);
ClassImp(HitInfo);
ClassImp(StageInfo);
ClassImp(SeedInfo);
ClassImp(TruthInfo);
ClassImp(RecoInfo);
ClassImp(WCSimOutputTree);
#endif

///////////////////////////////////////////////////////////////
//  METHODS FOR EVENTHEADER CLASS                            //
///////////////////////////////////////////////////////////////

EventHeader::EventHeader() :
	fUID(""),
	fInputFile(""),
	fInputEvent(-999),
	fFailed(false),
	fRecoType("")
{
	// Empty
}

EventHeader::EventHeader(std::string inputFile,
	    				 int inputEvent,
	    				 bool failed,
	    				 std::string recoType
						 ) :
	fInputFile(inputFile),
	fInputEvent(inputEvent),
	fFailed(failed),
	fRecoType(recoType)
{
	BuildUID();
}

EventHeader::EventHeader(const EventHeader& other) :
	fUID(other.fUID),
	fInputFile(other.fInputFile),
	fInputEvent(other.fInputEvent),
	fFailed(other.fFailed),
	fRecoType(other.fRecoType)
{
	// Empty
}

EventHeader& EventHeader::operator=(const EventHeader& rhs)
{
    if(&rhs != this)
    {
    	fUID = rhs.fUID;
    	fInputFile = rhs.fInputFile;
    	fInputEvent = rhs.fInputEvent;
    	fFailed = rhs.fFailed;
    	fRecoType = rhs.fRecoType;
    }
    return *this;
}

EventHeader::~EventHeader()
{
	// Empty
}

void EventHeader::Print()
{
    std::cout << "EventHeader::Print()..."       << std::endl
    		  << "UID -> " << fUID               << std::endl
    		  << "InputFile -> " << fInputFile   << std::endl
    		  << "InputEvent -> " << fInputEvent << std::endl
    		  << "Failed -> " << fFailed         << std::endl
    		  << "RecoType -> " << fRecoType     << std::endl;
    return;
}

/**
 * @brief Make unique identifier for a file/event combination.  Sets fUID to
 * md5("%s%d", md5(inputFile), eventNumber)
 *
**/
void EventHeader::BuildUID()
{
    TMD5 myMD5;

    // Build a char array to hold our filename string:
    TString inputStr = TString::Format("%s%06d", fInputFile.c_str(), fInputEvent);
    std::string input(inputStr.Data());
    const size_t sizeOfString(input.size());
    UChar_t charArray[sizeOfString];
    strcpy(reinterpret_cast<char*>(charArray), input.c_str());

    // Hash the filename
    UChar_t output[16];
    myMD5.Update(charArray, sizeOfString);
    myMD5.Final(output);
    fUID = myMD5.AsString();
    std::cout << "Hash of " << input.c_str() << " is " << fUID << std::endl;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR HITINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

HitInfo::HitInfo() : 
	fVeto(false),
	fNHits(0),
	fNHitsUpstream(0),
	fNHitsDownstream(0),
	fFracHitsUpstream(0.0),
	fFracHitsDownstream(0.0),
	fTotalQ(0.0),
	fTotalQUpstream(0.0),
	fTotalQDownstream(0.0),
	fFracQUpstream(0.0),
	fFracQDownstream(0.0)
{
	// Empty
}

HitInfo::HitInfo(bool veto,
				 int NHits,
				 int NHitsUpstream,
				 float totalQ,
				 float totalQUpstream
        		 ) :
	fVeto(veto),
	fNHits(NHits),
	fNHitsUpstream(NHitsUpstream),
	fTotalQ(totalQ),
	fTotalQUpstream(totalQUpstream)
{
    fNHitsDownstream = fNHits - fNHitsUpstream;
    fFracHitsUpstream = fNHitsUpstream/(static_cast<float>(fNHits));
    fFracHitsDownstream = 1 - fFracHitsUpstream;

    fTotalQDownstream = fTotalQ - fTotalQUpstream;
    fFracQUpstream = fTotalQUpstream / fTotalQ;
    fFracQDownstream = 1 - fFracQUpstream;
}

HitInfo::HitInfo(const HitInfo& other) : 
	fVeto(other.fVeto),
    fNHits(other.fNHits),
    fNHitsUpstream(other.fNHitsUpstream),
    fNHitsDownstream(other.fNHitsDownstream),
    fFracHitsUpstream(other.fFracHitsUpstream),
    fFracHitsDownstream(other.fFracHitsDownstream),
    fTotalQ(other.fTotalQ),
    fTotalQUpstream(other.fTotalQUpstream),
    fTotalQDownstream(other.fTotalQDownstream),
    fFracQUpstream(other.fFracQUpstream),
    fFracQDownstream(other.fFracQDownstream)
{
	// Empty
}

HitInfo& HitInfo::operator=(const HitInfo& rhs)
{
    if(&rhs != this)
    {
        fVeto = rhs.fVeto;
        fNHits = rhs.fNHits;
        fNHitsUpstream = rhs.fNHitsUpstream;
        fNHitsDownstream = rhs.fNHitsDownstream;
        fFracHitsUpstream = rhs.fFracHitsUpstream;
        fFracHitsDownstream = rhs.fFracHitsDownstream;
        fTotalQ = rhs.fTotalQ;
        fTotalQUpstream = rhs.fTotalQUpstream;
        fTotalQDownstream = rhs.fTotalQDownstream;
        fFracQUpstream = rhs.fFracQUpstream;
        fFracQDownstream = rhs.fFracQDownstream;
    }
    return *this;
}

HitInfo::~HitInfo()
{
    // Empty
}

void HitInfo::Print()
{
    std::cout << "HitInfo::Print()..."                           << std::endl
    		  << "Veto -> " << fVeto                             << std::endl
    		  << "NHits -> " << fNHits                           << std::endl
    		  << "NHitsUpstream -> " << fNHitsUpstream           << std::endl
    		  << "NHitsDownstream -> " << fNHitsDownstream       << std::endl
    		  << "FracHitsUpstream -> " << fFracHitsUpstream     << std::endl
    		  << "FracHitsDownstream -> " << fFracHitsDownstream << std::endl
    		  << "TotalQ -> " << fTotalQ                         << std::endl
    		  << "TotalQUpstream -> " << fTotalQUpstream         << std::endl
    		  << "TotalQDownstream -> " << fTotalQDownstream     << std::endl
    		  << "FracQUpstream -> " << fFracQUpstream           << std::endl
    		  << "FracQDownstream -> " << fFracQDownstream       << std::endl;
    return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR STAGEINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

StageInfo::StageInfo()
{
	fStageNcalls.clear();
	fStagePreds.clear();
	fStageTracks.clear();
}

StageInfo::StageInfo(const StageInfo& other) :
	fStageNcalls(other.fStageNcalls),
	fStagePreds(other.fStagePreds),
	fStageTracks(other.fStageTracks)
{
	// Empty
}

StageInfo& StageInfo::operator=(const StageInfo& rhs)
{
    if(this != &rhs)
    {
    	fStageNcalls = rhs.fStageNcalls;
    	fStagePreds = rhs.fStagePreds;
    	fStageTracks = rhs.fStageTracks;
    }
    return *this;
}

StageInfo::~StageInfo()
{
	// Empty
}

void StageInfo::Print()
{
    std::cout << "StageInfo::Print()..."       << std::endl;
	std::cout << "Number of Stages -> " << fStageNcalls.size() << std::endl;
	std::cout << "Number of Tracks -> " << fStageTracks[0].size() << std::endl << std::endl;

	for(int stage = 0; (size_t)stage<fStageNcalls.size(); stage++){
		if(stage==0){ std::cout << "(seed)..." << std::endl; }
		else{ std::cout << "(stage " << stage << ")..." << std::endl; }

		std::cout << "NCalls -> " << fStageNcalls[stage] << std::endl;

		// Print the track parameters at this stage...
		for(int t=0; (size_t)t<fStageTracks[stage].size(); t++ ){
			WCSimLikelihoodTrackBase* track = fStageTracks[stage][t];
			std::cout << "Track " << t << ", vtx = (" << track->GetX() << "," << track->GetY() << "," << track->GetZ() << "), dir = (" << track->GetTheta() << "," << track->GetPhi() << "), Energy = " << track->GetE() << std::endl;
		}

		// Print the likelihoods and charge at this stage...
		double totPredQ = 0;
		double totTotal2LnL = 0,  totTime2LnL = 0, totCharge2LnL = 0;
		for(int p=0; (size_t)p<fStagePreds[stage].size(); p++){
			totPredQ += fStagePreds[stage][p].GetPredictedCharge();
			totTotal2LnL += fStagePreds[stage][p].GetTotal2LnL();
			totTime2LnL += fStagePreds[stage][p].GetTime2LnL();
			totCharge2LnL += fStagePreds[stage][p].GetCharge2LnL();
		}
		std::cout << "NCalls -> " << fStageNcalls[stage] << ", TotalPredQ -> " << totPredQ << std::endl;
		std::cout << "TotTotal2LnL -> " << totTotal2LnL << "TotTime2LnL -> " << totTime2LnL << "totCharge2LnL -> " << totCharge2LnL << std::endl << std::endl;
	}
	return;
}

void StageInfo::SetStageInfo(int stageNCalls,
				  std::vector< WCSimHitPrediction > stagePreds,
				  std::vector<WCSimLikelihoodTrackBase*> stageTracks)
{
	fStageNcalls.push_back(stageNCalls);
	fStagePreds.push_back(stagePreds);
	fStageTracks.push_back(stageTracks);
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR SEEDINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

SeedInfo::SeedInfo() :
	fNSlices(-999)
{
	fTracks.clear();
	fRingHeight.clear();
	fRingTime.clear();
	fRingAngle.clear();
	fRingVtx.clear();
	fRingDir.clear();
}

SeedInfo::SeedInfo(std::vector<WCSimLikelihoodTrackBase*> tracks,
				   int slices,
				   std::vector<WCSimRecoRing*> ringVec,
				   std::vector<double> ringTime
				   ) :
	fTracks(tracks),
	fNSlices(slices),
	fRingTime(ringTime)
{
	fRingHeight.clear();
	fRingAngle.clear();
	fRingVtx.clear();
	fRingDir.clear();
	for(int r=0; (size_t)r<ringVec.size(); r++){
		fRingHeight.push_back(ringVec[r]->GetHeight());
    	fRingAngle.push_back(ringVec[r]->GetAngle());
    	fRingVtx.push_back(TVector3(ringVec[r]->GetVtxX(), ringVec[r]->GetVtxY(), ringVec[r]->GetVtxZ()));
    	fRingDir.push_back(TVector3(ringVec[r]->GetDirX(), ringVec[r]->GetDirY(), ringVec[r]->GetDirZ()));
	}
}

SeedInfo::SeedInfo(const SeedInfo& other):
	fTracks(other.fTracks),
	fNSlices(other.fNSlices),
	fRingHeight(other.fRingHeight),
	fRingTime(other.fRingTime),
	fRingAngle(other.fRingAngle),
	fRingVtx(other.fRingVtx),
	fRingDir(other.fRingDir)
{
	// Empty...
}

SeedInfo& SeedInfo::operator=(const SeedInfo& rhs)
{
    if(this != &rhs)
    {
    	fTracks = rhs.fTracks;
    	fNSlices = rhs.fNSlices;
    	fRingHeight = rhs.fRingHeight;
    	fRingTime = rhs.fRingTime;
    	fRingAngle = rhs.fRingAngle;
    	fRingVtx = rhs.fRingVtx;
    	fRingDir = rhs.fRingDir;
    }
    return *this;
}

SeedInfo::~SeedInfo()
{
	// Empty...
}

void SeedInfo::Print()
{
    std::cout << "SeedInfo::Print()..."       << std::endl;
	std::cout << "NSlices -> " << fNSlices << std::endl;
	std::cout << "NRings -> " << fRingTime.size() << std::endl;
	std::cout << "NTracks -> " << fTracks.size() << std::endl;

	// Print out the rings...
	for(int r=0; (size_t)r<fRingTime.size(); r++){
		std::cout << "Ring " << r << ": Height->" << fRingHeight[r] << ", Time->" << fRingTime[r] << ", Angle->" << fRingAngle[r] << std::endl;
		std::cout << "VtxPos->(" << fRingVtx[r].X() << "," << fRingVtx[r].Y() << "," << fRingVtx[r].Z() << ") Dir->(" << fRingDir[r].X() << "," << fRingDir[r].Y() << "," << fRingDir[r].Z() << ")" << std::endl;
	}

	// Print out the track given to the fitter...
	for(int t=0; (size_t)t<fTracks.size(); t++){
		WCSimLikelihoodTrackBase* track = fTracks[t];
		std::cout << "Track " << t << ", vtx = (" << track->GetX() << "," << track->GetY() << "," << track->GetZ() << "), dir = (" << track->GetTheta() << "," << track->GetPhi() << "), Energy = " << track->GetE() << std::endl;
	}
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR TRUTHINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

TruthInfo::TruthInfo() :
	fType(-999),
	fBeamPDG(-999),
	fBeamEnergy(0.0),
	fNPrimaries(-999),
	fIsCC(false),
	fIsNC(false),
	fIsRes(false),
	fIsDIS(false),
	fIsCoherent(false),
	fIsNueElectronElastic(false),
	fIsInverseMuonDecay(false),
	fIsOther(false),
	fVtxTime(0.0),
	fVtxX(0.0),
	fVtxY(0.0),
	fVtxZ(0.0),
	fBeamDirX(0.0),
	fBeamDirY(0.0),
	fBeamDirZ(0.0)
{
    fPrimaryPDGs.clear();
    fPrimaryEnergies.clear();
    fPrimaryDirs.clear(); 
}

TruthInfo::TruthInfo(int type,
					 int beamPDG,
					 float beamEnergy,
					 int nPrimaries,
					 std::vector<int> primaryPDGs,
                     std::vector<double> primaryEnergies,
					 std::vector<TVector3> primaryDirs
					 ) :
	fType(type),
	fBeamPDG(beamPDG),
	fBeamEnergy(beamEnergy),
	fNPrimaries(nPrimaries),
    fPrimaryPDGs(primaryPDGs),
	fPrimaryEnergies(primaryEnergies),
	fPrimaryDirs(primaryDirs)
{
    fIsCC = WCSimTruthSummary::TypeIsCCEvent(fType);
    fIsNC = WCSimTruthSummary::TypeIsNCEvent(fType);
    fIsQE = WCSimTruthSummary::TypeIsQEEvent(fType);
    fIsRes = WCSimTruthSummary::TypeIsResEvent(fType);
    fIsDIS = WCSimTruthSummary::TypeIsDISEvent(fType);
    fIsCoherent = WCSimTruthSummary::TypeIsCohEvent(fType);
    fIsNueElectronElastic = WCSimTruthSummary::TypeIsNueElectronElasticEvent(type);
    fIsInverseMuonDecay = WCSimTruthSummary::TypeIsInverseMuonDecayEvent(type);
    fIsOther = fIsCC && !(fIsQE || fIsRes || fIsDIS || fIsCoherent 
                          || fIsNueElectronElastic || fIsInverseMuonDecay);
    fVtxTime = 0.0;
    fVtxX = 0.0; 
    fVtxY = 0.0; 
    fVtxZ = 0.0;
    fBeamDirX = 0.0; 
    fBeamDirY = 0.0; 
    fBeamDirZ = 0.0;   
}

TruthInfo::TruthInfo(const TruthInfo& other):
    fType(other.fType),
    fBeamPDG(other.fBeamPDG),
    fBeamEnergy(other.fBeamEnergy),
    fNPrimaries(other.fNPrimaries),
    fPrimaryPDGs(other.fPrimaryPDGs),
    fPrimaryEnergies(other.fPrimaryEnergies),
    fPrimaryDirs(other.fPrimaryDirs),
    fIsCC(other.fIsCC),
    fIsNC(other.fIsNC),
    fIsQE(other.fIsQE),
    fIsRes(other.fIsRes),
    fIsDIS(other.fIsDIS),
    fIsCoherent(other.fIsCoherent),
    fIsNueElectronElastic(other.fIsNueElectronElastic),
    fIsInverseMuonDecay(other.fIsInverseMuonDecay),
    fIsOther(other.fIsOther),
    fVtxTime(other.fVtxTime),
    fVtxX(other.fVtxX),
    fVtxY(other.fVtxY),
    fVtxZ(other.fVtxZ),
    fBeamDirX(other.fBeamDirX),
    fBeamDirY(other.fBeamDirY),
    fBeamDirZ(other.fBeamDirZ)
{
    // Empty
}

TruthInfo& TruthInfo::operator=(const TruthInfo& rhs)
{
    if(this != &rhs)
    {
        fType = rhs.fType;
        fBeamPDG = rhs.fBeamPDG;
        fBeamEnergy = rhs.fBeamEnergy;
        fNPrimaries = rhs.fNPrimaries;
        fPrimaryPDGs = rhs.fPrimaryPDGs;
        fPrimaryEnergies = rhs.fPrimaryEnergies;
        fPrimaryDirs = rhs.fPrimaryDirs;
        fIsCC = rhs.fIsCC;
        fIsNC = rhs.fIsNC;
        fIsQE = rhs.fIsQE;
        fIsRes = rhs.fIsRes;
        fIsDIS = rhs.fIsDIS;
        fIsCoherent = rhs.fIsCoherent;
        fIsNueElectronElastic = rhs.fIsNueElectronElastic;
        fIsInverseMuonDecay = rhs.fIsInverseMuonDecay;
        fIsOther = rhs.fIsOther;
        fVtxTime = rhs.fVtxTime;
        fVtxX = rhs.fVtxX;
        fVtxY = rhs.fVtxY;
        fVtxZ = rhs.fVtxZ;
        fBeamDirX = rhs.fBeamDirX;
        fBeamDirY = rhs.fBeamDirY;
        fBeamDirZ = rhs.fBeamDirZ;
    }
    return *this;
}

TruthInfo::~TruthInfo()
{
    // Empty
}

void TruthInfo::Print()
{
    std::cout << "TruthInfo::Print()..."       << std::endl;
    std::cout << "Beam: PDG->" << fBeamPDG << ", Energy->" << fBeamEnergy << ", Dir->(" << fBeamDirX << "," << fBeamDirY << "," << fBeamDirZ << ")" << std::endl;
    std::cout << "Vertex: Position->(" << fVtxX << "," << fVtxY << "," << fVtxZ << "), Time->" << fVtxTime << std::endl;

    std::cout << "Interaction type code -> " << fType << std::endl;
    std::cout << "Interaction is";
    if(fIsCC){ std::cout << ", Charged-current"; }
    if(fIsNC){ std::cout << ", Neutral-current"; }
    if(fIsQE){ std::cout << ", Quasielastic"; }
    if(fIsRes){ std::cout << ", Resonant"; }
    if(fIsDIS){ std::cout << ", Deep Inelastic Scattering"; }
    if(fIsCoherent){ std::cout << ", Coherent"; }
    if(fIsNueElectronElastic){ std::cout << ", Electron elastic"; }
    if(fIsInverseMuonDecay){ std::cout << ", Inverse muon decay"; }
    if(fIsOther){ std::cout << ", Other"; }
    std::cout << std::endl;

    // Print info about the primaries...
    for(int p=0; p<fNPrimaries; p++){
    	std::cout << "Primary " << p << ": PDG->" << fPrimaryPDGs[p] << ", Energy->" << fPrimaryEnergies[p] << ", Dir->(" << fPrimaryDirs[p].X() << "," << fPrimaryDirs[p].Y() << "," << fPrimaryDirs[p].Z() << ")" << std::endl;
    }
	return;
}

void TruthInfo::SetVtxTime(float t)
{
	fVtxTime = t;
	return;
}

void TruthInfo::SetVtx(float x,
					   float y,
					   float z)
{
	fVtxX = x;
	fVtxY = y;
	fVtxZ = z;
	return;
}

void TruthInfo::SetBeamDir(float x,
						   float y,
						   float z)
{
	fBeamDirX = x;
	fBeamDirY = y;
	fBeamDirZ = z;
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR RECOINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

RecoInfo::RecoInfo() : 
    fTotal2LnL(0.0), 
    fHit2LnL(0.0), fCharge2LnL(0.0), fTime2LnL(0.0), fCutoff2LnL(0.0),
    fTotalQInRing(0.0), fTotalQOutsideRing(0.0), fTotalQInRingHole(0.0),
    fNHitsInRing(0), fNHitsOutsideRing(0), fNHitsInRingHole(0),
    fPredQInRing(0.0), fPredQOutsideRing(0.0), fPredQInRingHole(0.0),
    fFracTotalQInRing(0.0), fFracTotalQOutsideRing(0.0), fFracTotalQInRingHole(0.0),
    fFracNHitsInRing(0.0), fFracNHitsOutsideRing(0.0), fFracNHitsInRingHole(0.0),
    fFracPredQInRing(0.0), fFracPredQOutsideRing(0.0), fFracPredQInRingHole(0.0),
    fPredictedOverTotalCharge(0.0),
    fVtxTime(0.0),
	fEnergy(0.0),
    fVtxX(0.0), fVtxY(0.0), fVtxZ(0.0), fVtxRho(0.0),
    fEndX(0.0), fEndY(0.0), fEndZ(0.0), fEndRho(0.0),
    fDirX(0.0), fDirY(0.0), fDirZ(0.0),
    fEscapes(false)
{
    // Empty
}

RecoInfo::RecoInfo(const RecoInfo& other) : 
    fTotal2LnL(other.fTotal2LnL),
    fHit2LnL(other.fHit2LnL),
    fCharge2LnL(other.fCharge2LnL),
    fTime2LnL(other.fTime2LnL),
    fCutoff2LnL(other.fCutoff2LnL),
    fTotalQInRing(other.fTotalQInRing),
    fTotalQOutsideRing(other.fTotalQOutsideRing),
    fTotalQInRingHole(other.fTotalQInRingHole),
    fNHitsInRing(other.fNHitsInRing),
    fNHitsOutsideRing(other.fNHitsOutsideRing),
    fNHitsInRingHole(other.fNHitsInRingHole),
    fPredQInRing(other.fPredQInRing),
    fPredQOutsideRing(other.fPredQOutsideRing),
    fPredQInRingHole(other.fPredQInRingHole),
    fFracTotalQInRing(other.fFracTotalQInRing),
    fFracTotalQOutsideRing(other.fFracTotalQOutsideRing),
    fFracTotalQInRingHole(other.fFracTotalQInRingHole),
    fFracNHitsInRing(other.fFracNHitsInRing),
    fFracNHitsOutsideRing(other.fFracNHitsOutsideRing),
    fFracNHitsInRingHole(other.fFracNHitsInRingHole),
    fFracPredQInRing(other.fFracPredQInRing),
    fFracPredQOutsideRing(other.fFracPredQOutsideRing),
    fFracPredQInRingHole(other.fFracPredQInRingHole),
    fPredictedOverTotalCharge(other.fPredictedOverTotalCharge),
    fVtxTime(other.fVtxTime),
	fEnergy(other.fEnergy),
    fVtxX(other.fVtxX),
    fVtxY(other.fVtxY),
    fVtxZ(other.fVtxZ),
    fVtxRho(other.fVtxRho),
    fEndX(other.fEndX),
    fEndY(other.fEndY),
    fEndZ(other.fEndZ),
    fEndRho(other.fEndRho),
    fDirX(other.fDirX),
    fDirY(other.fDirY),
    fDirZ(other.fDirZ),
    fEscapes(other.fEscapes)
{
    // Empty
}

RecoInfo& RecoInfo::operator=(const RecoInfo& other)
{
    if(this == &other)
    {
        fTotal2LnL = other.fTotal2LnL;
        fHit2LnL = other.fHit2LnL;
        fCharge2LnL = other.fCharge2LnL;
        fTime2LnL = other.fTime2LnL;
        fCutoff2LnL = other.fCutoff2LnL;

        fTotalQInRing = other.fTotalQInRing;
        fTotalQOutsideRing = other.fTotalQOutsideRing;
        fTotalQInRingHole = other.fTotalQInRingHole;

        fNHitsInRing = other.fNHitsInRing;
        fNHitsOutsideRing = other.fNHitsOutsideRing;
        fNHitsInRingHole = other.fNHitsInRingHole;

        fPredQInRing = other.fPredQInRing;
        fPredQOutsideRing = other.fPredQOutsideRing;
        fPredQInRingHole = other.fPredQInRingHole;

        fFracTotalQInRing = other.fFracTotalQInRing;
        fFracTotalQOutsideRing = other.fFracTotalQOutsideRing;
        fFracTotalQInRingHole = other.fFracTotalQInRingHole;

        fFracNHitsInRing = other.fFracNHitsInRing;
        fFracNHitsOutsideRing = other.fFracNHitsOutsideRing;
        fFracNHitsInRingHole = other.fFracNHitsInRingHole;

        fFracPredQInRing = other.fFracPredQInRing;
        fFracPredQOutsideRing = other.fFracPredQOutsideRing;
        fFracPredQInRingHole = other.fFracPredQInRingHole;

        fPredictedOverTotalCharge = other.fPredictedOverTotalCharge;
	
        fVtxTime = other.fVtxTime;
        fEnergy = other.fEnergy;
        fVtxX = other.fVtxX;
        fVtxY = other.fVtxY;
        fVtxZ = other.fVtxZ;
        fVtxRho = other.fVtxRho;

        fEndX = other.fEndX;
        fEndY = other.fEndY;
        fEndZ = other.fEndZ;
        fEndRho = other.fEndRho;

        fDirX = other.fDirX;
        fDirY = other.fDirY;
        fDirZ = other.fDirZ;

        fEscapes = other.fEscapes;
    }
    return *this;
}

RecoInfo::~RecoInfo()
{
	// Empty
}

void RecoInfo::Print()
{
    std::cout << "RecoInfo::Print()..."                           << std::endl
    		  << "fTotal2LnL" << fTotal2LnL << std::endl
    		  << "fHit2LnL" << fHit2LnL << std::endl
    		  << "fCharge2LnL" << fCharge2LnL << std::endl
    		  << "fTime2LnL" << fTime2LnL << std::endl
    		  << "fCutoff2LnL" << fCutoff2LnL << std::endl
    		  << "fTotalQInRing" << fTotalQInRing << std::endl
    		  << "fTotalQOutsideRing" << fTotalQOutsideRing << std::endl
    		  << "fTotalQInRingHole" << fTotalQInRingHole << std::endl
    		  << "fNHitsInRing" << fNHitsInRing << std::endl
    		  << "fNHitsOutsideRing" << fNHitsOutsideRing << std::endl
    		  << "fNHitsInRingHole" << fNHitsInRingHole << std::endl
    		  << "fPredQInRing" << fPredQInRing << std::endl
    		  << "fPredQOutsideRing" << fPredQOutsideRing << std::endl
    		  << "fPredQInRingHole" << fPredQInRingHole << std::endl
    		  << "fFracTotalQInRing" << fFracTotalQInRing << std::endl
    		  << "fFracTotalQOutsideRing" << fFracTotalQOutsideRing << std::endl
    		  << "fFracTotalQInRingHole" << fFracTotalQInRingHole << std::endl
    		  << "fFracNHitsInRing" << fFracNHitsInRing << std::endl
    		  << "fFracNHitsOutsideRing" << fFracNHitsOutsideRing << std::endl
    		  << "fFracNHitsInRingHole" << fFracNHitsInRingHole << std::endl
    		  << "fFracPredQInRing" << fFracPredQInRing << std::endl
    		  << "fFracPredQOutsideRing" << fFracPredQOutsideRing << std::endl
    		  << "fFracPredQInRingHole" << fFracPredQInRingHole << std::endl
    		  << "fPredictedOverTotalCharge" << fPredictedOverTotalCharge << std::endl
    		  << "fVtxTime" << fVtxTime << std::endl
    		  << "fEnergy" << fEnergy << std::endl
    		  << "fVtxX" << fVtxX << std::endl
    		  << "fVtxY" << fVtxY << std::endl
    		  << "fVtxZ" << fVtxZ << std::endl
    		  << "fVtxRho" << fVtxRho << std::endl
    		  << "fEndX" << fEndX << std::endl
    		  << "fEndY" << fEndY << std::endl
    		  << "fEndZ" << fEndZ << std::endl
    		  << "fEndRho" << fEndRho << std::endl
    		  << "fDirX" << fDirX << std::endl
    		  << "fDirY" << fDirY << std::endl
    		  << "fDirZ" << fDirZ << std::endl
    		  << "fEscapes" << fEscapes << std::endl;
    return;
}

void RecoInfo::SetLikelihoods(float total,
							  float hit,
							  float charge,
							  float time)
{
    fTotal2LnL = total;
    fHit2LnL = hit;
    fCharge2LnL = charge;
    fTime2LnL = time;
    fCutoff2LnL = (hit + charge + time) - total;
	return;
}

void RecoInfo::SetTotalQ(float in,
						 float out,
						 float hole)
{
    fTotalQInRing = in;
    fTotalQOutsideRing = out;
    fTotalQInRingHole = hole;
    float totalQ = in + out + hole;

    if(totalQ > 0)
    {
        fFracTotalQInRing = fTotalQInRing/totalQ;
        fFracTotalQOutsideRing = fTotalQOutsideRing/totalQ;
        fFracTotalQInRingHole = fTotalQInRingHole/totalQ;
    }
    else
    {
        fFracTotalQInRing = 0.0;
        fFracTotalQOutsideRing = 0.0;
        fFracTotalQInRingHole = 0.0;
    }

    if(totalQ > 0)
    {
        fPredictedOverTotalCharge = 
            (fPredQInRing + fPredQOutsideRing + fPredQInRingHole) / totalQ;
    }
    else
    {
        fPredictedOverTotalCharge = 0.0;
    }
	return;
}

void RecoInfo::SetNHits(int in,
						int out,
						int hole)
{
    fNHitsInRing = in;
    fNHitsOutsideRing = out;
    fNHitsInRingHole = hole;
    float nHits = in + out + hole;

    if(nHits > 0)
    {
        fFracNHitsInRing = fNHitsInRing/nHits;
        fFracNHitsOutsideRing = fNHitsOutsideRing/nHits;
        fFracNHitsInRingHole = fNHitsInRingHole/nHits;
    }
    else
    {
        fFracNHitsInRing = 0.0;
        fFracNHitsOutsideRing = 0.0;
        fFracNHitsInRingHole = 0.0;
    }
	return;
}

void RecoInfo::SetPredQ(float in,
						float out,
						float hole)
{
    fPredQInRing = in;
    fPredQOutsideRing = out;
    fPredQInRingHole = hole;
    float predQ = in + out + hole;

    if(predQ > 0)
    {
        fFracPredQInRing = fPredQInRing/predQ;
        fFracPredQOutsideRing = fPredQOutsideRing/predQ;
        fFracPredQInRingHole = fPredQInRingHole/predQ;
    }
    else
    {
        fFracPredQInRing = 0.0;
        fFracPredQOutsideRing = 0.0;
        fFracPredQInRingHole = 0.0;
    }

    float totalQ = fTotalQInRing + fTotalQOutsideRing + fTotalQInRingHole;
    if(totalQ > 0)
    {
        fPredictedOverTotalCharge = predQ/totalQ;
    }
    else
    {
        fPredictedOverTotalCharge = 0.0;
    }
    return;
}


void RecoInfo::SetEnergy(float E)
{
	fEnergy = E;
	return;
}

void RecoInfo::SetVtxTime(float t)
{
	fVtxTime = t;
	return;
}

void RecoInfo::SetVtx(float x,
					  float y,
					  float z)
{
    fVtxX = x;
    fVtxY = y;
    fVtxZ = z;
    fVtxRho = TMath::Sqrt(x*x + y*y);
	return;
}

void RecoInfo::SetEnd(float x,
					  float y,
					  float z)
{
    fEndX = x;
    fEndY = y;
    fEndZ = z;
    fEndRho = TMath::Sqrt(x*x + y*y);
	return;
}

void RecoInfo::SetDir(float x,
					  float y,
					  float z)
{
    fDirX = x;
    fDirY = y;
    fDirZ = z;
	return;
}

void RecoInfo::SetEscapes(bool escapes)
{
    fEscapes = escapes;
	return;
}


///////////////////////////////////////////////////////////////
//  METHODS FOR WCSimOutputTree CLASS                        //
///////////////////////////////////////////////////////////////

WCSimOutputTree::WCSimOutputTree(const TString &saveFileName) : 
	fSaveFile(0x0),
	fSaveFileName(saveFileName),
	fResultsTree(0x0),
	fGeoTree(0x0),
	fGeometry(0x0),
    fEventHeader(0x0),
    fTruthInfo(0x0),
	fRecoSummary(0x0),
    fHitInfo(0x0),
    fRecoInfo(0x0),
    fSeedInfo(0x0),
    fStageInfo(0x0),
    fHitComparison(0x0)
{
	fSaveFileName.Form("%s_tree.root", saveFileName.Data());
    fSaveFile = 0x0;
}

WCSimOutputTree::~WCSimOutputTree()
{
    if(fEventHeader != 0x0)
    {
        std::cout << "delete Event Header" << std::endl;
        delete fEventHeader;
        fEventHeader = 0x0;
        std::cout << "done" << std::endl;
    }

    if(fTruthInfo != 0x0)
    {
        std::cout << "delete truth info" << std::endl;
        delete fTruthInfo;
        fTruthInfo = 0x0;
        std::cout << "done" << std::endl;
    }

    if(fRecoSummary != 0x0)
    {
        std::cout << "delete reco summary" << std::endl;
        delete fRecoSummary;
        std::cout << "done" << std::endl;
    }

    if(fHitInfo != 0x0)
    {
        std::cout << "delete hit info" << std::endl;
        delete fHitInfo;
        fHitInfo = 0x0;
        std::cout << "done" << std::endl;
    }

    if(fRecoInfo != 0x0)
    {
        std::cout << "delete reco info" << std::endl;
        delete fRecoInfo;
        fRecoInfo = 0x0;
        std::cout << "done" << std::endl;
    }

    if(fSeedInfo != 0x0)
    {
        std::cout << "delete seed info" << std::endl;
        delete fSeedInfo;
        fSeedInfo = 0x0;
        std::cout << "done" << std::endl;
    }

    if(fStageInfo != 0x0)
    {
        std::cout << "delete stage info" << std::endl;
        delete fStageInfo;
        fStageInfo = 0x0;
        std::cout << "done" << std::endl;
    }

    if(fHitComparison != 0x0)
    {
        std::cout << "delete hit comp" << std::endl;
        delete fHitComparison;
        fHitComparison = 0x0;
        std::cout << "done" << std::endl;
    }
}

void WCSimOutputTree::MakeTree()
{
    std::cout << "fSaveFileName = " << fSaveFileName << " and file is " << fSaveFile << std::endl;
	fSaveFile->cd();

    DeletePointersIfExist();

    fGeoTree = new TTree("fGeoTree","GeoTree");
    fGeometry = new WCSimRootGeom();
	fGeoTree->Branch("wcsimrootgeom", "WCSimRootGeom", &fGeometry, 64000,0);

    fResultsTree = new TTree("fResultsTree","ResultsTree");

    fEventHeader = new EventHeader();
    fTruthInfo = new TruthInfo();
    fRecoSummary = new WCSimRecoSummary();
    fHitInfo = new HitInfo();
    fRecoInfo = new RecoInfo();
    fSeedInfo = new SeedInfo();
    fStageInfo = new StageInfo();
    fHitComparison = new WCSimHitComparison();

    fResultsTree->Branch("EventHeader", "EventHeader", &fEventHeader, 64000, 1);
    fResultsTree->Branch("TruthInfo", "TruthInfo", &fTruthInfo, 64000, 1);
    fResultsTree->Branch("RecoSummary", "WCSimRecoSummary", &fRecoSummary, 64000, 0);
    fResultsTree->Branch("HitInfo", "HitInfo", &fHitInfo, 64000, 1);
    fResultsTree->Branch("RecoInfo", "RecoInfo", &fRecoInfo, 64000, 1);

    if(WCSimAnalysisConfig::Instance()->GetSaveSeedInfo()){
    	fResultsTree->Branch("SeedInfo", "SeedInfo", &fSeedInfo, 64000, 1);
    }
    if(WCSimAnalysisConfig::Instance()->GetSaveStageInfo()){
    fResultsTree->Branch("StageInfo", "StageInfo", &fStageInfo, 64000, 1);
    }
    if(WCSimAnalysisConfig::Instance()->GetSaveHitComparison()){
    fResultsTree->Branch("HitComparison", "HitComparison", &fHitComparison, 64000, 1);
    }
    return;
}

void WCSimOutputTree::SetSeedInfo(SeedInfo& seedInfo)
{
	if( fResultsTree == 0x0)
	{
		std::cerr << "Calling SetSeedInfo before the resultsTree has been made!!!" << std::endl;
		MakeTree();
	}
	delete fSeedInfo;
	fSeedInfo = new SeedInfo(seedInfo);
	return;
}

void WCSimOutputTree::SetStageInfo(int stageNCalls, std::vector< WCSimHitPrediction > stagePreds, std::vector<WCSimLikelihoodTrackBase*> stageTracks)
{
	if( fResultsTree == 0x0)
	{
		std::cerr << "Calling SetStageInfo before the resultsTree has been made!!!" << std::endl;
		MakeTree();
	}
	fStageInfo->SetStageInfo(stageNCalls, stagePreds, stageTracks);
	return;
}

void WCSimOutputTree::SetHitComparison(WCSimHitComparison& hitComp)
{
	if( fResultsTree == 0x0)
	{
		std::cerr << "Calling SetHitComparison before the resultsTree has been made!!!" << std::endl;
		MakeTree();
	}
	delete fHitComparison;
    fHitComparison = new WCSimHitComparison(hitComp);
	return;
}

void WCSimOutputTree::Fill(EventHeader& eventHead,
          	  	  	  	   TruthInfo& truthInfo,
						   WCSimRecoSummary& summ,
						   HitInfo& hitInfo,
						   RecoInfo& recoInfo)
{

	// Fill the truth and fit tree entries
	if( fResultsTree == 0x0)
	{
		std::cerr << "Calling Fill before the resultsTree has been made!!!" << std::endl;
		MakeTree();
	}

	delete fEventHeader;
	fEventHeader = new EventHeader(eventHead);

    delete fTruthInfo;
    fTruthInfo = new TruthInfo(truthInfo);

    delete fRecoSummary;
    fRecoSummary = new WCSimRecoSummary(summ);

    delete fHitInfo;
    fHitInfo = new HitInfo(hitInfo);

    delete fRecoInfo;
    fRecoInfo = new RecoInfo(recoInfo);

    fResultsTree->Fill();

	// Fill the geometry tree
	if(fGeoTree->GetEntries() == 0)
	{
		*fGeometry = *(WCSimGeometry::Instance()->GetWCSimGeometry());
		WCSimGeometry::Instance()->PrintGeometry();
		fGeoTree->Fill();
	}

    SaveTree();
    return;
}

void WCSimOutputTree::SaveTree() {

	TDirectory* tmpd = 0;
	tmpd = gDirectory;
	fSaveFile->cd();
	std::cout << " *** WCSimFitter::SaveTree() *** " << std::endl;
	std::cout << "Save file name = " << fSaveFileName << std::endl;

	assert(fSaveFile->IsOpen());
	fResultsTree->Write("",TObject::kOverwrite);
	fGeoTree->Write("",TObject::kOverwrite);

	tmpd->cd();
	return;
}

void WCSimOutputTree::SetSaveFileName(TString saveName) {
	fSaveFileName = saveName;
	return;
}

TString WCSimOutputTree::GetSaveFileName() const
{
	return fSaveFileName;
}

void WCSimOutputTree::MakeSaveFile()
{
    std::cout << " *** WCSimOutputTree::MakeSaveFile *** " << std::endl;
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
    return;
}

void WCSimOutputTree::SetInputFile(const TString &saveFileName)
{
    fInputFile = std::string(saveFileName.Data());
	return;
}

TString WCSimOutputTree::GetInputFileName() const
{
	return fInputFile;
}

void WCSimOutputTree::DeletePointersIfExist()
{
    if(fResultsTree != 0x0)
    {
        delete fResultsTree;
        fResultsTree = 0x0;
    }

    if(fGeoTree != 0x0)
    {
        delete fGeoTree;
        fGeoTree = 0x0;
    }

    if(fEventHeader != 0x0)
    {
        delete fEventHeader;
        fEventHeader = 0x0;
    }

    if(fTruthInfo != 0x0)
    {
        delete fTruthInfo;
        fTruthInfo = 0x0;
    }

    if(fRecoSummary != 0x0)
    {
        delete fRecoSummary;
        fRecoSummary = 0x0;
    }

    if(fHitInfo != 0x0)
    {
        delete fHitInfo;
        fHitInfo = 0x0;
    }

    if(fRecoInfo != 0x0)
    {
        delete fRecoInfo;
        fRecoInfo = 0x0;
    }

    if(fSeedInfo != 0x0)
    {
        delete fSeedInfo;
        fSeedInfo = 0x0;
    }

    if(fStageInfo != 0x0)
    {
        delete fStageInfo;
        fStageInfo = 0x0;
    }

    if(fHitComparison != 0x0)
    {
        delete fHitComparison;
        fHitComparison = 0x0;
    }
    return;
}

