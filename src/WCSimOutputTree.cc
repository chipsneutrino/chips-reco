/*
 * WCSimOutputTree.cc
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#include "WCSimParameters.hh"
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
ClassImp (EventHeader);
ClassImp (ParameterInfo);
ClassImp (TruthInfo);
ClassImp (PidInfo);
ClassImp (SeedInfo);
ClassImp (StageInfo);
ClassImp (WCSimOutputTree);
#endif

///////////////////////////////////////////////////////////////
//  METHODS FOR EVENTHEADER CLASS                            //
///////////////////////////////////////////////////////////////

EventHeader::EventHeader() :
		fUID(""), fInputFile(""), fInputEvent(-999), fFailed(false), fRecoType(
				"") {
	// Empty
}

EventHeader::EventHeader(std::string inputFile, int inputEvent, bool failed,
		std::string recoType) :
		fInputFile(inputFile), fInputEvent(inputEvent), fFailed(failed), fRecoType(
				recoType) {
	BuildUID();
}

EventHeader::EventHeader(const EventHeader& other) :
		fUID(other.fUID), fInputFile(other.fInputFile), fInputEvent(
				other.fInputEvent), fFailed(other.fFailed), fRecoType(
				other.fRecoType) {
	// Empty
}

EventHeader& EventHeader::operator=(const EventHeader& rhs) {
	if (&rhs != this) {
		fUID = rhs.fUID;
		fInputFile = rhs.fInputFile;
		fInputEvent = rhs.fInputEvent;
		fFailed = rhs.fFailed;
		fRecoType = rhs.fRecoType;
	}
	return *this;
}

EventHeader::~EventHeader() {
	// Empty
}

void EventHeader::Print() {
	std::cout << "EventHeader::Print()..." << std::endl << "UID -> " << fUID
			<< std::endl << "InputFile -> " << fInputFile << std::endl
			<< "InputEvent -> " << fInputEvent << std::endl << "Failed -> "
			<< fFailed << std::endl << "RecoType -> " << fRecoType << std::endl;
	return;
}

/**
 * @brief Make unique identifier for a file/event combination.  Sets fUID to
 * md5("%s%d", md5(inputFile), eventNumber)
 *
 **/
void EventHeader::BuildUID() {
	TMD5 myMD5;

	// Build a char array to hold our filename string:
	TString inputStr = TString::Format("%s%06d", fInputFile.c_str(),
			fInputEvent);
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
//  METHODS FOR PARAMETERINFO CLASS                          //
///////////////////////////////////////////////////////////////

ParameterInfo::ParameterInfo() {
	// Need to set everything from the WCSimParameters instance...
	WCSimParameters* params = WCSimParameters::Instance();

	// Slicer parameters
	fSlicerClusterDistance = params->GetSlicerClusterDistance();
	fSlicerMinSize = params->GetSlicerMinSize();
	fSlicerChargeCut = params->GetSlicerMinSize();
	fSlicerTimeCut = params->GetSlicerTimeCut();
	fIterateSlicing = params->IterateSlicing();

	// Veto slicing parameters
	fVetoClusterDistance = params->GetVetoClusterDistance();
	fVetoMinSize = params->GetVetoMinSize();
	fVetoMinChargeCut = params->GetVetoMinChargeCut();
	fVetoMaxChargeCut = params->GetVetoMaxChargeCut();
	fVetoTimeCut = params->GetVetoTimeCut();

	// Integral parameters
	fCalculateIntegrals = params->CalculateIntegrals();
	fTruncateIntegrals = params->TruncateIntegrals();
	fConstrainExtent = params->ConstrainExtent();

	// Likelihood tuner parameters
	fUseTransmission = params->UseTransmission();
	fUseAngularEfficiency = params->UseAngularEfficiency();
	fUseGlassCathodeReflection = params->UseGlassCathodeReflection();
	fUseScatteringTable = params->UseScatteringTable();
	fUseNewAngularEfficiency = params->UseNewAngularEfficiency();
	fUseTrackFit = params->UseTrackFit();

	// Fitter parameters
	fUseTime = params->UseTime();
	fUseCharge = params->UseCharge();
	fEqualiseChargeAndTime = params->EqualiseChargeAndTime();
	fSaveWCSimRootEvent = params->EqualiseChargeAndTime();
	fDigiType = params->GetDigiType();
	fSaveSeedInfo = params->SaveSeedInfo();
	fSaveStageInfo = params->SaveStageInfo();
	fSaveHitComparison = params->SaveHitComparison();
	fSaveParameters = params->SaveParameterInfo();

	// Speed of light parameters
	fUseCustomParticleSpeed = params->UseCustomParticleSpeed();
	fUseCustomSpeedOfLight = params->UseCustomSpeedOfLight();
	fUseFittedSpeedOfLight = params->UseFittedSpeedOfLight();
	fCustomParticleSpeed = params->GetCustomParticleSpeed();
	fCustomSpeedOfLight = params->GetCustomSpeedOfLight();
	fFittedSpeedOfLight = params->GetFittedSpeedOfLight();

	// Not currently used parameters
	fUseSimpleTimeResolution = params->GetUseSimpleTimeResolution();
	fUseSimpleTimeSlew = params->GetUseSimpleTimeSlew();
	fUseSimpleRefractiveIndex = params->GetUseSimpleRefractiveIndex();
}

ParameterInfo::~ParameterInfo() {
	// Empty
}

void ParameterInfo::Print() {
	std::cout << "ParameterInfo::Print()..." << std::endl << "SlicerClusterDistance -> " << fSlicerClusterDistance
			<< std::endl << "SlicerMinSize -> " << fSlicerMinSize << std::endl << "SlicerChargeCut -> "
			<< fSlicerChargeCut << std::endl << "SlicerTimeCut -> " << fSlicerTimeCut << std::endl
			<< "IterateSlicing -> " << fIterateSlicing << std::endl << "VetoClusterDistance -> " << fVetoClusterDistance
			<< std::endl << "VetoMinSize -> " << fVetoMinSize << std::endl << "VetoMinChargeCut -> "
			<< fVetoMinChargeCut << std::endl << "VetoMaxChargeCut -> " << fVetoMaxChargeCut << std::endl
			<< "VetoTimeCut -> " << fVetoTimeCut << std::endl << "CalculateIntegrals -> " << fCalculateIntegrals
			<< std::endl << "TruncateIntegrals -> " << fTruncateIntegrals << std::endl << "ConstrainExtent -> "
			<< fConstrainExtent << std::endl << "UseTransmission -> " << fUseTransmission << std::endl
			<< "UseAngularEfficiency -> " << fUseAngularEfficiency << std::endl << "UseGlassCathodeReflection -> "
			<< fUseGlassCathodeReflection << std::endl << "UseScatteringTable -> " << fUseScatteringTable << std::endl
			<< "UseNewAngularEfficiency -> " << fUseNewAngularEfficiency << std::endl << "UseTrackFit -> " << fUseTrackFit << std::endl
			<< "UseTime -> " << fUseTime << std::endl << "UseCharge -> " << fUseCharge << std::endl
			<< "EqualiseChargeAndTime -> " << fEqualiseChargeAndTime << std::endl << "SaveWCSimRootEvent -> "
			<< fSaveWCSimRootEvent << std::endl << "DigiType -> " << fDigiType << std::endl << "SaveSeedInfo -> "
			<< fSaveSeedInfo << std::endl << "SaveStageInfo -> " << fSaveStageInfo << std::endl
			<< "SaveHitComparison -> " << fSaveHitComparison << std::endl << "SaveParameters -> " << fSaveParameters << std::endl
			<< "UseCustomParticleSpeed -> " << fUseCustomParticleSpeed << std::endl << "UseCustomSpeedOfLight -> "
			<< fUseCustomSpeedOfLight << std::endl << "UseFittedSpeedOfLight -> " << fUseFittedSpeedOfLight << std::endl
			<< "CustomParticleSpeed -> " << fCustomParticleSpeed << std::endl << "CustomSpeedOfLight -> "
			<< fCustomSpeedOfLight << std::endl << "FittedSpeedOfLight -> " << fFittedSpeedOfLight << std::endl
			<< "UseSimpleTimeResolution -> " << fUseSimpleTimeResolution << std::endl << "UseSimpleTimeSlew -> "
			<< fUseSimpleTimeSlew << std::endl << "UseSimpleRefractiveIndex -> " << fUseSimpleRefractiveIndex
			<< std::endl;
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR TRUTHINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

TruthInfo::TruthInfo() :
		fType(-999), fBeamPDG(-999), fBeamEnergy(0.0), fNPrimaries(-999), fIsCC(
				false), fIsNC(false), fIsRes(false), fIsDIS(false), fIsCoherent(
				false), fIsNueElectronElastic(false), fIsInverseMuonDecay(
				false), fIsOther(false), fVtxTime(0.0), fVtxX(0.0), fVtxY(0.0), fVtxZ(
				0.0), fBeamDirX(0.0), fBeamDirY(0.0), fBeamDirZ(0.0) {
	fPrimaryPDGs.clear();
	fPrimaryEnergies.clear();
	fPrimaryDirs.clear();
}

TruthInfo::TruthInfo(int type, int beamPDG, float beamEnergy, int nPrimaries,
		std::vector<int> primaryPDGs, std::vector<double> primaryEnergies,
		std::vector<TVector3> primaryDirs) :
		fType(type), fBeamPDG(beamPDG), fBeamEnergy(beamEnergy), fNPrimaries(
				nPrimaries), fPrimaryPDGs(primaryPDGs), fPrimaryEnergies(
				primaryEnergies), fPrimaryDirs(primaryDirs) {
	fIsCC = WCSimTruthSummary::TypeIsCCEvent(fType);
	fIsNC = WCSimTruthSummary::TypeIsNCEvent(fType);
	fIsQE = WCSimTruthSummary::TypeIsQEEvent(fType);
	fIsRes = WCSimTruthSummary::TypeIsResEvent(fType);
	fIsDIS = WCSimTruthSummary::TypeIsDISEvent(fType);
	fIsCoherent = WCSimTruthSummary::TypeIsCohEvent(fType);
	fIsNueElectronElastic = WCSimTruthSummary::TypeIsNueElectronElasticEvent(
			type);
	fIsInverseMuonDecay = WCSimTruthSummary::TypeIsInverseMuonDecayEvent(type);
	fIsOther = fIsCC
			&& !(fIsQE || fIsRes || fIsDIS || fIsCoherent
					|| fIsNueElectronElastic || fIsInverseMuonDecay);
	fVtxTime = 0.0;
	fVtxX = 0.0;
	fVtxY = 0.0;
	fVtxZ = 0.0;
	fBeamDirX = 0.0;
	fBeamDirY = 0.0;
	fBeamDirZ = 0.0;
}

TruthInfo::TruthInfo(const TruthInfo& other) :
		fType(other.fType), fBeamPDG(other.fBeamPDG), fBeamEnergy(
				other.fBeamEnergy), fNPrimaries(other.fNPrimaries), fPrimaryPDGs(
				other.fPrimaryPDGs), fPrimaryEnergies(other.fPrimaryEnergies), fPrimaryDirs(
				other.fPrimaryDirs), fIsCC(other.fIsCC), fIsNC(other.fIsNC), fIsQE(
				other.fIsQE), fIsRes(other.fIsRes), fIsDIS(other.fIsDIS), fIsCoherent(
				other.fIsCoherent), fIsNueElectronElastic(
				other.fIsNueElectronElastic), fIsInverseMuonDecay(
				other.fIsInverseMuonDecay), fIsOther(other.fIsOther), fVtxTime(
				other.fVtxTime), fVtxX(other.fVtxX), fVtxY(other.fVtxY), fVtxZ(
				other.fVtxZ), fBeamDirX(other.fBeamDirX), fBeamDirY(
				other.fBeamDirY), fBeamDirZ(other.fBeamDirZ) {
	// Empty
}

TruthInfo& TruthInfo::operator=(const TruthInfo& rhs) {
	if (this != &rhs) {
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

TruthInfo::~TruthInfo() {
	// Empty
}

void TruthInfo::Print() {
	std::cout << "TruthInfo::Print()..." << std::endl;
	std::cout << "Beam: PDG->" << fBeamPDG << ", Energy->" << fBeamEnergy
			<< ", Dir->(" << fBeamDirX << "," << fBeamDirY << "," << fBeamDirZ
			<< ")" << std::endl;
	std::cout << "Vertex: Position->(" << fVtxX << "," << fVtxY << "," << fVtxZ
			<< "), Time->" << fVtxTime << std::endl;

	std::cout << "Interaction type code -> " << fType << std::endl;
	std::cout << "Interaction is";
	if (fIsCC) {
		std::cout << ", Charged-current";
	}
	if (fIsNC) {
		std::cout << ", Neutral-current";
	}
	if (fIsQE) {
		std::cout << ", Quasielastic";
	}
	if (fIsRes) {
		std::cout << ", Resonant";
	}
	if (fIsDIS) {
		std::cout << ", Deep Inelastic Scattering";
	}
	if (fIsCoherent) {
		std::cout << ", Coherent";
	}
	if (fIsNueElectronElastic) {
		std::cout << ", Electron elastic";
	}
	if (fIsInverseMuonDecay) {
		std::cout << ", Inverse muon decay";
	}
	if (fIsOther) {
		std::cout << ", Other";
	}
	std::cout << std::endl;

	// Print info about the primaries...
	for (int p = 0; p < fNPrimaries; p++) {
		std::cout << "Primary " << p << ": PDG->" << fPrimaryPDGs[p]
				<< ", Energy->" << fPrimaryEnergies[p] << ", Dir->("
				<< fPrimaryDirs[p].X() << "," << fPrimaryDirs[p].Y() << ","
				<< fPrimaryDirs[p].Z() << ")" << std::endl;
	}
	return;
}

void TruthInfo::SetVtxTime(float t) {
	fVtxTime = t;
	return;
}

void TruthInfo::SetVtx(float x, float y, float z) {
	fVtxX = x;
	fVtxY = y;
	fVtxZ = z;
	return;
}

void TruthInfo::SetBeamDir(float x, float y, float z) {
	fBeamDirX = x;
	fBeamDirY = y;
	fBeamDirZ = z;
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR PIDINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

PidInfo::PidInfo() :
		fVeto(false), fNVetoHits(0), fNHits(0), fNHitsUpstream(0), fNHitsDownstream(
				0), fNHitsInBottom(0), fNHitsInTop(0), fNHitsAboveMid(0), fNHitsBelowMid(
				0), fFracHitsUpstream(0.0), fFracHitsDownstream(0.0), fFracHitsInBottom(
				0.0), fFracHitsInTop(0.0), fFracHitsAboveMid(0.0), fFracHitsBelowMid(
				0.0), fTotalQ(0.0), fTotalQUpstream(0.0), fTotalQDownstream(
				0.0), fTotalQInBottom(0.0), fTotalQInTop(0.0), fTotalQAboveMid(
				0.0), fTotalQBelowMid(0.0), fFracQUpstream(0.0), fFracQDownstream(
				0.0), fFracQInBottom(0.0), fFracQInTop(0.0), fFracQAboveMid(
				0.0), fFracQBelowMid(0.0), fTotal2LnL(0.0), fHit2LnL(0.0), fCharge2LnL(
				0.0), fTime2LnL(0.0), fCutoff2LnL(0.0), fTotalQInRing(0.0), fTotalQOutsideRing(
				0.0), fTotalQInRingHole(0.0), fNHitsInRing(0), fNHitsOutsideRing(
				0), fNHitsInRingHole(0), fPredQInRing(0.0), fPredQOutsideRing(
				0.0), fPredQInRingHole(0.0), fFracTotalQInRing(0.0), fFracTotalQOutsideRing(
				0.0), fFracTotalQInRingHole(0.0), fFracNHitsInRing(0.0), fFracNHitsOutsideRing(
				0.0), fFracNHitsInRingHole(0.0), fFracPredQInRing(0.0), fFracPredQOutsideRing(
				0.0), fFracPredQInRingHole(0.0), fPredictedOverTotalCharge(0.0), fVtxTime(
				0.0), fEnergy(0.0), fVtxX(0.0), fVtxY(0.0), fVtxZ(0.0), fVtxRho(
				0.0), fEndX(0.0), fEndY(0.0), fEndZ(0.0), fEndRho(0.0), fDirX(
				0.0), fDirY(0.0), fDirZ(0.0), fEscapes(false) {
	// Empty
}

PidInfo::PidInfo(const PidInfo& other) :
		fVeto(other.fVeto), fNVetoHits(other.fNVetoHits), fNHits(other.fNHits), fNHitsUpstream(
				other.fNHitsUpstream), fNHitsDownstream(other.fNHitsDownstream), fNHitsInBottom(
				other.fNHitsInBottom), fNHitsInTop(other.fNHitsInTop), fNHitsAboveMid(
				other.fNHitsAboveMid), fNHitsBelowMid(other.fNHitsBelowMid), fFracHitsUpstream(
				other.fFracHitsUpstream), fFracHitsDownstream(
				other.fFracHitsDownstream), fFracHitsInBottom(
				other.fFracHitsInBottom), fFracHitsInTop(other.fFracHitsInTop), fFracHitsAboveMid(
				other.fFracHitsAboveMid), fFracHitsBelowMid(
				other.fFracHitsBelowMid), fTotalQ(other.fTotalQ), fTotalQUpstream(
				other.fTotalQUpstream), fTotalQDownstream(
				other.fTotalQDownstream), fTotalQInBottom(
				other.fTotalQInBottom), fTotalQInTop(other.fTotalQInTop), fTotalQAboveMid(
				other.fTotalQAboveMid), fTotalQBelowMid(other.fTotalQBelowMid), fFracQUpstream(
				other.fFracQUpstream), fFracQDownstream(other.fFracQDownstream), fFracQInBottom(
				other.fFracQInBottom), fFracQInTop(other.fFracQInTop), fFracQAboveMid(
				other.fFracQAboveMid), fFracQBelowMid(other.fFracQBelowMid), fTotal2LnL(
				other.fTotal2LnL), fHit2LnL(other.fHit2LnL), fCharge2LnL(
				other.fCharge2LnL), fTime2LnL(other.fTime2LnL), fCutoff2LnL(
				other.fCutoff2LnL), fTotalQInRing(other.fTotalQInRing), fTotalQOutsideRing(
				other.fTotalQOutsideRing), fTotalQInRingHole(
				other.fTotalQInRingHole), fNHitsInRing(other.fNHitsInRing), fNHitsOutsideRing(
				other.fNHitsOutsideRing), fNHitsInRingHole(
				other.fNHitsInRingHole), fPredQInRing(other.fPredQInRing), fPredQOutsideRing(
				other.fPredQOutsideRing), fPredQInRingHole(
				other.fPredQInRingHole), fFracTotalQInRing(
				other.fFracTotalQInRing), fFracTotalQOutsideRing(
				other.fFracTotalQOutsideRing), fFracTotalQInRingHole(
				other.fFracTotalQInRingHole), fFracNHitsInRing(
				other.fFracNHitsInRing), fFracNHitsOutsideRing(
				other.fFracNHitsOutsideRing), fFracNHitsInRingHole(
				other.fFracNHitsInRingHole), fFracPredQInRing(
				other.fFracPredQInRing), fFracPredQOutsideRing(
				other.fFracPredQOutsideRing), fFracPredQInRingHole(
				other.fFracPredQInRingHole), fPredictedOverTotalCharge(
				other.fPredictedOverTotalCharge), fVtxTime(other.fVtxTime), fEnergy(
				other.fEnergy), fVtxX(other.fVtxX), fVtxY(other.fVtxY), fVtxZ(
				other.fVtxZ), fVtxRho(other.fVtxRho), fEndX(other.fEndX), fEndY(
				other.fEndY), fEndZ(other.fEndZ), fEndRho(other.fEndRho), fDirX(
				other.fDirX), fDirY(other.fDirY), fDirZ(other.fDirZ), fEscapes(
				other.fEscapes) {
	// Empty
}

PidInfo& PidInfo::operator=(const PidInfo& rhs) {
	if (&rhs != this) {
		fVeto = rhs.fVeto;
		fNVetoHits = rhs.fNVetoHits;
		fNHits = rhs.fNHits;
		fNHitsUpstream = rhs.fNHitsUpstream;
		fNHitsDownstream = rhs.fNHitsDownstream;
		fNHitsInBottom = rhs.fNHitsInBottom;
		fNHitsInTop = rhs.fNHitsInTop;
		fNHitsAboveMid = rhs.fNHitsAboveMid;
		fNHitsBelowMid = rhs.fNHitsBelowMid;
		fFracHitsUpstream = rhs.fFracHitsUpstream;
		fFracHitsDownstream = rhs.fFracHitsDownstream;
		fFracHitsInBottom = rhs.fFracHitsInBottom;
		fFracHitsInTop = rhs.fFracHitsInTop;
		fFracHitsAboveMid = rhs.fFracHitsAboveMid;
		fFracHitsBelowMid = rhs.fFracHitsBelowMid;
		fTotalQ = rhs.fTotalQ;
		fTotalQUpstream = rhs.fTotalQUpstream;
		fTotalQDownstream = rhs.fTotalQDownstream;
		fTotalQInBottom = rhs.fTotalQInBottom;
		fTotalQInTop = rhs.fTotalQInTop;
		fTotalQAboveMid = rhs.fTotalQAboveMid;
		fTotalQBelowMid = rhs.fTotalQBelowMid;
		fFracQUpstream = rhs.fFracQUpstream;
		fFracQDownstream = rhs.fFracQDownstream;
		fFracQInBottom = rhs.fFracQInBottom;
		fFracQInTop = rhs.fFracQInTop;
		fFracQAboveMid = rhs.fFracQAboveMid;
		fFracQBelowMid = rhs.fFracQBelowMid;

		fTotal2LnL = rhs.fTotal2LnL;
		fHit2LnL = rhs.fHit2LnL;
		fCharge2LnL = rhs.fCharge2LnL;
		fTime2LnL = rhs.fTime2LnL;
		fCutoff2LnL = rhs.fCutoff2LnL;
		fTotalQInRing = rhs.fTotalQInRing;
		fTotalQOutsideRing = rhs.fTotalQOutsideRing;
		fTotalQInRingHole = rhs.fTotalQInRingHole;
		fNHitsInRing = rhs.fNHitsInRing;
		fNHitsOutsideRing = rhs.fNHitsOutsideRing;
		fNHitsInRingHole = rhs.fNHitsInRingHole;
		fPredQInRing = rhs.fPredQInRing;
		fPredQOutsideRing = rhs.fPredQOutsideRing;
		fPredQInRingHole = rhs.fPredQInRingHole;
		fFracTotalQInRing = rhs.fFracTotalQInRing;
		fFracTotalQOutsideRing = rhs.fFracTotalQOutsideRing;
		fFracTotalQInRingHole = rhs.fFracTotalQInRingHole;
		fFracNHitsInRing = rhs.fFracNHitsInRing;
		fFracNHitsOutsideRing = rhs.fFracNHitsOutsideRing;
		fFracNHitsInRingHole = rhs.fFracNHitsInRingHole;
		fFracPredQInRing = rhs.fFracPredQInRing;
		fFracPredQOutsideRing = rhs.fFracPredQOutsideRing;
		fFracPredQInRingHole = rhs.fFracPredQInRingHole;
		fPredictedOverTotalCharge = rhs.fPredictedOverTotalCharge;
		fVtxTime = rhs.fVtxTime;
		fEnergy = rhs.fEnergy;
		fVtxX = rhs.fVtxX;
		fVtxY = rhs.fVtxY;
		fVtxZ = rhs.fVtxZ;
		fVtxRho = rhs.fVtxRho;
		fEndX = rhs.fEndX;
		fEndY = rhs.fEndY;
		fEndZ = rhs.fEndZ;
		fEndRho = rhs.fEndRho;
		fDirX = rhs.fDirX;
		fDirY = rhs.fDirY;
		fDirZ = rhs.fDirZ;
		fEscapes = rhs.fEscapes;
	}
	return *this;
}

PidInfo::~PidInfo() {
	// Empty
}

void PidInfo::Print() {
	std::cout << "HitInfo::Print()..." << std::endl << "Veto -> " << fVeto
			<< std::endl << "fNVetoHits -> " << fNVetoHits << std::endl
			<< "NHits -> " << fNHits << std::endl << "NHitsUpstream -> "
			<< fNHitsUpstream << std::endl << "NHitsDownstream -> "
			<< fNHitsDownstream << std::endl << "NHitsInBottom -> "
			<< fNHitsInBottom << std::endl << "NHitsInTop -> " << fNHitsInTop
			<< std::endl << "fNHitsAboveMid -> " << fNHitsAboveMid << std::endl
			<< "fNHitsBelowMid -> " << fNHitsBelowMid << std::endl
			<< "FracHitsUpstream -> " << fFracHitsUpstream << std::endl
			<< "FracHitsDownstream -> " << fFracHitsDownstream << std::endl
			<< "fFracHitsInBottom -> " << fFracHitsInBottom << std::endl
			<< "fFracHitsInTop -> " << fFracHitsInTop << std::endl
			<< "fFracHitsAboveMid -> " << fFracHitsAboveMid << std::endl
			<< "fFracHitsBelowMid -> " << fFracHitsBelowMid << std::endl
			<< "TotalQ -> " << fTotalQ << std::endl << "TotalQUpstream -> "
			<< fTotalQUpstream << std::endl << "TotalQDownstream -> "
			<< fTotalQDownstream << std::endl << "TotalQInBottom -> "
			<< fTotalQInBottom << std::endl << "TotalQInTop -> " << fTotalQInTop
			<< std::endl << "fTotalQAboveMid -> " << fTotalQAboveMid
			<< std::endl << "fTotalQBelowMid -> " << fTotalQBelowMid
			<< std::endl << "FracQUpstream -> " << fFracQUpstream << std::endl
			<< "FracQDownstream -> " << fFracQDownstream << std::endl
			<< "FracQInBottom -> " << fFracQInBottom << std::endl
			<< "FracQInTop -> " << fFracQInTop << std::endl
			<< "fFracQAboveMid -> " << fFracQAboveMid << std::endl
			<< "fFracQBelowMid -> " << fFracQBelowMid << std::endl;

	std::cout << "RecoInfo::Print()..." << std::endl << "fTotal2LnL"
			<< fTotal2LnL << std::endl << "fHit2LnL" << fHit2LnL << std::endl
			<< "fCharge2LnL" << fCharge2LnL << std::endl << "fTime2LnL"
			<< fTime2LnL << std::endl << "fCutoff2LnL" << fCutoff2LnL
			<< std::endl << "fTotalQInRing" << fTotalQInRing << std::endl
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
			<< "fPredictedOverTotalCharge" << fPredictedOverTotalCharge
			<< std::endl << "fVtxTime" << fVtxTime << std::endl << "fEnergy"
			<< fEnergy << std::endl << "fVtxX" << fVtxX << std::endl << "fVtxY"
			<< fVtxY << std::endl << "fVtxZ" << fVtxZ << std::endl << "fVtxRho"
			<< fVtxRho << std::endl << "fEndX" << fEndX << std::endl << "fEndY"
			<< fEndY << std::endl << "fEndZ" << fEndZ << std::endl << "fEndRho"
			<< fEndRho << std::endl << "fDirX" << fDirX << std::endl << "fDirY"
			<< fDirY << std::endl << "fDirZ" << fDirZ << std::endl << "fEscapes"
			<< fEscapes << std::endl;
	return;
}

void PidInfo::SetHitInfo(bool veto, int NVetoHits, int NHits, int NHitsUpstream,
		int NHitsInBottom, int NHitsInTop, int NHitsAboveMid, float totalQ,
		float totalQUpstream, float totalQInBottom, float totalQInTop,
		float totalQAboveMid) {
	fVeto = veto;
	fNHits = NHits;
	fNVetoHits = NVetoHits;
	fNHitsUpstream = NHitsUpstream;
	fNHitsInBottom = NHitsInBottom;
	fNHitsInTop = NHitsInTop;
	fNHitsAboveMid = NHitsAboveMid;
	fTotalQ = totalQ;
	fTotalQUpstream = totalQUpstream;
	fTotalQInBottom = totalQInBottom;
	fTotalQInTop = totalQInTop;
	fTotalQAboveMid = totalQAboveMid;

	fNHitsDownstream = fNHits - fNHitsUpstream;
	fFracHitsUpstream = fNHitsUpstream / (static_cast<float>(fNHits));
	fFracHitsDownstream = 1 - fFracHitsUpstream;
	fFracHitsInBottom = fNHitsInBottom / (static_cast<float>(fNHits));
	fFracHitsInTop = fNHitsInTop / (static_cast<float>(fNHits));

	fNHitsBelowMid = fNHits - fNHitsAboveMid;
	fFracHitsAboveMid = fNHitsAboveMid / (static_cast<float>(fNHits));
	fFracHitsBelowMid = 1 - fFracHitsAboveMid;

	fTotalQDownstream = fTotalQ - fTotalQUpstream;
	fFracQUpstream = fTotalQUpstream / fTotalQ;
	fFracQDownstream = 1 - fFracQUpstream;
	fFracQInBottom = fTotalQInBottom / fTotalQ;
	fFracQInTop = fTotalQInTop / fTotalQ;

	fTotalQBelowMid = fTotalQ - fTotalQAboveMid;
	fFracQAboveMid = fTotalQAboveMid / fTotalQ;
	fFracQBelowMid = 1 - fFracQAboveMid;

	return;
}

void PidInfo::SetLikelihoods(float total, float hit, float charge, float time) {
	fTotal2LnL = total;
	fHit2LnL = hit;
	fCharge2LnL = charge;
	fTime2LnL = time;
	fCutoff2LnL = (hit + charge + time) - total;
	return;
}

void PidInfo::SetTotalQ(float in, float out, float hole) {
	fTotalQInRing = in;
	fTotalQOutsideRing = out;
	fTotalQInRingHole = hole;
	float totalQ = in + out + hole;

	if (totalQ > 0) {
		fFracTotalQInRing = fTotalQInRing / totalQ;
		fFracTotalQOutsideRing = fTotalQOutsideRing / totalQ;
		fFracTotalQInRingHole = fTotalQInRingHole / totalQ;
	} else {
		fFracTotalQInRing = 0.0;
		fFracTotalQOutsideRing = 0.0;
		fFracTotalQInRingHole = 0.0;
	}

	if (totalQ > 0) {
		fPredictedOverTotalCharge = (fPredQInRing + fPredQOutsideRing
				+ fPredQInRingHole) / totalQ;
	} else {
		fPredictedOverTotalCharge = 0.0;
	}
	return;
}

void PidInfo::SetNHits(int in, int out, int hole) {
	fNHitsInRing = in;
	fNHitsOutsideRing = out;
	fNHitsInRingHole = hole;
	float nHits = in + out + hole;

	if (nHits > 0) {
		fFracNHitsInRing = fNHitsInRing / nHits;
		fFracNHitsOutsideRing = fNHitsOutsideRing / nHits;
		fFracNHitsInRingHole = fNHitsInRingHole / nHits;
	} else {
		fFracNHitsInRing = 0.0;
		fFracNHitsOutsideRing = 0.0;
		fFracNHitsInRingHole = 0.0;
	}
	return;
}

void PidInfo::SetPredQ(float in, float out, float hole) {
	fPredQInRing = in;
	fPredQOutsideRing = out;
	fPredQInRingHole = hole;
	float predQ = in + out + hole;

	if (predQ > 0) {
		fFracPredQInRing = fPredQInRing / predQ;
		fFracPredQOutsideRing = fPredQOutsideRing / predQ;
		fFracPredQInRingHole = fPredQInRingHole / predQ;
	} else {
		fFracPredQInRing = 0.0;
		fFracPredQOutsideRing = 0.0;
		fFracPredQInRingHole = 0.0;
	}

	float totalQ = fTotalQInRing + fTotalQOutsideRing + fTotalQInRingHole;
	if (totalQ > 0) {
		fPredictedOverTotalCharge = predQ / totalQ;
	} else {
		fPredictedOverTotalCharge = 0.0;
	}
	return;
}

void PidInfo::SetEnergy(float E) {
	fEnergy = E;
	return;
}

void PidInfo::SetVtxTime(float t) {
	fVtxTime = t;
	return;
}

void PidInfo::SetVtx(float x, float y, float z) {
	fVtxX = x;
	fVtxY = y;
	fVtxZ = z;
	fVtxRho = TMath::Sqrt(x * x + y * y);
	return;
}

void PidInfo::SetEnd(float x, float y, float z) {
	fEndX = x;
	fEndY = y;
	fEndZ = z;
	fEndRho = TMath::Sqrt(x * x + y * y);
	return;
}

void PidInfo::SetDir(float x, float y, float z) {
	fDirX = x;
	fDirY = y;
	fDirZ = z;
	return;
}

void PidInfo::SetEscapes(bool escapes) {
	fEscapes = escapes;
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR SEEDINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

SeedInfo::SeedInfo() :
		fNSlices(-999) {
	fTracks.clear();
	fRingHeight.clear();
	fRingTime.clear();
	fRingAngle.clear();
	fRingVtx.clear();
	fRingDir.clear();
}

SeedInfo::SeedInfo(std::vector<WCSimLikelihoodTrackBase*> tracks, int slices,
		std::vector<WCSimRecoRing*> ringVec, std::vector<double> ringTime) :
		fTracks(tracks), fNSlices(slices), fRingTime(ringTime) {
	fRingHeight.clear();
	fRingAngle.clear();
	fRingVtx.clear();
	fRingDir.clear();
	for (int r = 0; (size_t) r < ringVec.size(); r++) {
		fRingHeight.push_back(ringVec[r]->GetHeight());
		fRingAngle.push_back(ringVec[r]->GetAngle());
		fRingVtx.push_back(
				TVector3(ringVec[r]->GetVtxX(), ringVec[r]->GetVtxY(),
						ringVec[r]->GetVtxZ()));
		fRingDir.push_back(
				TVector3(ringVec[r]->GetDirX(), ringVec[r]->GetDirY(),
						ringVec[r]->GetDirZ()));
	}
}

SeedInfo::SeedInfo(const SeedInfo& other) :
		fTracks(other.fTracks), fNSlices(other.fNSlices), fRingHeight(
				other.fRingHeight), fRingTime(other.fRingTime), fRingAngle(
				other.fRingAngle), fRingVtx(other.fRingVtx), fRingDir(
				other.fRingDir) {
	// Empty...
}

SeedInfo& SeedInfo::operator=(const SeedInfo& rhs) {
	if (this != &rhs) {
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

SeedInfo::~SeedInfo() {
	// Empty...
}

void SeedInfo::Print() {
	std::cout << "SeedInfo::Print()..." << std::endl;
	std::cout << "NSlices -> " << fNSlices << std::endl;
	std::cout << "NRings -> " << fRingTime.size() << std::endl;
	std::cout << "NTracks -> " << fTracks.size() << std::endl;

	// Print out the rings...
	for (int r = 0; (size_t) r < fRingTime.size(); r++) {
		std::cout << "Ring " << r << ": Height->" << fRingHeight[r]
				<< ", Time->" << fRingTime[r] << ", Angle->" << fRingAngle[r]
				<< std::endl;
		std::cout << "VtxPos->(" << fRingVtx[r].X() << "," << fRingVtx[r].Y()
				<< "," << fRingVtx[r].Z() << ") Dir->(" << fRingDir[r].X()
				<< "," << fRingDir[r].Y() << "," << fRingDir[r].Z() << ")"
				<< std::endl;
	}

	// Print out the track given to the fitter...
	for (int t = 0; (size_t) t < fTracks.size(); t++) {
		WCSimLikelihoodTrackBase* track = fTracks[t];
		std::cout << "Track " << t << ", vtx = (" << track->GetX() << ","
				<< track->GetY() << "," << track->GetZ() << "), dir = ("
				<< track->GetTheta() << "," << track->GetPhi() << "), Energy = "
				<< track->GetE() << std::endl;
	}
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR STAGEINFO CLASS                            	 //
///////////////////////////////////////////////////////////////

StageInfo::StageInfo() {
	fStageNcalls.clear();
	fStagePreds.clear();
	fStageTracks.clear();
}

StageInfo::StageInfo(const StageInfo& other) :
		fStageNcalls(other.fStageNcalls), fStagePreds(other.fStagePreds), fStageTracks(
				other.fStageTracks) {
	// Empty
}

StageInfo& StageInfo::operator=(const StageInfo& rhs) {
	if (this != &rhs) {
		fStageNcalls = rhs.fStageNcalls;
		fStagePreds = rhs.fStagePreds;
		fStageTracks = rhs.fStageTracks;
	}
	return *this;
}

StageInfo::~StageInfo() {
	// Empty
}

void StageInfo::Print() {
	std::cout << "StageInfo::Print()..." << std::endl;
	std::cout << "Number of Stages -> " << fStageNcalls.size() << std::endl;
	std::cout << "Number of Tracks -> " << fStageTracks[0].size() << std::endl
			<< std::endl;

	for (int stage = 0; (size_t) stage < fStageNcalls.size(); stage++) {
		if (stage == 0) {
			std::cout << "(seed)..." << std::endl;
		} else {
			std::cout << "(stage " << stage << ")..." << std::endl;
		}

		std::cout << "NCalls -> " << fStageNcalls[stage] << std::endl;

		// Print the track parameters at this stage...
		for (int t = 0; (size_t) t < fStageTracks[stage].size(); t++) {
			WCSimLikelihoodTrackBase* track = fStageTracks[stage][t];
			std::cout << "Track " << t << ", vtx = (" << track->GetX() << ","
					<< track->GetY() << "," << track->GetZ() << "), dir = ("
					<< track->GetTheta() << "," << track->GetPhi()
					<< "), Energy = " << track->GetE() << std::endl;
		}

		// Print the likelihoods and charge at this stage...
		double totPredQ = 0;
		double totTotal2LnL = 0, totTime2LnL = 0, totCharge2LnL = 0;
		for (int p = 0; (size_t) p < fStagePreds[stage].size(); p++) {
			totPredQ += fStagePreds[stage][p].GetPredictedCharge();
			totTotal2LnL += fStagePreds[stage][p].GetTotal2LnL();
			totTime2LnL += fStagePreds[stage][p].GetTime2LnL();
			totCharge2LnL += fStagePreds[stage][p].GetCharge2LnL();
		}
		std::cout << "NCalls -> " << fStageNcalls[stage] << ", TotalPredQ -> "
				<< totPredQ << std::endl;
		std::cout << "TotTotal2LnL -> " << totTotal2LnL << "TotTime2LnL -> "
				<< totTime2LnL << "totCharge2LnL -> " << totCharge2LnL
				<< std::endl << std::endl;
	}
	return;
}

void StageInfo::SetStageInfo(int stageNCalls,
		std::vector<WCSimHitPrediction> stagePreds,
		std::vector<WCSimLikelihoodTrackBase*> stageTracks) {
	fStageNcalls.push_back(stageNCalls);
	fStagePreds.push_back(stagePreds);
	fStageTracks.push_back(stageTracks);
	return;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR WCSimOutputTree CLASS                        //
///////////////////////////////////////////////////////////////

WCSimOutputTree::WCSimOutputTree() :
		fSaveFile(0x0), fResultsTree(0x0), fGeoTree(0x0), fGeometry(0x0), fParameterTree(0x0),
		fParameterInfo(0x0), fEventHeader(0x0), fTruthInfo(0x0), fRecoSummary(0x0), fPidInfo(0x0),
		fSeedInfo(0x0), fStageInfo(0x0), fHitComparison(0x0) {
	//fSaveFileName.Form("%s_tree.root", saveFileName.Data());
	fSaveFileName = "";
	fInputFile = "";
	fSaveFile = 0x0;
	fModifyFile = false;

	fEventHeader_branch = 0x0;
	fTruthInfo_branch = 0x0;
	fRecoSummary_branch = 0x0;
	fPidInfo_branch = 0x0;
	fSeedInfo_branch = 0x0;
	fStageInfo_branch = 0x0;
	fHitComparison_branch = 0x0;
}

WCSimOutputTree::~WCSimOutputTree() {
	if (fEventHeader != 0x0) {
		std::cout << "delete Event Header" << std::endl;
		delete fEventHeader;
		fEventHeader = 0x0;
		std::cout << "done" << std::endl;
	}

	if (fTruthInfo != 0x0) {
		std::cout << "delete truth info" << std::endl;
		delete fTruthInfo;
		fTruthInfo = 0x0;
		std::cout << "done" << std::endl;
	}

	if (fRecoSummary != 0x0) {
		std::cout << "delete reco summary" << std::endl;
		delete fRecoSummary;
		std::cout << "done" << std::endl;
	}

	if (fPidInfo != 0x0) {
		std::cout << "delete pid info" << std::endl;
		delete fPidInfo;
		fPidInfo = 0x0;
		std::cout << "done" << std::endl;
	}

	if (fSeedInfo != 0x0) {
		std::cout << "delete seed info" << std::endl;
		delete fSeedInfo;
		fSeedInfo = 0x0;
		std::cout << "done" << std::endl;
	}

	if (fStageInfo != 0x0) {
		std::cout << "delete stage info" << std::endl;
		delete fStageInfo;
		fStageInfo = 0x0;
		std::cout << "done" << std::endl;
	}

	if (fHitComparison != 0x0) {
		std::cout << "delete hit comp" << std::endl;
		delete fHitComparison;
		fHitComparison = 0x0;
		std::cout << "done" << std::endl;
	}
}

TString WCSimOutputTree::FormTime() {
	TTimeStamp ts;
	unsigned int year, month, day, hour, minute, second;
	ts.GetDate(true, 0, &year, &month, &day);
	ts.GetTime(true, 0, &hour, &minute, &second);
	TString time = Form("%02d%02d%02d", hour, minute, second);
	return time;
}

std::string WCSimOutputTree::FitType(WCSimFitterConfig * config) {
	if (config->GetIsCosmicFit()) {
		return "cosmic";
	} else if (config->GetIsPiZeroFit()) {
		return "piZero";
	} else {
		return TrackType::AsString(config->GetTrackType(0));
	}
}

TString WCSimOutputTree::InitOutputTree(TString fileName, bool modifyInput,
		WCSimFitterConfig * fitterConfig) {
	std::cout << "*** WCSimOutputTree::InitOutputTree() ***" << std::endl;
	// This will require checks to be made that the events are matching as well as the EventHeader and TruthInfo.
	// TODO should be able to determine what type of file it is automatically
	fInputFile = fileName;
	fModifyFile = modifyInput;

	assert(fSaveFile == NULL);
	DeletePointersIfExist();

	// Initialise the input and output names, open the file and get tree...
	TString time = FormTime();
	if (fModifyFile) {
		// TODO We need to update the time component of the name
		std::cout << "Modifying existing file..." << std::endl;
		fSaveFileName = fileName;

		TDirectory * tmpd = gDirectory;
		fSaveFile = new TFile(fSaveFileName.Data(), "update");
		std::cout << "Output to be saved in: " << fSaveFileName.Data()
				<< std::endl;
		tmpd->cd();

		GetExistingTree();

	} else {
		// The input file is a WCSim file and a new output needs to be created
		std::cout << "Creating new file..." << std::endl;
		TString basename = gSystem->BaseName(fileName.Data());
		basename.ReplaceAll(".root", "");

		TString saveNameTree = Form("fit_%s_%04d_to_%04d_%s_tree.root",
				basename.Data(), fitterConfig->GetFirstEventToFit(),
				fitterConfig->GetNumEventsToFit()
						+ fitterConfig->GetFirstEventToFit(), time.Data());

		fSaveFileName = saveNameTree;

		TDirectory * tmpd = gDirectory;
		fSaveFile = new TFile(fSaveFileName.Data(), "CREATE");
		std::cout << "Output to be saved in: " << fSaveFileName.Data()
				<< std::endl;
		tmpd->cd();

		MakeNewTree();
	}

	MakeNewFit(FitType(fitterConfig));
	return fSaveFileName;
}

void WCSimOutputTree::GetExistingTree() {
	std::cout << "*** WCSimOutputTree::GetExistingTree() ***" << std::endl;
	fSaveFile->cd();
	assert(fSaveFile->GetListOfKeys()->Contains("fResultsTree")); // Check it already has a results tree

	fResultsTree = (TTree*) fSaveFile->Get("fResultsTree");
}

void WCSimOutputTree::MakeNewTree() {
	std::cout << "*** WCSimOutputTree::MakeNewTree() ***" << std::endl;
	fSaveFile->cd();

	// Set up the Geometry tree
	fGeoTree = new TTree("fGeoTree", "GeoTree");
	fGeometry = new WCSimRootGeom();
	fGeoTree->Branch("wcsimrootgeom", "WCSimRootGeom", &fGeometry, 64000, 0);

	// Set up the Parameter tree
	if (WCSimParameters::Instance()->SaveParameterInfo()) {
		fParameterTree = new TTree("fParameterTree", "ParameterTree");
		fParameterInfo = new ParameterInfo();
		fParameterTree->Branch("parameterTree", "ParameterInfo", &fParameterInfo, 64000, 0);
	}

	// Set up the results tree
	fResultsTree = new TTree("fResultsTree", "ResultsTree");

	if (fEventHeader == 0x0) {
		fEventHeader = new EventHeader();
	}
	if (fTruthInfo == 0x0) {
		fTruthInfo = new TruthInfo();
	}

	fEventHeader_branch = fResultsTree->Branch("EventHeader", "EventHeader",
			&fEventHeader, 64000, 1);
	fTruthInfo_branch = fResultsTree->Branch("TruthInfo", "TruthInfo",
			&fTruthInfo, 64000, 1);
}

void WCSimOutputTree::MakeNewFit(std::string type) {
	std::cout << "*** WCSimOutputTree::MakeNewFit() ***" << std::endl;
	fSaveFile->cd();

	if (fRecoSummary == 0x0) {
		fRecoSummary = new WCSimRecoSummary();
	}
	if (fPidInfo == 0x0) {
		fPidInfo = new PidInfo();
	}
	if (fSeedInfo == 0x0) {
		fSeedInfo = new SeedInfo();
	}
	if (fStageInfo == 0x0) {
		fStageInfo = new StageInfo();
	}
	if (fHitComparison == 0x0) {
		fHitComparison = new WCSimHitComparison();
	}

	fRecoSummary_branch = fResultsTree->Branch(("RecoSummary_" + type).c_str(),
			"WCSimRecoSummary", &fRecoSummary, 64000, 0);
	fPidInfo_branch = fResultsTree->Branch(("PidInfo_" + type).c_str(),
			"PidInfo", &fPidInfo, 64000, 1);

	// TODO need to add something even if it's empty!
	if (WCSimParameters::Instance()->SaveSeedInfo()) {
		std::cout << "WCSimOutputTree::MakeTree()...Saving Seed Info"
				<< std::endl;
		fSeedInfo_branch = fResultsTree->Branch(("SeedInfo_" + type).c_str(),
				"SeedInfo", &fSeedInfo, 64000, 1);
	}
	if (WCSimParameters::Instance()->SaveStageInfo()) {
		std::cout << "WCSimOutputTree::MakeTree()...Saving Stage Info"
				<< std::endl;
		fStageInfo_branch = fResultsTree->Branch(("StageInfo_" + type).c_str(),
				"StageInfo", &fStageInfo, 64000, 1);
	}
	if (WCSimParameters::Instance()->SaveHitComparison()) {
		std::cout << "WCSimOutputTree::MakeTree()...Saving Hit Comparison"
				<< std::endl;
		fHitComparison_branch = fResultsTree->Branch(
				("HitComparison_" + type).c_str(), "WCSimHitComparison",
				&fHitComparison, 64000, 1);
	}
}

void WCSimOutputTree::FillSeedInfo(SeedInfo& seedInfo) {
	assert(fResultsTree != 0x0);

	delete fSeedInfo;
	fSeedInfo = new SeedInfo(seedInfo);
	fSeedInfo_branch->Fill();
}

void WCSimOutputTree::FillStageInfo(int stageNCalls,
		std::vector<WCSimHitPrediction> stagePreds,
		std::vector<WCSimLikelihoodTrackBase*> stageTracks) {
	assert(fResultsTree != 0x0);

	fStageInfo->SetStageInfo(stageNCalls, stagePreds, stageTracks);
	fStageInfo_branch->Fill();
}

void WCSimOutputTree::FillHitComparison(WCSimHitComparison& hitComp) {
	assert(fResultsTree != 0x0);

	delete fHitComparison;
	fHitComparison = new WCSimHitComparison(hitComp);
	fHitComparison_branch->Fill();
}

void WCSimOutputTree::Fill(EventHeader& eventHead, TruthInfo& truthInfo,
		WCSimRecoSummary& summ, PidInfo& pidInfo) {
	std::cout << " *** WCSimFitter::Fill() *** " << std::endl;
	// Fill the truth and fit tree entries
	assert(fResultsTree != 0x0);

	if (!fModifyFile) {
		delete fEventHeader;
		fEventHeader = new EventHeader(eventHead);
		fEventHeader_branch->Fill();

		delete fTruthInfo;
		fTruthInfo = new TruthInfo(truthInfo);
		fTruthInfo_branch->Fill();

		// Fill the geometry tree
		if (fGeoTree->GetEntries() == 0) {
			*fGeometry = *(WCSimGeometry::Instance()->GetWCSimGeometry());
			WCSimGeometry::Instance()->PrintGeometry();
			fGeoTree->Fill();
		}

		// Fill the parameter tree
		if (fParameterTree->GetEntries() == 0 && WCSimParameters::Instance()->SaveParameterInfo()) {
			if (fParameterInfo == 0x0) {
				fParameterInfo = new ParameterInfo();
			}
			fParameterInfo->Print();
			fParameterTree->Fill();
		}
	}

	delete fRecoSummary;
	fRecoSummary = new WCSimRecoSummary(summ);
	fRecoSummary_branch->Fill();

	delete fPidInfo;
	fPidInfo = new PidInfo(pidInfo);
	fPidInfo_branch->Fill();
}

void WCSimOutputTree::SaveTree() {
	std::cout << " *** WCSimOutputTree::SaveTree() *** " << std::endl;
	fResultsTree->SetEntries(fRecoSummary_branch->GetEntries());

	TDirectory* tmpd = 0;
	tmpd = gDirectory;
	fSaveFile->cd();

	assert(fSaveFile->IsOpen());

	fResultsTree->Write("", TObject::kOverwrite);

	if (!fModifyFile) {
		fGeoTree->Write("", TObject::kOverwrite);
		if (WCSimParameters::Instance()->SaveParameterInfo()) {
			fParameterTree->Write("", TObject::kOverwrite);
		}
	}

	fSaveFile->Close();
	delete fSaveFile;

	tmpd->cd();
}

void WCSimOutputTree::SaveFile() {
	TDirectory* tmpd = 0;
	tmpd = gDirectory;
	fSaveFile->cd();
	std::cout << " *** WCSimFitter::SaveFile() *** " << std::endl;
	std::cout << "Save file name = " << fSaveFileName << std::endl;

	fSaveFile->Close();
	delete fSaveFile;
}

void WCSimOutputTree::SetSaveFileName(TString saveName) {
	fSaveFileName = saveName;
	return;
}

TString WCSimOutputTree::GetSaveFileName() const {
	return fSaveFileName;
}

void WCSimOutputTree::SetModifyFile(bool modify) {
	fModifyFile = modify;
	return;
}

void WCSimOutputTree::SetInputFile(TString inputFile) {
	fInputFile = inputFile;
	return;
}

TString WCSimOutputTree::GetInputFileName() const {
	return fInputFile;
}

void WCSimOutputTree::DeletePointersIfExist() {
	if (fResultsTree != 0x0) {
		delete fResultsTree;
		fResultsTree = 0x0;
	}

	if (fGeoTree != 0x0) {
		delete fGeoTree;
		fGeoTree = 0x0;
	}

	if (fParameterTree != 0x0) {
		delete fParameterTree;
		fParameterTree = 0x0;
	}

	if (fEventHeader != 0x0) {
		delete fEventHeader;
		fEventHeader = 0x0;
	}

	if (fTruthInfo != 0x0) {
		delete fTruthInfo;
		fTruthInfo = 0x0;
	}

	if (fRecoSummary != 0x0) {
		delete fRecoSummary;
		fRecoSummary = 0x0;
	}

	if (fPidInfo != 0x0) {
		delete fPidInfo;
		fPidInfo = 0x0;
	}

	if (fSeedInfo != 0x0) {
		delete fSeedInfo;
		fSeedInfo = 0x0;
	}

	if (fStageInfo != 0x0) {
		delete fStageInfo;
		fStageInfo = 0x0;
	}

	if (fHitComparison != 0x0) {
		delete fHitComparison;
		fHitComparison = 0x0;
	}
	return;
}
