#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class WCSimGeometry+;
#pragma link C++ class WCSimInterface+;
#pragma link C++ class WCSimParameters+;

#pragma link C++ class WCSimTrueEvent+;
#pragma link C++ class WCSimTrueTrack+;

#pragma link C++ class WCSimRecoObjectTable+;
#pragma link C++ class WCSimRecoSlicer+;
#pragma link C++ class WCSimRecoSeed+;
#pragma link C++ class WCSimRecoEvent+;
#pragma link C++ class WCSimRecoDigit+;
#pragma link C++ class WCSimRecoCluster+;
#pragma link C++ class WCSimRecoClusterDigit+;
#pragma link C++ class WCSimRecoClusteringUtil+;
#pragma link C++ class WCSimRecoRing+;
#pragma link C++ class WCSimRecoVertex+;

#pragma link C++ class WCSimCosmicSeed+;
#pragma link C++ class WCSimCosmicFitter+;

#pragma link C++ class WCSimVertexFinder+;
#pragma link C++ class WCSimRingFinder+;
#pragma link C++ class WCSimDataCleaner+;

#pragma link C++ class WCSimVertexGeometry+;
#pragma link C++ class WCSimHoughTransform+;
#pragma link C++ class WCSimHoughTransformArray+;

#pragma link C++ class WCSimMsg+;

#pragma link C++ class WCSimTotalLikelihood+;
#pragma link C++ class WCSimLikelihoodFitter+;
#pragma link C++ class WCSimLikelihoodTuner+;
#pragma link C++ class WCSimLikelihoodTunerCache+;
#pragma link C++ class WCSimChargePredictor+;
#pragma link C++ class WCSimLikelihoodDigit+;
#pragma link C++ class vector<WCSimLikelihoodDigit>+;
#pragma link C++ class vector<vector<WCSimLikelihoodDigit> >;
#pragma link C++ class WCSimLikelihoodDigitArray+;
#pragma link C++ class WCSimLikelihoodTrack+;
#pragma link C++ class TrackType+;
#pragma link C++ enum TrackType::Type;
#pragma link C++ class FitterParameterType+;
#pragma link C++ enum FitterParameterType::Type;
#pragma link C++ class WCSimDigitizerLikelihood+;
#pragma link C++ class vector<WCSimLikelihoodTrack>;
#pragma link C++ class vector<vector<WCSimLikelihoodTrack> >;
#pragma link C++ class vector<WCSimLikelihoodTrack*>;
#pragma link C++ class vector<WCSimChargePredictor>;

#pragma link C++ class WCSimFitterConfig+;
#pragma link C++ class WCSimFitterPlots+;
#pragma link C++ class WCSimFitterInterface+;
#pragma link C++ class WCSimFitterParameter+;
#pragma link C++ class WCSimFitterSingleTrackParameters+;
#pragma link C++ class WCSimFitterParameters+;
#pragma link C++ class vector<WCSimFitterSingleTrackParameters>+;
#pragma link C++ class WCSimEmissionProfiles+;
#pragma link C++ enum EmissionProfile_t::Type;
#pragma link C++ class WCSimEmissionProfileManager+;
#pragma link C++ class WCSimRecoSummary+;
#pragma link C++ class WCSimRecoEvDisplay+;
#pragma link C++ class WCSimIntegralLookupMaker3D+;
#pragma link C++ class WCSimIntegralLookupHistArray+;
#pragma link C++ class WCSimIntegralLookupReader+;
#pragma link C++ class WCSimFitterTrackParMap+;
#pragma link C++ class WCSimLikelihoodTrackBase+;
#pragma link C++ class vector<WCSimLikelihoodTrackBase*>;
#pragma link C++ class WCSimLikelihoodPhotonTrack+;
#pragma link C++ class vector<WCSimLikelihoodPhotonTrack>;
#pragma link C++ class vector<vector<WCSimLikelihoodPhotonTrack> >;
#pragma link C++ class vector<WCSimLikelihoodPhotonTrack*>;
#pragma link C++ class WCSimLikelihoodUnknownTrack+;
#pragma link C++ class vector<WCSimLikelihoodUnknownTrack>;
#pragma link C++ class vector<vector<WCSimLikelihoodUnknownTrack> >;
#pragma link C++ class vector<WCSimLikelihoodUnknownTrack*>;
#pragma link C++ class WCSimLikelihoodTrackFactory+;
#pragma link C++ class WCSimTimePredictor+;
#pragma link C++ class WCSimDetectorParameters+;
#pragma link C++ class WCSimPiZeroFitter+;
#pragma link C++ class WCSimPiZeroSeeder+;
#pragma link C++ class WCSimPiZeroSeedGenerator+;
#pragma link C++ class WCSimPiZeroSeed+;
#pragma link C++ class WCSimPiZeroSingleElectronSeeder+;
#pragma link C++ class WCSimPiZeroHoughSeeder+;
#pragma link C++ class WCSimPiZeroElectronAdjuster+;
#pragma link C++ class WCSimTransmissionFunctionLookupTable+;
#pragma link C++ class vector<WCSimTransmissionFunctionLookupTable>+;
#pragma link C++ class WCSimTransmissionFunctionLookup+;
#pragma link C++ class WCSimTimeLikelihood3+;
#pragma link C++ class WCSimIntegralLookup3D+;

#pragma link C++ class WCSimHitPrediction+;
#pragma link C++ class vector<WCSimHitPrediction>+;
#pragma link C++ class vector<vector<WCSimHitPrediction> >;
#pragma link C++ class WCSimSingleHitComparison+;
#pragma link C++ class vector<WCSimSingleHitComparison>+;
#pragma link C++ class EventHeader+;
#pragma link C++ class ParameterInfo+;
#pragma link C++ class TruthInfo+;
#pragma link C++ class PidInfo+;
#pragma link C++ class SeedInfo+;
#pragma link C++ class StageInfo+;
#pragma link C++ class WCSimOutputTree+;

#pragma link C++ class WCSimPIDTree+;
#pragma link C++ class WCSimPID+;

#pragma link C++ class WCSimDigitizerPDFMaker+;

#endif
