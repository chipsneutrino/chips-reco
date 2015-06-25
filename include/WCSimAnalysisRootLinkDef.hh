#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class WCSimDisplay+;
#pragma link C++ class WCSimDisplayFactory+;
#pragma link C++ class WCSimDisplayViewer+;

#pragma link C++ class WCSimDisplayAB+;
#pragma link C++ class WCSimEveDisplay+;
#pragma link C++ class WCSimVertexViewer+;
#pragma link C++ class WCSimRingViewer+;

#pragma link C++ class WCSimEventWriter+;
#pragma link C++ class WCSimGeometry+;
#pragma link C++ class WCSimInterface+;
#pragma link C++ class WCSimParameters+;

#pragma link C++ class WCSimTrueEvent+;
#pragma link C++ class WCSimTrueTrack+;

#pragma link C++ class WCSimRecoObjectTable+;
#pragma link C++ class WCSimRecoFactory+;
#pragma link C++ class WCSimReco+;
#pragma link C++ class WCSimRecoAB+;
#pragma link C++ class WCSimRecoEvent+;
#pragma link C++ class WCSimRecoDigit+;
#pragma link C++ class WCSimRecoCluster+;
#pragma link C++ class WCSimRecoClusterDigit+;
#pragma link C++ class WCSimRecoRing+;
#pragma link C++ class WCSimRecoVertex+;

#pragma link C++ class WCSimVertexFinder+;
#pragma link C++ class WCSimRingFinder+;
#pragma link C++ class WCSimDataCleaner+;

#pragma link C++ class WCSimVertexGeometry+;
#pragma link C++ class WCSimHoughTransform+;
#pragma link C++ class WCSimHoughTransformArray+;

#pragma link C++ class WCSimNtupleFactory+;
#pragma link C++ class WCSimNtupleWriter+;

#pragma link C++ class WCSimNtuple+;
#pragma link C++ class WCSimRecoNtuple+;
#pragma link C++ class WCSimVertexNtuple+;
#pragma link C++ class WCSimVertexSeedNtuple+;

#pragma link C++ class WCSimMsg+;

#pragma link C++ class WCSimTotalLikelihood+;
#pragma link C++ class WCSimLikelihoodFitter+;
#pragma link C++ class WCSimLikelihoodTuner+;
#pragma link C++ class WCSimChargePredictor+;
#pragma link C++ class WCSimLikelihoodDigit+;
#pragma link C++ class WCSimLikelihoodDigitArray+;
#pragma link C++ class WCSimLikelihoodTrack+;
#pragma link C++ class TrackType+;
#pragma link C++ enum TrackType::Type;
#pragma link C++ class FitterParameterType+;
#pragma link C++ enum FitterParameterType::Type;
#pragma link C++ class WCSimTimeLikelihood+;
#pragma link C++ class WCSimDigitizerLikelihood+;
#pragma link C++ class vector<WCSimLikelihoodTrack>;
#pragma link C++ class vector<vector<WCSimLikelihoodTrack> >;
#pragma link C++ class vector<WCSimLikelihoodTrack*>;
#pragma link C++ class vector<WCSimChargePredictor>;

#pragma link C++ class WCSimAnalysisConfig+;
#pragma link C++ class WCSimFitterConfig+;
#pragma link C++ class WCSimFitterPlots+;
#pragma link C++ class WCSimFitterTree+;
#pragma link C++ class WCSimFitterInterface+;
#pragma link C++ class WCSimFitterParameter+;
#pragma link C++ class WCSimFitterSingleTrackParameters+;
#pragma link C++ class WCSimFitterParameters+;
#pragma link C++ class vector<WCSimFitterSingleTrackParameters>+;
#pragma link C++ class WCSimEmissionProfiles+;
#pragma link C++ enum EmissionProfile_t::Type;
#pragma link C++ class WCSimRecoSummary+;
#pragma link C++ class WCSimRecoEvDisplay+;
#pragma link C++ class WCSimIntegralLookupMaker+;
#pragma link C++ class WCSimIntegralLookupMaker3D+;
#pragma link C++ class WCSimIntegralLookupReader+;
#pragma link C++ class WCSimFitterTrackParMap+;
#pragma link C++ class WCSimLikelihoodTrackBase+;
#pragma link C++ class WCSimLikelihoodPhotonTrack+;
#pragma link C++ class WCSimLikelihoodTrackFactory+;

#endif
