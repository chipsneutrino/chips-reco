/********************************************************************
* ./src/WCSimAnalysisRootDict.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error ./src/WCSimAnalysisRootDict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtableWCSimAnalysisRootDict();
extern void G__cpp_setup_inheritanceWCSimAnalysisRootDict();
extern void G__cpp_setup_typetableWCSimAnalysisRootDict();
extern void G__cpp_setup_memvarWCSimAnalysisRootDict();
extern void G__cpp_setup_globalWCSimAnalysisRootDict();
extern void G__cpp_setup_memfuncWCSimAnalysisRootDict();
extern void G__cpp_setup_funcWCSimAnalysisRootDict();
extern void G__set_cpp_environmentWCSimAnalysisRootDict();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "WCSimTimeLikelihood.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimChargeLikelihood.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimDisplayViewer.hh"
#include "WCSimDisplayFactory.hh"
#include "WCSimDisplay.hh"
#include "WCSimDisplayAB.hh"
#include "WCSimEveDisplay.hh"
#include "WCSimEventWriter.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimParameters.hh"
#include "WCSimRecoObjectTable.hh"
#include "WCSimRecoFactory.hh"
#include "WCSimReco.hh"
#include "WCSimRecoAB.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoCluster.hh"
#include "WCSimRecoClusterDigit.hh"
#include "WCSimRecoRing.hh"
#include "WCSimRecoVertex.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimTrueEvent.hh"
#include "WCSimTrueTrack.hh"
#include "WCSimHoughTransform.hh"
#include "WCSimHoughTransformArray.hh"
#include "WCSimDataCleaner.hh"
#include "WCSimVertexFinder.hh"
#include "WCSimVertexGeometry.hh"
#include "WCSimVertexViewer.hh"
#include "WCSimRingFinder.hh"
#include "WCSimRingViewer.hh"
#include "WCSimNtupleFactory.hh"
#include "WCSimNtuple.hh"
#include "WCSimRecoNtuple.hh"
#include "WCSimVertexNtuple.hh"
#include "WCSimVertexSeedNtuple.hh"
#include "WCSimNtupleWriter.hh"
#include "WCSimMsg.hh"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TClass;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TBuffer;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMemberInspector;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TObject;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TString;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEfloatcOallocatorlEfloatgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_string;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TObjArray;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TClonesArray;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoDigit;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRootCherenkovDigiHit;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRootTrigger;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEintcOallocatorlEintgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRootEvent;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimLikelihoodDigit;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimLikelihoodDigitArray;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TElementActionTlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TElementPosActionTlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTRow_constlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTRowlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTDiag_constlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTColumn_constlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTFlat_constlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTSub_constlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTSparseRow_constlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTSparseDiag_constlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTColumnlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTDiaglEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTFlatlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTSublEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTSparseRowlEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMatrixTSparseDiaglEfloatgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TVector3;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimLikelihoodTrack;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimLikelihoodTrackcLcLTrackType;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoVertex;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoRing;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoEvent;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimRecoDigitmUcOallocatorlEWCSimRecoDigitmUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimRecoDigitmUcOallocatorlEWCSimRecoDigitmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimRecoRingmUcOallocatorlEWCSimRecoRingmUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimRecoRingmUcOallocatorlEWCSimRecoRingmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimTimeLikelihood;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TFile;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimLikelihoodFitter;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRootGeom;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimLikelihoodTuner;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimChargeLikelihood;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TCanvas;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TBox;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TLegend;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TButton;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TTree;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TChain;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TStyle;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimDisplay;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimReco;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimTrueEvent;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimDisplayViewer;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEconstsPcharmUcOallocatorlEconstsPcharmUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEconstsPcharmUcOallocatorlEconstsPcharmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimDisplaymUcOallocatorlEWCSimDisplaymUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimDisplaymUcOallocatorlEWCSimDisplaymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimDisplayFactory;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimTrueTrack;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimTrueTrackmUcOallocatorlEWCSimTrueTrackmUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimTrueTrackmUcOallocatorlEWCSimTrueTrackmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TLatex;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TH1D;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TH2D;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TLine;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TEllipse;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMarker;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimDisplayAB;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlETMarkermUcOallocatorlETMarkermUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlETMarkermUcOallocatorlETMarkermUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlETPolyMarkermUcOallocatorlETPolyMarkermUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlETPolyMarkermUcOallocatorlETPolyMarkermUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlETLinemUcOallocatorlETLinemUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlETLinemUcOallocatorlETLinemUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimEveDisplay;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimEventWriter;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimGeometry;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimGeometrycLcLEGeoConfiguration;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimGeometrycLcLEGeoType;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimGeometrycLcLEGeoRegion;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimInterface;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimParameters;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoObjectTable;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoFactory;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoVertexcLcLEFitStatus;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoAB;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoCluster;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoClusterDigit;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimRecoClusterDigitmUcOallocatorlEWCSimRecoClusterDigitmUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimRecoClusterDigitmUcOallocatorlEWCSimRecoClusterDigitmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimHoughTransform;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimHoughTransformArray;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEvectorlEfloatcOallocatorlEfloatgRsPgRcOallocatorlEvectorlEfloatcOallocatorlEfloatgRsPgRsPgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEvectorlEfloatcOallocatorlEfloatgRsPgRcOallocatorlEvectorlEfloatcOallocatorlEfloatgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimHoughTransformmUcOallocatorlEWCSimHoughTransformmUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimHoughTransformmUcOallocatorlEWCSimHoughTransformmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimDataCleaner;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimDataCleanercLcLEFilterConfig;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimRecoClustermUcOallocatorlEWCSimRecoClustermUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimRecoClustermUcOallocatorlEWCSimRecoClustermUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TMinuit;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimVertexFinder;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_vectorlEWCSimRecoVertexmUcOallocatorlEWCSimRecoVertexmUgRsPgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_reverse_iteratorlEvectorlEWCSimRecoVertexmUcOallocatorlEWCSimRecoVertexmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimVertexGeometry;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_TF1;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimVertexViewer;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRingFinder;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRingViewer;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimNtuple;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimNtupleFactory;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimRecoNtuple;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimVertexNtuple;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimVertexSeedNtuple;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimNtupleWriter;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimMsg;
extern G__linked_taginfo G__WCSimAnalysisRootDictLN_WCSimMsgcLcLEMsgLevel;

/* STUB derived class for protected member access */
