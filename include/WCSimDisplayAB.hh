#ifndef WCSIMDISPLAYAB_HH
#define WCSIMDISPLAYAB_HH

#include "WCSimDisplay.hh"

#include <vector>

#include "TString.h"
#include "TLatex.h"

class TCanvas;
class TPad;
class TH1D;
class TH2D;

class TLine;
class TBox;
class TEllipse;
class TMarker;
class TPolyMarker;
class TPolyLine;
class TLegend;
class WCSimRootEvent;

class WCSimDisplayAB : public WCSimDisplay {

 public:
  WCSimDisplayAB();
  ~WCSimDisplayAB();
  
  void DrawDisplay(WCSimRecoEvent* recoevent);
  void DrawCleanDisplay(WCSimRecoEvent* recoevent);

  void DrawRecoEvent(WCSimRecoEvent* recoevent);
  void DrawTrueEvent(WCSimRootEvent* rootevent);

  void ResetDisplay();
  void PrintDisplay();

 private:

  void Initialize();
  void Reset();

  void BuildGeometry();
  
  void DrawNewDisplay();
  void UpdateDisplay();

  void PrintDisplayGIF();
  void PrintDisplayEPS();

  Int_t QPEtoCOL(Double_t Q);
  Int_t TIMEtoCOL(Double_t T);

  void MakeLegendQPE();
  void MakeLegendTIME();
  void MakeLegendRECO();

  void DrawVertex(Double_t vx, Double_t vy, Double_t vz,
                  Double_t px, Double_t py, Double_t pz,
                  Int_t colourcode);

  void DrawRing(Double_t vx, Double_t vy, Double_t vz,
                Double_t px, Double_t py, Double_t pz,
                Double_t angle, Int_t colourcode);

  // Geometry
  Int_t fGeoType;

  Double_t fCylRadius;
  Double_t fCylLength;
  Double_t fCylDiagonal;

  Double_t fMailBoxX;
  Double_t fMailBoxY;
  Double_t fMailBoxZ;
  Double_t fMailBoxDiagonal;

  Double_t fScale;


  // Event Info
  Int_t fRunNumber;
  Int_t fEventNumber;
  Int_t fTriggerNumber;


  // Spatial Display
  TCanvas* wcCanvas;
  TH2D* wcDisplay;
  TLegend* wcLegend;

  Double_t fU;
  Double_t fV;

  Double_t canvasWidth;
  Double_t canvasHeight;
  Int_t canvasU;
  Int_t canvasV;

  Double_t binsWidth;
  Int_t binsU;
  Int_t binsV;

  TBox* wcCylEdgeSide;
  TEllipse* wcCylEdgeTop;
  TEllipse* wcCylEdgeBottom;  

  TBox* wcBoxEdgeSide;
  TBox* wcBoxEdgeTop;
  TBox* wcBoxEdgeBottom;
  TLine* wcBoxEdgeLine1;
  TLine* wcBoxEdgeLine2;
  TLine* wcBoxEdgeLine3;
  TLine* wcBoxEdgeLine4;

  TString wcTitleString;
  TLatex wcTitleLatex;

  Double_t titleU;
  Double_t titleV;


  // Timing Display
  Bool_t fMakeTimeDisplay;

  TCanvas* wcTimeCanvas;
  TLegend* wcTimeLegend;


  // Reco/True
  TLegend* wcRecoLegend;

  // Stores for markers
  std::vector<TMarker*> wcDigits;
  std::vector<TMarker*> wcTimeDigits;
  std::vector<TMarker*> wcCommonDigits;
  std::vector<TPolyMarker*> wcRings;

  // Markers for legends (better way?)
  std::vector<TMarker*> legDigits;

  // Lines for legends (better way?)
  std::vector<TLine*> legRings;

  ClassDef(WCSimDisplayAB,0)

};

#endif
