#ifndef WCSIMRINGVIEWER_HH
#define WCSIMRINGVIEWER_HH

#include "TObject.h"

class TCanvas;
class TPad;
class TH1D;
class TH2D;
class TMarker;
class TLegend;

class WCSimRecoEvent;

class WCSimRingViewer : public TObject {

 public:
  WCSimRingViewer();
  ~WCSimRingViewer();

  void Display(WCSimRecoEvent* evt);

 private:

  void Initialize();

  TCanvas* fCanvas;

  ClassDef(WCSimRingViewer,0)

};

#endif

