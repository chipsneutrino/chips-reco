#ifndef WCSIMDISPLAYVIEWER_HH
#define WCSIMDISPLAYVIEWER_HH

#include "TObject.h"
#include "TCanvas.h"
#include "TButton.h"
#include "TChain.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"

#include <vector>

class WCSimDisplay;
class WCSimReco;

class WCSimRecoEvent;
class WCSimTrueEvent;

class WCSimRootTrigger;
class WCSimRootEvent;
class WCSimRootGeom;

class WCSimDisplayViewer : public TObject
{
 public:

  WCSimDisplayViewer();
  virtual ~WCSimDisplayViewer();

  void RootStyle();

  void Initialize();

  void controlStart();
  void controlNext();
  void controlPrevious();
  void controlGoTo();
  void controlPrint();

  void BuildDisplay();
  void UpdateDisplay(Int_t entry);

  void DisplayEvent(Int_t entry);

  void DisplayRecoClusters(Bool_t yesno = 1){ fDisplayRecoClusters = yesno; }
  void DisplayRecoEvent(Bool_t yesno = 1){ fDisplayRecoEvent = yesno; }
  void DisplayTrueEvent(Bool_t yesno = 1){ fDisplayTrueEvent = yesno; }

  void SetPulseHeightCut(Double_t cut);

  void UseDisplay(const char* display);

  void PrintGIF(Bool_t yesno = 1);
  void PrintEPS(Bool_t yesno = 1);

  void PrintDisplayList();

 private:

  TStyle* fRootStyle;

  TCanvas* fControlCanvas;
  TButton* buttonStart;
  TButton* buttonNext;
  TButton* buttonPrevious;
  TButton* buttonGoTo;
  TButton* buttonPrint;
  TButton* buttonQuit;

  Int_t fCounter;

  WCSimReco* fReco;
  WCSimRecoEvent* fRecoEvent;
  WCSimTrueEvent* fTrueEvent;

  Int_t fRunNum;
  Int_t fEventNum;
  Int_t fTriggerNum;

  Bool_t fDisplayRecoClusters;
  Bool_t fDisplayRecoEvent;
  Bool_t fDisplayTrueEvent;

  Double_t fPulseHeightCut;

  Bool_t fPrintGIF;
  Bool_t fPrintEPS;

  Bool_t fDisplayReady;

  std::vector<const char*> fDisplayNameList;
  std::vector<WCSimDisplay*> fDisplayList;

  ClassDef(WCSimDisplayViewer,0)
};

#endif 

