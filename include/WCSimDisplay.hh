#ifndef WCSIMDISPLAY_HH
#define WCSIMDISPLAY_HH

#include "TObject.h"

#include "WCSimRecoEvent.hh"
#include "WCSimTrueEvent.hh"
#include "WCSimRootEvent.hh"

class WCSimDisplay : public TObject {

 public:  
  WCSimDisplay();
  virtual ~WCSimDisplay();

  virtual void DrawDisplay(WCSimRecoEvent* recoevent) { }
  virtual void DrawCleanDisplay(WCSimRecoEvent* recoevent) { }

  virtual void DrawRecoEvent(WCSimRecoEvent* recoevent) { }
  virtual void DrawTrueEvent(WCSimRootEvent* rootevent) { }

  virtual void ResetDisplay() { }
  virtual void PrintDisplay() { } 

  void SetPulseHeightCut(Double_t cut) { fPulseHeightCut = cut; }
  Double_t GetPulseHeightCut() { return fPulseHeightCut; }

  void PrintGIF(Bool_t yesno=1) { fPrintGIF = yesno; }
  void PrintEPS(Bool_t yesno=1) { fPrintEPS = yesno; }

 protected:

  Double_t fPulseHeightCut;

  Bool_t fPrintGIF;
  Bool_t fPrintEPS;

 private:

  ClassDef(WCSimDisplay,0)

};

#endif
