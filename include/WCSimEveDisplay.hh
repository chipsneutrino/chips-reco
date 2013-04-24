#ifndef WCSIMEVEDISPLAY_HH
#define WCSIMEVEDISPLAY_HH

#include "WCSimDisplay.hh"

class WCSimEveDisplay : public WCSimDisplay {

 public:
  WCSimEveDisplay();
  ~WCSimEveDisplay();

  void DrawDisplay(WCSimRecoEvent* recoevent);
  void DrawCleanDisplay(WCSimRecoEvent* recoevent);

  // nothing else implemented
  void DrawRecoEvent(WCSimRecoEvent*) ;
  void DrawTrueEvent(WCSimTrueEvent*);
  void ResetDisplay();        
  void PrintDisplay();

private:

  void Initialize();
  void BuildGeometry();
 

  ClassDef(WCSimEveDisplay,0)

};

#endif

