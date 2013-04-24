#ifndef WCSIMDISPLAYFACTORY_HH
#define WCSIMDISPLAYFACTORY_HH

#include "TObject.h"

class WCSimDisplay;

class WCSimDisplayFactory : public TObject
{
 public:

  static WCSimDisplayFactory* Instance();

  WCSimDisplay* MakeDisplay(const char* name = "default");

 private:

  WCSimDisplayFactory();
  ~WCSimDisplayFactory();

  WCSimDisplay* fDisplayAB;
  WCSimDisplay* fDisplayEve;
  WCSimDisplay* fDisplayVtx;

  ClassDef(WCSimDisplayFactory,0)
};

#endif 
