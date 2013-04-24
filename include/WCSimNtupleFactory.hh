#ifndef WCSIMNTUPLEFACTORY_HH
#define WCSIMNTUPLEFACTORY_HH

#include "TObject.h"

class WCSimNtuple;

class WCSimNtupleFactory : public TObject
{
 public:

  static WCSimNtupleFactory* Instance();

  WCSimNtuple* MakeNtuple(const char* name = "default");

 private:

  WCSimNtupleFactory();
  ~WCSimNtupleFactory();

  WCSimNtuple* fRecoNtuple;
  WCSimNtuple* fVertexNtuple;
  WCSimNtuple* fVertexSeedNtuple;

  ClassDef(WCSimNtupleFactory,0)
};

#endif 
