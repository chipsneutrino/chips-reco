#ifndef WCSIMNTUPLEWRITER_HH
#define WCSIMNTUPLEWRITER_HH

#include "TObject.h"
#include "TString.h"

class WCSimNtuple;
class WCSimReco;

class WCSimRecoEvent;
class WCSimTrueEvent;

class WCSimNtupleWriter : public TObject {

 public:
  WCSimNtupleWriter();
  ~WCSimNtupleWriter();

  void Initialize();

  void UseNtuple( const char* ntuple ) { fNtupleName = ntuple; }
  void SetFileName( const char* filename ) { fNtpFileName = filename; }
 
  void Run( Int_t entries = -1 );
  void Run( Int_t i1, Int_t i2 );

 private:

  void BuildNtuple();

  WCSimReco* fReco;
  WCSimRecoEvent* fRecoEvent;
  WCSimTrueEvent* fTrueEvent;

  WCSimNtuple* fNtuple;

  TString fNtupleName;
  TString fNtpFileName;

  ClassDef(WCSimNtupleWriter,0)
};

#endif
