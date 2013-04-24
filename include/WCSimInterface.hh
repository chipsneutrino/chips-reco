#ifndef WCSIMINTERFACE_HH
#define WCSIMINTERFACE_HH

#include "TObject.h"
#include "TChain.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"

#include <vector>

class WCSimRecoDigit;
class WCSimRecoEvent;

class WCSimTrueTrack;
class WCSimTrueEvent;

class WCSimInterface : public TObject {

 public: 
  static WCSimInterface* Instance();

  static Bool_t TouchData();
  static Bool_t TouchEvent();

  static void LoadData(const char* file);
  static void LoadEvent(Int_t ievent);
  static Int_t GetNumEvents();

  static WCSimRootTrigger* FilterTrigger(WCSimRootEvent* event);

  static WCSimTrueEvent* TrueEvent();
  static WCSimRecoEvent* RecoEvent();
  
  static WCSimRootEvent* WCSimEvent();
  static WCSimRootTrigger* WCSimTrigger();

  static Int_t GetRunNumber();
  static Int_t GetEventNumber();
  static Int_t GetTriggerNumber();

  static void Reset();

  static void SetEnergyThreshold(Double_t input_mev);
  static void SetRangeThreshold(Double_t input_cm);

  void SetTrueEnergyThreshold(Double_t ke) { fEnergyThreshold = ke; }
  void SetTrueRangeThreshold(Double_t ds) { fRangeThreshold = ds; }

  WCSimTrueEvent* GetTrueEvent(){ return fTrueEvent; }
  WCSimRecoEvent* GetRecoEvent(){ return fRecoEvent; }

  WCSimRootEvent* GetWCSimEvent(){ return fEvent; }
  WCSimRootTrigger* GetWCSimTrigger(){ return fTrigger; }

  WCSimRootEvent* GetWCSimEvent(Int_t ievent);
  WCSimRootTrigger* GetWCSimTrigger(Int_t ievent);

  void AddFile(const char* file);
  void BuildEvent(Int_t ievent);
  Bool_t CheckEvent();
  Int_t GetEntries();

  void ResetForNewSample();

 private:
  WCSimInterface();
  ~WCSimInterface();

  void BuildEvent(WCSimRootTrigger* trigger);

  void BuildTrueEvent(WCSimRootTrigger* trigger);
  void BuildRecoEvent(WCSimRootTrigger* trigger);

  void ResetTrueEvent();
  void ResetRecoEvent();

  WCSimTrueEvent* fTrueEvent;
  WCSimRecoEvent* fRecoEvent;

  std::vector<WCSimRecoDigit*>* fDigitList;
  std::vector<WCSimTrueTrack*>* fTrackList;

  TChain* fChain;
  TChain* fChainGeom;

  WCSimRootTrigger* fTrigger;
  WCSimRootEvent* fEvent;
  WCSimRootGeom* fGeometry;

  Double_t fEnergyThreshold;
  Double_t fRangeThreshold;

  ClassDef(WCSimInterface,0)

};

#endif
