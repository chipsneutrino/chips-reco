#ifndef WCSIMRECOSEED_HH
#define WCSIMRECOSEED_HH

#include "WCSimReco.hh"

#include <vector>

class WCSimRecoDigit;
class WCSimRecoVertex;
class WCSimRecoRing;

class WCSimRecoSeed : public WCSimReco {

 public:
  WCSimRecoSeed();
  ~WCSimRecoSeed();

  // Reconstruction Methods
  // ======================
  void Run();
  void Run(WCSimRecoEvent* evt);
  std::vector<WCSimRecoEvent*> RunSeed(WCSimRecoEvent* evt);
  void RunFilter(WCSimRecoEvent* evt);
  std::vector<WCSimRecoEvent*> RunSlicer(WCSimRecoEvent* evt);
  void RunRecoVertex(WCSimRecoEvent* evt);
  void RunRecoRings(WCSimRecoEvent* evt);


 private:

  ClassDef(WCSimRecoSeed,0)
};

#endif
