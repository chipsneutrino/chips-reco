#ifndef WCSIM_RECOSLICER_HH
#define WCSIM_RECOSLICER_HH

#include <vector>

class WCSimRecoEvent;
class WCSimRecoDigit;
class WCSimCosmicSeed;
class WCSimRecoClusteringUtil;

// A class to slice up the main event for the purpose of running the
// basic reconstruction for seeding. The output events are basically
// transient and will not exist beyond the end of the seeding.
class WCSimRecoSlicer {

  public:

  WCSimRecoSlicer();
  ~WCSimRecoSlicer();

  // Function to perform the slicing
  void Run();

  void SetInputEvent(WCSimRecoEvent* evt);
  std::vector<WCSimRecoEvent*> GetSlicedEvents(){return fSlicedEvents;};

  void SetDistanceCut(double val){fDistanceCut = val;};
  double GetDistanceCut() const {return fDistanceCut;};

  void SetTimeCut(double val){fTimeCut = val;};
  double GetTimeCut() const {return fTimeCut;};

  void SetMinSliceSize(unsigned int val){fMinSliceSize = val;};
  unsigned int GetMinSliceSize(unsigned int val) const {return fMinSliceSize;};

  void SetChargeCut(double val){fChargeCut = val;};
  double GetChargeCut() const {return fChargeCut;};

  void SetIterateSlicing(bool val){fIterateSlicing = val;};
  bool GetIterateSlicing() const{return fIterateSlicing;};

  void SetNumberOfTracks(unsigned int val){fNTracks = val;};
  unsigned int GetNumberOfTracks() const {return fNTracks;};

  void SetCosmicSeed(WCSimCosmicSeed* seed){fCosmicSeed = seed;};

  private:

  // Actually do the slicing
  void SliceTheEvent();
  
  void UpdateParameters();
  void ResetSlices();

  // Build events from the filtered digits.
  void BuildEvents();

  WCSimRecoEvent* fInputEvent;
  std::vector<WCSimRecoEvent*> fSlicedEvents;

  std::vector< std::vector<WCSimRecoDigit*> > fSlicedDigits;

  double fDistanceCut;
  double fTimeCut;
  unsigned int fMinSliceSize;
  double fChargeCut;

  bool fIterateSlicing; ///< Keep trying to slice, increasing the charge cut, in order to get the number of slices
                        ///< to match the number of tracks.

  unsigned int fNTracks;

  WCSimCosmicSeed* fCosmicSeed;

  WCSimRecoClusteringUtil* fClusteringUtil;

};

#endif
