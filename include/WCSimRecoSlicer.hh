#ifndef WCSIM_RECOSLICER_HH
#define WCSIM_RECOSLICER_HH

#include <vector>

class WCSimRecoEvent;
class WCSimRecoDigit;

// A class to slice up the main event for the purpose of running the
// basic reconstruction for seeding. The output events are basically
// transient and will not exist beyond the end of the seeding.
class WCSimRecoSlicer {

  public:

  WCSimRecoSlicer();
  ~WCSimRecoSlicer();

  // Function to perform the slicing
  void SliceTheEvent();

  void SetInputEvent(WCSimRecoEvent* evt);
  std::vector<WCSimRecoEvent*> GetSlicedEvents(){return fSlicedEvents;};

  private:
  
  void UpdateParameters();

  bool AddDigitIfClose(WCSimRecoDigit*,std::vector<WCSimRecoDigit*>&);
  double GetDistanceBetweenDigits(WCSimRecoDigit*, WCSimRecoDigit*);

  unsigned int GetFirstUnclusteredDigit();
  bool AreAllDigitsClustered();

  // Build events from the filtered digits.
  void BuildEvents();

  WCSimRecoEvent* fInputEvent;
  std::vector<WCSimRecoEvent*> fSlicedEvents;

  std::vector<bool> fIsDigitClustered;

  std::vector< std::vector<WCSimRecoDigit*> > fSlicedDigits;

  double fDistanceCut;
  unsigned int fMinSliceSize;
  double fChargeCut;
};

#endif
