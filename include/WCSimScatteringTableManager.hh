#ifndef WCSIMSCATTERINGTABLEMANAGER_HH
#define WCSIMSCATTERINGTABLEMANAGER_HH

#include <vector>
#include "THnSparse.h"

class TFile;

class WCSimScatteringTableManager{

public:

  static WCSimScatteringTableManager* Instance();
  ~WCSimScatteringTableManager();

  float GetScatteringValue(std::vector<float> &values, int pmtLocation);
  float GetScatteringValue(int *bin, int pmtLocation);

private:
  WCSimScatteringTableManager();

  void LoadScatteringTables();

  int CalculateBin(int dimension, float value, int pmtLocation);

  // This class is a singleton, so we need an instance pointer.
  static WCSimScatteringTableManager* fManager;

  // Since THnSparses don't have a SetDirectory function, also need the file.
  TFile* fScatteringFile;

  // Need three THnSparseFs to store the tables.
  THnSparseF* fScatteringTableTop;
  THnSparseF* fScatteringTableBottom;
  THnSparseF* fScatteringTableBarrel;

};
#endif

