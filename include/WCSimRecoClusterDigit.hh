#ifndef WCSIMRECOCLUSTERDIGIT_HH
#define WCSIMRECOCLUSTERDIGIT_HH

#include "TObject.h"

#include "WCSimRecoDigit.hh"

#include <vector>

class WCSimRecoClusterDigit : public TObject {

 public:
  WCSimRecoClusterDigit( WCSimRecoDigit* recoDigit );
  ~WCSimRecoClusterDigit();

  void AddClusterDigit(WCSimRecoClusterDigit* clusterDigit);
  
  Int_t GetNClusterDigits(); 
  WCSimRecoClusterDigit* GetClusterDigit(Int_t idigit); 
  std::vector<WCSimRecoClusterDigit*>* GetClusterDigitList();

  Double_t GetX(){ return fRecoDigit->GetX(); }
  Double_t GetY(){ return fRecoDigit->GetY(); }
  Double_t GetZ(){ return fRecoDigit->GetZ(); }
  Double_t GetTime(){ return fRecoDigit->GetTime(); }

  void SetClustered( Bool_t yesno = 1 ){ fIsClustered = yesno; }
  Bool_t IsClustered(){ return fIsClustered; }

  Bool_t IsAllClustered();

  WCSimRecoDigit* GetRecoDigit(){ return fRecoDigit; }

 private:

  Bool_t fIsClustered;
  Bool_t fIsAllClustered;

  WCSimRecoDigit* fRecoDigit;

  std::vector<WCSimRecoClusterDigit*>* fClusterDigitList;

  ClassDef(WCSimRecoClusterDigit,0)
};

#endif
