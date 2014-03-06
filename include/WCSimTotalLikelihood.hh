#ifndef WCSIMTOTALLIKELIHOOD_H
#define WCSIMTOTALLIKELIHOOD_H

#include "WCSimChargeLikelihood.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"

#include "TObject.h"
#include <vector>

class WCSimTotalLikelihood : public TObject
{
  public:
      WCSimTotalLikelihood( WCSimLikelihoodDigitArray * myLikelihoodDigitArray );
      virtual ~WCSimTotalLikelihood();
      void SetTracks(std::vector<WCSimLikelihoodTrack>);
      void ResetTracks();
      Double_t Calc2LnL();
 
  protected:
  private:
      WCSimLikelihoodDigitArray * fLikelihoodDigitArray;
      WCSimChargeLikelihood fChargeLikelihood;
      //WCSimTimeLikelihood fTimeLikelihood;
      std::vector<WCSimLikelihoodTrack> fTracks;

	ClassDef(WCSimTotalLikelihood,1)
		
};

#endif
