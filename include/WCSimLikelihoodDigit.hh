#ifndef WCSIMLIKELIHOODDIGIT_H
#define WCSIMLIKELIHOODDIGIT_H

#include <vector>
#include "WCSimRootEvent.hh"
#include "TClonesArray.h"
#include "TObject.h"

class WCSimLikelihoodDigit : public TObject
{
  public:
      WCSimLikelihoodDigit( Double_t x, Double_t y, Double_t z, Double_t t, Double_t Q, Int_t tubeId, Double_t faceX, Double_t faceY, Double_t faceZ );
      WCSimLikelihoodDigit( WCSimRootCherenkovDigiHit * myDigiHit );
      virtual ~WCSimLikelihoodDigit();
	
      int GetTubeId();
      double GetQ();
      double GetT();
      double GetX();
      double GetY();
      double GetZ();
      double GetFaceX();
      double GetFaceY();
      double GetFaceZ();
      void Print();
 
  protected:
  private:
	Int_t fTubeId;
	Double_t fQ;
	Double_t fT;
	Double_t fPos[3];
	Double_t fFace[3];

	ClassDef(WCSimLikelihoodDigit,1)
		
};

#endif
