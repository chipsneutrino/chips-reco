#ifndef WCSIMLIKELIHOODDIGIT_H
#define WCSIMLIKELIHOODDIGIT_H

#include <vector>
#include "WCSimRootEvent.hh"
#include "TClonesArray.h"
#include "TObject.h"
#include "TVector3.h"

class WCSimLikelihoodDigit : public TObject
{
  public:
      WCSimLikelihoodDigit( Double_t x, Double_t y, Double_t z, Double_t t, Double_t Q, Int_t tubeId, Double_t faceX, Double_t faceY, Double_t faceZ );
      WCSimLikelihoodDigit( WCSimRootCherenkovDigiHit * myDigiHit );
      virtual ~WCSimLikelihoodDigit();
	
      int GetTubeId() const;
      double GetQ() const;
      double GetT() const;

      TVector3 GetPos() const;
      double GetX() const;
      double GetY() const;
      double GetZ() const;

      TVector3 GetFace() const;
      double GetFaceX() const;
      double GetFaceY() const;
      double GetFaceZ() const;

      void Print() const;
 
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
