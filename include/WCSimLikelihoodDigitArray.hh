#ifndef WCSIMLIKELIHOODDIGITARRAY_H
#define WCSIMLIKELIHOODDIGITARRAY_H

#include <vector>
#include "WCSimLikelihoodDigit.hh"
#include "WCSimRootEvent.hh"
#include "TClonesArray.h"
#include "TObject.h"

class WCSimLikelihoodDigitArray : public TObject
{
    public:
        WCSimLikelihoodDigitArray();
        WCSimLikelihoodDigitArray(WCSimRootEvent * myRootEvent);
        WCSimLikelihoodDigitArray( WCSimRootEvent * myRootEvent, Bool_t useUndigitized);
        virtual ~WCSimLikelihoodDigitArray();

        WCSimLikelihoodDigit * GetDigit( Int_t digit);
        Int_t GetNDigits();
        Bool_t IsCylinder();
        Bool_t IsMailBox();
        Double_t GetExtent(Int_t i);
 
    protected:
    private:
				TClonesArray * fLikelihoodDigitArray;
				Int_t fNLikelihoodDigits;
        Int_t fGeomType;  // Cylinder = 0, Mailbox = 1, Unkown = -1
        Double_t fExtent[3];  // Maximum x, y, and z in the detector
                              // uses volume, not fiducial volume

		ClassDef(WCSimLikelihoodDigitArray,1)
};

#endif // WCSIMCHARGELIKELIHOOD_H
