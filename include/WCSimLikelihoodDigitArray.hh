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
        WCSimLikelihoodDigitArray(WCSimRootEvent *);
        virtual ~WCSimLikelihoodDigitArray();

        WCSimLikelihoodDigit * GetDigit( Int_t digit);
        Int_t GetNDigits();
 
    protected:
    private:
				TClonesArray * fLikelihoodDigitArray;
				Int_t fNLikelihoodDigits;

		ClassDef(WCSimLikelihoodDigitArray,1)
};

#endif // WCSIMCHARGELIKELIHOOD_H
