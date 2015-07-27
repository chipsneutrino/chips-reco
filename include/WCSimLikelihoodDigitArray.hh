/**
 * \class WCSimLikelihoodDigitArray
 * A class to hold one WCSimLikelihoodDigit for every
 * PMT in the detector, hit or unhit.
 * It contains all the detector information used to fit
 * a single event
 */
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
    
        /**
         * Geometry type enum
         */
        enum GeomType_t
        {
          kUnknown  = 0,//!< Unknown type
          kCylinder = 1,//!< Detector is a cylinder
          kMailBox  = 2 //!< Detector is a cuboid
        };
    
        //TODO: we have a pointer as class member, so there should be
        //      a proper copy constructor and proper destructor...

        WCSimLikelihoodDigitArray();
        /**
         * Constructor
         * @param myRootEvent Event object from WCSim to build hit list from
         */
        WCSimLikelihoodDigitArray(WCSimRootEvent * myRootEvent);

        /**
         * Constructor using undigitized hits for debugging
         * @param myRootEvent Event object from WCSim to build the hit list from
         * @param useUndigitized A hack to call the undigitized constructor.  Has to be true.
         */
        WCSimLikelihoodDigitArray( WCSimRootEvent * myRootEvent, Bool_t useUndigitized);

        virtual ~WCSimLikelihoodDigitArray();

        /**
         * Get a single WCSimLikelihoodDigit object
         * @param digit The digit's position in the TClonesArray (NOT its WCSim tubeID)
         * @return The WCSimLikelihoodDigit at this place in the array
         */
        WCSimLikelihoodDigit * GetDigit( Int_t digit);

        /**
         * @return Total number of digits in the array (i.e. PMTs in the detector)
         */
        Int_t GetNDigits();


        /** 
         * @return Total number of digits in the array with a non-zero hit (that
         * passed the filtering step)
         * */
        Int_t GetNHits();

        /**
         * @return True if detector is a cylinder
         */
        Bool_t IsCylinder();

        /**
         * @return True if detector is a mailbox (cuboid)
         */
        Bool_t IsMailBox();

        /**
         * @return The geometry type enum for this detector
         */
        GeomType_t GetGeomType();

        /**
         * @param i (0,1,2) for (x,y,z) if a mailbox and (r,r,z) if a cylinder
         * @return The co-ordinate at which a particle will exit the detector volume
         */
        Double_t GetExtent(Int_t i);
 
    protected:
    private:
        TClonesArray * fLikelihoodDigitArray; ///< Array with every PMT and its hit information
        Int_t fNLikelihoodDigits; ///< Total number of PMTs
        Int_t fNumHitPMTs; ///< Total number of hit PMTs, ie. ones with Q>0 that weren't filtered out
        GeomType_t fGeomType;  ///< Detector geometry type
        Double_t fExtent[3];  ///< Maximum (x, y, z) or (r,r,z)  in the detector

		ClassDef(WCSimLikelihoodDigitArray,1)
};

#endif // WCSIMCHARGELIKELIHOOD_H
