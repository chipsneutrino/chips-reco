/**
 * \class
 * This class holds hit information for a single PMT
 * such as its position, orientation, recorded charge
 * and hit time.  WCSimAnalysis interfaces with hits
 * from WCSim using this class.
 */
#ifndef WCSIMLIKELIHOODDIGIT_H
#define WCSIMLIKELIHOODDIGIT_H

#include <vector>
#include "WCSimDetectorParameters.hh"
#include "WCSimRootEvent.hh"
#include "TClonesArray.h"
#include "TObject.h"
#include "TVector3.h"
#include "TString.h"

class WCSimLikelihoodDigit: public TObject
{
public:

    /**
     * Constructor
     * @param x PMT x co-ordinate
     * @param y PMT y co-ordinate
     * @param z PMT x co-ordinate
     * @param t Time of hit
     * @param Q Recorded charge (P.E.)
     * @param tubeId Unique PMT ID number from WCSim
     * @param faceX x co-ordinate of PMT normal
     * @param faceY y co-ordinate of PMT normal
     * @param faceZ z co-ordinate of PMT normal
     */
    WCSimLikelihoodDigit(Double_t x, Double_t y, Double_t z, Double_t t,
            Double_t Q, Int_t tubeId, Double_t faceX, Double_t faceY,
            Double_t faceZ, TString name, TGraph * wlWeightedQE, Double_t wlWeightedRefIndex, double exposeHeight);

    /**
     * Constructor
     * @param myDigiHit Digitized Cherenkov hit object from WCSim
     */
    WCSimLikelihoodDigit(WCSimRootCherenkovDigiHit * myDigiHit);

    /**
     * Copy Constructor
     */
    WCSimLikelihoodDigit(const WCSimLikelihoodDigit &otherLikelihoodDigit);

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

    double GetAverageQE(const double &distanceToPMT) const;
    double GetAverageRefIndex() const;
    double GetExposeHeight() const;

    TString GetPMTName() const;

    void Print() const;

protected:
private:
    Int_t fTubeId;     ///< Unique PMT ID number from WCSim
    Double_t fQ;       ///< Digitized charge (P.E.)
    Double_t fT;       ///< Time of hit
    Double_t fPos[3];  ///< (x,y,z) co-ordinates of the PMT location
    Double_t fFace[3]; ///< (x,y,z) components of the direction normal to the PMT
    TString fPMTName; ///< Name of PMT type, e.g. 3_inch_HQE

    TGraph * fAverageQE; ///< Average QE of the PMT, from weighting QE(wavelength) by the average Cherenkov spectrum
    double fAverageRefIndex; ///< Weight WCSim's refractive index by (wavelength * PMT QE(wavelength))
    double fExposeHeight; ///< Heigh of PMT dome expose through the detector liner


ClassDef(WCSimLikelihoodDigit,1)

};

#endif
