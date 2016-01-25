/*
 * WCSimLikelihoodTrackBase.hh
 *
 *  Created on: 10 Jun 2015
 *      Author: andy
 */

#ifndef WCSIMLIKELIHOODTRACKBASE_HH_
#define WCSIMLIKELIHOODTRACKBASE_HH_
#include "WCSimTrackParameterEnums.hh"
class WCSimTrueTrack;
#include "TVector3.h"
#include "TMath.h"

class WCSimLikelihoodTrackBase {
public:
      
      /*
       * Constructor from true WCSim track object
      * @param trueTrack Track to convert into a WCSimLikelihoodTrack
      */
     WCSimLikelihoodTrackBase( WCSimTrueTrack * trueTrack);

     virtual ~WCSimLikelihoodTrackBase();
     virtual void Print();

     // Setters
     void SetX(double x);
     void SetY(double y);
     void SetZ(double z);
     void SetT(double t);
     void SetTheta(double th);
     void SetPhi(double phi);
     void SetE(double E);
     virtual void SetType(const TrackType::Type &type) = 0;

     // Getters
     double GetX() const;
     double GetY() const;
     double GetZ() const;
     TVector3 GetVtx() const;
     virtual TVector3 GetFirstEmissionVtx() const = 0;
     double GetT() const;
     double GetTheta() const;
     double GetPhi() const;
     double GetDirX() const;
     double GetDirY() const;
     double GetDirZ() const;
     double GetConversionDistance() const;
     TVector3 GetDir() const;
     double GetE() const;
     TrackType GetType() const;
     virtual double GetTrackParameter(const FitterParameterType::Type &type) const = 0;

     static bool EnergyGreaterThanOrEqual(const WCSimLikelihoodTrackBase &a, const WCSimLikelihoodTrackBase &b);
     static bool EnergyGreaterThanOrEqualPtrs(WCSimLikelihoodTrackBase *a, WCSimLikelihoodTrackBase *b);
     /**
      * Function to calculate the position a distance s away from the track vertex in
      * the direction of travel of the track
      * @param s Distance to propagate from the vertex
      * @return The position of the particle after it propagates a distance s from
      *         the vertex
      */
     virtual TVector3 GetPropagatedPos(const Double_t &s) const = 0;
     Int_t GetPDG() const;
     bool operator == (const WCSimLikelihoodTrackBase &b) const;
     bool IsSameTrack(WCSimLikelihoodTrackBase * b) const;
     double GetPropagationSpeedFrac() const;

 protected:
   void LoadEmissionProfile(WCSimLikelihoodTrackBase * myTrack);
   double fVtx[3];	    ///< Vertex position
   double fDir[3];		///< (x,y,z) components of direction - calculate it once every time theta and phi are set
   double fT0;			///< Vertex time
   double fTheta0;	    ///< Polar angle to z axis of direction
   double fPhi0;		///< Azimuthal angle of direction
   double fE0;			///< Kinetic energy
   double fConversionDistance; ///< Distance from the vertex to the point where Cherenkov light is first emitted (zero unless a photon)

   TrackType::Type fType; ///< Particle type (e, mu, pi...)
   WCSimLikelihoodTrackBase();
};

#endif /* WCSIMLIKELIHOODTRACKBASE_HH_ */
