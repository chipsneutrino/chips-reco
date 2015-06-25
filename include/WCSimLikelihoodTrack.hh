/**
 * \class WCSimLikelihoodTrack
 * This class represents a single charged particle
 * track.  It holds its initial position and direction
 * of propagation, as well as its energy and the particle's
 * type.
 */

#ifndef WCSIMLIKELIHOODTRACK_H
#define WCSIMLIKELIHOODTRACK_H

#include "TObject.h"
#include "TVector3.h"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include <string>
#include <cstdlib>
class WCSimTrueTrack;

class WCSimLikelihoodTrack : public WCSimLikelihoodTrackBase
{
    
    public:
		WCSimLikelihoodTrack();
		WCSimLikelihoodTrack(double x, double y, double z, double t, double theta, double phi, double E, TrackType::Type myType);
		void Print();
		double GetTrackParameter(const FitterParameterType::Type &type) const;
		TVector3 GetPropagatedPos(const Double_t &s) const;
		void SetType(const TrackType::Type &type);
		ClassDef(WCSimLikelihoodTrack,1)
};



#endif

