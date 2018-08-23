/**
 * \class WCSimLikelihoodTrack
 * This class represents a single charged particle
 * track.  It holds its initial position and direction
 * of propagation, as well as its energy and the particle's
 * type.
 */

#ifndef WCSIMLIKELIHOODTRACK_HH
#define WCSIMLIKELIHOODTRACK_HH

#include "TObject.h"
#include "TVector3.h"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include <string>
#include <cstdlib>
class WCSimTrueTrack;

class WCSimLikelihoodTrack: public WCSimLikelihoodTrackBase {
	public:
		WCSimLikelihoodTrack();
		WCSimLikelihoodTrack(double x, double y, double z, double t, double theta, double phi, double E,
				TrackType::Type myType);
		virtual ~WCSimLikelihoodTrack();
		double GetTrackParameter(const FitterParameterType::Type &type) const;
		TVector3 GetPropagatedPos(const Double_t &s) const;
		void SetType(const TrackType::Type &type);
		TVector3 GetFirstEmissionVtx() const;ClassDef(WCSimLikelihoodTrack,1)
};

class WCSimLikelihoodPhotonTrack: public WCSimLikelihoodTrackBase {
	public:
		WCSimLikelihoodPhotonTrack();
		WCSimLikelihoodPhotonTrack(double x, double y, double z, double t, double theta, double phi, double E,
				double convDistance);
		virtual ~WCSimLikelihoodPhotonTrack();
		void SetConversionDistance(const double &convDist);
		double GetTrackParameter(const FitterParameterType::Type &type) const;
		TVector3 GetPropagatedPos(const Double_t &s) const;
		void SetType(const TrackType::Type &type);
		TVector3 GetFirstEmissionVtx() const;
	private:
		TVector3 fFirstEmissionVtx;ClassDef(WCSimLikelihoodPhotonTrack,1)
};

class WCSimLikelihoodUnknownTrack: public WCSimLikelihoodTrackBase {
	public:
		WCSimLikelihoodUnknownTrack();
		WCSimLikelihoodUnknownTrack(double x, double y, double z, double t, double theta, double phi, double E,
				double convDistance);
		virtual ~WCSimLikelihoodUnknownTrack();
		void SetConversionDistance(const double &convDist);
		double GetTrackParameter(const FitterParameterType::Type &type) const;
		TVector3 GetPropagatedPos(const Double_t &s) const;
		void SetType(const TrackType::Type &type);
		TVector3 GetFirstEmissionVtx() const;
	private:
		TVector3 fFirstEmissionVtx;ClassDef(WCSimLikelihoodUnknownTrack,1)
};

#endif /* WCSIMLIKELIHOODTRACK_HH */

