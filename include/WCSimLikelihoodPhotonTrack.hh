/*
 * WCSimLikelihoodPhotonTrack.hh
 *
 *  Created on: 9 Jun 2015
 *      Author: andy
 */

#ifndef WCSIMLIKELIHOODPHOTONTRACK_HH_
#define WCSIMLIKELIHOODPHOTONTRACK_HH_

#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrackParameterEnums.hh"

class WCSimLikelihoodPhotonTrack: public WCSimLikelihoodTrackBase {
public:
	WCSimLikelihoodPhotonTrack();
    WCSimLikelihoodPhotonTrack( double x, double y, double z, double t,
                          	    double theta, double phi, double E,
                          	    double convDistance);
	virtual ~WCSimLikelihoodPhotonTrack();
	void SetConversionDistance(const double &convDist);
	double GetTrackParameter(const FitterParameterType::Type &type) const;
	TVector3 GetPropagatedPos(const Double_t &s) const;
	void SetType(const TrackType::Type &type);
  TVector3 GetFirstEmissionVtx() const;


private:
  TVector3 fFirstEmissionVtx;

};

#endif /* WCSIMLIKELIHOODPHOTONTRACK_HH_ */
