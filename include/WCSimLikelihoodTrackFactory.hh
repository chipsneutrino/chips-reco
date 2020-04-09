/*
 * WCSimLikelihoodTrackFactory.hh
 *
 *  Created on: 10 Jun 2015
 *      Author: andy
 */

#ifndef WCSIMLIKELIHOODTRACKFACTORY_HH_
#define WCSIMLIKELIHOODTRACKFACTORY_HH_

#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimTrackParameterEnums.hh"
#include <map>

class WCSimLikelihoodTrackFactory
{
private:
	WCSimLikelihoodTrackFactory();
	virtual ~WCSimLikelihoodTrackFactory();

public:
	static WCSimLikelihoodTrackBase *MakeTrack(const char *type);
	static WCSimLikelihoodTrackBase *MakeTrack(const TrackType &type);

	static WCSimLikelihoodTrackBase *MakeTrack(const TrackType &type, const double &x, const double &y,
											   const double &z, const double &t, const double &theta, const double &phi, const double &energy,
											   std::map<FitterParameterType::Type, double> &extraPars);
	static WCSimLikelihoodTrackBase *MakeTrack(const TrackType &type, const double &x, const double &y,
											   const double &z, const double &t, const double &theta, const double &phi, const double &energy);
	static WCSimLikelihoodTrackBase *MakeTrack(const TrackType &type, const double &x, const double &y,
											   const double &z, const double &t, const double &theta, const double &phi, const double &energy,
											   const double &conv);
};

#endif /* WCSIMLIKELIHOODTRACKFACTORY_HH_ */
