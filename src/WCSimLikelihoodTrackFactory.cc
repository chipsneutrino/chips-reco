/*
 * WCSimLikelihoodTrackFactory.cc
 *
 *  Created on: 10 Jun 2015
 *      Author: andy
 */

#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodPhotonTrack.hh"
#include "WCSimTrackParameterEnums.hh"
#include <map>

WCSimLikelihoodTrackFactory::WCSimLikelihoodTrackFactory() {
	// TODO Auto-generated constructor stub

}

WCSimLikelihoodTrackFactory::~WCSimLikelihoodTrackFactory() {
	// TODO Auto-generated destructor stub
}

WCSimLikelihoodTrackBase* WCSimLikelihoodTrackFactory::MakeTrack(const char* type) {
	return MakeTrack(TrackType::FromName(type));
}

WCSimLikelihoodTrackBase* WCSimLikelihoodTrackFactory::MakeTrack(
		const TrackType &type) {
  WCSimLikelihoodTrackBase * myTrack = 0x0;
	if(	   type == TrackType::ElectronLike
		|| type == TrackType::MuonLike 	)
	{
		myTrack = new WCSimLikelihoodTrack();
    myTrack->SetType(type);
	}
	else if( type == TrackType::PhotonLike)
	{
    myTrack = new WCSimLikelihoodPhotonTrack();
    myTrack->SetType(type);
	}
  else
  {
	  std::cerr << "Error in WCSimLikelihoodTrackFactory::MakeTrack" << std::endl
	  		<< "Don't know how to make a track of type " << type
		  	<< " (" << TrackType::AsString(type) << ")" << std::endl;
	  assert(false);
	}
  return myTrack;
}


WCSimLikelihoodTrackBase* WCSimLikelihoodTrackFactory::MakeTrack(
		const TrackType& type, const double& x,
		const double& y, const double& z, const double& t, const double& theta,
		const double& phi, const double& energy,
		std::map<FitterParameterType, double> &extraPars) {

	if(	   type == TrackType::ElectronLike
		|| type == TrackType::MuonLike 	)
	{
		return new WCSimLikelihoodTrack(x, y, z, t, theta, phi, energy, type);
	}
	else if( type == TrackType::PhotonLike)
	{
		std::map<FitterParameterType, double>::iterator mapItr = extraPars.find(FitterParameterType::kConversionDistance);
		double conversionDistance = 0.0;
		if( mapItr != extraPars.end() )
		{
			conversionDistance = mapItr->second;
		}
		return new WCSimLikelihoodPhotonTrack(x, y, z, t, theta, phi, energy, conversionDistance);
	}
	std::cerr << "Error in WCSimLikelihoodTrackFactory::MakeTrack" << std::endl
			<< "Don't know how to make a track of type " << type
			<< " (" << TrackType::AsString(type) << ")" << std::endl;
	assert(false);
	return 0x0;

}

WCSimLikelihoodTrackBase* WCSimLikelihoodTrackFactory::MakeTrack(
		const TrackType& type, const double& x,
		const double& y, const double& z, const double& t, const double& theta,
		const double& phi, const double& energy) {
	std::map<FitterParameterType, double> emptyMap;
	return MakeTrack(type, x, y, z, t, theta, phi, energy, emptyMap);
}