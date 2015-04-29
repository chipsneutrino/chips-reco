/*
 * WCSimIntegralLookupReader.hh
 *
 *  Created on: 12 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUPREADER_HH_
#define WCSIMINTEGRALLOOKUPREADER_HH_
#include "WCSimLikelihoodTrack.hh"
#include "WCSimIntegralLookup3D.hh"
#include <TObject.h>
#include <map>


class WCSimIntegralLookupReader : public TObject{
public:
	void LoadIntegrals(WCSimLikelihoodTrack * myTrack);
	void LoadIntegrals(const WCSimLikelihoodTrack::TrackType &type);
	static WCSimIntegralLookupReader * Instance();
	WCSimIntegralLookup3D * GetIntegralLookup3D(WCSimLikelihoodTrack::TrackType type);
  
  double GetTrackLengthForPercentile(const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &percentile);

	double GetRhoIntegral( const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &s);
	double GetRhoSIntegral( const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &s);
	double GetRhoSSIntegral( const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &s);

	double GetRhoGIntegral( const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSIntegral( const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSSIntegral( const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &s, const double &R0, const double &cosTh0);
	virtual ~WCSimIntegralLookupReader();
private:
	WCSimIntegralLookupReader();
	std::map<WCSimLikelihoodTrack::TrackType, WCSimIntegralLookup3D*> fLookupMap;
	TString GetLookupFilename(WCSimLikelihoodTrack * track);
	TString GetLookupFilename(const WCSimLikelihoodTrack::TrackType &type);
      
  WCSimLikelihoodTrack::TrackType fLastPercentileType;
  double fLastPercentileEnergy;
  double fLastPercentileLength;
      

	ClassDef(WCSimIntegralLookupReader,1)
};

#endif /* WCSIMINTEGRALLOOKUPREADER_HH_ */
