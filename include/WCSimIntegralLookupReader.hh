/*
 * WCSimIntegralLookupReader.hh
 *
 *  Created on: 12 Mar 2015
 *      Author: ajperch
 */

#ifndef WCSIMINTEGRALLOOKUPREADER_HH_
#define WCSIMINTEGRALLOOKUPREADER_HH_
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimIntegralLookup3D.hh"
#include <TObject.h>
#include <map>


class WCSimIntegralLookupReader : public TObject{
public:
	void LoadIntegrals(WCSimLikelihoodTrackBase * myTrack);
	void LoadIntegrals(const TrackType::Type &type);
	static WCSimIntegralLookupReader * Instance();
	WCSimIntegralLookup3D * GetIntegralLookup3D(TrackType::Type type);
  
	double GetRhoIntegral( const TrackType::Type &type, const double &E, const double &s);
	double GetRhoSIntegral( const TrackType::Type &type, const double &E, const double &s);
	double GetRhoSSIntegral( const TrackType::Type &type, const double &E, const double &s);

	double GetRhoGIntegral( const TrackType::Type &type, const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSIntegral( const TrackType::Type &type, const double &E, const double &s, const double &R0, const double &cosTh0);
	double GetRhoGSSIntegral( const TrackType::Type &type, const double &E, const double &s, const double &R0, const double &cosTh0);
	virtual ~WCSimIntegralLookupReader();


    void SaveIntegrals( const TrackType::Type &type, const double E, const double s, const double R0, const double cosTh0);
private:
	WCSimIntegralLookupReader();
	std::map<TrackType::Type, WCSimIntegralLookup3D*> fLookupMap;
	TString GetLookupFilename(WCSimLikelihoodTrackBase * track);
	TString GetLookupFilename(const TrackType::Type &type);
      
      

	ClassDef(WCSimIntegralLookupReader,1)
};

#endif /* WCSIMINTEGRALLOOKUPREADER_HH_ */
