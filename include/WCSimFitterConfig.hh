/*
 * WCSimFitterConfig.hh
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#ifndef WCSIMFITTERCONFIG_HH_
#define WCSIMFITTERCONFIG_HH_
#include "WCSimFitterParameters.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrackParameterEnums.hh"

class WCSimFitterConfig{
public:
	WCSimFitterConfig();
	void Print();

	virtual ~WCSimFitterConfig();
	void SetNumTracks(int nTracks);
  unsigned int GetNumTracks() const;
	unsigned int GetNumParameters();
  unsigned int GetNumIndependentParameters();

	void SetTrackType(unsigned int numTrack, const char * typeName);
	TrackType::Type GetTrackType(const unsigned int &numTrack) const;
	WCSimFitterParameters GetFitterParameters(){ return fFitterParameters; }

	void FixTrackParameter(int  numTrack, const char * name, bool doIt = true);
	void FreeTrackParameter(int numTrack, const char * name, bool doIt = true);
	bool GetIsFixedParameter(int numTrack, const char * name);

  void SetParameter(unsigned int numTrack, const char * name, double min, double max, double start, bool fixIt);

	void SetParMin( int numTrack, const char * name, double min);
	double GetParMin( unsigned int numTrack, const char * name );

	void SetParMax( int numTrack, const char * name, double max);
	double GetParMax(unsigned int numTrack, const char * name);

	void SetParStart( int numTrack, const char * name, double start);
	double GetParStart( unsigned int numTrack, const char * name );

	void SetParRange( int numTrack, const char * name, double min, double max);
	std::pair<double, double> GetParRange( unsigned int numTrack, const char * name );

	void SetNumEventsToFit(int numEvents);
	int GetNumEventsToFit();

	void SetFirstEventToFit(unsigned int iEvt);
	int GetFirstEventToFit() const;

	void SetJoinParametersTogether(unsigned int numTrack1, unsigned int numTrack2, const char * name);
	bool GetJoinParametersTogether(unsigned int numTrack1, unsigned int numTrack2, const char * name);
	bool GetIsParameterJoined(unsigned int numTrack, const char * name);
	unsigned int GetTrackIsJoinedWith(unsigned int numTrack, const char * name);
	bool GetIsParameterJoined(unsigned int numTrack, FitterParameterType::Type type);
	unsigned int GetTrackIsJoinedWith(unsigned int numTrack, FitterParameterType::Type type);

	bool GetMakeFits() const;
	void SetMakeFits(const bool &makeFits = true);

	void SetForcePiZeroMass(const bool &doIt = true);
	bool GetForcePiZeroMass();

	bool GetIsPiZeroFit() const;
	void SetIsPiZeroFit(bool isPiZero);
private:



    bool fMakeFits;
	WCSimFitterParameters fFitterParameters;
    int fNumTracks;
	int fNumEventsToFit;
	int fNumParameters;
	int fFirstEventToFit;
	bool fIsPiZeroFit;
	bool fForcePiZeroMass;
};

#endif /* WCSIMFITTERCONFIG_HH_ */
