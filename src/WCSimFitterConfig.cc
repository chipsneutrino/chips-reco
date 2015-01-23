/*
 * WCSimFitterConfig.cc
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */
#include "WCSimLikelihoodTrack.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimInterface.hh"

static WCSimFitterConfig * fgFitterConfig = 0;


WCSimFitterConfig::WCSimFitterConfig() : fNumTracks(1), fNumEventsToFit(0){
	// TODO Auto-generated constructor stub

}

WCSimFitterConfig::~WCSimFitterConfig() {
	// TODO Auto-generated destructor stub
}

WCSimFitterConfig* WCSimFitterConfig::Instance() {
	if( fgFitterConfig == 0 ){ fgFitterConfig = new WCSimFitterConfig(); }
	return fgFitterConfig;
}


void WCSimFitterConfig::SetNumTracks(int nTracks) {
  fFitterParameters.SetNumTracks(nTracks);
	fNumTracks = nTracks;
}

unsigned int WCSimFitterConfig::GetNumTracks() const{
  return fNumTracks;
}

unsigned int WCSimFitterConfig::GetNumParameters()
{
	return fFitterParameters.GetNumParameters();
}

unsigned int WCSimFitterConfig::GetNumIndependentParameters()
{
	return fFitterParameters.GetNumIndependentParameters();
}

void WCSimFitterConfig::SetParameter( unsigned int numTrack, const char * name, double min, double max, double start, bool fixed)
{
  std::cout << "In WCSimFitterConfig::SetParameter" << std::endl;
  FitterParameterType::Type type = FitterParameterType::FromName(name);
  WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
  trackPars->SetParMin(type, min);
  trackPars->SetParMax(type, max);
  trackPars->SetParStart(type, start);
  trackPars->SetParIsFixed(type, fixed);
  return;
}

void WCSimFitterConfig::FixTrackParameter(int numTrack, const char* name,
		bool doIt) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);

	trackPars->SetParIsFixed(type, true);
}

void WCSimFitterConfig::FreeTrackParameter(int numTrack, const char* name,
		bool doIt) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParIsFixed(type, false);
}

bool WCSimFitterConfig::GetIsFixedParameter(int numTrack, const char* name) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParIsFixed(type);
}

void WCSimFitterConfig::SetParMin(int numTrack, const char* name,
		double min) {
  std::cout << "In WCSimFitterConfig::SetParMin" << std::endl;
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParMin(type, min);
}

double WCSimFitterConfig::GetParMin(unsigned int numTrack, const char* name) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParMin(type);
}

void WCSimFitterConfig::SetParMax(int numTrack, const char* name,
		double max) {
  std::cout << "In WCSimFitterConfig::SetParMax" << std::endl;
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParMax(type, max);
}

double WCSimFitterConfig::GetParMax(unsigned int numTrack, const char * name) {

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParMax(type);
}

void WCSimFitterConfig::SetParStart(int numTrack, const char* name,
		double start) {
  std::cout << "In WCSimFitterConfig::SetParStart" << std::endl;

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParStart(type, start);
}

double WCSimFitterConfig::GetParStart(unsigned int numTrack, const char* name) {

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParStart(type);
}

void WCSimFitterConfig::SetParRange(int numTrack, const char* name,
		double min, double max) {
  std::cout << "In WCSimFitterConfig::SetParRange" << std::endl;
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParRange(type, min, max);
}

std::pair<double, double> WCSimFitterConfig::GetParRange(unsigned int numTrack,
		const char* name) {

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters * trackPars = fFitterParameters.GetTrackParameters(numTrack);
	std::pair<double, double> minMax;
	minMax = std::make_pair(trackPars->GetParMin(type),trackPars->GetParMax(type));
	return minMax;
}

void WCSimFitterConfig::SetNumEventsToFit(unsigned int numEvents) {
  int maxNumEvents = WCSimInterface::GetNumEvents();
  if(numEvents > maxNumEvents)
  {
    std::cerr << "Warning: Requested to fit " << numEvents << " events, but the file only contains " << maxNumEvents << std::endl;
    std::cerr << "         Fitting " << maxNumEvents << " events instead" << std::endl;
  }
	fNumEventsToFit = numEvents;
}

unsigned int WCSimFitterConfig::GetNumEventsToFit() {
	return fNumEventsToFit;
}


void WCSimFitterConfig::Print() {
	std::cout << " *** WCSimFitterConfig::Print() *** " << std::endl;
	std::cout << "Should add some statements in here..." << std::endl;
}

void WCSimFitterConfig::SetJoinParametersTogether(unsigned int numTrack1,
		unsigned int numTrack2, const char* name) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	fFitterParameters.JoinParametersTogether(numTrack1, numTrack2, type);
}

bool WCSimFitterConfig::GetJoinParametersTogether(unsigned int numTrack1,
		unsigned int numTrack2, const char* name) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	return fFitterParameters.GetJoinParametersTogether(numTrack1, numTrack2, type);
}

bool WCSimFitterConfig::GetIsParameterJoined(unsigned int numTrack, const char* name) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	return fFitterParameters.GetIsParameterJoined(numTrack, type);
}

unsigned int WCSimFitterConfig::GetTrackIsJoinedWith(unsigned int numTrack, const char* name) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	return fFitterParameters.GetTrackIsJoinedWith(numTrack, type);
}

void WCSimFitterConfig::SetTrackType(unsigned int numTrack,
		const char* typeName)
{
  std::cout << "WCSimFitterConfig::SetTrackType(" << numTrack << ", " << typeName << ")" << std::endl;
	WCSimLikelihoodTrack::TrackType type = WCSimLikelihoodTrack::GetTypeFromName(typeName);
	fFitterParameters.SetTrackType(numTrack, type);
}

WCSimLikelihoodTrack::TrackType WCSimFitterConfig::GetTrackType(const unsigned int &numTrack) const
{
	return fFitterParameters.GetTrackType(numTrack);
}

