/*
 * WCSimFitterConfig.cc
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */
#include "WCSimLikelihoodTrack.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimInterface.hh"
#include "WCSimTrackParameterEnums.hh"

WCSimFitterConfig::WCSimFitterConfig() : fNumTracks(1), fNumEventsToFit(0), fFirstEventToFit(0)
{
	fIsPiZeroFit = 0;
	fIsCosmicFit = 0;
	fNumParameters = 0;
	fForcePiZeroMass = 0;
	fMakeFits = 1;

	// TODO Auto-generated constructor stub
}

WCSimFitterConfig::~WCSimFitterConfig()
{
	// TODO Auto-generated destructor stub
}

void WCSimFitterConfig::SetNumTracks(int nTracks)
{
	fFitterParameters.SetNumTracks(nTracks);
	fNumTracks = nTracks;
}

unsigned int WCSimFitterConfig::GetNumTracks() const
{
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

void WCSimFitterConfig::SetParameter(unsigned int numTrack, const char *name,
									 double min, double max, double start, double step, bool isFixed)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	std::cout << "In WCSimFitterConfig::SetParameter: " << type << std::endl;
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParMin(type, min);
	trackPars->SetParMax(type, max);
	trackPars->SetParStart(type, start);
	trackPars->SetParStep(type, step);
	trackPars->SetParIsFixed(type, isFixed);
	return;
}

void WCSimFitterConfig::FixTrackParameter(int numTrack, const char *name,
										  bool doIt)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);

	trackPars->SetParIsFixed(type, true);
}

void WCSimFitterConfig::FreeTrackParameter(int numTrack, const char *name,
										   bool doIt)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParIsFixed(type, false);
}

bool WCSimFitterConfig::GetIsFixedParameter(int numTrack, const char *name)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParIsFixed(type);
}

void WCSimFitterConfig::SetParMin(int numTrack, const char *name, double min)
{
	//std::cout << "In WCSimFitterConfig::SetParMin" << std::endl;
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParMin(type, min);
}

double WCSimFitterConfig::GetParMin(unsigned int numTrack, const char *name)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParMin(type);
}

void WCSimFitterConfig::SetParMax(int numTrack, const char *name, double max)
{
	//std::cout << "In WCSimFitterConfig::SetParMax" << std::endl;
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParMax(type, max);
}

double WCSimFitterConfig::GetParMax(unsigned int numTrack, const char *name)
{

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParMax(type);
}

void WCSimFitterConfig::SetParStart(int numTrack, const char *name,
									double start)
{
	//std::cout << "In WCSimFitterConfig::SetParStart" << std::endl;

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParStart(type, start);
}

double WCSimFitterConfig::GetParStart(unsigned int numTrack, const char *name)
{

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParStart(type);
}

void WCSimFitterConfig::SetParStep(int numTrack, const char *name,
								   double step)
{
	//std::cout << "In WCSimFitterConfig::SetParStep" << std::endl;

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParStep(type, step);
}

double WCSimFitterConfig::GetParStep(unsigned int numTrack, const char *name)
{

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	return trackPars->GetParStep(type);
}

void WCSimFitterConfig::SetParRange(int numTrack, const char *name, double min,
									double max)
{
	//std::cout << "In WCSimFitterConfig::SetParRange" << std::endl;
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	trackPars->SetParRange(type, min, max);
}

std::pair<double, double> WCSimFitterConfig::GetParRange(unsigned int numTrack,
														 const char *name)
{

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	WCSimFitterSingleTrackParameters *trackPars =
		fFitterParameters.GetTrackParameters(numTrack);
	std::pair<double, double> minMax;
	minMax = std::make_pair(trackPars->GetParMin(type),
							trackPars->GetParMax(type));
	return minMax;
}

void WCSimFitterConfig::SetNumEventsToFit(int numEvents)
{
	int maxNumEvents = WCSimInterface::GetNumEvents();
	if (numEvents > maxNumEvents)
	{
		std::cerr << "Warning: Requested to fit " << numEvents
				  << " events, but the file only contains " << maxNumEvents
				  << std::endl;
		std::cerr << "         Fitting " << maxNumEvents << " events instead"
				  << std::endl;
		fNumEventsToFit = maxNumEvents;
	}
	else if ((numEvents + fNumEventsToFit) > maxNumEvents)
	{
		std::cerr << "Warning: Requested to fit " << numEvents
				  << " starting at event " << fFirstEventToFit
				  << " but there are only " << maxNumEvents << " in total"
				  << std::endl;
		std::cerr << "         Will start at " << fFirstEventToFit
				  << " and fit up to the end of the file" << std::endl;
		fNumEventsToFit = maxNumEvents - fFirstEventToFit;
	}
	else if (numEvents == 0)
	{
		std::cerr
			<< "Warning: Requested to fit 0 events.  Will default to fitting all "
			<< maxNumEvents << " events" << std::endl;
		fNumEventsToFit = maxNumEvents;
	}
	else
	{
		fNumEventsToFit = numEvents;
	}
}

int WCSimFitterConfig::GetNumEventsToFit()
{
	return fNumEventsToFit;
}

void WCSimFitterConfig::Print()
{
	std::cout << " *** WCSimFitterConfig::Print() *** " << std::endl;
	std::cout << "fNumTracks = " << fNumTracks << std::endl;
	std::cout << "fNumEventsToFit = " << fNumEventsToFit << std::endl;
	std::cout << "fNumParameters = " << fNumParameters << std::endl;
	std::cout << "fFirstEventToFit = " << fFirstEventToFit << std::endl;
	std::cout << "fIsPiZeroFit = " << fIsPiZeroFit << std::endl;
	std::cout << "fForcePiZeroMass = " << fForcePiZeroMass << std::endl;
	std::cout << "**************************************" << std::endl;
	//std::cout << "Should add some statements in here..." << std::endl;
}

void WCSimFitterConfig::SetJoinParametersTogether(unsigned int numTrack1,
												  unsigned int numTrack2, const char *name)
{
	assert(numTrack1 != numTrack2);
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	if (numTrack1 < numTrack2)
	{
		fFitterParameters.JoinParametersTogether(numTrack1, numTrack2, type);
	}
	else
	{
		fFitterParameters.JoinParametersTogether(numTrack2, numTrack1, type);
	}
}

bool WCSimFitterConfig::GetJoinParametersTogether(unsigned int numTrack1,
												  unsigned int numTrack2, const char *name)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	return fFitterParameters.GetJoinParametersTogether(numTrack1, numTrack2,
													   type);
}

bool WCSimFitterConfig::GetIsParameterJoined(unsigned int numTrack,
											 FitterParameterType::Type type)
{
	return fFitterParameters.GetIsParameterJoined(numTrack, type);
}

unsigned int WCSimFitterConfig::GetTrackIsJoinedWith(unsigned int numTrack,
													 FitterParameterType::Type type)
{
	return fFitterParameters.GetTrackIsJoinedWith(numTrack, type);
}

bool WCSimFitterConfig::GetIsParameterJoined(unsigned int numTrack,
											 const char *name)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	return fFitterParameters.GetIsParameterJoined(numTrack, type);
}

unsigned int WCSimFitterConfig::GetTrackIsJoinedWith(unsigned int numTrack,
													 const char *name)
{
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	return fFitterParameters.GetTrackIsJoinedWith(numTrack, type);
}

void WCSimFitterConfig::SetTrackType(unsigned int numTrack,
									 const char *typeName)
{
	// std::cout << "WCSimFitterConfig::SetTrackType(" << numTrack << ", " << typeName << ")" << std::endl;
	TrackType::Type type = TrackType::FromName(typeName);
	fFitterParameters.SetTrackType(numTrack, type);
}

TrackType::Type WCSimFitterConfig::GetTrackType(
	const unsigned int &numTrack) const
{
	return fFitterParameters.GetTrackType(numTrack);
}

void WCSimFitterConfig::SetFirstEventToFit(unsigned int iEvt)
{
	unsigned int maxNumEvents = WCSimInterface::GetNumEvents();
	if (iEvt > maxNumEvents)
	{
		std::cerr << "Warning: Requested to start fitting at event " << iEvt
				  << ", but the file only contains " << maxNumEvents << std::endl;
		std::cerr << "         Starting at event 0 instead" << std::endl;
		fFirstEventToFit = 0;
	}
	else if (iEvt + fNumEventsToFit > maxNumEvents)
	{
		std::cerr << "Warning: Requested to fit " << fNumEventsToFit
				  << " events, starting at " << iEvt << std::endl;
		std::cerr << "         But there are only " << maxNumEvents
				  << " - will start at " << iEvt
				  << " and fit up to the end of the file" << std::endl;
		fFirstEventToFit = iEvt;
		SetNumEventsToFit(maxNumEvents - iEvt);
	}
	else
	{
		fFirstEventToFit = iEvt;
	}
}

int WCSimFitterConfig::GetFirstEventToFit() const
{
	return fFirstEventToFit;
}

bool WCSimFitterConfig::GetIsPiZeroFit() const
{
	return fIsPiZeroFit;
}

bool WCSimFitterConfig::GetIsCosmicFit() const
{
	return fIsCosmicFit;
}

bool WCSimFitterConfig::GetMakeFits() const
{
	return fMakeFits;
}

void WCSimFitterConfig::SetMakeFits(const bool &makeFits)
{
	fMakeFits = makeFits;
}

void WCSimFitterConfig::SetIsPiZeroFit(bool isPiZero)
{
	fIsPiZeroFit = isPiZero;
}

void WCSimFitterConfig::SetIsCosmicFit(bool val)
{
	fIsCosmicFit = val;
}

void WCSimFitterConfig::SetForcePiZeroMass(const bool &doIt)
{
	if (!fIsPiZeroFit)
	{
		std::cerr
			<< "Can't force the mass to that of a pi zero without setting WCSimFitterConfig::fIsPiZeroFit first"
			<< std::endl;
		assert(fIsPiZeroFit);
	}
	fForcePiZeroMass = doIt;
}

bool WCSimFitterConfig::GetForcePiZeroMass()
{
	return fForcePiZeroMass;
}
