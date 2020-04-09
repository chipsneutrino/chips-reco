/*
 * WCSimIntegralLookupReader.cc
 *
 *  Created on: 12 Mar 2015
 *      Author: ajperch
 */

#include "WCSimIntegralLookupReader.hh"
#include "WCSimParameters.hh"
#include "WCSimTrackParameterEnums.hh"
#include "TH1F.h"
#include "TH2D.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimIntegralLookupReader)
#endif

static WCSimIntegralLookupReader *fgIntegralLookupReader = 0x0;

WCSimIntegralLookupReader *WCSimIntegralLookupReader::Instance()
{
	if (!fgIntegralLookupReader)
	{
		fgIntegralLookupReader = new WCSimIntegralLookupReader();
	}

	// die if finder hasn't actually been created
	if (!fgIntegralLookupReader)
	{
		assert(fgIntegralLookupReader);
	}

	return fgIntegralLookupReader;
}

void WCSimIntegralLookupReader::LoadIntegrals(WCSimLikelihoodTrackBase *myTrack)
{
	LoadIntegrals(myTrack->GetType());
}

WCSimIntegralLookup3D *WCSimIntegralLookupReader::GetIntegralLookup3D(TrackType::Type type)
{
	if (fLookupMap.find(type) != fLookupMap.end())
	{
		return fLookupMap[type];
	}
	else
	{
		LoadIntegrals(type);
		if (fLookupMap.find(type) != fLookupMap.end())
		{
			return fLookupMap[type];
		}
		else
		{
			std::cerr << "Could not open an integral file for track type " << TrackType::AsString(type) << std::endl;
			assert(0);
		}
	}

	return 0x0;
}

WCSimIntegralLookupReader::WCSimIntegralLookupReader()
{
	// TODO Auto-generated constructor stub
}

WCSimIntegralLookupReader::~WCSimIntegralLookupReader()
{
	// TODO Auto-generated destructor stub
	std::map<TrackType::Type, WCSimIntegralLookup3D *>::iterator mapItr;
	for (mapItr = fLookupMap.begin(); mapItr != fLookupMap.end(); ++mapItr)
	{
		delete mapItr->second;
		mapItr->second = 0x0;
	}
	fLookupMap.clear();
}

TString WCSimIntegralLookupReader::GetLookupFilename(WCSimLikelihoodTrackBase *track)
{
	return GetLookupFilename(track->GetType());
}

void WCSimIntegralLookupReader::LoadIntegrals(const TrackType::Type &type)
{
	std::cout << "Loading integrals" << std::endl;
	if (fLookupMap.find(type) == fLookupMap.end())
	{
		std::cout << "Making new truncated integral lookup" << std::endl;
		WCSimIntegralLookup3D *myLookup = new WCSimIntegralLookup3D(GetLookupFilename(type));
		fLookupMap.insert(std::make_pair(type, myLookup));
	}
	std::cout << "Loaded integrals" << std::endl;
}

double WCSimIntegralLookupReader::GetRhoIntegral(const TrackType::Type &type, const double &E, const double &s)
{
	if (fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if (fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << TrackType::AsString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoIntegral(E, s);
}

double WCSimIntegralLookupReader::GetRhoSIntegral(const TrackType::Type &type, const double &E, const double &s)
{
	if (fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if (fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << TrackType::AsString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoSIntegral(E, s);
}

double WCSimIntegralLookupReader::GetRhoSSIntegral(const TrackType::Type &type, const double &E, const double &s)
{
	if (fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if (fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << TrackType::AsString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoSSIntegral(E, s);
}

double WCSimIntegralLookupReader::GetRhoGIntegral(const TrackType::Type &type, const double &E, const double &s,
												  const double &R0, const double &cosTh0)
{
	if (fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if (fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << TrackType::AsString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoGIntegral(E, s, R0, cosTh0);
}

double WCSimIntegralLookupReader::GetRhoGSIntegral(const TrackType::Type &type, const double &E, const double &s,
												   const double &R0, const double &cosTh0)
{
	// std::cout << "Looking for type " << type << std::endl;
	if (fLookupMap.find(type) == fLookupMap.end())
	{
		// std::cout << "Need to load them" << std::endl;
		LoadIntegrals(type);
		// std::cout << "Loaded" << std::endl;
	}
	else
	{
		// std::cout << "Already exists" << std::endl;
	}

	if (fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << TrackType::AsString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	// std::cout << fLookupMap.size() << std::endl;
	// std::cout << fLookupMap[type] << std::endl;
	// std::cout << fLookupMap[type]->GetRhoGIntegral(E, s, R0, cosTh0) << std::endl;
	// std::cout << "Trying to return the thing" << std::endl;
	return fLookupMap[type]->GetRhoGSIntegral(E, s, R0, cosTh0);
}

double WCSimIntegralLookupReader::GetRhoGSSIntegral(const TrackType::Type &type, const double &E, const double &s,
													const double &R0, const double &cosTh0)
{
	if (fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if (fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << TrackType::AsString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoGSSIntegral(E, s, R0, cosTh0);
}

TString WCSimIntegralLookupReader::GetLookupFilename(const TrackType::Type &type)
{
	TString str(getenv("CHIPSRECO"));
	if (type == TrackType::MuonLike)
	{
		if (WCSimParameters::Instance()->TruncateIntegrals())
		{
			str += "/config/muonIntegrals.root";
		}
		else
		{
			//str += "/config/electronIntegralsSmall.root";
			str += "/config/muonIntegralsSmall.root";
		}
	}
	else if (type == TrackType::ElectronLike || type == TrackType::PhotonLike)
	{
		if (WCSimParameters::Instance()->TruncateIntegrals())
		{
			str += "/config/electronIntegrals.root";
		}
		else
		{
			str += "/config/electronIntegralsSmall.root";
		}
	}
	else
	{
		std::cerr << "I don't know the filename for a track of type " << TrackType::AsString(type) << std::endl;
		assert(type == TrackType::MuonLike || type == TrackType::ElectronLike);
	}
	std::cout << "Lookup file name = " << str << std::endl;
	return str;
}

void WCSimIntegralLookupReader::SaveIntegrals(const TrackType::Type &type, const double E, const double s,
											  const double R0, const double cosTh0)
{
	fLookupMap[type]->SaveIntegrals(E, s, R0, cosTh0);
}
