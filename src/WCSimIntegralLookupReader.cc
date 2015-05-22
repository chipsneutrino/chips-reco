/*
 * WCSimIntegralLookupReader.cc
 *
 *  Created on: 12 Mar 2015
 *      Author: ajperch
 */

#include "WCSimIntegralLookupReader.hh"
#include "WCSimAnalysisConfig.hh"
#include "TH1F.h"
#include "TH2D.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimIntegralLookupReader)
#endif

static WCSimIntegralLookupReader* fgIntegralLookupReader = 0x0;

WCSimIntegralLookupReader* WCSimIntegralLookupReader::Instance() {
	  if(!fgIntegralLookupReader){
		  fgIntegralLookupReader = new WCSimIntegralLookupReader();
	  }

	  // die if finder hasn't actually been created
	  if(!fgIntegralLookupReader){
	    assert(fgIntegralLookupReader);
	  }

	  return fgIntegralLookupReader;
}

void WCSimIntegralLookupReader::LoadIntegrals(WCSimLikelihoodTrack* myTrack) {
	LoadIntegrals(myTrack->GetType());
}

WCSimIntegralLookup3D* WCSimIntegralLookupReader::GetIntegralLookup3D(WCSimLikelihoodTrack::TrackType type) {
	if( fLookupMap.find( type ) != fLookupMap.end())
	{
		return fLookupMap[type];
	}
	else
	{
		LoadIntegrals(type);
		if(fLookupMap.find(type) != fLookupMap.end())
		{
			return fLookupMap[type];
		}
		else
		{
			std::cerr << "Could not open an integral file for track type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
			assert(0);
		}
	}

	return 0x0;
}

WCSimIntegralLookupReader::WCSimIntegralLookupReader() {
	// TODO Auto-generated constructor stub
  fLastPercentileLength = -999.9;
  fLastPercentileType = WCSimLikelihoodTrack::Unknown;
  fLastPercentileEnergy = -999.9;
}

WCSimIntegralLookupReader::~WCSimIntegralLookupReader() {
	// TODO Auto-generated destructor stub
	std::map<WCSimLikelihoodTrack::TrackType, WCSimIntegralLookup3D*>:: iterator mapItr;
	for( mapItr = fLookupMap.begin(); mapItr != fLookupMap.end(); ++mapItr)
	{
		delete mapItr->second;
		mapItr->second = 0x0;
	}
	fLookupMap.clear();
}

TString WCSimIntegralLookupReader::GetLookupFilename(WCSimLikelihoodTrack* track) {
	return GetLookupFilename(track->GetType());
}

void WCSimIntegralLookupReader::LoadIntegrals(const WCSimLikelihoodTrack::TrackType& type) {
	if( fLookupMap.find( type ) == fLookupMap.end())
	{
		if(WCSimAnalysisConfig::Instance()->GetTruncateIntegrals())
		{
			WCSimIntegralLookup * myLookup = new WCSimIntegralLookup(GetLookupFilename(type));
			fLookupMap.insert(std::make_pair(type, (WCSimIntegralLookup3D*)myLookup));
		}
		else
		{
			WCSimIntegralLookup3D * myLookup = new WCSimIntegralLookup3D(GetLookupFilename(type));
			fLookupMap.insert(std::make_pair(type, myLookup));
		}
	}
}

double WCSimIntegralLookupReader::GetRhoIntegral(const WCSimLikelihoodTrack::TrackType& type, const double& E,
		const double& s) {
	if(fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if(fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	else
	{
		return fLookupMap[type]->GetRhoIntegral(E, s);
	}
}

double WCSimIntegralLookupReader::GetRhoSIntegral(const WCSimLikelihoodTrack::TrackType& type, const double& E,
		const double& s) {
	if(fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if(fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoSIntegral(E, s);
}

double WCSimIntegralLookupReader::GetRhoSSIntegral(const WCSimLikelihoodTrack::TrackType& type, const double& E,
		const double& s) {
	if(fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if(fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoSSIntegral(E, s);
}

double WCSimIntegralLookupReader::GetRhoGIntegral(const WCSimLikelihoodTrack::TrackType& type, const double& E,
		const double& s, const double& R0, const double& cosTh0) {
	if(fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if(fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoGIntegral(E, s, R0, cosTh0);
}

double WCSimIntegralLookupReader::GetRhoGSIntegral(const WCSimLikelihoodTrack::TrackType& type, const double& E,
		const double& s, const double& R0, const double& cosTh0) {
	if(fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if(fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoGSIntegral(E, s, R0, cosTh0);
}

double WCSimIntegralLookupReader::GetRhoGSSIntegral(const WCSimLikelihoodTrack::TrackType& type, const double& E,
		const double& s, const double& R0, const double& cosTh0) {
	if(fLookupMap.find(type) == fLookupMap.end())
	{
		LoadIntegrals(type);
	}

	if(fLookupMap.find(type) == fLookupMap.end())
	{
		std::cerr << "Could not find a table of integrals for type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
		assert(fLookupMap.find(type) != fLookupMap.end());
	}
	return fLookupMap[type]->GetRhoGSSIntegral(E, s, R0, cosTh0);
}

TString WCSimIntegralLookupReader::GetLookupFilename(const WCSimLikelihoodTrack::TrackType& type) {
	TString str;
	if( type == WCSimLikelihoodTrack::MuonLike)
	{
		if(WCSimAnalysisConfig::Instance()->GetTruncateIntegrals())
		{
			str = "config/muonIntegrals.root";
		}
		else
		{
			str = "config/muonIntegralsSmall.root";
		}
	}
	else if( type == WCSimLikelihoodTrack::ElectronLike)
	{
		if(WCSimAnalysisConfig::Instance()->GetTruncateIntegrals())
		{
			str = "config/electronIntegrals.root";
		}
		else
		{
			str = "config/electronIntegralsSmall.root";
		}
	}
	else
	{
		std::cerr << "I don't know the filename for a track of type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
		assert( type == WCSimLikelihoodTrack::MuonLike || type == WCSimLikelihoodTrack::ElectronLike);
	}
	std::cout << "Lookup file name = " << str << std::endl;
	return str;
}

double WCSimIntegralLookupReader::GetTrackLengthForPercentile(const WCSimLikelihoodTrack::TrackType &type, const double &E, const double &percentile)
{
  double dist = fLastPercentileLength;
  if(type != fLastPercentileType || E != fLastPercentileEnergy)
  {

	  if(fLookupMap.find(type) == fLookupMap.end())
	  {
	  	LoadIntegrals(type);
	  }

	  if(fLookupMap.find(type) == fLookupMap.end())
	  {
	  	std::cerr << "Could not find a table of integrals for type " << WCSimLikelihoodTrack::TrackTypeToString(type) << std::endl;
	  	assert(fLookupMap.find(type) != fLookupMap.end());
	  }

    TH2D * tempRhoInt = ((TH2D*)fLookupMap[type]->GetRhoIntegralHist()->Projection(0,1));
    int bin = tempRhoInt->GetYaxis()->FindBin(E);
    TH1D * tempRhoInt1D = tempRhoInt->ProjectionX("tempRhoInt1D",bin, bin);
    dist = tempRhoInt1D->GetBinCenter(tempRhoInt1D->FindFirstBinAbove(percentile));
    delete tempRhoInt1D;
    delete tempRhoInt;
    std::cout << "Length (last) = " << dist << " (" << fLastPercentileLength << ")" << std::endl;
    fLastPercentileType = type;
    fLastPercentileLength = dist;
    fLastPercentileEnergy = E;
  }
  return dist;
}
 
