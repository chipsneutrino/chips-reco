#include "TAxis.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "WCSimEmissionProfileManager.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimEmissionProfiles.hh"
#include "WCSimFastMath.hh"
#include <cassert>
#include <map>
#include <vector>


#ifndef REFLEX_DICTIONARY
ClassImp(WCSimEmissionProfileManager)
#endif


WCSimEmissionProfileManager::WCSimEmissionProfileManager()
{
  fNumTracksToCache = 10;
  fLastNearbyEnergyBin = -999;
  fLastNearbyType = TrackType::Unknown;

  fLastLengthType = TrackType::Unknown;
  fLastLengthEnergy = -999.9;
  fLastLength = -999.9;
}

WCSimEmissionProfileManager::~WCSimEmissionProfileManager()
{
  std::map<TrackType::Type, std::map<double, WCSimEmissionProfiles*> >::iterator outerItr = fProfileMap.begin();
  while(outerItr != fProfileMap.end())
  {
    std::map<double, WCSimEmissionProfiles*>::iterator innerItr = (outerItr->second).begin();
    while(innerItr != outerItr->second.end())
    {
      std::cout << "Deleting " << innerItr->second << std::endl;
      delete innerItr->second;
      innerItr->second = 0x0;
      ++innerItr;
    }
    (outerItr->second).clear();
    ++outerItr;
  }
  fProfileMap.clear();
  for(size_t i = 0; i < fSCosThetaForTimeHists.size(); ++i)
  {
      std::cout << "Deleting fSCosThetaForTimeHists[" << i << "] = " << fSCosThetaForTimeHists[i] << std::endl;
      std::cout << fSCosThetaForTimeHists[i]->GetEntries() << std::endl;
      delete fSCosThetaForTimeHists[i];
      fSCosThetaForTimeHists[i] = 0x0;
  }
  fSCosThetaForTimeHists.clear();
  for(size_t i = 0; i < fSForTimeHists.size(); ++i)
  {
      std::cout << "Deleting fSForTimeHists[" << i << "] = " << fSForTimeHists[i] << std::endl;
      delete fSForTimeHists[i];
      fSForTimeHists[i] = 0x0;
  }
  fSForTimeHists.clear();

}

WCSimEmissionProfiles * WCSimEmissionProfileManager::GetEmissionProfile(TrackType::Type type, double energy)
{
  //  std::cout << "Getting profile for energy " << energy << std::endl;

  if( fProfileMap.find(type) == fProfileMap.end() || fProfileMap[type].find(energy) == fProfileMap[type].end())
  {

    AddNewProfile(type, energy);
  }
  UpdateRecentlyUsed(type, energy);
    
  return fProfileMap[type][energy];
}

void WCSimEmissionProfileManager::UpdateRecentlyUsed(const TrackType::Type &type, const double &energy)
{
  // First make an empty vector if we haven't used any tracks of this type yet
  // std::cout << "Recently used... " << energy << std::endl;


  if(fOrderOfUse.find(type) == fOrderOfUse.end())
  {
    std::vector<double> empty;
    fOrderOfUse[type] = empty;
  }

  // std::cout << "Vector has size " << fOrderOfUse[type].size() << " and contains" << std::endl;
  // if(fOrderOfUse[type].size() > 0)
  // {
  //   for(unsigned int iUsed = 0; iUsed < fOrderOfUse[type].size(); ++iUsed)
  //   {
  //     std::cout << fOrderOfUse[type].at(iUsed) << std::endl;
  //   }
  // }
  // See if this energy has been used before
  std::vector<double>::iterator energyItr = std::find(fOrderOfUse[type].begin(), fOrderOfUse[type].end(), energy);
  if(energyItr == fOrderOfUse[type].end())
  {
    // It hasn't, so add it to the vector
    fOrderOfUse[type].push_back(energy);
  }
  else
  {
    // It already exists: we need to remove it from its place and then put it
    // back as the last element
    energyItr = fOrderOfUse[type].erase(energyItr);
    fOrderOfUse[type].push_back(energy);
  }

  // std::cout << "Now vector has size " << fOrderOfUse[type].size() << " and contains" << std::endl;

  // for(unsigned int iUsed = 0; iUsed < fOrderOfUse[type].size(); ++iUsed)
  // {
  //   std::cout << fOrderOfUse[type].at(iUsed) << std::endl;

  // }
  return;
}

void WCSimEmissionProfileManager::RemoveFromRecentEnergies(const TrackType::Type &type, const double &energy)
{
  std::map<TrackType::Type, std::vector<double> >::iterator orderItr;
  orderItr = fOrderOfUse.find(type);
  assert(orderItr != fOrderOfUse.end());

  std::vector<double>::iterator energySearch = std::find((orderItr->second).begin(), (orderItr->second).end(), energy);
  assert(energySearch != (orderItr->second).end());
  (orderItr->second).erase(energySearch);

  return;
}
  


WCSimEmissionProfiles * WCSimEmissionProfileManager::GetEmissionProfile(WCSimLikelihoodTrackBase * myTrack)
{
  return GetEmissionProfile(myTrack->GetType(), myTrack->GetE());
}

void WCSimEmissionProfileManager::SetNumTracksToCache(const unsigned int numTracks)
{
  if(numTracks >= 4)
  {
    fNumTracksToCache = numTracks;
  }
  else
  {
    std::cerr << "Error: trying to cache " << numTracks << " tracks - you need to cache at least four for the time likelihood" << std::endl;
    fNumTracksToCache = 4;
  }
  std::map<TrackType::Type, std::map<double, WCSimEmissionProfiles*> >::iterator mapItr = fProfileMap.begin();
  while( (mapItr->second).size() > fNumTracksToCache)
  {
    double energy = GetLeastRecentlyUsedEnergy(mapItr->first);
    WCSimEmissionProfiles * temp = (mapItr->second)[energy];
    mapItr->second.erase(energy);
    delete temp;
    RemoveFromRecentEnergies(mapItr->first, energy);
  }
}

void WCSimEmissionProfileManager::CacheAtLeast(const unsigned int numTracks)
{
  if( numTracks > fNumTracksToCache )
  {
    SetNumTracksToCache(numTracks);
  }
}

unsigned int WCSimEmissionProfileManager::GetNumTracksToCache() const
{
  return fNumTracksToCache;
}

void WCSimEmissionProfileManager::AddNewProfile(const TrackType::Type &type, const double &energy)
{
  // How many sets of profiles do we have cached for this track type?
  unsigned int numThisType = GetNumCached(type);
 
  // If zero, make a new one
  if(numThisType == 0)
  {
    fProfileMap[type][energy] = new WCSimEmissionProfiles(type, energy);
  }

  // If > 0 but < number to cache, make a new one **TODO** would be nice to copy
  // an old one instead of opening a bunch of new files here
  if(numThisType > 0 && numThisType < fNumTracksToCache)
  {
    // std::cout << "Adding 2 " << std::endl;
    fProfileMap[type][energy] = new WCSimEmissionProfiles(type, energy);
  }
  
  // If  == number to cache, adapt the energy in one of them (the last used)
  if(numThisType == fNumTracksToCache)
  {

    double replaceEnergy = GetLeastRecentlyUsedEnergy(type);
    WCSimEmissionProfiles * changeMe = fProfileMap[type][replaceEnergy];

    fProfileMap[type].erase(replaceEnergy);
    changeMe->SetEnergy(energy);
    fProfileMap[type][energy] = changeMe;
    fOrderOfUse[type].erase(fOrderOfUse[type].begin());
    UpdateRecentlyUsed(type, energy);
  }
}

double WCSimEmissionProfileManager::GetLeastRecentlyUsedEnergy(const TrackType::Type &type)
{
  assert(fOrderOfUse.find(type) != fOrderOfUse.end());
  return fOrderOfUse[type].at(0);
}

double WCSimEmissionProfileManager::GetMostRecentlyUsedEnergy(const TrackType::Type &type)
{
  assert(fOrderOfUse.find(type) != fOrderOfUse.end());
  assert(fOrderOfUse[type].size() > 0);
  return fOrderOfUse[type].at(fOrderOfUse[type].size() - 1);
}

std::vector<double> WCSimEmissionProfileManager::GetFourNearestEnergies(
		WCSimLikelihoodTrackBase* myTrack) 
{
    return GetFourNearestEnergies(myTrack->GetType(), myTrack->GetE());
}


std::vector<double> WCSimEmissionProfileManager::GetFourNearestEnergies(
        const TrackType::Type& type, const double& energy)
{

	// Construct a vector of energies so that the energy we want lies between
	// index 1 and 2:
    // std::cout << "Get fEnergies" << std::endl;
	std::vector<double> energies;
	TAxis * energyAxis = GetEnergyHist(type)->GetXaxis();
	int bin = energyAxis->FindBin(energy);
	if(bin != fLastNearbyEnergyBin || fLastNearbyType != type)
	{
		ResetSCosThetaForTimeHists();
		ResetSForTimeHists();
		fEnergies.clear();
		fEnergies.resize(4, 0.0);

        if(bin == 0)
        {
            double lowVal = energyAxis->GetBinLowEdge(1);
            double width = energyAxis->GetBinLowEdge(1);
            while( lowVal >= (energy - width) )
            {
                lowVal = energy - width;
            }
            fEnergies.at(0) = lowVal;
            fEnergies.at(1) = lowVal +   width;
            fEnergies.at(2) = lowVal + 2*width;
            fEnergies.at(3) = lowVal + 3*width;

        }
		else if(bin == energyAxis->GetNbins())
		{
            // We're in the final energy bin, which goes from our highest simulated energy to 30GeV
            // Just want to return four equally-spaced energies inside this bin such that the 
            // track energy is between the middle two
            double binWidth = energyAxis->GetBinWidth(1);
            fEnergies.at(1) = floor(energy / binWidth) * binWidth;
            fEnergies.at(0) = fEnergies.at(1) - binWidth;
            fEnergies.at(2) = fEnergies.at(1) + binWidth;
            fEnergies.at(3) = fEnergies.at(2) + binWidth;
		}
        else 
        {
            if(bin == 1)
		    {
		    	fEnergies.at(0) = energyAxis->GetBinLowEdge(bin) - energyAxis->GetBinWidth(bin);
		        fEnergies.at(1) = energyAxis->GetBinLowEdge(bin);
		    }
		    else
		    {
		    	fEnergies.at(0) = energyAxis->GetBinLowEdge(bin-1);
		        fEnergies.at(1) = energyAxis->GetBinLowEdge(bin);
		    }

		    if(bin == energyAxis->GetNbins()-1)
		    {
		    	fEnergies.at(2) = energyAxis->GetBinLowEdge(bin+1);
		    	fEnergies.at(3) = 2*fEnergies.at(2) - fEnergies.at(1);
		    }
		    else
		    {
		    	fEnergies.at(2) = energyAxis->GetBinLowEdge(bin+1);
		    	fEnergies.at(3) = energyAxis->GetBinLowEdge(bin+2);
		    }
        }
		fLastNearbyEnergyBin = bin;
		fLastNearbyType = type;
	}
	return fEnergies;
}

std::vector<TH1F*> WCSimEmissionProfileManager::GetFourNearestSForTimeHists(
		WCSimLikelihoodTrackBase* myTrack) {
    // std::cout << "Get nearest s for time" << std::endl;

	std::vector<double> energies = GetFourNearestEnergies(myTrack);
	if(fSForTimeHists.size() == 0)
	{
		fSForTimeHists.resize(4, NULL);
		for(size_t i = 0; i < energies.size(); ++i)
		{
			fSForTimeHists[i] = (TH1F*)GetEmissionProfile(myTrack->GetType(), energies[i])->GetSForTime()->Clone();
            fSForTimeHists[i]->SetDirectory(0);
		}
	}
	return fSForTimeHists;
}

std::vector<TH2F*> WCSimEmissionProfileManager::GetFourNearestSCosThetaForTimeHists(
		WCSimLikelihoodTrackBase* myTrack) {

   	std::vector<double> energies = GetFourNearestEnergies(myTrack);
    if(fSCosThetaForTimeHists.size() == 0)
    {
    	ResetSCosThetaForTimeHists();
    	fSCosThetaForTimeHists.resize(energies.size(), NULL);

        for(size_t i = 0; i < energies.size(); ++i)
        {
            fSCosThetaForTimeHists[i] = (TH2F*)GetEmissionProfile(myTrack->GetType(), energies[i])->GetSCosThetaForTime()->Clone();
            fSCosThetaForTimeHists[i]->SetDirectory(0);
        }
    }
    assert(fSCosThetaForTimeHists.size() > 0);
    return fSCosThetaForTimeHists;
}

void WCSimEmissionProfileManager::ResetSForTimeHists() {
	for(size_t i = 0; i < fSForTimeHists.size(); ++i)
	{
		delete fSForTimeHists[i];
	}
	fSForTimeHists.clear();
}

void WCSimEmissionProfileManager::ResetSCosThetaForTimeHists() {
	for(size_t i = 0 ; i < fSCosThetaForTimeHists.size(); ++i)
	{
		delete fSCosThetaForTimeHists[i];
	}
	fSCosThetaForTimeHists.clear();
}

unsigned int WCSimEmissionProfileManager::GetNumCached(const TrackType::Type &type)
{
  unsigned int size = 0;
  if(fProfileMap.find(type) != fProfileMap.end())
  {
    size = fProfileMap[type].size(); 
  }
  return size;
}

double WCSimEmissionProfileManager::GetStoppingDistance(WCSimLikelihoodTrackBase * myTrack)
{
    double dist = (4.24226 + 0.516017*myTrack->GetE() - 5.08037e-6 * myTrack->GetE() * myTrack->GetE());
    if(dist < 50) { dist = 50; }
    return dist;
  //return GetEmissionProfile(myTrack)->GetStoppingDistance();
}

TH1F * WCSimEmissionProfileManager::GetEnergyHist(const TrackType::Type &type)
{
  std::map<TrackType::Type, std::map<double, WCSimEmissionProfiles*> >::iterator mapItr;
  mapItr = fProfileMap.find(type);
  if(mapItr != fProfileMap.end())
  {
    double energy = GetMostRecentlyUsedEnergy(type);
    return GetEmissionProfile(type, energy)->GetEnergyHist();
  }
  
  return GetEmissionProfile(type, 1000.0)->GetEnergyHist();
}

TH1F * WCSimEmissionProfileManager::GetRho(const TrackType::Type &type, const double &energy)
{
  return GetEmissionProfile(type, energy)->GetRho();
}

std::pair<TH2F *, TH2F *>  WCSimEmissionProfileManager::GetG(const TrackType::Type &type, const double &energy)
{
  return GetEmissionProfile(type, energy)->GetG();
}

double WCSimEmissionProfileManager::GetLightFlux(WCSimLikelihoodTrackBase * myTrack)
{
  return GetEmissionProfile(myTrack)->GetLightFlux(myTrack->GetType(), myTrack->GetE());
}

double WCSimEmissionProfileManager::GetTrackLengthForPercentile(WCSimLikelihoodTrackBase * myTrack, double percentile)
{
  double dist =  (-26.0305 + 0.326113*myTrack->GetE() - 3.86161e-6 * myTrack->GetE() * myTrack->GetE()) * percentile/0.70;
  if(dist < 50){ dist = 50; }
  return dist;
  return GetEmissionProfile(myTrack)->GetTrackLengthForPercentile(percentile);
}

std::vector<Double_t> WCSimEmissionProfileManager::GetRhoIntegrals(std::vector<Int_t> sPowers, WCSimLikelihoodTrackBase* myTrack, Double_t startS, Double_t endS, Bool_t multiplyByWidth)
{
  return GetEmissionProfile(myTrack)->GetRhoIntegrals(myTrack, sPowers, startS, endS, multiplyByWidth);
}

std::vector<double> WCSimEmissionProfileManager::GetRhoGIntegrals(WCSimLikelihoodTrackBase* myTrack, WCSimLikelihoodDigit* myDigit, std::vector<int> sPowers, double cutoffS, bool multiplyByWidth)
{
  return GetEmissionProfile(myTrack)->GetRhoGIntegrals(myTrack, myDigit, sPowers, cutoffS, multiplyByWidth);
}

bool WCSimEmissionProfileManager::EnergiesInSameBin(const TrackType::Type& type, const double energy1, const double energy2)
{
    TH1F * energyHist = GetEnergyHist(type);
    return (energyHist->GetXaxis()->FindBin(energy1) == energyHist->GetXaxis()->FindBin(energy2));
}

double WCSimEmissionProfileManager::GetTrackLength(const TrackType::Type& type, const double& energy)
{
    if(!(fLastLengthType == type && fLastLengthEnergy == energy))
    {
        std::vector<double> energies = GetFourNearestEnergies(type, energy);
        std::vector<double> lengths(energies.size());
        for(size_t i = 0; i < energies.size(); ++i)
        {
            lengths[i] = GetEmissionProfile(type, energy)->GetTrackLengthForPercentile(99);
        }
        fLastLengthType = type;
        fLastLengthEnergy = energy;
        fLastLength = WCSimFastMath::CatmullRomSpline(&(energies[0]), &(lengths[0]), energy);
    }
    return fLastLength;

}
