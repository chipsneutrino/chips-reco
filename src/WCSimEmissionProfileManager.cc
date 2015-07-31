#include "WCSimEmissionProfileManager.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimEmissionProfiles.hh"
#include <cassert>
#include <map>
#include <vector>


#ifndef REFLEX_DICTIONARY
ClassImp(WCSimEmissionProfileManager)
#endif


WCSimEmissionProfileManager::WCSimEmissionProfileManager()
{
  fNumTracksToCache = 4;
}

WCSimEmissionProfileManager::~WCSimEmissionProfileManager()
{
  std::map<TrackType::Type, std::map<double, WCSimEmissionProfiles*> >::iterator outerItr = fProfileMap.begin();
  while(outerItr != fProfileMap.end())
  {
    std::map<double, WCSimEmissionProfiles*>::iterator innerItr = (outerItr->second).begin();
    while(innerItr != outerItr->second.end())
    {
      delete innerItr->second;
      innerItr->second = 0x0;
      ++innerItr;
    }
    (outerItr->second).clear();
    ++outerItr;
  }
  fProfileMap.clear();
  fProfileMap.clear();
}

WCSimEmissionProfiles * WCSimEmissionProfileManager::GetEmissionProfile(TrackType::Type type, double energy)
{

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
  if(numTracks >= 1)
  {
    fNumTracksToCache = numTracks;
  }
  else
  {
    std::cerr << "Error: trying to cache " << numTracks << " tracks - you need to cache at least one" << std::endl;
    fNumTracksToCache = 1;
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
  // std::cout << "Looking for type = " << type << " and energy = " << energy << std::endl;
  // std::cout << "We have..." << std::endl;
  // std::map<TrackType::Type, std::map<double, WCSimEmissionProfiles*> >::iterator outerItr = fProfileMap.begin();
  // for(; outerItr != fProfileMap.end(); ++outerItr)
  // {
    // std::map<double, WCSimEmissionProfiles*>::iterator innerItr = outerItr->second.begin();
    // for ( ; innerItr != outerItr->second.end(); ++innerItr)
    // {
      // std::cout << outerItr->first << "   "  << innerItr->first << std::endl;
    // }

  // }
  // How many sets of profiles do we have cached for this track type?
  // std::cout << "Type = " << type << " and energy = " << energy << std::endl;
  unsigned int numThisType = GetNumCached(type);
  // std::cout << "Num this type = " << numThisType << std::endl;
 
  // If zero, make a new one
  if(numThisType == 0)
  {
    // std::cout << "Adding" << std::endl;
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

    // std::cout << "Replacing" << std::endl;
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
  return GetEmissionProfile(myTrack)->GetStoppingDistance();
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
  return GetEmissionProfile(myTrack)->GetTrackLengthForPercentile(percentile);
}

std::vector<Double_t> WCSimEmissionProfileManager::GetRhoIntegrals(std::vector<Int_t> sPowers, WCSimLikelihoodTrackBase* myTrack, Double_t startS, Double_t endS, Bool_t multiplyByWidth)
{
  return GetEmissionProfile(myTrack)->GetRhoIntegrals(sPowers, startS, endS, multiplyByWidth);
}

std::vector<double> WCSimEmissionProfileManager::GetRhoGIntegrals(WCSimLikelihoodTrackBase* myTrack, WCSimLikelihoodDigit* myDigit, std::vector<int> sPowers, double cutoffS, bool multiplyByWidth)
{
  return GetEmissionProfile(myTrack)->GetRhoGIntegrals(myTrack, myDigit, sPowers, cutoffS, multiplyByWidth);
}
