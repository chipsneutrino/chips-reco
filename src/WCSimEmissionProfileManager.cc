#include "TAxis.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
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
  for(size_t i = 0; i < fSCosThetaForTimeHists.size(); ++i)
  {
      delete fSCosThetaForTimeHists[i];
      fSCosThetaForTimeHists[i] = 0x0;
  }
  fSCosThetaForTimeHists.clear();
  for(size_t i = 0; i < fSForTimeHists.size(); ++i)
  {
      delete fSForTimeHists[i];
      fSForTimeHists[i] = 0x0;
  }
  fSForTimeHists.clear();

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

std::vector<double> WCSimEmissionProfileManager::GetNearestEnergies(
		WCSimLikelihoodTrackBase* myTrack) {

	// Construct a vector of energies so that the energy we want lies between
	// index 1 and 2:
    // std::cout << "Get fEnergies" << std::endl;
	if(fEnergies.size() == 0)
    {
	    TAxis * energyAxis = GetEnergyHist(myTrack->GetType())->GetXaxis();
	    assert(energyAxis->GetNbins() >= 99);
        for(int bin = 0; bin < energyAxis->GetNbins(); ++bin)
        {
            fEnergies.push_back(energyAxis->GetBinLowEdge(bin+1));
        }
    }
    return fEnergies;
/*
	if(bin == 1)
	{
		energies.at(0) = energyAxis->GetBinLowEdge(bin) - energyAxis->GetBinWidth(bin);
	}
	else
	{
		energies.at(0) = energyAxis->GetBinLowEdge(bin-1);
	}
	energies.at(1) = energyAxis->GetBinLowEdge(bin);

	if(bin == energyAxis->GetNbins())
	{
		energies.at(2) = energyAxis->GetBinLowEdge(bin) + energyAxis->GetBinWidth(bin);
		energies.at(3) = 2*energies.at(2) - energies.at(1);
	}
	else if(bin == energyAxis->GetNbins()-1)
	{
		energies.at(2) = energyAxis->GetBinLowEdge(bin+1);
		energies.at(3) = 2*energies.at(2) - energies.at(1);
	}
	else
	{
		energies.at(2) = energyAxis->GetBinLowEdge(bin+1);
		energies.at(3) = energyAxis->GetBinLowEdge(bin+2);
	}
	return energies;
    */
}

std::vector<TH1F*> WCSimEmissionProfileManager::GetNearestSForTimeHists(
		WCSimLikelihoodTrackBase* myTrack) {
    // std::cout << "Get nearest s for time" << std::endl;

    if(fSForTimeHists.size() == 0)
    {
	    TAxis * energyAxis = GetEnergyHist(myTrack->GetType())->GetXaxis();
	    assert(fNumTracksToCache >= 4);
	    assert(energyAxis->GetNbins() >= 4);
	    std::vector<TH1F*> hists(energyAxis->GetNbins(), NULL);

	    std::vector<double> energies;
        energies.reserve(energyAxis->GetNbins());
        for(int bin = 1; bin <= energyAxis->GetNbins(); ++bin)
        {
            // std::cout << "Energy bin " << bin << "/" << energyAxis->GetNbins() << std::endl;
            energies.push_back(energyAxis->GetBinLowEdge(bin));
        }
        for(size_t i = 0; i < energies.size(); ++i)
        {
            fSForTimeHists.push_back((TH1F*)GetEmissionProfile(myTrack->GetType(), energies.at(i))->GetSForTime()->Clone());
        }
    }
    assert(fSForTimeHists.size() > 0);
    return fSForTimeHists;

    /*
	TAxis * energyAxis = GetEnergyHist(myTrack->GetType())->GetXaxis();
	assert(fNumTracksToCache >= energyAxis->GetNbins());
	assert(energyAxis->GetNbins() >= 4);
	//int bin = energyAxis->FindBin(myTrack->GetE());

	std::vector<TH1F*> hists(energyAxis->GetNbins(), NULL);

	std::vector<double> energies(energyAxis->GetNbins());
    for(int bin = 1; bin <= energyAxis->GetNbins(); ++bin)
    {
        energies.push_back(energyAxis->GetBinLowEdge(bin));
    }
    for(size_t i = 0; i < energies.size(); ++i)
    {
        hists[i] = GetEmissionProfile(myTrack->GetType(), energies.at(i))->GetSForTime();
    }


	if(bin == 1)
	{
		energies.push_back(energyAxis->GetBinLowEdge(1));
		energies.push_back(energyAxis->GetBinLowEdge(1));
		energies.push_back(energyAxis->GetBinLowEdge(2));
		energies.push_back(energyAxis->GetBinLowEdge(3));
	}
	else if(bin == energyAxis->GetNbins()-1)
	{
		energies.push_back(energyAxis->GetBinLowEdge(bin-1));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin+1));
		energies.push_back(energyAxis->GetBinLowEdge(bin+1));
	}
	else if(bin == energyAxis->GetNbins())
	{
		energies.push_back(energyAxis->GetBinLowEdge(bin-1));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
	}
	else
	{
		energies.push_back(energyAxis->GetBinLowEdge(bin-1));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin+1));
		energies.push_back(energyAxis->GetBinLowEdge(bin+2));
	}
	for(size_t i = 0; i < 4; ++i)
	{
		hists[i] = GetEmissionProfile(myTrack->GetType(), energies.at(i))->GetSForTime();
	}
	return hists;
*/
}

std::vector<TH2F*> WCSimEmissionProfileManager::GetNearestSCosThetaForTimeHists(
		WCSimLikelihoodTrackBase* myTrack) {
    // std::cout << "Get nearest sCosTheta for time" << std::endl;

    if(fSCosThetaForTimeHists.size() == 0)
    {
	    TAxis * energyAxis = GetEnergyHist(myTrack->GetType())->GetXaxis();
	    assert(fNumTracksToCache >= 4);
	    assert(energyAxis->GetNbins() >= 4);
	    std::vector<TH2F*> hists(energyAxis->GetNbins(), NULL);

	    std::vector<double> energies;
        energies.reserve(energyAxis->GetNbins());
        for(int bin = 1; bin <= energyAxis->GetNbins(); ++bin)
        {
            energies.push_back(energyAxis->GetBinLowEdge(bin));
        }
        for(size_t i = 0; i < energies.size(); ++i)
        {
            // std::cout << "Making fSCosThetaForTimeHists, bin " << i << "/" << energies.size() << " total size " << fSCosThetaForTimeHists.size() << std::endl;
            fSCosThetaForTimeHists.push_back((TH2F*)GetEmissionProfile(myTrack->GetType(), energies.at(i))->GetSCosThetaForTime()->Clone());
        }
    }
    assert(fSCosThetaForTimeHists.size() > 0);
    return fSCosThetaForTimeHists;



    /*
	int bin = energyAxis->FindBin(myTrack->GetE());

	std::vector<TH2F*> hists(4, NULL);

	std::vector<double> energies;
    energies.reserve(4);
	if(bin == 1)
	{
		energies.push_back(energyAxis->GetBinLowEdge(1));
		energies.push_back(energyAxis->GetBinLowEdge(1));
		energies.push_back(energyAxis->GetBinLowEdge(2));
		energies.push_back(energyAxis->GetBinLowEdge(3));
	}
	else if(bin == energyAxis->GetNbins()-1)
	{
		energies.push_back(energyAxis->GetBinLowEdge(bin-1));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin+1));
		energies.push_back(energyAxis->GetBinLowEdge(bin+1));
	}
	else if(bin == energyAxis->GetNbins())
	{
		energies.push_back(energyAxis->GetBinLowEdge(bin-1));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
	}
	else
	{
		energies.push_back(energyAxis->GetBinLowEdge(bin-1));
		energies.push_back(energyAxis->GetBinLowEdge(bin));
		energies.push_back(energyAxis->GetBinLowEdge(bin+1));
		energies.push_back(energyAxis->GetBinLowEdge(bin+2));
	}

	for(size_t i = 0; i < 4; ++i)
	{
		hists[i] = GetEmissionProfile(myTrack->GetType(), energies.at(i))->GetSCosThetaForTime();
        // TCanvas * can = new TCanvas(TString::Format("bin_%d", (int)i).Data(), "", 1200,600);
        // can->Divide(2,1);
        // can->cd(1);
        // hists[i]->Draw("COLZ");
        // can->cd(2);
        // std::cout << "Energy is " << energies.at(i) << std::endl;
        // GetEmissionProfile(myTrack->GetType(), energies.at(i))->GetRho()->Draw();
        // can->SaveAs(TString::Format("bin_%d.png", (int)i).Data());
	}
	return hists;
    */
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
