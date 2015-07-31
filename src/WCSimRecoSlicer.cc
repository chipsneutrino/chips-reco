#include <iostream>
#include <vector>
#include <algorithm> 
#include <cassert>
#include <cmath>

#include "WCSimRecoSlicer.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimParameters.hh"

// A little sort function for events
bool SortTheSlices(WCSimRecoEvent* a, WCSimRecoEvent *b);

// A little sort function for digits
bool SortDigitsByCharge(WCSimRecoDigit* a, WCSimRecoDigit *b);

WCSimRecoSlicer::WCSimRecoSlicer(){
  fInputEvent = 0x0;
  this->UpdateParameters();
}

WCSimRecoSlicer::~WCSimRecoSlicer(){
  // Created events must be deleted by the user at a later point
}

void WCSimRecoSlicer::UpdateParameters(){
  WCSimParameters* pars=WCSimParameters::Instance();
  fDistanceCut = pars->GetSlicerClusterDistance();
  fMinSliceSize = pars->GetSlicerMinSize();
  fChargeCut = pars->GetSlicerChargeCut();
}

void WCSimRecoSlicer::SetInputEvent(WCSimRecoEvent* evt){
  fInputEvent = evt;

}

void WCSimRecoSlicer::SliceTheEvent(){

  if(fInputEvent == 0x0){
    std::cerr << "WCSimRecoSlicer: Can't slice an empty event!" << std::endl;
    assert(0);
  }

  if(!fInputEvent->IsFilterDone()){
    std::cerr << "WCSimRecoSlicer: Need to apply filter to digits before slicing." << std::endl;
    assert(0);
  }

  std::vector<WCSimRecoDigit*> digitVector;
  // Make a charge sorted vector of those digits above the charge threshold
  for(unsigned int v = 0; v < fInputEvent->GetFilterDigitList()->size(); ++v){
    if((*(fInputEvent->GetFilterDigitList()))[v]->GetRawQPEs() >= fChargeCut){
      digitVector.push_back((*(fInputEvent->GetFilterDigitList()))[v]);
    }
  }
  std::sort(digitVector.begin(),digitVector.end(),SortDigitsByCharge);
//  std::cout << " Slicing " << digitVector.size() << " digits from " << digitVector[0]->GetRawQPEs() 
//            << " to " << digitVector[digitVector.size()-1]->GetRawQPEs() << " PE." << std::endl;

  fIsDigitClustered.clear();
  for(unsigned int d = 0; d < digitVector.size(); ++d){
    fIsDigitClustered.push_back(0);
  }

  std::vector<std::vector<WCSimRecoDigit*> > slicedDigits;
  // Don't do any slicing if we don't have enough digits
  if(digitVector.size() < 3){
    std::cerr << "Not performing slicing, too few digits" << std::endl;
  }
  else{
    unsigned int clusterNum = 0;
    std::vector<WCSimRecoDigit*> tempVec;

    while(!AreAllDigitsClustered()){
      slicedDigits.push_back(tempVec);
      unsigned int nextDigitSeed = GetFirstUnclusteredDigit();
      // Seed the cluster with the next unclustered digit.
      slicedDigits[clusterNum].push_back(digitVector[nextDigitSeed]);
      fIsDigitClustered[nextDigitSeed] = true;
  
      // Loop over the filtered digits in the event
      bool addedHit = true;
      while(addedHit){
        int nAdded = 0;
        for(unsigned int d = 0; d < digitVector.size(); ++d){
          if(fIsDigitClustered[d]) continue;

          // Check if we should add the hit.
          fIsDigitClustered[d] = AddDigitIfClose(digitVector[d],slicedDigits[clusterNum]);  
          if(fIsDigitClustered[d]) ++nAdded;
        } // End the for loop.

        // If we added any hits then we should iterate over the list again.
        if(nAdded > 0){
          addedHit = true;
        }
        else{
          addedHit = false;
        }
      } // End clustering while loop

      ++clusterNum;
    } // End the while loop
  }

  // Now we should have all of our clusters.
  // Want to copy those that are big enough into the main slice vector.
  for(unsigned int v = 0; v < slicedDigits.size(); ++v){
    if(slicedDigits[v].size() >= fMinSliceSize){
      fSlicedDigits.push_back(slicedDigits[v]);
    }
  }

  // Now we need to turn our slices back into WCSimRecoEvents
  BuildEvents();

}

bool WCSimRecoSlicer::AddDigitIfClose(WCSimRecoDigit* digit, std::vector<WCSimRecoDigit*>& cluster){

  for(unsigned int d = 0; d < cluster.size(); ++d){
    if(GetDistanceBetweenDigits(digit,cluster[d]) < fDistanceCut){
      cluster.push_back(digit);
      return true;
    }
  }
  return false;
}

double WCSimRecoSlicer::GetDistanceBetweenDigits(WCSimRecoDigit* digit1, WCSimRecoDigit* digit2){
  double dist = (digit1->GetX()-digit2->GetX())*(digit1->GetX()-digit2->GetX());
  dist += (digit1->GetY()-digit2->GetY())*(digit1->GetY()-digit2->GetY());
  dist += (digit1->GetZ()-digit2->GetZ())*(digit1->GetZ()-digit2->GetZ());
  return sqrt(dist);
}

unsigned int WCSimRecoSlicer::GetFirstUnclusteredDigit(){
  for(unsigned int i = 0; i < fIsDigitClustered.size(); ++i){
    if(fIsDigitClustered[i] == false){
      return i;
    }
  }
  return 999999;
}

bool WCSimRecoSlicer::AreAllDigitsClustered(){

  if(GetFirstUnclusteredDigit() == 999999) return true;
  else return false;

}

void WCSimRecoSlicer::BuildEvents(){

//  std::vector<std::pair<unsigned int, WCSimRecoEvent*> > sortVec;
  for(unsigned int v=0; v < fSlicedDigits.size(); ++v){
    WCSimRecoEvent* newEvt = new WCSimRecoEvent();
    // Initialise the event from the input event
    newEvt->SetHeader(fInputEvent->GetRun(), fInputEvent->GetEvent(), fInputEvent->GetTrigger());
    newEvt->SetVertex(fInputEvent->GetVtxX(), fInputEvent->GetVtxY(), fInputEvent->GetVtxZ(), fInputEvent->GetVtxTime());
    newEvt->SetDirection(fInputEvent->GetVtxX(), fInputEvent->GetVtxY(), fInputEvent->GetVtxZ());
    newEvt->SetConeAngle(fInputEvent->GetConeAngle());
    newEvt->SetTrackLength(fInputEvent->GetTrackLength());
    newEvt->SetVtxFOM(fInputEvent->GetVtxFOM(), fInputEvent->GetVtxIterations(), fInputEvent->GetVtxPass());
    newEvt->SetVtxStatus(fInputEvent->GetVtxStatus());

    if(fInputEvent->IsFilterDone()){ newEvt->SetFilterDone();}
    if(fInputEvent->IsVertexFinderDone()){ newEvt->SetVertexFinderDone();}
    if(fInputEvent->IsRingFinderDone()){ newEvt->SetRingFinderDone();}

    // Add the digits
    for(unsigned int d = 0; d < fSlicedDigits[v].size(); ++d){
      newEvt->AddDigit(fSlicedDigits[v][d]);
      newEvt->AddFilterDigit(fSlicedDigits[v][d]);
    }
    fSlicedEvents.push_back(newEvt);
//    sortVec.push_back(std::make_pair<unsigned int, WCSimRecoEvent*>(fSlicedDigits[v].size(),newEvt));
  }

  // Now we want to sort the vector by the number of digits.
  std::sort(fSlicedEvents.begin(),fSlicedEvents.end(),SortTheSlices);

}

// A quick sort function for the reco events.
bool SortTheSlices(WCSimRecoEvent* a, WCSimRecoEvent *b){
  unsigned int aSize = a->GetFilterDigitList()->size();
  unsigned int bSize = b->GetFilterDigitList()->size();
  return aSize > bSize;
}

// Sort the digits time
bool SortDigitsByCharge(WCSimRecoDigit *a, WCSimRecoDigit *b){
  return a->GetRawQPEs() > b->GetRawQPEs();
}

