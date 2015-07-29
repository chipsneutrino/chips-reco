#include <iostream>
#include <vector>
#include <algorithm> 
#include <cassert>
#include <cmath>

#include "WCSimRecoSlicer.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoDigit.hh"

// A little sort function for events
bool SortTheSlices(WCSimRecoEvent* a, WCSimRecoEvent *b);

WCSimRecoSlicer::WCSimRecoSlicer(){
  fInputEvent = 0x0;

  fDistanceCut = 500.0;
  fMinSliceSize = 10;
}

WCSimRecoSlicer::~WCSimRecoSlicer(){
  // Delete the created events?
}

void WCSimRecoSlicer::SetInputEvent(WCSimRecoEvent* evt){
  fInputEvent = evt;

  fIsDigitClustered.clear();
  for(unsigned int d = 0; d < fInputEvent->GetFilterDigitList()->size(); ++d){
    fIsDigitClustered.push_back(0);
  }
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

  std::vector<std::vector<WCSimRecoDigit*> > slicedDigits;
  // Don't do any slicing if we don't have enough digits
  if(fInputEvent->GetFilterDigitList()->size() < 3){
    std::cerr << "Not performing slicing, too few digits" << std::endl;
  }
  else{
    unsigned int clusterNum = 0;
    std::vector<WCSimRecoDigit*> tempVec;

    while(!AreAllDigitsClustered()){
      slicedDigits.push_back(tempVec);
      unsigned int nextDigitSeed = GetFirstUnclusteredDigit();
      // Seed the cluster with the next unclustered digit.
      slicedDigits[clusterNum].push_back((*(fInputEvent->GetFilterDigitList()))[nextDigitSeed]);
      fIsDigitClustered[nextDigitSeed] = true;
  
      // Loop over the filtered digits in the event
      for(unsigned int d = nextDigitSeed+1; d < fInputEvent->GetFilterDigitList()->size(); ++d){
        if(fIsDigitClustered[d]) continue;
        fIsDigitClustered[d] = AddDigitIfClose((*(fInputEvent->GetFilterDigitList()))[d],slicedDigits[clusterNum]);  
  
      } // End the for loop

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
    // Add the digits
    for(unsigned int d = 0; d < fSlicedDigits[v].size(); ++d){
      newEvt->AddDigit(fSlicedDigits[v][d]);
      newEvt->AddFilterDigit(fSlicedDigits[v][d]);
    }
    fSlicedEvents.push_back(newEvt);
//    sortVec.push_back(std::make_pair<unsigned int, WCSimRecoEvent*>(fSlicedDigits[v].size(),newEvt));
  }

  // Now we want to sort the vector by the number of digits.
//  std::sort(sortVec.begin(),sortVec.end());
  std::sort(fSlicedEvents.begin(),fSlicedEvents.end(),SortTheSlices);
//  for(unsigned int i = 0; i < sortVec.size(); ++i){ 
  for(unsigned int i = 0; i < fSlicedEvents.size(); ++i){
    std::cout << "Built sliced event with " << fSlicedEvents[i]->GetFilterDigitList()->size() << " digits." << std::endl;
  }

}

// A quick sort function for the reco events.
bool SortTheSlices(WCSimRecoEvent* a, WCSimRecoEvent *b){
  unsigned int aSize = a->GetFilterDigitList()->size();
  unsigned int bSize = b->GetFilterDigitList()->size();
  return aSize > bSize;
}


