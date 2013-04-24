#include "WCSimRecoAB.hh"

#include "WCSimInterface.hh"

#include "WCSimRecoDigit.hh"
#include "WCSimRecoRing.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoVertex.hh"

#include "WCSimGeometry.hh"
#include "WCSimHoughTransform.hh"

#include "WCSimRingFinder.hh"
#include "WCSimVertexFinder.hh"
#include "WCSimDataCleaner.hh"

#include "TMath.h"

#include <iostream>
#include <cmath>

ClassImp(WCSimRecoAB)

WCSimRecoAB::WCSimRecoAB()
{

}

WCSimRecoAB::~WCSimRecoAB()
{

}

void WCSimRecoAB::Run()
{
  WCSimRecoEvent* myEvent = WCSimInterface::RecoEvent();

  this->Run(myEvent);

  return;
}
  
void WCSimRecoAB::Run(WCSimRecoEvent* myEvent)
{
  // filter digits
  // =============
  this->RunFilter(myEvent);

  // reconstruct vertex
  // ==================
  this->RunRecoVertex(myEvent);

  // reconstruct rings
  // =================
  this->RunRecoRings(myEvent);

  return;
}

void WCSimRecoAB::RunFilter(WCSimRecoEvent* myEvent)
{  
  // Get Digit List
  // ==============
  std::vector<WCSimRecoDigit*>* myDigitList = (std::vector<WCSimRecoDigit*>*)(myEvent->GetDigitList()); 

  // Get Data Cleaner
  // ================
  WCSimDataCleaner* myDataCleaner = WCSimDataCleaner::Instance();  
  
  // Run Data Cleaner
  // ================
  std::vector<WCSimRecoDigit*>* myFilterDigitList = (std::vector<WCSimRecoDigit*>*)(myDataCleaner->Run(myDigitList));

  // Add Filtered Digits
  // ===================
  for( UInt_t n=0; n<myFilterDigitList->size(); n++ ){
    WCSimRecoDigit* recoDigit = (WCSimRecoDigit*)(myFilterDigitList->at(n));
    myEvent->AddFilterDigit(recoDigit);
  }

  // Done!
  // =====
  myEvent->SetFilterDone();

  return;
}

void WCSimRecoAB::RunRecoVertex(WCSimRecoEvent* myEvent)
{
  // Run Filter (if necessary)
  // =========================
  if( myEvent->IsFilterDone()==0 ){
    this->RunFilter(myEvent);
  }

  // Check Filter Digits (bail out if necessary)
  // ===========================================
  if( myEvent->GetNFilterDigits()==0 ){
    myEvent->SetVertexFinderDone();
    return;
  }

  // Get Vertex Finder
  // =================
  WCSimVertexFinder* myVertexFinder = WCSimVertexFinder::Instance();

  // Run Vertex Finder
  // =================
  WCSimRecoVertex* myVertex = (WCSimRecoVertex*)(myVertexFinder->Run(myEvent));
    
  // Set Vertex
  // ==========
  if( myVertex->FoundVertex() ){
    myEvent->SetVertex(myVertex->GetX(),
                       myVertex->GetY(),
                       myVertex->GetZ(),
                       myVertex->GetTime());
  }

  // Set Direction
  // =============
  if( myVertex->FoundDirection() ){
    myEvent->SetDirection(myVertex->GetDirX(),
                          myVertex->GetDirY(),
                          myVertex->GetDirZ());
  }

  // Set FoM
  // =======
  myEvent->SetConeAngle(myVertex->GetConeAngle());
  myEvent->SetTrackLength(myVertex->GetTrackLength());
  myEvent->SetVtxFOM(myVertex->GetFOM(),
                     myVertex->GetIterations(),
                     myVertex->GetPass() );
  myEvent->SetVtxStatus(myVertex->GetStatus());

  // Done!
  // =====
  myEvent->SetVertexFinderDone();

  return;
}

void WCSimRecoAB::RunRecoRings(WCSimRecoEvent* myEvent)
{
  // Find Vertex (if necessary)
  // ==========================
  if( myEvent->IsVertexFinderDone()==0 ){
    this->RunRecoVertex(myEvent);
  }

  // Check Vertex (bail out if necessary)
  // ====================================
  if( myEvent->FoundVertex()==0 ){
    myEvent->SetRingFinderDone();
    return;
  }

  // Get Vertex
  // ==========
  WCSimRecoVertex* myVertex = (WCSimRecoVertex*)(myEvent->GetVertex());

  // Get Ring Finder
  // ===============
  WCSimRingFinder* myRingFinder = WCSimRingFinder::Instance();

  // Run Ring Finder
  // ===============
  std::vector<WCSimRecoRing*>* myRingList = (std::vector<WCSimRecoRing*>*)(myRingFinder->Run(myEvent,myVertex));

  // Add Rings
  // =========
  for( UInt_t n=0; n<myRingList->size(); n++ ){
	 std::cout << "Adding ring" << std::endl;
    WCSimRecoRing* recoRing = (WCSimRecoRing*)(myRingList->at(n));
    myEvent->AddRing(recoRing);
  }

  // Done!
  // =====
  myEvent->SetRingFinderDone();

  return;
}

