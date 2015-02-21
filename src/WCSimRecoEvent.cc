#include "WCSimRecoEvent.hh"
#include "WCSimRecoRing.hh"
#include "WCSimRecoVertex.hh"
#include "WCSimRecoDigit.hh"

#include "WCSimRecoObjectTable.hh"
#include "WCSimRootEvent.hh"

#include <iostream>
#include <fstream>

ClassImp(WCSimRecoEvent)

WCSimRecoEvent::WCSimRecoEvent()
{
  fVertex = new WCSimRecoVertex();

  fDigitList = new std::vector<WCSimRecoDigit*>;
  fFilterDigitList =  new std::vector<WCSimRecoDigit*>;
  fRingList = new std::vector<WCSimRecoRing*>;

  fRunNum = -1;
  fEventNum = -1;
  fTriggerNum = -1;

  fIsFilterDone = 0;
  fIsVertexFinderDone = 0;
  fIsRingFinderDone = 0;

  WCSimRecoObjectTable::Instance()->NewEvent();
}

WCSimRecoEvent::~WCSimRecoEvent()
{
  this->Reset();

  delete fVertex;

  delete fDigitList;
  delete fFilterDigitList;
  delete fRingList;

  WCSimRecoObjectTable::Instance()->DeleteEvent();
}

void WCSimRecoEvent::Reset()
{  
  fVertex->Reset();

  this->ClearDigits();
  this->ClearFilterDigits();
  this->ClearRings();

  fRunNum = -1;
  fEventNum = -1;
  fTriggerNum = -1;

  fIsFilterDone = 0;
  fIsVertexFinderDone = 0;
  fIsRingFinderDone = 0;
}

void WCSimRecoEvent::SetHeader(Int_t run, Int_t event, Int_t trigger)
{
  fRunNum = run;
  fEventNum = event;
  fTriggerNum = trigger;
}

void WCSimRecoEvent::AddDigit(WCSimRecoDigit* digit)
{
  fDigitList->push_back(digit);
}

void WCSimRecoEvent::AddFilterDigit(WCSimRecoDigit* digit)
{
  fFilterDigitList->push_back(digit);
}

void WCSimRecoEvent::AddRing(WCSimRecoRing* ring)
{
  fRingList->push_back(ring);
}

void WCSimRecoEvent::ClearDigits()
{
  fDigitList->clear();
}

void WCSimRecoEvent::ClearFilterDigits()
{
  fFilterDigitList->clear();
}

void WCSimRecoEvent::ClearRings()
{
  fRingList->clear();
}

WCSimRecoDigit* WCSimRecoEvent::GetDigit(Int_t n)
{
  return (WCSimRecoDigit*)(fDigitList->at(n));
}
  
Int_t WCSimRecoEvent::GetNDigits()
{
  return fDigitList->size();
}

WCSimRecoDigit* WCSimRecoEvent::GetFilterDigit(Int_t n)
{
  return (WCSimRecoDigit*)(fFilterDigitList->at(n));
}
  
Int_t WCSimRecoEvent::GetNFilterDigits()
{
  return fFilterDigitList->size();
}

WCSimRecoRing* WCSimRecoEvent::GetRing(Int_t n)
{
  return (WCSimRecoRing*)(fRingList->at(n));
}
  
Int_t WCSimRecoEvent::GetNRings()
{
  return fRingList->size();
}

WCSimRecoRing* WCSimRecoEvent::GetPrimaryRing()
{
  if( fRingList->size()>0 ){
    return (WCSimRecoRing*)(fRingList->at(0));
  }
  else{
    return 0;
  }
}

void WCSimRecoEvent::SetVertex( Double_t x, Double_t y, Double_t z, Double_t t )
{
  fVertex->SetVertex(x,y,z,t);
}

void WCSimRecoEvent::SetDirection( Double_t px, Double_t py, Double_t pz )
{
  fVertex->SetDirection(px,py,pz);  
}  

void WCSimRecoEvent::SetConeAngle( Double_t angle )
{
  fVertex->SetConeAngle(angle);
}

void WCSimRecoEvent::SetTrackLength( Double_t length )
{
  fVertex->SetTrackLength(length);
}

void WCSimRecoEvent::SetVtxFOM( Double_t fom, Int_t nsteps, Bool_t pass)
{
  fVertex->SetFOM(fom,nsteps,pass);
}

void WCSimRecoEvent::SetVtxStatus( Int_t status )
{
  fVertex->SetStatus(status);
}

WCSimRecoVertex* WCSimRecoEvent::GetVertex()
{ 
  return fVertex; 
}

Double_t WCSimRecoEvent::GetVtxX() 
{ 
  return fVertex->GetX(); 
}
  
Double_t WCSimRecoEvent::GetVtxY() 
{ 
  return fVertex->GetY(); 
}
  
Double_t WCSimRecoEvent::GetVtxZ() 
{ 
  return fVertex->GetZ(); 
}
  
Double_t WCSimRecoEvent::GetVtxTime()
{ 
  return fVertex->GetTime(); 
}
  
Double_t WCSimRecoEvent::GetDirX()
{
  return fVertex->GetDirX();
}
  
Double_t WCSimRecoEvent::GetDirY()
{
  return fVertex->GetDirY();
}
  
Double_t WCSimRecoEvent::GetDirZ()
{
  return fVertex->GetDirZ();
}
  
Double_t WCSimRecoEvent::GetConeAngle()
{
  return fVertex->GetConeAngle();
}
 
Double_t WCSimRecoEvent::GetTrackLength()
{
  return fVertex->GetTrackLength();
}
 
Double_t WCSimRecoEvent::GetVtxFOM()
{ 
  return fVertex->GetFOM(); 
}

Int_t WCSimRecoEvent::GetVtxIterations()
{
  return fVertex->GetIterations();
}

Bool_t WCSimRecoEvent::GetVtxPass()
{
  return fVertex->GetPass();
}

Int_t WCSimRecoEvent::GetVtxStatus()
{
  return fVertex->GetStatus();
}
  
Bool_t WCSimRecoEvent::FoundVertex()
{
  return fVertex->FoundVertex();
}
  
Bool_t WCSimRecoEvent::FoundDirection()
{
  return fVertex->FoundDirection();
}

Bool_t WCSimRecoEvent::FoundRings()
{
  if( fRingList->size()>0 ) return 1;
  else return 0;
}

void WCSimRecoEvent::PrintDigitList(const char* filename)
{
  std::ofstream output(filename);

  for( Int_t i=0; i<this->GetNDigits(); i++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(this->GetDigit(i));
    output << myDigit->GetX() << " " << myDigit->GetY() << " " << myDigit->GetZ() << " " << myDigit->GetTime() << " " << myDigit->GetQPEs() << std::endl; 
  }

  output.close();
}

void WCSimRecoEvent::PrintFilterDigitList(const char* filename)
{
  std::ofstream output(filename);

  for( Int_t i=0; i<this->GetNFilterDigits(); i++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(this->GetFilterDigit(i));
    output << myDigit->GetX() << " " << myDigit->GetY() << " " << myDigit->GetZ() << " " << myDigit->GetTime() << " " << myDigit->GetQPEs() << std::endl;
  }

  output.close();
}

void WCSimRecoEvent::PrintEvent()
{
  std::cout << " *** WCSimRecoEvent::PrintEvent() *** " << std::endl
            << " * VtxX = " << fVertex->GetX() << std::endl
	    << " * VtxY = " << fVertex->GetY() << std::endl
            << " * VtxZ = " << fVertex->GetZ() << std::endl
            << " * VtxTime = = " << fVertex->GetTime() << std::endl
            << " * DirX = " << fVertex->GetDirX() << std::endl
            << " * DirY = " << fVertex->GetDirY() << std::endl
            << " * DirZ = " << fVertex->GetDirZ() << std::endl
            << " * VtxFoM = " << fVertex->GetFOM() << std::endl
            << " ************************************ " << std::endl;

  return;
}
