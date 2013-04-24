#include "WCSimEveDisplay.hh"

#include "WCSimGeometry.hh"

#include "WCSimRecoDigit.hh"

#include "TGeometry.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"

#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveViewer.h"
#include "TEveGeoNode.h"
#include "TEvePointSet.h"
#include "TEveStraightLineSet.h"
#include "TEveArrow.h"
#include "TEveText.h"

#include <iostream>
#include <vector>

ClassImp(WCSimEveDisplay)

WCSimEveDisplay::WCSimEveDisplay() : WCSimDisplay()
{
  this->Initialize();
}

WCSimEveDisplay::~WCSimEveDisplay()
{
  
}

void WCSimEveDisplay::Initialize()
{
  std::cout << " *** WCSimEveDisplay::Initialize() *** " << std::endl;

  this->BuildGeometry();
}

void WCSimEveDisplay::BuildGeometry()
{
  std::cout << " *** WCSimEveDisplay::BuildGeometry() *** " << std::endl;

  // Event Display Manager
  // =====================
  TEveManager::Create();

  // TGeo Geometry
  // =============
  TGeoManager *geom = new TGeoManager("DetectorGeometry", "detector geometry");

  // materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum",0.0,0.0,0.0);
  TGeoMaterial *matWater = new TGeoMaterial("Water",18.0,8.0,1.0);

  // media
  TGeoMedium* Vacuum = new TGeoMedium("Vacuum",1,matVacuum);
  TGeoMedium* Water = new TGeoMedium("Water",2,matWater);

  // top volume
  TGeoVolume *top = geom->MakeBox("Detector", Vacuum, 10000.0, 10000.0, 10000.0);
  geom->SetTopVolume(top);

  // cylindrical detector
  if( WCSimGeometry::Instance()->GetGeoType()==WCSimGeometry::kCylinder ){
    Double_t fCylRadius = WCSimGeometry::Instance()->GetCylRadius();
    Double_t fCylLength = WCSimGeometry::Instance()->GetCylLength();

    TGeoVolume* myCylinder = geom->MakeTube("Cylinder",Water,0.0,fCylRadius,0.5*fCylLength);
    myCylinder->SetLineColor(kCyan);
    myCylinder->SetTransparency(70);  // percentage transparency [0-100]
    myCylinder->SetVisibility(1);     
    top->AddNode( myCylinder, 0, new TGeoTranslation(0, 0, 0));
  }

  // mailbox detector
  if( WCSimGeometry::Instance()->GetGeoType()==WCSimGeometry::kMailBox ){
    Double_t fMailBoxX = WCSimGeometry::Instance()->GetMailBoxX();
    Double_t fMailBoxY = WCSimGeometry::Instance()->GetMailBoxY();
    Double_t fMailBoxZ = WCSimGeometry::Instance()->GetMailBoxZ();

    TGeoVolume* myMailBox = geom->MakeBox("MailBox",Water,0.5*fMailBoxX,0.5*fMailBoxY,0.5*fMailBoxZ);
    myMailBox->SetLineColor(kCyan);
    myMailBox->SetTransparency(70);  // percentage transparency [0-100]
    myMailBox->SetVisibility(1);     
    top->AddNode( myMailBox, 0, new TGeoTranslation(0, 0, 0));
  }

  // close geometry
  geom->CloseGeometry();

  // Create Geometry in Event Display
  // ================================
  TGeoNode* node = gGeoManager->GetTopNode();
  TEveGeoTopNode* eveNode = new TEveGeoTopNode(gGeoManager, node);
  eveNode->SetVisLevel(1);
  gEve->AddGlobalElement(eveNode);

  // Draw Display
  // ============
  gEve->Redraw3D(kTRUE);

}

void WCSimEveDisplay::DrawDisplay(WCSimRecoEvent* myRecoEvent)
{  
  // Reset Event Display
  // ===================
  TEveEventManager* currEvent = gEve->GetCurrentEvent();
  if( currEvent ) currEvent->DestroyElements();

  // Check for Event
  // ===============
  if( myRecoEvent==0 ) return;

  // Containers for Hits
  // ===================
  Int_t markerStyle = 4;   
  Double_t markerSize = 0.75;

  /*
  Int_t colourCode1 = kYellow-7;
  Int_t colourCode2 = kCyan-7;
  Int_t colourCode3 = kCyan+1;
  Int_t colourCode4 = kBlue-4;
  Int_t colourCode5 = kBlue+1;
  Int_t colourCode6 = kMagenta+1;
  Int_t colourCode7 = kMagenta+2;
  Int_t colourCode8 = kRed-4;
  Int_t colourCode9 = kRed; 
  */

  Int_t colourCode1 = kBlue+1;
  Int_t colourCode2 = kCyan+1;
  Int_t colourCode3 = kGreen;
  Int_t colourCode4 = kYellow;
  Int_t colourCode5 = kOrange;
  Int_t colourCode6 = kOrange+7;
  Int_t colourCode7 = kRed;
  Int_t colourCode8 = kRed;
  Int_t colourCode9 = kRed;

  TEvePointSet* fHitList1 = new TEvePointSet(); 
  fHitList1->SetMarkerColor(colourCode1);
  fHitList1->SetMarkerSize(markerSize);
  fHitList1->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList1);

  TEvePointSet* fHitList2 = new TEvePointSet(); 
  fHitList2->SetMarkerColor(colourCode2);
  fHitList2->SetMarkerSize(markerSize);
  fHitList2->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList2);

  TEvePointSet* fHitList3 = new TEvePointSet();  
  fHitList3->SetMarkerColor(colourCode3);
  fHitList3->SetMarkerSize(markerSize);
  fHitList3->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList3);

  TEvePointSet* fHitList4 = new TEvePointSet();  
  fHitList4->SetMarkerColor(colourCode4);
  fHitList4->SetMarkerSize(markerSize);
  fHitList4->SetMarkerStyle(markerStyle);  
  gEve->AddElement(fHitList4);

  TEvePointSet* fHitList5 = new TEvePointSet();  
  fHitList5->SetMarkerColor(colourCode5);
  fHitList5->SetMarkerSize(markerSize);
  fHitList5->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList5);

  TEvePointSet* fHitList6 = new TEvePointSet();  
  fHitList6->SetMarkerColor(colourCode6);
  fHitList6->SetMarkerSize(markerSize);
  fHitList6->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList6);

  TEvePointSet* fHitList7 = new TEvePointSet();  
  fHitList7->SetMarkerColor(colourCode7);
  fHitList7->SetMarkerSize(markerSize);
  fHitList7->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList7);

  TEvePointSet* fHitList8 = new TEvePointSet(); 
  fHitList8->SetMarkerColor(colourCode8);
  fHitList8->SetMarkerSize(markerSize);
  fHitList8->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList8);

  TEvePointSet* fHitList9 = new TEvePointSet();
  fHitList9->SetMarkerColor(colourCode9);
  fHitList9->SetMarkerSize(markerSize);
  fHitList9->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList9);


  // Loop over digits
  // ================
  for( Int_t nDigit=0; nDigit<myRecoEvent->GetNDigits(); nDigit++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(myRecoEvent->GetDigit(nDigit));

    Double_t Q = myDigit->GetQPEs();
    // Double_t T = myDigit->GetTime();

    if( Q<GetPulseHeightCut() ) continue;

    Double_t x = myDigit->GetX();
    Double_t y = myDigit->GetY();
    Double_t z = myDigit->GetZ();

    Int_t listNumber = 0;

    if( Q<0.8 )             listNumber = 1;
    if( Q>=0.8 && Q<1.5 )   listNumber = 2;
    if( Q>=1.5 && Q<2.5 )   listNumber = 3;
    if( Q>=2.5 && Q<5.0 )   listNumber = 4;
    if( Q>=5.0 && Q<10.0 )  listNumber = 5;
    if( Q>=10.0 && Q<15.0 ) listNumber = 6;
    if( Q>=15.0 && Q<20.0 ) listNumber = 7;
    if( Q>=20.0 && Q<30.0 ) listNumber = 8;
    if( Q>=30.0 )           listNumber = 9;

    switch( listNumber ){
      case 1: fHitList1->SetNextPoint(x,y,z); break;
      case 2: fHitList2->SetNextPoint(x,y,z); break;
      case 3: fHitList3->SetNextPoint(x,y,z); break;
      case 4: fHitList4->SetNextPoint(x,y,z); break;
      case 5: fHitList5->SetNextPoint(x,y,z); break;
      case 6: fHitList6->SetNextPoint(x,y,z); break;
      case 7: fHitList7->SetNextPoint(x,y,z); break;
      case 8: fHitList8->SetNextPoint(x,y,z); break;
      case 9: fHitList9->SetNextPoint(x,y,z); break;
      default: break;
    }
  }

  // Re-draw Event Display
  // =====================
  gEve->Redraw3D();

  return;
}

void WCSimEveDisplay::DrawCleanDisplay(WCSimRecoEvent* myRecoEvent)
{  
  // Reset Event Display
  // ===================
  TEveEventManager* currEvent = gEve->GetCurrentEvent();
  if( currEvent ) currEvent->DestroyElements();

  // Check for Event
  // ===============
  if( myRecoEvent==0 ) return;

  // Containers for Hits
  // ===================
  Int_t markerStyle = 4;   
  Double_t markerSize = 0.75;

  /*
  Int_t colourCode1 = kYellow-7;
  Int_t colourCode2 = kCyan-7;
  Int_t colourCode3 = kCyan+1;
  Int_t colourCode4 = kBlue-4;
  Int_t colourCode5 = kBlue+1;
  Int_t colourCode6 = kMagenta+1;
  Int_t colourCode7 = kMagenta+2;
  Int_t colourCode8 = kRed-4;
  Int_t colourCode9 = kRed; 
  */

  Int_t colourCode1 = kBlue+1;
  Int_t colourCode2 = kCyan+1;
  Int_t colourCode3 = kGreen;
  Int_t colourCode4 = kYellow;
  Int_t colourCode5 = kOrange;
  Int_t colourCode6 = kOrange+7;
  Int_t colourCode7 = kRed;
  Int_t colourCode8 = kRed;
  Int_t colourCode9 = kRed;

  TEvePointSet* fHitList1 = new TEvePointSet(); 
  fHitList1->SetMarkerColor(colourCode1);
  fHitList1->SetMarkerSize(markerSize);
  fHitList1->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList1);

  TEvePointSet* fHitList2 = new TEvePointSet(); 
  fHitList2->SetMarkerColor(colourCode2);
  fHitList2->SetMarkerSize(markerSize);
  fHitList2->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList2);

  TEvePointSet* fHitList3 = new TEvePointSet();  
  fHitList3->SetMarkerColor(colourCode3);
  fHitList3->SetMarkerSize(markerSize);
  fHitList3->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList3);

  TEvePointSet* fHitList4 = new TEvePointSet();  
  fHitList4->SetMarkerColor(colourCode4);
  fHitList4->SetMarkerSize(markerSize);
  fHitList4->SetMarkerStyle(markerStyle);  
  gEve->AddElement(fHitList4);

  TEvePointSet* fHitList5 = new TEvePointSet();  
  fHitList5->SetMarkerColor(colourCode5);
  fHitList5->SetMarkerSize(markerSize);
  fHitList5->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList5);

  TEvePointSet* fHitList6 = new TEvePointSet();  
  fHitList6->SetMarkerColor(colourCode6);
  fHitList6->SetMarkerSize(markerSize);
  fHitList6->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList6);

  TEvePointSet* fHitList7 = new TEvePointSet();  
  fHitList7->SetMarkerColor(colourCode7);
  fHitList7->SetMarkerSize(markerSize);
  fHitList7->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList7);

  TEvePointSet* fHitList8 = new TEvePointSet(); 
  fHitList8->SetMarkerColor(colourCode8);
  fHitList8->SetMarkerSize(markerSize);
  fHitList8->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList8);

  TEvePointSet* fHitList9 = new TEvePointSet();
  fHitList9->SetMarkerColor(colourCode9);
  fHitList9->SetMarkerSize(markerSize);
  fHitList9->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList9);


  // Loop over digits
  // ================
  for( Int_t nDigit=0; nDigit<myRecoEvent->GetNFilterDigits(); nDigit++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(myRecoEvent->GetFilterDigit(nDigit));

    Double_t Q = myDigit->GetQPEs();
    // Double_t T = myDigit->GetTime();

    if( Q<GetPulseHeightCut() ) continue;

    Double_t x = myDigit->GetX();
    Double_t y = myDigit->GetY();
    Double_t z = myDigit->GetZ();

    Int_t listNumber = 0;

    if( Q<0.8 )             listNumber = 1;
    if( Q>=0.8 && Q<1.5 )   listNumber = 2;
    if( Q>=1.5 && Q<2.5 )   listNumber = 3;
    if( Q>=2.5 && Q<5.0 )   listNumber = 4;
    if( Q>=5.0 && Q<10.0 )  listNumber = 5;
    if( Q>=10.0 && Q<15.0 ) listNumber = 6;
    if( Q>=15.0 && Q<20.0 ) listNumber = 7;
    if( Q>=20.0 && Q<30.0 ) listNumber = 8;
    if( Q>=30.0 )           listNumber = 9;

    switch( listNumber ){
      case 1: fHitList1->SetNextPoint(x,y,z); break;
      case 2: fHitList2->SetNextPoint(x,y,z); break;
      case 3: fHitList3->SetNextPoint(x,y,z); break;
      case 4: fHitList4->SetNextPoint(x,y,z); break;
      case 5: fHitList5->SetNextPoint(x,y,z); break;
      case 6: fHitList6->SetNextPoint(x,y,z); break;
      case 7: fHitList7->SetNextPoint(x,y,z); break;
      case 8: fHitList8->SetNextPoint(x,y,z); break;
      case 9: fHitList9->SetNextPoint(x,y,z); break;
      default: break;
    }
  }

  // Re-draw Event Display
  // =====================
  gEve->Redraw3D();

  return;
}

void WCSimEveDisplay::DrawTrueEvent(WCSimTrueEvent*) 
{ 
  std::cout << " --- WCSimEveDisplay::DrawTrueEvent(...) [NOT IMPLEMENTED YET] --- " << std::endl;  

  //
  // note: draw true tracks as a TEveStraightLineSet here
  //

  return; 
}
    
void WCSimEveDisplay::DrawRecoEvent(WCSimRecoEvent*)
{ 
  std::cout << " --- WCSimEveDisplay::DrawRecoEvent(...) [NOT IMPLEMENTED YET] --- " << std::endl;  

  //
  // note: draw reco tracks as a TEveStraightLineSet here
  //

  return; 
}
  
void WCSimEveDisplay::ResetDisplay()                   
{ 
  std::cout << " --- WCSimEveDisplay::ResetDisplay() [NOT IMPLEMENTED YET] --- " << std::endl;
  return; 
}
  
void WCSimEveDisplay::PrintDisplay()
{ 
  std::cout << " --- WCSimEveDisplay::PrintDisplay() [NOT IMPLEMENTED YET] --- " << std::endl;
  return; 
}
