#include "WCSimEveDisplay.hh"

#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimTruthSummary.hh"

#include "TDatabasePDG.h"
#include "TGeometry.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGLViewer.h"
#include "TParticlePDG.h"
#include "TVector3.h"

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
  fEventNum = -999;
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
    myCylinder->SetLineColor(kWhite);
    myCylinder->SetTransparency(50);  // percentage transparency [0-100]
    myCylinder->SetVisibility(1);     
    top->AddNode( myCylinder, 0, new TGeoTranslation(0, 0, 0));
  }

  // mailbox detector
  if( WCSimGeometry::Instance()->GetGeoType()==WCSimGeometry::kMailBox ){
    Double_t fMailBoxX = WCSimGeometry::Instance()->GetMailBoxX();
    Double_t fMailBoxY = WCSimGeometry::Instance()->GetMailBoxY();
    Double_t fMailBoxZ = WCSimGeometry::Instance()->GetMailBoxZ();

    TGeoVolume* myMailBox = geom->MakeBox("MailBox",Water,0.5*fMailBoxX,0.5*fMailBoxY,0.5*fMailBoxZ);
    myMailBox->SetLineColor(kWhite);
    myMailBox->SetTransparency(80);  // percentage transparency [0-100]
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
  fEventNum = myRecoEvent->GetEvent();
  
  // Re-make the Eve point sets
  // =========================
  this->BuildEvePointSets();


  // Loop over digits
  // ================
  for( Int_t nDigit=0; nDigit<myRecoEvent->GetNDigits(); nDigit++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(myRecoEvent->GetDigit(nDigit));

    Double_t Q = myDigit->GetQPEs();
    Double_t T = myDigit->GetTime() - myRecoEvent->GetVtxTime();

    if( Q<GetPulseHeightCut() ) continue;

    Double_t x = myDigit->GetX();
    Double_t y = myDigit->GetY();
    Double_t z = myDigit->GetZ();

    Int_t listNumberQ = 0, listNumberT = 0;
    if( Q<0.8 )              listNumberQ = 1;
    if( Q>=0.8 && Q<1.5 )    listNumberQ = 2;
    if( Q>=1.5 && Q<2.5 )    listNumberQ = 3;
    if( Q>=2.5 && Q<5.0 )    listNumberQ = 4;
    if( Q>=5.0 && Q<10.0 )   listNumberQ = 5;
    if( Q>=10.0 && Q<15.0 )  listNumberQ = 6;
    if( Q>=15.0 && Q<20.0 )  listNumberQ = 7;
    if( Q>=20.0 && Q<30.0 )  listNumberQ = 8;
    if( Q>=30.0 )            listNumberQ = 9;

    if( T < 5.0 )            listNumberT = 0;
    if( T >= 5.0 && T<10.0)  listNumberT = 1;
    if( T >= 10.0 && T<15.0) listNumberT = 2;
    if( T >= 15.0 && T<20.0) listNumberT = 3;
    if( T >= 20.0 && T<25.0) listNumberT = 4;
    if( T >= 25.0 && T<30.0) listNumberT = 5;
    if( T >= 30.0 && T<35.0) listNumberT = 6;
    if( T >= 35.0 && T<40.0) listNumberT = 7;
    if( T >= 40.0 && T<45.0) listNumberT = 8;
    if( T >= 50.0)           listNumberT = 9;

   
    fEvePointSets[listNumberQ][listNumberT]->SetNextPoint(x,y,z);

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
  fEventNum = myRecoEvent->GetEvent();

  // Containers for Hits
  // ===================

  // Re-make the Eve point sets
  // =========================
  this->BuildEvePointSets();



  // Loop over digits
  // ================
  for( Int_t nDigit=0; nDigit<myRecoEvent->GetNFilterDigits(); nDigit++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(myRecoEvent->GetFilterDigit(nDigit));

    Double_t Q = myDigit->GetQPEs();
    Double_t T = myDigit->GetTime() - myRecoEvent->GetVtxTime();

    if( Q<GetPulseHeightCut() ) continue;

    Double_t x = myDigit->GetX();
    Double_t y = myDigit->GetY();
    Double_t z = myDigit->GetZ();

    Int_t listNumberQ = 0, listNumberT = 0;

    if( Q<0.8 )              listNumberQ = 1;
    if( Q>=0.8 && Q<1.5 )    listNumberQ = 2;
    if( Q>=1.5 && Q<2.5 )    listNumberQ = 3;
    if( Q>=2.5 && Q<5.0 )    listNumberQ = 4;
    if( Q>=5.0 && Q<10.0 )   listNumberQ = 5;
    if( Q>=10.0 && Q<15.0 )  listNumberQ = 6;
    if( Q>=15.0 && Q<20.0 )  listNumberQ = 7;
    if( Q>=20.0 && Q<30.0 )  listNumberQ = 8;
    if( Q>=30.0 )            listNumberQ = 9;

    if( T < 5.0 )            listNumberT = 0;
    if( T >= 5.0 && T<10.0)  listNumberT = 1;
    if( T >= 10.0 && T<15.0) listNumberT = 2;
    if( T >= 15.0 && T<20.0) listNumberT = 3;
    if( T >= 20.0 && T<25.0) listNumberT = 4;
    if( T >= 25.0 && T<30.0) listNumberT = 5;
    if( T >= 30.0 && T<35.0) listNumberT = 6;
    if( T >= 35.0 && T<40.0) listNumberT = 7;
    if( T >= 40.0 && T<45.0) listNumberT = 8;
    if( T >= 50.0)           listNumberT = 9;

    fEvePointSets[listNumberQ][listNumberT]->SetNextPoint(x,y,z);
  }

  // Re-draw Event Display
  // =====================
  gEve->Redraw3D();

  return;
}

void WCSimEveDisplay::DrawTrueEvent(WCSimRootEvent* myEvent) 
{ 
  WCSimTruthSummary * mySummary = new WCSimTruthSummary(WCSimInterface::Instance()->GetTruthSummary());
  mySummary->Print();
  std::cout << " --- WCSimEveDisplay::DrawTrueEvent(...) --- " << std::endl;  

  BuildTrueTracks(mySummary);
  delete mySummary;

  gEve->Redraw3D();
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
  std::cout << " --- WCSimEveDisplay::PrintDisplay() --- " << std::endl;
  TString filename = Form("event_%d.png", fEventNum);
  gEve->GetDefaultGLViewer()->SavePicture(filename.Data());
  std::cout << " --- Saved to " << filename << " --- " << std::endl;
  return; 
}

void WCSimEveDisplay::BuildEvePointSets()
{
  // Containers for Hits
  // ===================
  Int_t markerStyle = 21;

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

  Int_t colourCode0 = kMagenta+3;
  Int_t colourCode1 = kBlue+1;
  Int_t colourCode2 = kCyan+1;
  Int_t colourCode3 = kGreen+1;
  Int_t colourCode4 = kGreen-9;
  Int_t colourCode5 = kYellow;
  Int_t colourCode6 = kOrange;
  Int_t colourCode7 = kOrange+7;
  Int_t colourCode8 = kRed-3;
  Int_t colourCode9 = kRed;

  Int_t colours[fNumBins] = { colourCode0, colourCode1, colourCode2, colourCode3, colourCode4,
                              colourCode5, colourCode6, colourCode7, colourCode8, colourCode9 };
  Float_t markerSize[fNumBins] = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

  ClearEvePointSets();
  for(int iCharge = 0; iCharge < fNumBins; ++iCharge)
  {
    std::vector<TEvePointSet*> blank;
    fEvePointSets.push_back(blank);
    for(int iTime = 0; iTime < fNumBins; ++iTime)
    {
      TEvePointSet * pointSet = new TEvePointSet();
      pointSet->SetMarkerColor(colours[iTime]);
      pointSet->SetMarkerSize(markerSize[iCharge]);
      pointSet->SetMarkerStyle(markerStyle);
      fEvePointSets[iCharge].push_back(pointSet);
      gEve->AddElement(fEvePointSets[iCharge][iTime]);
    }
  }
}

void WCSimEveDisplay::ClearEvePointSets()
{
  // Think Eve handles deletion of the object it knows about
  /*
  std::vector<std::vector<TEvePointSet*> >::iterator setQItr = fEvePointSets.begin();
  while(setQItr != fEvePointSets.end())
  {
    std::vector<TEvePointSet*>::iterator setTItr = setQItr->begin();
    while(setTItr != setQItr->end())
    {
      delete (*setTItr);
      (*setTItr) = NULL;
      ++setTItr;
    }
    setQItr++;
  }
  */
  fEvePointSets.clear();
  fEvePointSets.resize(0);
}

void WCSimEveDisplay::BuildTrueTracks(WCSimTruthSummary *mySummary)
{
  // Make a TEveStraightLineSet
  TEveStraightLineSet* trueCherenkovTracks = new TEveStraightLineSet("True Cherenkov tracks");
  trueCherenkovTracks->SetLineColor(7);
  trueCherenkovTracks->SetLineStyle(kSolid);
  trueCherenkovTracks->SetLineWidth(2);
  TEveStraightLineSet* leadingTrack = new TEveStraightLineSet("Leading track");
  leadingTrack->SetLineColor(kGreen);
  leadingTrack->SetLineStyle(kSolid);
  leadingTrack->SetLineWidth(2);
  TEveStraightLineSet* neutrinoTrack = new TEveStraightLineSet("Neutrino track");
  neutrinoTrack->SetLineColor(kMagenta);
  neutrinoTrack->SetLineStyle(kSolid);
  neutrinoTrack->SetLineWidth(4);
  TDatabasePDG myDB;
	
  if(mySummary->IsParticleGunEvent())
	{
    
    // If we have a pi-zero gun, make sure we treat it properly.
    if (mySummary->GetBeamPDG() == 111){
    }
    else{
		  double mm_to_cm = 0.1;
		  TVector3 vtx(mySummary->GetVertexX() * mm_to_cm,
		  						 mySummary->GetVertexY() * mm_to_cm,
		  						 mySummary->GetVertexZ() * mm_to_cm);
      TVector3 dir(mySummary->GetBeamDir().X(), mySummary->GetBeamDir().Y(), mySummary->GetBeamDir().Z());
      TVector3 end = vtx + WCSimGeometry::Instance()->ForwardProjectionToEdge(vtx.X(), vtx.Y(), vtx.Z(), dir.X(), dir.Y(), dir.Z()) * dir;
      trueCherenkovTracks->AddLine(vtx.X(), vtx.Y(), vtx.Z(), end.X(), end.Y(), end.Z());
    }
	}
	else
	{
    for(size_t iTrack = 0; iTrack < mySummary->GetNPrimaries(); ++iTrack)
    {
		  double mm_to_cm = 0.1;
      TParticlePDG * particle = myDB.GetParticle(mySummary->GetPrimaryPDG(iTrack));
		  for(unsigned int i = 0; i < mySummary->GetNPrimaries(); ++i)
		  {
        if( abs(particle->PdgCode()) == 12 || abs(particle->PdgCode()) == 14)
        {   
            TVector3 vtx = mySummary->GetVertex() * mm_to_cm;
            TVector3 dir = mySummary->GetPrimaryDir(iTrack);
            TVector3 end = vtx + WCSimGeometry::Instance()->BackwardProjectionToEdge(vtx.X(), vtx.Y(), vtx.Z(), dir.X(), dir.Y(), dir.Z()) * dir;
            neutrinoTrack->AddLine(vtx.X(), vtx.Y(), vtx.Z(), end.X(), end.Y(), end.Z());
        }

        if( (particle->Charge() == 0 && particle->PdgCode() != 22) || mySummary->GetPrimaryEnergy(iTrack) < 100)
        {
          continue;
        }

        if( (particle->Charge() != 0  || particle->PdgCode() == 22) && mySummary->GetPrimaryEnergy(iTrack) > 100)
        {
          double gamma = mySummary->GetPrimaryEnergy(iTrack) / (1000.0 * particle->Mass());
          double beta = sqrt( (gamma * gamma -1) / (gamma * gamma) );
          if( beta > 3/4. ) // Assume refractive index of 4/3 for the purposes of drawing
          {
            TVector3 vtx = mySummary->GetVertex() * mm_to_cm;
            TVector3 dir = mySummary->GetPrimaryDir(iTrack);
            TVector3 end = vtx + WCSimGeometry::Instance()->ForwardProjectionToEdge(vtx.X(), vtx.Y(), vtx.Z(), dir.X(), dir.Y(), dir.Z()) * dir;

            trueCherenkovTracks->AddLine(vtx.X(), vtx.Y(), vtx.Z(), end.X(), end.Y(), end.Z());
          }
        }
      }
    }
	}
  gEve->AddElement(trueCherenkovTracks);
  gEve->AddElement(neutrinoTrack);
  gEve->AddElement(leadingTrack);
}
