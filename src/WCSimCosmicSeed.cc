#include <iostream>
#include <cmath>
#include <algorithm> 

#include "WCSimCosmicSeed.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimGeometry.hh"
#include "WCSimRecoClusteringUtil.hh"
#include "WCSimParameters.hh"

#include <TVector2.h>
#include <TVector3.h>
#include <TMath.h>

ClassImp(WCSimCosmicSeed)

bool SortVetoClusters(std::vector<WCSimRecoDigit*> a, std::vector<WCSimRecoDigit*> b);

WCSimCosmicSeed::WCSimCosmicSeed(){
  fSuccess = false;
  fVetoDigits = 0x0;
  fInnerDigits = 0x0;
  fMaxDistBetweenRingPoints = -999.;
  fClusteringUtil = 0x0;
}

WCSimCosmicSeed::WCSimCosmicSeed(std::vector<WCSimRecoDigit*>* digits){
  fSuccess = false;
  fVetoDigits = digits;
  fInnerDigits = 0x0;
  fMaxDistBetweenRingPoints = -999.;
  fClusteringUtil = 0x0;
}

WCSimCosmicSeed::~WCSimCosmicSeed(){

  if(fClusteringUtil!=0x0){
    delete fClusteringUtil;
  }

}

void WCSimCosmicSeed::CalcSeedVtxAndDir(){

  if(fVetoDigits==0x0){
    std::cout << " WCSimCosmicSeed has a null veto digit vector, exiting." << std::endl;
    return;
  }
  if(fVetoDigits->size()==0){
    std::cout << " WCSimCosmicSeed has no veto digits, exiting." << std::endl;
    return;
  }

  this->RunSimpleEntrySeed();
  this->RunClusteringSeed();

  // Generate the projected ring points in case we want to shadow part
  // of the detector in the slicer.
  if(fSuccess){
    this->GenerateRingPoints();
  }
  else{
    this->RunSimpleCrossingSeed();
    if(fSuccess){
      this->GenerateRingPoints();
    }
  }
//  else{
    // Use the inner hits to try to get a direction to go with the simple seed.
    // Get the highest charge hit on the bottom for now.
//    this->GetSimpleDirectionFromInner();
//    if(fSuccess){
//      this->GenerateRingPoints();
//    }
//  }

}

void WCSimCosmicSeed::RunSimpleEntrySeed(){

  // Find the highest charge hit.
  double maxQ = -1;
  unsigned int maxIndex = 999999;
  for(unsigned int i = 0; i < fVetoDigits->size(); ++i){

    WCSimRecoDigit* digit = fVetoDigits->at(i);    

//    std::cout << "- Veto hit: " << digit->GetX() << ", " << digit->GetY() << ", " << digit->GetZ() << std::endl;

    // Top half only.
    if(digit->GetZ() < 0.) continue;

    double thisQ = digit->GetRawQPEs();
    if(thisQ > maxQ){
      maxQ = thisQ;
      maxIndex = i;
    } 
  }

  TVector3 vtx(0.,0.,0.);
  double vtxT = 0.0;
  if(maxIndex != 999999){
    WCSimRecoDigit* qDigit = fVetoDigits->at(maxIndex);
    vtx = TVector3(qDigit->GetX(),qDigit->GetY(),qDigit->GetZ());
    vtxT = qDigit->GetRawTime();
  }
  std::cout << " Best guess simple cosmic vertex = " << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << std::endl;
  fFittedVtx = vtx;
  fFittedVtxT = vtxT;
}

void WCSimCosmicSeed::RunSimpleCrossingSeed(){

  fSuccess = false;

  if(!this->HasVetoDigits()) return;

  bool isTopHit = false;
  bool isSideHit = false;
  bool isSideUpperHit = false;
  bool isSideLowerHit = false;
  bool isBottomHit = false;

  double qMaxTop = -999.;
  unsigned int indexMaxTop = 999999;
  double qMaxSide = -999.;
  unsigned int indexMaxSide = 999999;
  double qMaxSideUpper = -999.;
  unsigned int indexMaxSideUpper = 999999;
  double qMaxSideLower = -999.;
  unsigned int indexMaxSideLower = 999999;
  double qMaxBottom = -999.;
  unsigned int indexMaxBottom = 999999;

  WCSimGeometry* myGeo = WCSimGeometry::Instance();
  double h = myGeo->GetCylLength() - 10.0; // Give a little leeway.

  std::cout << "- Total of " << fVetoDigits->size() << " veto hits." << std::endl;

  for(unsigned int i = 0; i < fVetoDigits->size(); ++i){

    WCSimRecoDigit* digit = fVetoDigits->at(i);
    double thisQ = digit->GetRawQPEs();

    if(fabs(digit->GetZ()) < 0.5*h){
      isSideHit = true;
      if(thisQ > qMaxSide){
        qMaxSide = thisQ;
        indexMaxSide = i;
      }
      if(digit->GetZ() > 0){
        isSideUpperHit = true;
        if(thisQ > qMaxSideUpper){
          qMaxSideUpper = thisQ;
          indexMaxSideUpper = i;
        }
      }
      else{
        isSideLowerHit = true;
        if(thisQ > qMaxSideLower){
          qMaxSideLower = thisQ;
          indexMaxSideLower = i;
        }
      }
    }
    else{
      if(digit->GetZ() < 0){
        isBottomHit = true;
        if(thisQ > qMaxBottom){
          qMaxBottom = thisQ;
          indexMaxBottom = i;
        }
      }
      else{
        isTopHit = true;
        if(thisQ > qMaxTop){
          qMaxTop = thisQ;
          indexMaxTop = i;
        }
      }
    }
  }

  // Now have the highest charge hit in each of the regions.
  double minSeedQ = 4.0;
  WCSimRecoDigit* upperHit = 0x0;
  WCSimRecoDigit* lowerHit = 0x0;
  //int upperRegion = -1; // 0 = top, 1 = upper side, 2 = lower side

  // Update the booleans to see if we have sizable hits
  isTopHit = isTopHit && (qMaxTop > minSeedQ);
  isSideHit = isSideHit && (qMaxSide > minSeedQ);
  isSideUpperHit = isSideUpperHit && (qMaxSideUpper > minSeedQ);
  isSideLowerHit = isSideLowerHit && (qMaxSideLower > minSeedQ);
  isBottomHit = isBottomHit && (qMaxBottom > minSeedQ);

  // Look to see if we want to treat the barrel as a single object or two.
  bool isBarrelSeparate = false;
  WCSimRecoDigit* tempDig = 0x0;
  if(isSideUpperHit && isSideLowerHit){
    tempDig = fVetoDigits->at(indexMaxSideUpper); 
    TVector2 vSideUpper2D(tempDig->GetX(),tempDig->GetY());
    TVector3 vSideUpper3D(tempDig->GetX(),tempDig->GetY(),tempDig->GetZ()); 

    tempDig = fVetoDigits->at(indexMaxSideLower); 
    TVector2 vSideLower2D(tempDig->GetX(),tempDig->GetY());
    TVector3 vSideLower3D(tempDig->GetX(),tempDig->GetY(),tempDig->GetZ());

    double dPhi = vSideLower2D.DeltaPhi(vSideUpper2D);
    if((dPhi > TMath::Pi()/2.0) || ((vSideLower3D-vSideUpper3D).Mag() > 500)){
      isBarrelSeparate = true;
    }
  }

  // Treat the barrel as a single region
  if(!isBarrelSeparate){
    // If all three are hit we need to work out which one we should ignore.
    if(isTopHit && isSideHit && isBottomHit){
      if((qMaxTop > qMaxSide) && (qMaxSide > qMaxBottom)){
        isBottomHit = false;
      }
      else if((qMaxBottom > qMaxSide) && (qMaxSide > qMaxTop)){
        isTopHit = false;
      }
      else{
        isSideHit = false;
      }
    }

    // We are guaranteed to only have two of these true now, so make the seeds.
    if(isTopHit && isSideHit && !isBottomHit){
      upperHit = fVetoDigits->at(indexMaxTop);
      lowerHit = fVetoDigits->at(indexMaxSide);
    }
    else if(isTopHit && !isSideHit && isBottomHit){
      upperHit = fVetoDigits->at(indexMaxTop);
      lowerHit = fVetoDigits->at(indexMaxBottom);
    }
    else if(!isTopHit && isSideHit && isBottomHit){
      upperHit = fVetoDigits->at(indexMaxSide);
      lowerHit = fVetoDigits->at(indexMaxBottom);
    }
  }
  else{
    // Otherwise the upper and lower barrel regions could contain useful information
    bool mergeTopAndUpperSide = false;
    bool mergeBottomAndLowerSide = false;

    if(isTopHit && isSideUpperHit){
      tempDig = fVetoDigits->at(indexMaxTop);
      TVector3 topPos(tempDig->GetX(),tempDig->GetY(),tempDig->GetZ());
      tempDig = fVetoDigits->at(indexMaxSideUpper);
      TVector3 sidePos(tempDig->GetX(),tempDig->GetY(),tempDig->GetZ());
      if((topPos - sidePos).Mag() < 500){
        mergeTopAndUpperSide = true;
      }
    }
    if(isBottomHit && isSideLowerHit){
      tempDig = fVetoDigits->at(indexMaxBottom);
      TVector3 bottomPos(tempDig->GetX(),tempDig->GetY(),tempDig->GetZ());
      tempDig = fVetoDigits->at(indexMaxSideLower);
      TVector3 sidePos(tempDig->GetX(),tempDig->GetY(),tempDig->GetZ());
      if((bottomPos - sidePos).Mag() < 500){
        mergeBottomAndLowerSide = true;
      }
    }

    if(mergeTopAndUpperSide && mergeBottomAndLowerSide){
      // We just have two positions as expected.
      int upperIndex = indexMaxSideUpper;
      if(qMaxTop > qMaxSideUpper){
        upperIndex = indexMaxTop;
      }
      int lowerIndex = indexMaxSideLower;
      if(qMaxBottom > qMaxSideLower){
        lowerIndex = indexMaxBottom;
      }
      upperHit = fVetoDigits->at(upperIndex);
      lowerHit = fVetoDigits->at(lowerIndex);
    }
    else if(mergeTopAndUpperSide && !mergeBottomAndLowerSide){
      int upperIndex = indexMaxSideUpper;
      if(qMaxTop > qMaxSideUpper){
        upperIndex = indexMaxTop;
      } 
      upperHit = fVetoDigits->at(upperIndex);
      if(isBottomHit) lowerHit = fVetoDigits->at(indexMaxBottom);
      else lowerHit = fVetoDigits->at(indexMaxSideLower);
    }
    else if(!mergeTopAndUpperSide && mergeBottomAndLowerSide){
      int lowerIndex = indexMaxSideLower;
      if(qMaxBottom > qMaxSideLower){
        lowerIndex = indexMaxBottom;
      }
      if(isTopHit) upperHit = fVetoDigits->at(indexMaxTop);
      else upperHit = fVetoDigits->at(indexMaxSideUpper);
      lowerHit = fVetoDigits->at(lowerIndex); 
    }
    else{
      // No merging...
      std::cout << "Unable to merge top and top side and bottom and bottom side..." << std::endl;
      // Let's just try choosing the highest charge of top or top side, and bottom or lowerside.
      int topIndex = -999;
      int bottomIndex = -999;
      if(isTopHit || isSideUpperHit){
        double topCharge = -999.; 
        if(isTopHit){
          topIndex = indexMaxTop;
          topCharge = qMaxTop;
        }
        if(isSideUpperHit){
          if(qMaxSideUpper > topCharge){
            topIndex = indexMaxSideUpper;
          }
        }
      }
      if(isBottomHit || isSideLowerHit){
        double bottomCharge = -999.; 
        if(isBottomHit){
          bottomIndex = indexMaxBottom;
          bottomCharge = qMaxBottom;
        }
        if(isSideLowerHit){
          if(qMaxSideLower > bottomCharge){
            bottomIndex = indexMaxSideLower;
          }
        }
      }

      // If we have no top index, try just using the bottom regions.
      if((topIndex == -999) && (isBottomHit && isSideLowerHit)){
        topIndex = indexMaxSideLower;
        bottomIndex = indexMaxBottom;
      }
      // If we have no bottom index, try just using the top regions.
      if((bottomIndex == -999) && (isTopHit && isSideUpperHit)){
        topIndex = indexMaxTop;
        bottomIndex = indexMaxSideUpper;
      }

      if(topIndex != -999){
        upperHit = fVetoDigits->at(topIndex);
      }
      if(bottomIndex != -999){
        lowerHit = fVetoDigits->at(bottomIndex);
      }
    }

  }


  // Give up if there isn't a suitable seed.
  if(upperHit == 0x0){
    std::cout << "Cosmic seed failed: No upper hit" << std::endl;
    return;
  }
  // Check if we have a lower hit.
  if(lowerHit == 0x0){
    std::cout << "Cosmic seed failed: No lower hit" << std::endl;
    return;
  }

  // Now draw a line between our hits.
  TVector3 upperPos(upperHit->GetX(),upperHit->GetY(),upperHit->GetZ());
  TVector3 lowerPos(lowerHit->GetX(),lowerHit->GetY(),lowerHit->GetZ());
  TVector3 dir = lowerPos - upperPos;

  // Check there is some reasonable distance between them
  if(dir.Mag() < 1000){ 
    std::cout << "Cosmic seed failed: Distance between entry and exit too low" << std::endl;
    std::cout << "Upper Hit = " << upperPos.X() << ", " << upperPos.Y() << ", " << upperPos.Z() << std::endl;
    std::cout << "Lower Hit = " << lowerPos.X() << ", " << lowerPos.Y() << ", " << lowerPos.Z() << std::endl;
    std::cout << "Distance  = " << dir.Mag() << std::endl;
    return;
  }

  fFittedVtx = upperPos;
  fFittedDir = dir.Unit();
  fSuccess = true;
  
  std::cout << "Cosmic seed: " << std::endl;
  std::cout << "Upper Hit = " << upperPos.X() << ", " << upperPos.Y() << ", " << upperPos.Z() << std::endl;
  std::cout << "Lower Hit = " << lowerPos.X() << ", " << lowerPos.Y() << ", " << lowerPos.Z() << std::endl;
  std::cout << "Direction = " << fFittedDir.X() << ", " << fFittedDir.Y() << ", " << fFittedDir.Z() << std::endl;

  return;

}

void WCSimCosmicSeed::RunClusteringSeed(){

  // Load parameters from WCSimParameters
  double chargeCut = WCSimParameters::Instance()->GetVetoMinChargeCut();
  double maxCharge = WCSimParameters::Instance()->GetVetoMaxChargeCut();
  double timeCut = WCSimParameters::Instance()->GetVetoTimeCut();
  double distCut = WCSimParameters::Instance()->GetVetoClusterDistance();
  unsigned int minSize = WCSimParameters::Instance()->GetVetoMinSize();

  std::vector<std::vector<WCSimRecoDigit*> > vetoClusters;

  if(fClusteringUtil != 0x0){
    delete fClusteringUtil;
  }
  fClusteringUtil = new WCSimRecoClusteringUtil(*fVetoDigits,chargeCut,timeCut,distCut,minSize);
  fClusteringUtil->PerformClustering();
  vetoClusters = fClusteringUtil->GetFullSlicedDigits();

 // Do we need to increase the charge cut and try again?
  if(vetoClusters.size() < 2){
    while((vetoClusters.size() < 2) && (chargeCut < (maxCharge - 0.99))){
      ++chargeCut;
      std::cout << "WCSimCosmicSeed: Could not form two clusters. Raising charge cut to " << chargeCut << std::endl;
      fClusteringUtil->SetChargeCut(chargeCut);
      fClusteringUtil->PerformClustering();
      vetoClusters = fClusteringUtil->GetFullSlicedDigits();
    }
    // Do we have the right number of slices?
    if(vetoClusters.size() < 2){
      // It wasn't possible. Go back to default values.
      std::cout << "Not possible to find two veto clusters." << std::endl;
      chargeCut = 4.;
      fClusteringUtil->SetChargeCut(chargeCut);
      fClusteringUtil->PerformClustering();
      vetoClusters = fClusteringUtil->GetFullSlicedDigits();
    }
  }
  
  // Sort the clusters by number of hits to get the biggest first
  std::sort(vetoClusters.begin(),vetoClusters.end(),SortVetoClusters);

  if(vetoClusters.size() > 1){
    WCSimRecoDigit* hit0 = this->GetBiggestClusterHit(vetoClusters.at(0));
    WCSimRecoDigit* hit1 = this->GetBiggestClusterHit(vetoClusters.at(1));

    // Swap hits if hit1 is at a higher Z position
    if(hit1->GetZ() > hit0->GetZ()){
      WCSimRecoDigit* temp = hit1;
      hit1 = hit0;
      hit0 = temp;
    }

    TVector3 upper(hit0->GetX(),hit0->GetY(),hit0->GetZ());
    TVector3 lower(hit1->GetX(),hit1->GetY(),hit1->GetZ());
    TVector3 diff = lower - upper;

    fFittedVtx = upper;
    fFittedDir = diff.Unit();
    fSuccess = true;

    std::cout << "Cosmic seed: " << std::endl;
    std::cout << "Upper Hit = " << upper.X() << ", " << upper.Y() << ", " << upper.Z() << std::endl;
    std::cout << "Lower Hit = " << lower.X() << ", " << lower.Y() << ", " << lower.Z() << std::endl;
    std::cout << "Direction = " << fFittedDir.X() << ", " << fFittedDir.Y() << ", " << fFittedDir.Z() << std::endl;
  } 

}

WCSimRecoDigit* WCSimCosmicSeed::GetBiggestClusterHit(std::vector<WCSimRecoDigit*> clust){

  double maxQ = -999.;
  unsigned int maxIndex = 999999;
  for(unsigned int i = 0; i < clust.size(); ++i){
    WCSimRecoDigit* dig = clust.at(i);
    if(dig->GetRawQPEs() > maxQ){
      maxQ = dig->GetRawQPEs();
      maxIndex = i;
    }
  }
  return clust.at(maxIndex);
}

void WCSimCosmicSeed::GetSimpleDirectionFromInner(){

  if(fInnerDigits == 0x0){
    std::cerr << "No inner detector digits passed to the cosmic seeding." << std::endl;
    return;
  } 

  // Find the highest charge inner detector hit on the bottom of the detector.
  // Use this to form a cosmic direction with the veto simple seed.
  int bottomMaxIndex = -999;
  double bottomMaxQ = -999.;
  for(unsigned int i = 0; i < fInnerDigits->size(); ++i){
    WCSimRecoDigit* digit = fInnerDigits->at(i);
    if(digit->GetRegion() != 0) continue; // Top
    if(digit->GetRegion() == 1 && digit->GetX() > 0) continue; // Front barrel
    
    double thisQ = digit->GetRawQPEs();
    if(thisQ > bottomMaxQ){
      bottomMaxQ = thisQ;
      bottomMaxIndex = i; 
    }
  }

  WCSimRecoDigit* digit = fInnerDigits->at(bottomMaxIndex);
  TVector3 bottomPos(digit->GetX(),digit->GetY(),digit->GetZ());

  TVector3 dir = bottomPos - fFittedVtx;
  fFittedDir = dir.Unit();

  std::cout << "Using cosmic seed and direction: " << std::endl;
  std::cout << " - " << fFittedVtx.X() << "," << fFittedVtx.Y() << "," << fFittedVtx.Z() << std::endl;
  std::cout << " - " << bottomPos.X() << "," << bottomPos.Y() << "," << bottomPos.Z() << std::endl;
  std::cout << " - " << fFittedDir.X() << "," << fFittedDir.Y() << "," << fFittedDir.Z() << std::endl;

  fSuccess = true;
}

bool WCSimCosmicSeed::HasVetoDigits(){

  if(fVetoDigits->size() == 0){
    return false;
  }
  else{
    return true;
  }

}

void WCSimCosmicSeed::GenerateRingPoints(){

  if(!fSuccess) return;

  WCSimGeometry *myGeo = WCSimGeometry::Instance();

  // Projected vertex on the detector wall
  Double_t xproj = 0.0;
  Double_t yproj = 0.0;
  Double_t zproj = 0.0;

  // Points on the ring
  Double_t xring = 0.0;
  Double_t yring = 0.0;
  Double_t zring = 0.0;

  // Ring points projected on to the wall
  Double_t x = 0.0;
  Double_t y = 0.0;
  Double_t z = 0.0;

  // Direction to the projected ring points
  Double_t nx = 0.0;
  Double_t ny = 0.0;
  Double_t nz = 0.0;

  // Ring radius - we don't care about this here
  Double_t radius = 0.0;

  // Region of the detector - not needed here
  Int_t region = -1;

  // project vertex onto edge of detector
  myGeo->ProjectToFarEdge(fFittedVtx.X(),fFittedVtx.Y(),fFittedVtx.Z(),fFittedDir.X(),fFittedDir.Y(),fFittedDir.Z(), xproj, yproj, zproj, region );

  // draw ring from many markers...
  Int_t nPoints = 360;

  Double_t cone_angle = 42.0;        // both angles
  Double_t delta_omega = 360.0/double(nPoints); // in degrees

//  std::cout << "Cosmic event ring" << std::endl;
//  std::cout << "x/D:y:z" << std::endl;
  for( Int_t n=0; n<nPoints; n++ ){

    myGeo->FindCircle( xproj, yproj, zproj,
                               fFittedVtx.X(), fFittedVtx.Y(), fFittedVtx.Z(),
                               cone_angle, n*delta_omega,
                               xring, yring, zring,
                               nx, ny, nz,
                               radius );

    myGeo->ProjectToFarEdge(fFittedVtx.X(),fFittedVtx.Y(),fFittedVtx.Z(),nx,ny,nz,x,y,z,region);

    fRingPosX.push_back(x);
    fRingPosY.push_back(y);
    fRingPosZ.push_back(z);

//    std::cout << x << " " << y << " " << z << std::endl;
  }

}

bool WCSimCosmicSeed::IsHitShadowed(double hx, double hy, double hz){

  bool isShadowed = false;

  TVector3 hitPos(hx,hy,hz);

  // Best thing to do is to check if the hit is within a 42.0 cone angle of the cosmic seed
  TVector3 hitToVtx = hitPos - TVector3(fFittedVtx.X(),fFittedVtx.Y(),fFittedVtx.Z());
  hitToVtx = hitToVtx.Unit();
  TVector3 fittedDir(fFittedDir.X(),fFittedDir.Y(),fFittedDir.Z());

  double angle = (180.0 / TMath::Pi()) * hitToVtx.Angle(fittedDir);

//  std::cout << "Hit to Vtx   = " << hitToVtx.X() << "," << hitToVtx.Y() << "," << hitToVtx.Z() << std::endl;
//  std::cout << "Direction    = " << fittedDir.X() << "," << fittedDir.Y() << "," << fittedDir.Z() << std::endl;
//  std::cout << "Angle to hit = " << angle << std::endl;

  if(angle <= 42.0){
    isShadowed = true;
  }

  return isShadowed;

}

// Sort the veto clusters by size (biggest first)
bool SortVetoClusters(std::vector<WCSimRecoDigit*> a, std::vector<WCSimRecoDigit*> b){
  unsigned int aSize = a.size();
  unsigned int bSize = b.size();
  return aSize > bSize;
}
