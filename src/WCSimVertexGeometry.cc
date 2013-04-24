#include "WCSimVertexGeometry.hh"

#include "WCSimRecoEvent.hh"
#include "WCSimRecoVertex.hh"
#include "WCSimRecoDigit.hh"

#include "WCSimGeometry.hh"
#include "WCSimParameters.hh"

#include "TMath.h"
#include "TRandom.h"

#include <vector>

#include <cmath>
#include <iostream>
#include <cassert>

ClassImp(WCSimVertexGeometry)

static WCSimVertexGeometry* fgVertexGeometry = 0;

WCSimVertexGeometry* WCSimVertexGeometry::Instance()
{
  if( !fgVertexGeometry ){
    fgVertexGeometry = new WCSimVertexGeometry();
  }

  if( !fgVertexGeometry ){
    assert(fgVertexGeometry);
  }

  if( fgVertexGeometry ){

  }

  return fgVertexGeometry;
}

WCSimVertexGeometry::WCSimVertexGeometry()
{
  if( WCSimGeometry::TouchGeometry() ){
    fPMTs = WCSimGeometry::Instance()->GetNumPMTs();
  }
  else{
    fPMTs = 100000; // maximum number of digits
  }
  
  fRunNum = -1;
  fEventNum = -1;
  fTriggerNum = -1;
  
  fNDigits = 0;
  fNFilterDigits = 0;

  fThisDigit = 0;
  fLastEntry = 0; 
  fCounter = 0;

  fMeanQ = 0.0;
  fTotalQ = 0.0;
  fMeanFilteredQ = 0.0;
  fTotalFilteredQ = 0.0;

  fMinTime = 0.0;
  fMaxTime = 0.0;
  
  fVtxX1 = 0.0;
  fVtxY1 = 0.0;
  fVtxZ1 = 0.0;
  fVtxTime1 = 0.0;

  fVtxX2 = 0.0;
  fVtxY2 = 0.0;
  fVtxZ2 = 0.0;
  fVtxTime2 = 0.0;

  fIsFiltered = new Bool_t[fPMTs];

  fDigitX = new Double_t[fPMTs];
  fDigitY = new Double_t[fPMTs];
  fDigitZ = new Double_t[fPMTs];
  fDigitT = new Double_t[fPMTs];
  fDigitQ = new Double_t[fPMTs];

  fConeAngle = new Double_t[fPMTs];
  fZenith = new Double_t[fPMTs];
  fAzimuth = new Double_t[fPMTs];
  fSolidAngle = new Double_t[fPMTs];

  fDistPoint = new Double_t[fPMTs];
  fDistTrack = new Double_t[fPMTs];
  fDistPhoton = new Double_t[fPMTs];
  fDistScatter = new Double_t[fPMTs];

  fDeltaTime = new Double_t[fPMTs];
  fDeltaSigma = new Double_t[fPMTs];

  fDeltaAngle = new Double_t[fPMTs];
  fDeltaPoint = new Double_t[fPMTs];
  fDeltaTrack = new Double_t[fPMTs];
  fDeltaPhoton = new Double_t[fPMTs];
  fDeltaScatter = new Double_t[fPMTs];

  fPointPath = new Double_t[fPMTs];
  fExtendedPath = new Double_t[fPMTs];

  fPointResidual = new Double_t[fPMTs];
  fExtendedResidual = new Double_t[fPMTs];

  fDelta = new Double_t[fPMTs];

  for( Int_t n=0; n<fPMTs; n++ ){
    fIsFiltered[n] = 0.0;

    fDigitX[n] = 0.0;
    fDigitY[n] = 0.0;
    fDigitZ[n] = 0.0;
    fDigitT[n] = 0.0;
    fDigitQ[n] = 0.0;

    fConeAngle[n] = 0.0;
    fZenith[n] = 0.0;
    fAzimuth[n] = 0.0;
    fSolidAngle[n] = 0.0;

    fDistPoint[n] = 0.0;
    fDistTrack[n] = 0.0;
    fDistPhoton[n] = 0.0; 
    fDistScatter[n] = 0.0; 

    fDeltaTime[n] = 0.0;
    fDeltaSigma[n] = 0.0;

    fDeltaAngle[n] = 0.0;    
    fDeltaPoint[n] = 0.0;
    fDeltaTrack[n] = 0.0;
    fDeltaPhoton[n] = 0.0;
    fDeltaScatter[n] = 0.0;

    fPointPath[n] = 0.0;
    fExtendedPath[n] = 0.0;
    fPointResidual[n] = 0.0;
    fExtendedResidual[n] = 0.0;

    fDelta[n] = 0.0;
  }
}

WCSimVertexGeometry::~WCSimVertexGeometry()
{
  if( fIsFiltered ) delete [] fIsFiltered;

  if( fDigitX ) delete [] fDigitX;
  if( fDigitY ) delete [] fDigitY;
  if( fDigitZ ) delete [] fDigitZ;
  if( fDigitT ) delete [] fDigitT;
  if( fDigitQ ) delete [] fDigitQ; 

  if( fConeAngle ) delete [] fConeAngle;
  if( fZenith ) delete [] fZenith;
  if( fAzimuth ) delete [] fAzimuth;
  if( fSolidAngle ) delete [] fSolidAngle;

  if( fDistPoint )  delete [] fDistPoint;
  if( fDistTrack )  delete [] fDistTrack;
  if( fDistPhoton ) delete [] fDistPhoton;   
  if( fDistScatter ) delete [] fDistScatter;

  if( fDeltaTime )   delete [] fDeltaTime;
  if( fDeltaSigma ) delete [] fDeltaSigma;

  if( fDeltaAngle )  delete [] fDeltaAngle;  
  if( fDeltaPoint )  delete [] fDeltaPoint;
  if( fDeltaTrack )  delete [] fDeltaTrack;
  if( fDeltaPhoton ) delete [] fDeltaPhoton;
  if( fDeltaScatter ) delete [] fDeltaScatter;

  if( fPointPath ) delete [] fPointPath;
  if( fExtendedPath ) delete [] fExtendedPath;
  if( fPointResidual )  delete [] fPointResidual;
  if( fExtendedResidual ) delete [] fExtendedResidual;

  if( fDelta ) delete [] fDelta;  
}

void WCSimVertexGeometry::LoadEvent(WCSimRecoEvent* myEvent)
{
  // resetting arrays
  // ================
  if( myEvent==0 ) return;

  Int_t thisRunNum = myEvent->GetRun();
  Int_t thisEventNum = myEvent->GetEvent();
  Int_t thisTriggerNum = myEvent->GetTrigger();

  if( fRunNum==thisRunNum
   && fEventNum==thisEventNum
   && fTriggerNum==thisTriggerNum ){
    return;
  }

  fRunNum = thisRunNum;
  fEventNum = thisEventNum;
  fTriggerNum = thisTriggerNum;

  fNDigits = 0;
  fNFilterDigits = 0;
  
  fThisDigit = 0;
  fLastEntry = 0;
  fCounter = 0;

  // clear vertices
  // ==============
  for(UInt_t i=0; i<vVertexList.size(); i++ ){
    delete (WCSimRecoVertex*)(vVertexList.at(i));
  }
  vVertexList.clear();

  // clear seed vertices
  // ===================
  vSeedVtxX.clear();
  vSeedVtxY.clear();
  vSeedVtxZ.clear();
  vSeedVtxTime.clear();
  vSeedDigitList.clear();

  // get lists of digits
  // ==================
  std::vector<WCSimRecoDigit*>* myDigitList = (std::vector<WCSimRecoDigit*>*)(myEvent->GetDigitList());
  fNDigits = myDigitList->size();

  std::vector<WCSimRecoDigit*>* myFilterDigitList = (std::vector<WCSimRecoDigit*>*)(myEvent->GetFilterDigitList());
  fNFilterDigits = myFilterDigitList->size();

  // sanity checks
  // =============
  if( fNDigits<=0 ){
    std::cout << " *** WCSimVertexGeometry::LoadEvent(...) *** " << std::endl;
    std::cout << "   <warning> event has no digits! " << std::endl;
    return;
  }

  if( fNDigits>fPMTs ){
    std::cout << " *** WCSimVertexGeometry::LoadEvent(...) *** " << std::endl;
    std::cout << "   <warning> event will be truncated! " << std::endl;
    fNDigits = fPMTs;
  }

  if( fNFilterDigits>fNDigits ){
    fNFilterDigits = fNDigits;
  }

  // load digits
  // ===========
  fTotalQ = 0.0;
  fMeanQ = 0.0;

  fTotalFilteredQ = 0.0;
  fMeanFilteredQ = 0.0;

  fMinTime = -999999.9;
  fMaxTime = -999999.9;

  Double_t Swx = 0.0;
  Double_t Sw = 0.0;

  Double_t SFwx = 0.0;
  Double_t SFw = 0.0;
  
  for( Int_t idigit=0; idigit<fNDigits; idigit++ ){
    WCSimRecoDigit* recoDigit = (WCSimRecoDigit*)(myDigitList->at(idigit));    

    fDigitX[idigit] = recoDigit->GetX();
    fDigitY[idigit] = recoDigit->GetY();
    fDigitZ[idigit] = recoDigit->GetZ();
    fDigitT[idigit] = recoDigit->GetTime();
    fDigitQ[idigit] = recoDigit->GetQPEs();

    fIsFiltered[idigit] = recoDigit->IsFiltered();

    fConeAngle[idigit] = 0.0;
    fZenith[idigit] = 0.0;
    fAzimuth[idigit] = 0.0;
    fSolidAngle[idigit] = 0.0;

    fDistPoint[idigit] = 0.0;
    fDistTrack[idigit] = 0.0;
    fDistPhoton[idigit] = 0.0;
    fDistScatter[idigit] = 0.0;

    fDeltaTime[idigit] = 0.0;
    fDeltaSigma[idigit] = 0.0;
    
    fDeltaAngle[idigit] = 0.0;
    fDeltaPoint[idigit] = 0.0;
    fDeltaTrack[idigit] = 0.0;
    fDeltaPhoton[idigit] = 0.0;
    fDeltaScatter[idigit] = 0.0;

    fPointPath[idigit] = 0.0;
    fExtendedPath[idigit] = 0.0;
    fPointResidual[idigit] = 0.0;
    fExtendedResidual[idigit] = 0.0;

    fDelta[idigit] = 0.0;

    Swx += recoDigit->GetQPEs();
    Sw += 1.0;

    if( fMinTime<0 
     || recoDigit->GetTime()<fMinTime ){
      fMinTime = recoDigit->GetTime();
    }

    if( fMaxTime<0 
     || recoDigit->GetTime()>fMaxTime ){
      fMaxTime = recoDigit->GetTime();
    }

    if( recoDigit->IsFiltered() ){
      SFwx += recoDigit->GetQPEs();
      SFw += 1.0;
    }
  }

  if( Sw>0.0 ){
    fTotalQ = Swx;
    fMeanQ = Swx/Sw;
  }

  if( SFw>0.0 ){
    fTotalFilteredQ = SFwx;
    fMeanFilteredQ = SFwx/SFw;
  }

  if( fMinTime<0 ) fMinTime = 0.0;
  if( fMaxTime<0 ) fMaxTime = 0.0;

  return;
}
  
WCSimRecoVertex* WCSimVertexGeometry::CalcSimpleVertex(WCSimRecoEvent* myEvent)
{
  // load event
  // ==========
  this->LoadEvent(myEvent);

  // calculate simple vertex
  // =======================
  return this->CalcSimpleVertex(myEvent);
}

WCSimRecoVertex* WCSimVertexGeometry::CalcSimpleVertex()
{
  // simple vertex
  // =============
  Double_t vtxX = 0.0;
  Double_t vtxY = 0.0;
  Double_t vtxZ = 0.0;
  Double_t vtxTime = 0.0;

  Int_t itr = 1;
  Bool_t pass = 1; 
  Double_t fom = 1.0;

  // create new vertex
  // =================
  WCSimRecoVertex* newVertex = new WCSimRecoVertex();
  vVertexList.push_back(newVertex);

  // calculate vertex
  // ================
  this->CalcSimpleVertex(vtxX,vtxY,vtxZ,vtxTime);

  // set vertex
  // ==========
  newVertex->SetVertex(vtxX,vtxY,vtxZ,vtxTime);
  newVertex->SetFOM(fom,itr,pass);  

  // set status
  // ==========
  newVertex->SetStatus(WCSimRecoVertex::kOK);
    
  // return vertex
  // =============
  return newVertex;
}

void WCSimVertexGeometry::CalcSimpleVertex(Double_t& vtxX, Double_t& vtxY, Double_t& vtxZ, Double_t& vtxTime)
{
  // simple vertex
  // =============
  // just calculate average position of digits

  // default vertex
  // ==============
  vtxX = 0.0;
  vtxY = 0.0;
  vtxZ = 0.0;
  vtxTime = 0.0;  

  // loop over digits
  // ================
  Double_t Swx = 0.0;
  Double_t Swy = 0.0;
  Double_t Swz = 0.0;
  Double_t Swt = 0.0;
  Double_t Sw = 0.0;

  for( Int_t idigit=0; idigit<fNDigits; idigit++ ){
    if( fIsFiltered[idigit] ){
      Swx += fDigitQ[idigit]*fDigitX[idigit];
      Swy += fDigitQ[idigit]*fDigitY[idigit];
      Swz += fDigitQ[idigit]*fDigitZ[idigit];
      Swt += fDigitQ[idigit]*fDigitT[idigit];
      Sw  += fDigitQ[idigit];
    }
  }

  // average position
  // ================
  if( Sw>0.0 ){
    vtxX = Swx/Sw;
    vtxY = Swy/Sw;
    vtxZ = Swz/Sw;
    vtxTime = Swt/Sw;
  }   

  return;
}
  
WCSimRecoVertex* WCSimVertexGeometry::CalcSimpleDirection(WCSimRecoEvent* myEvent, WCSimRecoVertex* myVertex)
{
  // load event
  // ==========
  this->LoadEvent(myEvent);

  // calculate simple direction
  // ==========================
  return this->CalcSimpleDirection(myVertex);
}

WCSimRecoVertex* WCSimVertexGeometry::CalcSimpleDirection(WCSimRecoVertex* myVertex)
{
  // load vertex
  // ===========
  Double_t vtxX = myVertex->GetX();
  Double_t vtxY = myVertex->GetY();
  Double_t vtxZ = myVertex->GetZ();
  Double_t vtxTime = myVertex->GetTime();
    
  // current status
  // ==============
  Int_t status = myVertex->GetStatus();

  // create new vertex
  // =================
  WCSimRecoVertex* newVertex = new WCSimRecoVertex();
  vVertexList.push_back(newVertex);

  // loop over digits
  // ================
  Double_t Swx = 0.0;
  Double_t Swy = 0.0;
  Double_t Swz = 0.0;
  Double_t Sw = 0.0;

  for( Int_t idigit=0; idigit<fNDigits; idigit++ ){
    if( fIsFiltered[idigit] ){
      Double_t q = fDigitQ[idigit];

      Double_t dx = fDigitX[idigit] - vtxX;
      Double_t dy = fDigitY[idigit] - vtxY;
      Double_t dz = fDigitZ[idigit] - vtxZ;
      Double_t ds = sqrt(dx*dx+dy*dy+dz*dz);

      Double_t px = dx/ds;
      Double_t py = dy/ds;
      Double_t pz = dz/ds;

      Swx += q*px;
      Swy += q*py;
      Swz += q*pz;
      Sw  += q;
    }
  }

  // average direction
  // =================
  Double_t dirX = 0.0;
  Double_t dirY = 0.0;
  Double_t dirZ = 0.0;
  
  Int_t itr = 0;
  Bool_t pass = 0; 
  Double_t fom = 0.0;

  if( Sw>0.0 ){
    Double_t qx = Swx/Sw;
    Double_t qy = Swy/Sw;
    Double_t qz = Swz/Sw;
    Double_t qs = sqrt(qx*qx+qy*qy+qz*qz);

    dirX = qx/qs;
    dirY = qy/qs;
    dirZ = qz/qs;

    fom = 1.0;
    itr = 1;
    pass = 1; 
  }

  // set vertex and direction
  // ========================
  if( pass ){
    newVertex->SetVertex(vtxX,vtxY,vtxZ,vtxTime);
    newVertex->SetDirection(dirX,dirY,dirZ);
    newVertex->SetFOM(fom,itr,pass);
  }

  // set status
  // ==========
  if( !pass ) status |= WCSimRecoVertex::kFailSimpleDirection;
  newVertex->SetStatus(status);

  // return vertex
  // =============
  return newVertex;
}

void WCSimVertexGeometry::CalcResiduals(WCSimRecoEvent* myEvent, WCSimRecoVertex* myVertex)
{
  // load event
  // ==========
  this->LoadEvent(myEvent);

  // calculate residuals
  // ===================
  return this->CalcResiduals(myVertex);
}

void WCSimVertexGeometry::CalcResiduals(WCSimRecoVertex* myVertex)
{
  // sanity check
  // ============
  if( myVertex==0 ) return;

  // get vertex information
  // ======================
  Double_t vtxX = myVertex->GetX();
  Double_t vtxY = myVertex->GetY();
  Double_t vtxZ = myVertex->GetZ();
  Double_t vtxTime = myVertex->GetTime();

  Double_t dirX = myVertex->GetDirX();
  Double_t dirY = myVertex->GetDirY();
  Double_t dirZ = myVertex->GetDirZ();

  // calculate residuals
  // ===================
  return this->CalcResiduals( vtxX, vtxY, vtxZ, vtxTime,
                              dirX, dirY, dirZ );
}

void WCSimVertexGeometry::CalcPointResiduals(Double_t vtxX, Double_t vtxY, Double_t vtxZ, Double_t vtxTime, Double_t dirX, Double_t dirY, Double_t dirZ)
{
  this->CalcResiduals( vtxX, vtxY, vtxZ, vtxTime,
                       dirX, dirY, dirZ );

  for( Int_t idigit=0; idigit<fNDigits; idigit++ ){
    fDelta[idigit] = fPointResidual[idigit];
  }

  return;
}

void WCSimVertexGeometry::CalcExtendedResiduals(Double_t vtxX, Double_t vtxY, Double_t vtxZ, Double_t vtxTime, Double_t dirX, Double_t dirY, Double_t dirZ )
{
  this->CalcResiduals( vtxX, vtxY, vtxZ, vtxTime,
                       dirX, dirY, dirZ );

  for( Int_t idigit=0; idigit<fNDigits; idigit++ ){
    fDelta[idigit] = fExtendedResidual[idigit];
  }

  return;
}

void WCSimVertexGeometry::CalcResiduals(Double_t vtxX, Double_t vtxY, Double_t vtxZ, Double_t vtxTime, Double_t dirX, Double_t dirY, Double_t dirZ )
{
  // reset arrays
  // ============
  for( Int_t idigit=0; idigit<fNDigits; idigit++ ){
    fConeAngle[idigit] = 0.0;
    fZenith[idigit] = 0.0;
    fAzimuth[idigit] = 0.0;
    fSolidAngle[idigit] = 0.0;
    fDistPoint[idigit] = 0.0;
    fDistTrack[idigit] = 0.0;
    fDistPhoton[idigit] = 0.0;  
    fDistScatter[idigit] = 0.0; 
    fDeltaTime[idigit] = 0.0; 
    fDeltaSigma[idigit] = 0.0;
    fDeltaAngle[idigit] = 0.0;
    fDeltaPoint[idigit] = 0.0;
    fDeltaTrack[idigit] = 0.0;
    fDeltaPhoton[idigit] = 0.0;
    fDeltaScatter[idigit] = 0.0;
    fPointPath[idigit] = 0.0;
    fExtendedPath[idigit] = 0.0;
    fPointResidual[idigit] = 0.0;
    fExtendedResidual[idigit] = 0.0;
    fDelta[idigit] = 0.0;
  }

  // cone angle
  // ==========
  Double_t thetadeg = WCSimParameters::CherenkovAngle(); // degrees
  Double_t theta = thetadeg*(TMath::Pi()/180.0); // degrees->radians

  // loop over digits
  // ================
  for( Int_t idigit=0; idigit<fNDigits; idigit++ ){
    Double_t dx = fDigitX[idigit]-vtxX;
    Double_t dy = fDigitY[idigit]-vtxY;
    Double_t dz = fDigitZ[idigit]-vtxZ;
    Double_t ds = sqrt(dx*dx+dy*dy+dz*dz);

    Double_t px = dx/ds;
    Double_t py = dy/ds;
    Double_t pz = dz/ds;

    Double_t cosphi = 1.0;
    Double_t sinphi = 1.0;
    Double_t phi = 0.0; 
    Double_t phideg = 0.0;

    Double_t ax = 0.0;
    Double_t ay = 0.0;
    Double_t az = 0.0;
    Double_t azideg = 0.0;

    // calculate angles if direction is known
    if( dirX*dirX + dirY*dirY + dirZ*dirZ>0.0 ){

      // zenith angle
      cosphi = px*dirX+py*dirY+pz*dirZ;
      phi = acos(cosphi); // radians
      phideg = phi/(TMath::Pi()/180.0); // radians->degrees
      sinphi = sqrt(1.0-cosphi*cosphi);
      sinphi += 0.24*exp(-sinphi/0.24);
      sinphi /= 0.684;  // sin(phideg)/sin(thetadeg)

      // azimuthal angle
      if( dirX*dirX+dirY*dirY>0.0 ){
        ax = (px*dirZ-pz*dirX) - (py*dirX-px*dirY)*(1.0-dirZ)*dirY/sqrt(dirX*dirX+dirY*dirY);
        ay = (py*dirZ-pz*dirY) - (px*dirY-py*dirX)*(1.0-dirZ)*dirX/sqrt(dirX*dirX+dirY*dirY);
        az = pz*dirZ + py*dirY + px*dirX;
      }
      else{
        ax = px;
        ay = py;
        az = pz;
      }

      azideg = atan2(ay,ax)/(TMath::Pi()/180.0); // radians->degrees
    }

    Double_t Lpoint = ds;
    Double_t Ltrack = 0.0;
    Double_t Lphoton = 0.0;
    Double_t Lscatter = 0.0;
 
    if( phi<theta ){
      Ltrack = Lpoint*sin(theta-phi)/sin(theta);
      Lphoton = Lpoint*sin(phi)/sin(theta);
      Lscatter = 0.0;
    }
    else{
      Ltrack = 0.0;
      Lphoton = Lpoint;
      Lscatter = Lpoint*(phi-theta);
    }

    Double_t fC = WCSimParameters::SpeedOfLight();
    Double_t fN = WCSimParameters::RefractiveIndex(Lphoton);

    Double_t dt = fDigitT[idigit] - vtxTime;
    Double_t qpes = fDigitQ[idigit];
    Double_t res = WCSimParameters::TimeResolution(qpes);

    fConeAngle[idigit] = thetadeg; // degrees
    fZenith[idigit] = phideg;      // degrees
    fAzimuth[idigit] = azideg;     // degrees
    fSolidAngle[idigit] = sinphi;

    fDistPoint[idigit] = Lpoint;       
    fDistTrack[idigit] = Ltrack;
    fDistPhoton[idigit] = Lphoton;
    fDistScatter[idigit] = Lscatter;

    fDeltaTime[idigit] = dt;
    fDeltaSigma[idigit] = res;   
    
    fDeltaAngle[idigit] = phideg-thetadeg; // degrees
    fDeltaPoint[idigit] = Lpoint/(fC/fN);
    fDeltaTrack[idigit] = Ltrack/fC;
    fDeltaPhoton[idigit] = Lphoton/(fC/fN);
    fDeltaScatter[idigit] = Lscatter/(fC/fN);
 
    fPointPath[idigit] = fN*Lpoint;
    fExtendedPath[idigit] = Ltrack + fN*Lphoton;

    fPointResidual[idigit] = dt - Lpoint/(fC/fN);
    fExtendedResidual[idigit] = dt - Ltrack/fC - Lphoton/(fC/fN);

    fDelta[idigit] = fExtendedResidual[idigit]; // default
  }

  return;
}   

Double_t WCSimVertexGeometry::GetDeltaCorrection(Int_t idigit, Double_t Length)
{
  if( Length<=0.0
   || GetDistTrack(idigit)<Length ){
    return 0.0;
  }

  else{
    Double_t Lpoint = GetDistPoint(idigit);
    Double_t Ltrack = GetDistTrack(idigit);
    Double_t Lphoton = GetDistPhoton(idigit);
    Double_t AngleRad = (TMath::Pi()/180.0)*GetAngle(idigit);    
    Double_t ConeAngleRad = (TMath::Pi()/180.0)*GetConeAngle(idigit);  

    Double_t LtrackNew = Length;
    Double_t LphotonNew = sqrt( Lpoint*Lpoint + Length*Length
                                -2.0*Lpoint*Length*cos(AngleRad) );

    Double_t theta = ConeAngleRad;
    Double_t sinphi = (Lpoint/LphotonNew)*sin(AngleRad);
    Double_t phi = asin(sinphi);
    Double_t alpha = theta-phi;
    Double_t LphotonNewCorrected = LphotonNew*alpha/sin(alpha);

    Double_t fC = WCSimParameters::SpeedOfLight();
    Double_t fN = WCSimParameters::RefractiveIndex(Lpoint);

    return ( Ltrack/fC + Lphoton/(fC/fN) )
         - ( LtrackNew/fC + LphotonNewCorrected/(fC/fN) );
  }
}

void WCSimVertexGeometry::CalcVertexSeeds(WCSimRecoEvent* myEvent, Int_t NSeeds)
{
  // load event
  // ==========
  this->LoadEvent(myEvent);

  // calculate vertex seeds
  // ======================
  return this->CalcVertexSeeds(NSeeds);
}
  
void WCSimVertexGeometry::CalcVertexSeeds(Int_t NSeeds)
{
  // reset list of seeds
  // ===================
  vSeedVtxX.clear();
  vSeedVtxY.clear();
  vSeedVtxZ.clear();
  vSeedVtxTime.clear();

  // always calculate the simple vertex
  // ==================================
  this->CalcSimpleVertex(fVtxX1,fVtxY1,fVtxZ1,fVtxTime1);

  // add this vertex
  vSeedVtxX.push_back(fVtxX1); 
  vSeedVtxY.push_back(fVtxY1);
  vSeedVtxZ.push_back(fVtxZ1);
  vSeedVtxTime.push_back(fVtxTime1);

  // check limit
  if( NSeeds<=1 ) return;


  // form list of golden digits
  // ==========================
  vSeedDigitList.clear();
  
  for( fThisDigit=0; fThisDigit<fNDigits; fThisDigit++ ){
    if( fIsFiltered[fThisDigit] ){
      if( fDigitT[fThisDigit] - fMinTime>=0
       && fDigitT[fThisDigit] - fMinTime<300
       && fDigitQ[fThisDigit]>=5.0 ){
        vSeedDigitList.push_back(fThisDigit);
      }
    }
  }

  // check for enough digits
  if( vSeedDigitList.size()<=4 ) return;


  // generate new list of seeds
  // ==========================
 
  Double_t x0 = 0.0;
  Double_t y0 = 0.0;
  Double_t z0 = 0.0;
  Double_t t0 = 0.0;

  Double_t x1 = 0.0;
  Double_t y1 = 0.0;
  Double_t z1 = 0.0;
  Double_t t1 = 0.0; 

  Double_t x2 = 0.0;
  Double_t y2 = 0.0;
  Double_t z2 = 0.0;
  Double_t t2 = 0.0;

  Double_t x3 = 0.0;
  Double_t y3 = 0.0;
  Double_t z3 = 0.0;
  Double_t t3 = 0.0;

  UInt_t counter = 0;
  UInt_t NSeedsTarget = NSeeds;

  while( GetNSeeds()<NSeedsTarget && counter<100*NSeedsTarget ){
    counter++;

    // choose next four digits
    this->ChooseNextQuadruple(x0,y0,z0,t0,
                              x1,y1,z1,t1,
                              x2,y2,z2,t2,
                              x3,y3,z3,t3);
       
    //
    //std::cout << "   digit0: (x,y,z,t)=(" << x0 << "," << y0 << "," << z0 << "," << t0 << ") " << std::endl;
    //std::cout << "   digit1: (x,y,z,t)=(" << x1 << "," << y1 << "," << z1 << "," << t1 << ") " << std::endl;
    //std::cout << "   digit2: (x,y,z,t)=(" << x2 << "," << y2 << "," << z2 << "," << t2 << ") " << std::endl;
    //std::cout << "   digit3: (x,y,z,t)=(" << x3 << "," << y3 << "," << z3 << "," << t3 << ") " << std::endl;
    //

    // find common vertex
    WCSimGeometry::FindVertex(x0,y0,z0,t0,
                              x1,y1,z1,t1,
                              x2,y2,z2,t2,
                              x3,y3,z3,t3,
                              fVtxX1,fVtxY1,fVtxZ1,fVtxTime1,
                              fVtxX2,fVtxY2,fVtxZ2,fVtxTime2);

    //
    // std::cout << "   result: (x,y,z,t)=(" << fVtxX1 << "," << fVtxY1 << "," << fVtxZ1 << "," << fVtxTime1 << ") " << std::endl
    //           << "   result: (x,y,z,t)=(" << fVtxX2 << "," << fVtxY2 << "," << fVtxZ2 << "," << fVtxTime2 << ") " << std::endl;
    //

    // add first digit
    if( WCSimGeometry::Instance()->InsideDetector(fVtxX1,fVtxY1,fVtxZ1) ){
      vSeedVtxX.push_back(fVtxX1); 
      vSeedVtxY.push_back(fVtxY1);
      vSeedVtxZ.push_back(fVtxZ1);
      vSeedVtxTime.push_back(fVtxTime1);
    }

    // add second digit
    if( WCSimGeometry::Instance()->InsideDetector(fVtxX2,fVtxY2,fVtxZ2) ){
      vSeedVtxX.push_back(fVtxX2); 
      vSeedVtxY.push_back(fVtxY2);
      vSeedVtxZ.push_back(fVtxZ2);
      vSeedVtxTime.push_back(fVtxTime2);
    }

  }

  return;
}

void WCSimVertexGeometry::ChooseNextQuadruple(Double_t& x0, Double_t& y0, Double_t& z0, Double_t& t0, Double_t& x1, Double_t& y1, Double_t& z1, Double_t& t1, Double_t& x2, Double_t& y2, Double_t& z2, Double_t& t2, Double_t& x3, Double_t& y3, Double_t& z3, Double_t& t3)
{
  this->ChooseNextDigit(x0,y0,z0,t0);
  this->ChooseNextDigit(x1,y1,z1,t1);
  this->ChooseNextDigit(x2,y2,z2,t2);
  this->ChooseNextDigit(x3,y3,z3,t3);

  return;
}

void WCSimVertexGeometry::ChooseNextDigit(Double_t& xpos, Double_t& ypos, Double_t& zpos, Double_t& time)
{
  // default
  xpos=0; ypos=0; zpos=0; time=0;

  // ROOT random number generator
  // Double_t r = gRandom->Rndm();

  // pseudo-random number generator
  Int_t numEntries = vSeedDigitList.size();

  fCounter++;
  if( fCounter>=fNDigits ) fCounter = 0;
  fThisDigit = vSeedDigitList.at(fLastEntry);

  Double_t t0 = 0.5 + fDigitT[fCounter] - fMinTime;
  Double_t q0 = 0.5 + fDigitQ[fCounter];

  Double_t t1 = 0.5 + fDigitT[fThisDigit] - fMinTime;
  Double_t q1 = 0.5 + fDigitQ[fThisDigit];

  Double_t tq = 100.0*(t0*q0+t1*q1);
  Double_t r = tq - TMath::Floor(tq);

  fLastEntry = (Int_t)(r*numEntries);

  // return the new digit
  fThisDigit = vSeedDigitList.at(fLastEntry);
  xpos = fDigitX[fThisDigit];
  ypos = fDigitY[fThisDigit];
  zpos = fDigitZ[fThisDigit];
  time = fDigitT[fThisDigit];

  return;
}



