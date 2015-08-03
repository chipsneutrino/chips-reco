#include "WCSimRingFinder.hh"

#include "WCSimRecoDigit.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoVertex.hh"
#include "WCSimRecoRing.hh"
#include "WCSimHoughTransform.hh"
#include "WCSimHoughTransformArray.hh"

#include "WCSimGeometry.hh"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

ClassImp(WCSimRingFinder)

static WCSimRingFinder* fgRingFinder = 0;

WCSimRingFinder* WCSimRingFinder::Instance()
{
  if( !fgRingFinder ){
    fgRingFinder = new WCSimRingFinder();
  }

  if( !fgRingFinder ){
    assert(fgRingFinder);
  }

  if( fgRingFinder ){

  }

  return fgRingFinder;
}
  
void WCSimRingFinder::UseRecoVertex()
{
  WCSimRingFinder::Instance()->SetUsingRecoVertex();
}

// Set the x co-ordinate bin in Hough space 
void WCSimRingFinder::HoughX( Int_t x )
{
  WCSimRingFinder::Instance()->SetHoughX(x);
}
// Set the y co-ordinate bin in Hough space  
void WCSimRingFinder::HoughY( Int_t y )
{
  WCSimRingFinder::Instance()->SetHoughY(y);
}
 
// Set the number of Hough points
void WCSimRingFinder::HoughPoints( Int_t n )
{
  WCSimRingFinder::Instance()->SetHoughPoints(n);
}

// Set number of bins for the Cherenkov cone angle  
void WCSimRingFinder::ConeAngleBins( Int_t bins )
{
  WCSimRingFinder::Instance()->SetConeAngleBins(bins);
}

// Set range of Cherenkov cone angle
void WCSimRingFinder::ConeAngleMin( Double_t min )
{
  WCSimRingFinder::Instance()->SetConeAngleMin(min);
}

void WCSimRingFinder::ConeAngleMax( Double_t max )
{
  WCSimRingFinder::Instance()->SetConeAngleMax(max);
}

void WCSimRingFinder::PrintParameters()
{
  WCSimRingFinder::Instance()->RunPrintParameters();
}

WCSimRingFinder::WCSimRingFinder()
{
  // use reco vertex
  fUseVertex = 0;
	fUseTSpectrum2 = 0;

  // default hough transform parameters
  fHoughX = 240;          // bins in phi
  fHoughY = 120;           // bins in cos(theta)
  fHoughPoints = 720;     // hough points

  fConeAngleBins = 16;   // hough angle bins
  fConeAngleMin = 34.5;  // hough cone (min)
  fConeAngleMax = 50.5;  // hough cone (max)
  
  // hough transform object
  fHoughTransform = 0;

  // hough transform array
  fHoughTransformArray = 0;
}

// Delete the Hough transform, array and ell entries in the ringlist
WCSimRingFinder::~WCSimRingFinder()
{
  if( fHoughTransform ){
    delete fHoughTransform;
  }

  if( fHoughTransformArray ){
    delete fHoughTransformArray;
  }

}


void WCSimRingFinder::RunPrintParameters()
{
  std::cout << " *** WCSimRingFinder::PrintParameters() *** " << std::endl;

  std::cout << "  Ring Finder Parameters: " << std::endl
            << "   UseVertex = " << fUseVertex << std::endl
            << "   HoughX = " << fHoughX << std::endl
            << "   HoughY = " << fHoughY << std::endl
            << "   ConeAngleBins = " << fConeAngleBins << std::endl
            << "   ConeAngleMin = " << fConeAngleMin << std::endl
            << "   ConeAngleMax = " << fConeAngleMax << std::endl;

  return;
}

void WCSimRingFinder::Reset()
{
	// blank?
  return;
}

std::vector<WCSimRecoRing*>* WCSimRingFinder::Run(WCSimRecoVertex* myVertex)
{
  std::cout << " *** WCSimRingFinder::Run(Vertex) *** " << std::endl;

  // Make new ring, using vertex
  // ===========================
  WCSimRecoRing* myRing = new WCSimRecoRing(myVertex->GetX(),
                                            myVertex->GetY(),
                                            myVertex->GetZ(),   
                                            myVertex->GetDirX(),
                                            myVertex->GetDirY(),
					    					myVertex->GetDirZ(),
					    					myVertex->GetConeAngle(), 
                                            0.0); // height of hough peak

  std::vector<WCSimRecoRing*>* newRings = new std::vector<WCSimRecoRing*>();
  newRings->push_back(myRing);

  // Return Ring List
  // ================
  return newRings;
}

std::vector<WCSimRecoRing*>* WCSimRingFinder::Run(WCSimRecoEvent* myEvent, WCSimRecoVertex* myVertex)
{
  // get filter digit list
  // =====================
  std::vector<WCSimRecoDigit*>* myFilterDigitList = (std::vector<WCSimRecoDigit*>*)(myEvent->GetFilterDigitList());

  // run vertex finder
  // =================
  return (std::vector<WCSimRecoRing*>*)(this->Run(myFilterDigitList,myVertex));
}

std::vector<WCSimRecoRing*>* WCSimRingFinder::Run(std::vector<WCSimRecoDigit*>* myDigitList, WCSimRecoVertex* myVertex)
{
  std::cout << " *** WCSimRingFinder::Run(...) *** " << std::endl;

  // override 
  // ========
  if( fUseVertex ){
    std::cout << "  --- reconstruct ring from vertex --- " << std::endl;
    return this->Run(myVertex);
  }

  // apply Hough Transform
  // =====================
  WCSimHoughTransformArray* myHoughTransformArray = (WCSimHoughTransformArray*)(this->HoughTransformArray(myDigitList,myVertex));
  
  // Find Hough Peak
  // ===============
  Double_t houghVtxX = myVertex->GetX();
  Double_t houghVtxY = myVertex->GetY();
  Double_t houghVtxZ = myVertex->GetZ();

  std::vector<Double_t> houghDirX;
  std::vector<Double_t> houghDirY;
  std::vector<Double_t> houghDirZ;
  std::vector<Double_t> houghAngle;
  std::vector<Double_t> houghHeight;
 
  if(this->UsingTSpectrum2()){ 
    //myHoughTransformArray->FitTSpectrum2(houghDirX,houghDirY,houghDirZ,houghAngle,houghHeight);
    myHoughTransformArray->FitMultiPeaksSmooth(houghDirX,houghDirY,houghDirZ,houghAngle,houghHeight);
  }
	else myHoughTransformArray->FindPeak(houghDirX, houghDirY, houghDirZ, houghAngle, houghHeight);
 	std::cout << "The number of rings passed into the loop is..." << houghDirX.size() << std::endl;

  std::vector<WCSimRecoRing*>* newRings = new std::vector<WCSimRecoRing*>();
	for( int iRings = 0; iRings < (int)houghDirX.size(); iRings++ )
	{
  // Make New Ring
  // =============
		std::cout << "  found new ring: (vx,vy,vz)=(" << houghVtxX << "," << houghVtxY << "," << houghVtxZ << ") " << std::endl
		          << "                  (px,py,pz)=(" << houghDirX[iRings] << "," << houghDirY[iRings] << "," << houghDirZ[iRings] << ") " << std::endl
		          << "                   angle=" << houghAngle[iRings] << " height=" << houghHeight[iRings] << std::endl;
		
		WCSimRecoRing* myRing = new WCSimRecoRing(houghVtxX,houghVtxY,houghVtxZ,
		                                          houghDirX[iRings],houghDirY[iRings],houghDirZ[iRings],
		                                          houghAngle[iRings],houghHeight[iRings]);
		newRings->push_back(myRing);
	}
  return newRings;
}

WCSimHoughTransform* WCSimRingFinder::HoughTransform(WCSimRecoEvent* myEvent, WCSimRecoVertex* myVertex, Double_t myAngle)
{
  // get filter digit list
  // =====================
  std::vector<WCSimRecoDigit*>* myFilterDigitList = (std::vector<WCSimRecoDigit*>*)(myEvent->GetFilterDigitList());

  // run Hough Transform
  // ===================
  return (WCSimHoughTransform*)(this->HoughTransform(myFilterDigitList,myVertex,myAngle));
}

WCSimHoughTransform* WCSimRingFinder::HoughTransform(std::vector<WCSimRecoDigit*>* myFilterDigitList, WCSimRecoVertex* myVertex, Double_t myAngle)
{
  std::cout << " *** WCSimRingFinder::HoughTransform(...) *** " << std::endl;
  std::cout << "  calculate Hough Transform for cone angle: " << myAngle << std::endl;  
 
  // make new Hough Transform (if necessary)
  // =======================================
  if( fHoughTransform==0 ){
    fHoughTransform = new WCSimHoughTransform(fHoughX,fHoughY);
  }

  // reset Hough Transform
  // =====================
  fHoughTransform->Reset();

  // perform Hough Transform
  // =======================
  for( UInt_t idigit=0; idigit<myFilterDigitList->size(); idigit++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(myFilterDigitList->at(idigit));

    Double_t x = myDigit->GetX();
    Double_t y = myDigit->GetY();
    Double_t z = myDigit->GetZ();

    Double_t vx = myVertex->GetX();
    Double_t vy = myVertex->GetY();
    Double_t vz = myVertex->GetZ();

    Double_t omega = 0.0;
    Double_t rx = 0.0;
    Double_t ry = 0.0;
    Double_t rz = 0.0;
    Double_t nx = 0.0;
    Double_t ny = 0.0;
    Double_t nz = 0.0;
    Double_t r = 0.0;

    Double_t weight = 1.0;

    for( Int_t n=0; n<fHoughPoints; n++ ){
      omega = 360.0*((double)n/(double)fHoughPoints);

      WCSimGeometry::FindCircle(x, y, z,
                                vx, vy, vz,
                                myAngle, omega,
                                rx, ry, rz,
                                nx, ny, nz, r);

      fHoughTransform->Fill(nx,ny,nz,weight);
    }
  }

  // return result of Hough Transform
  // ================================
  return fHoughTransform;
}

WCSimHoughTransformArray* WCSimRingFinder::HoughTransformArray(WCSimRecoEvent* myEvent, WCSimRecoVertex* myVertex)
{
  // get filter digit list
  // =====================
  std::vector<WCSimRecoDigit*>* myFilterDigitList = (std::vector<WCSimRecoDigit*>*)(myEvent->GetFilterDigitList());

  // run Hough Transform
  // ===================
  return (WCSimHoughTransformArray*)(this->HoughTransformArray(myFilterDigitList,myVertex));
}

WCSimHoughTransformArray* WCSimRingFinder::HoughTransformArray(std::vector<WCSimRecoDigit*>* myFilterDigitList, WCSimRecoVertex* myVertex)
{
  std::cout << " *** WCSimRingFinder::HoughTransformArray(...) *** " << std::endl;
  
  // make new Hough Transform array (if necessary)
  // =============================================
  if( fHoughTransformArray==0 ){
    fHoughTransformArray = new WCSimHoughTransformArray( fConeAngleBins,fConeAngleMin,fConeAngleMax,
                                                         fHoughX,fHoughY);
  }

  // reset Hough Transform array
  // ===========================
  fHoughTransformArray->Reset();

  // perform Hough Transform
  // =======================
  for( UInt_t idigit=0; idigit<myFilterDigitList->size(); idigit++ ){
    WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(myFilterDigitList->at(idigit));

    Double_t x = myDigit->GetX();
    Double_t y = myDigit->GetY();
    Double_t z = myDigit->GetZ();

    Double_t vx = myVertex->GetX();
    Double_t vy = myVertex->GetY();
    Double_t vz = myVertex->GetZ();

    Double_t omega = 0.0;
    Double_t rx = 0.0;
    Double_t ry = 0.0;
    Double_t rz = 0.0;
    Double_t nx = 0.0;
    Double_t ny = 0.0;
    Double_t nz = 0.0;
    Double_t r = 0.0;

    Double_t weight = 1.0;

	// Loops over omega to step round the rim of cones, and over nAngle to repeat over cone angles
    for( Int_t n=0; n<fHoughPoints; n++ ){
      omega = 360.0*((double)n/(double)fHoughPoints);

      for( Int_t nAngle=0; nAngle<fHoughTransformArray->GetConeAngleBins(); nAngle++ ){
        Double_t myAngle = fHoughTransformArray->GetConeAngle(nAngle);
	WCSimHoughTransform* myHoughTransform = (WCSimHoughTransform*)(fHoughTransformArray->GetHoughTransform(nAngle));
        
        WCSimGeometry::FindCircle(x, y, z,
                                  vx, vy, vz,
                                  myAngle, omega,
                                  rx, ry, rz,
                                  nx, ny, nz, r);

        myHoughTransform->Fill(nx,ny,nz,weight);
      }
    }
  }

  // return result of Hough Transform
  // ================================
  return fHoughTransformArray;
}




