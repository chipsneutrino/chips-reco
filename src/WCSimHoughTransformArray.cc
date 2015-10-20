#include "TCanvas.h"
#include "TStyle.h"

#include "WCSimHoughTransformArray.hh"
#include "WCSimHoughTransform.hh"

#include "TMath.h"
#include "TH2D.h"
#include "TVector2.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TSpectrum3.h"

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

ClassImp(WCSimHoughTransformArray)

bool PairSort(const std::pair<double,TVector2> &a, const std::pair<double,TVector2> &b);

// Constructor
//============
WCSimHoughTransformArray::WCSimHoughTransformArray(Int_t coneAngleBins, Double_t coneAngleMin, Double_t coneAngleMax, Int_t phiBins, Int_t cosThetaBins)
{
  this->BuildArray(coneAngleBins,coneAngleMin,coneAngleMax,
                   phiBins,cosThetaBins);
}

// Destructor
//============
WCSimHoughTransformArray::~WCSimHoughTransformArray()
{
  this->DeleteArray();
}

// Steps through possible Cherenkov cone angles, performs a Hough transform assuming each one and puts the resulting transform into an array 
//============================================================================================================================================
void WCSimHoughTransformArray::BuildArray(Int_t coneAngleBins, Double_t coneAngleMin, Double_t coneAngleMax, Int_t phiBins, Int_t cosThetaBins)
{
  std::cout << " *** WCSimHoughTransformArray::BuildArray() *** " << std::endl;

  // set parameters
  fHoughX = phiBins;
  fHoughY = cosThetaBins;
  fConeAngleBins = coneAngleBins;
  fConeAngleMin = coneAngleMin;
  fConeAngleMax = coneAngleMax;

  std::cout << "  building Hough Transform Array: Cone angle (bins,min,max)=(" << fConeAngleBins << "," << fConeAngleMin << "," << fConeAngleMax << ") " << std::endl
						<< "  phi      (bins, min, max) =(" << fHoughX << ",-180,180)" << std::endl
						<< "  costheta (bins, min, max) =(" << fHoughY << "-1,1)" << std::endl; 

  // clear current array
  this->DeleteArray();

  // create new array
  for( Int_t n=0; n<fConeAngleBins; n++ ){
    WCSimHoughTransform* myHoughTransform = new WCSimHoughTransform(fHoughX,fHoughY);
    fHoughArray.push_back(myHoughTransform);
  }

  return;
}

void WCSimHoughTransformArray::FindPeak(			std::vector<Double_t> &houghDirX, std::vector<Double_t> &houghDirY,
																							std::vector<Double_t> &houghDirZ, std::vector<Double_t> &houghAngle, 
																							std::vector<Double_t> &houghHeight )
{
	Double_t x,y,z,angle,height;
	this->FindPeak(x,y,z,angle,height);
	houghDirX.push_back(x);
	houghDirY.push_back(y);
	houghDirZ.push_back(z);
	houghAngle.push_back(angle);
	houghHeight.push_back(height);

	return;

}

// Runs the Hough peak finder for each Cherenkov angle and returns the biggest
//=============================================================================
void WCSimHoughTransformArray::FindPeak(Int_t& bin, Double_t& height)
{
  Int_t bestBin = -1;  
  Double_t bestHeight = 0.0;

  Double_t myHeight = 0.0;  
  
  for( Int_t iBin=0; iBin<this->GetBins(); iBin++ ){
    WCSimHoughTransform* myHoughTransform = this->GetHoughTransform(iBin);
    myHoughTransform->FindPeak(myHeight);

    if( myHeight>bestHeight ){
      bestBin = iBin;
      bestHeight = myHeight;
    }
  }

  bin = bestBin;
  height = bestHeight;

  return;
}

// Runs the Hough peak finder for each cone angle and returns the best one's height, phi, costheta and cone angle
//===============================================================================================================
void WCSimHoughTransformArray::FindPeak(Double_t& phi, Double_t& costheta, Double_t& angle, Double_t& height)
{
  Int_t bestBin = -1;  
  Double_t bestPhi = 0.0;
  Double_t bestCosTheta = 0.0;
  Double_t bestHeight = 0.0;
    
  Double_t myPhi = 0.0;
  Double_t myCosTheta = 0.0;
  Double_t myHeight = 0.0;  

  for( Int_t iBin=0; iBin<this->GetBins(); iBin++ ){
    WCSimHoughTransform* myHoughTransform = this->GetHoughTransform(iBin);
    myHoughTransform->FindPeak(myPhi,myCosTheta,myHeight);

    if( myHeight>bestHeight ){
      bestBin = iBin;
      bestPhi = myPhi;
      bestCosTheta = myCosTheta;
      bestHeight = myHeight;
    }
  }  

  angle = this->GetAngle(bestBin);
  phi = bestPhi;
  costheta = bestCosTheta;
  height = bestHeight;

  return;
}

 // Runs the Hough peak finder for each cone angle and returns the best one's height, cone angle and x,y,z direction
//===================================================================================================================
void WCSimHoughTransformArray::FindPeak(Double_t& hx, Double_t& hy, Double_t& hz, Double_t& angle, Double_t& height)
{

  Int_t bestBin = -1;  
  Double_t bestDirX = 0.0;
  Double_t bestDirY = 0.0;
  Double_t bestDirZ = 0.0;
  Double_t bestHeight = 0.0;
    
  Double_t myDirX = 0.0;
  Double_t myDirY = 0.0;
  Double_t myDirZ = 0.0;
  Double_t myHeight = 0.0;

  // Loop to find the maximum peak height
  for( Int_t iBin=0; iBin<this->GetBins(); iBin++ ){
    WCSimHoughTransform* myHoughTransform = this->GetHoughTransform(iBin);
    myHoughTransform->FindPeak(myDirX,myDirY,myDirZ,myHeight);

    if( myHeight>bestHeight ){
      bestBin = iBin;
      bestDirX = myDirX;
      bestDirY = myDirY;
      bestDirZ = myDirZ;
      bestHeight = myHeight;
    }
  }

  angle = this->GetAngle(bestBin);
  hx = bestDirX;
  hy = bestDirY;
  hz = bestDirZ;
  height = bestHeight;

  return;
}

// If myAngle is an allowed cone angle, returns the cone angle bin containing it
//===============================================================================
Int_t WCSimHoughTransformArray::FindBin(Double_t myAngle)
{
  if( myAngle>=fConeAngleMin && myAngle<fConeAngleMax ){
    return (Int_t)(fConeAngleBins*(myAngle-fConeAngleMin)/(fConeAngleMax-fConeAngleMin));
  }
  else{
    return -1;
  }
}

// Returns the cone angle at the centre of a given bin
//=====================================================
Double_t WCSimHoughTransformArray::GetAngle(Int_t myBin)
{
  return this->GetConeAngle(myBin);
}

// Returns the cone angle at the centre of a given bin, unless we don't have any bins
// in which case, returns -45.0
//====================================================================================
Double_t WCSimHoughTransformArray::GetConeAngle(Int_t myBin)
{
  if( fConeAngleBins>0 ){
    return fConeAngleMin + ((double)(myBin+0.5)/(double)fConeAngleBins)*(fConeAngleMax-fConeAngleMin); 
	// Note: doesn't check myBin <= fConeAngleBins
  }
  else{
    return -45.0;
  }
}

// Return a pointer to the Hough transform at entry nAngle in the array
//======================================================================
WCSimHoughTransform* WCSimHoughTransformArray::GetHoughTransform(Int_t nAngle)
{
  if( nAngle>=0 && nAngle<fConeAngleBins ){
    return (WCSimHoughTransform*)(fHoughArray.at(nAngle));
  }
  else return 0;
}

// Delete all the members of the array
//=====================================
void WCSimHoughTransformArray::DeleteArray()
{
  for( UInt_t n=0; n<fHoughArray.size(); n++ ){
    delete (WCSimHoughTransform*)(fHoughArray.at(n));
  }

  fHoughArray.clear();

  return;
}

// Write out the Cherenkov angle bin numbers and central values
//==============================================================
void WCSimHoughTransformArray::PrintArray()
{
  std::cout << " *** WCSimHoughTransformArray::PrintArray() *** " << std::endl;

  for( Int_t n=0; n<fConeAngleBins; n++ ){
    std::cout << "  [" << n << "] coneAngle=" << this->GetConeAngle(n) << std::endl; 
  }

  return;
}

// Reset each hough transform
//============================
void WCSimHoughTransformArray::Reset()
{
  for( UInt_t n=0; n<fHoughArray.size(); n++ ){
    WCSimHoughTransform* myHoughTransform = (WCSimHoughTransform*)(fHoughArray.at(n));
    myHoughTransform->Reset();
  }
}

// Get maximum in each angle bin as a function of cone angle
//===========================================================
void WCSimHoughTransformArray::AngleMaxima( std::vector<std::vector<Float_t> > &fHoughArrayMaxAngle )
{
	
	std::cout << " *** WCSimHoughTransformArray::AngleMaxima() *** " << std::endl;

	std::cout << "Setting fHoughArrayMaxAngle" << std::endl;	

	// Find the maximum for each (x,y) bin over all cone angles and fill the array with it	
	Double_t maximum;
	Int_t maxAngleBin;
	std::vector<Float_t> heightAngle;
	for(int xbin = 0; xbin < fHoughX; xbin++){
		for(int ybin = 0; ybin < fHoughY; ybin++){
			heightAngle.clear();
			maximum = 0.0;
			maxAngleBin = 0;
			for( int conebin = 0; conebin < fConeAngleBins; conebin++){
				WCSimHoughTransform * myHoughTransform = (WCSimHoughTransform *)fHoughArray[conebin];
				
				if( (myHoughTransform->GetEntry(xbin, ybin)) > maximum ){
					maximum = (myHoughTransform->GetEntry(xbin, ybin));
					maxAngleBin = conebin;
				}
			}
			
			heightAngle.push_back(maximum);
			heightAngle.push_back(this->GetAngle(maxAngleBin));
			fHoughArrayMaxAngle.push_back(heightAngle);
		}
	}
	std::cout << "Returning angle maxima" << std::endl;
	return;
}

// Use a TSpectrum3 to find multiple peaks in the output of AngleMaxima
//======================================================================
void WCSimHoughTransformArray::FitTSpectrum3( std::vector<Double_t> &houghDirX, std::vector<Double_t> &houghDirY,
																							std::vector<Double_t> &houghDirZ, std::vector<Double_t> &houghAngle, 
																							std::vector<Double_t> &houghHeight )
{

  std::cout << " *** WCSimHoughTransformArray::FitTSpectrum3() *** " << std::endl;

	// Convert fHoughArrayMaxAngle into the right form to pass to TSpectrum2
	std::cout << "Converting to triple pointer" << std::endl;
	float *** source = new float **[fHoughX];
	float *** dest = new float **[fHoughX];     

	for(int xbin = 0; xbin < fHoughX; xbin++ ){
		source[xbin] = new float * [fHoughY];
		dest[xbin] = new float * [fHoughY];
		for(int ybin=0; ybin < fHoughY; ybin++ ){
	  	source[xbin][ybin]=new float [fConeAngleBins];
	  	dest[xbin][ybin]=new float [fConeAngleBins];
		}
	}          
	std::cout << "Initializing TSpectrum2" << std::endl;	
	TSpectrum3 *search = new TSpectrum3();
	
	std::cout << "Populating source" << std::endl;
	for( int conebin = 0; conebin < fConeAngleBins; conebin++){
		WCSimHoughTransform * myHoughTransform = (WCSimHoughTransform *)fHoughArray[conebin];
		for (int xbin = 0; xbin < fHoughX; xbin++){
		  for (int ybin = 0; ybin < fHoughY; ybin++){
		       source[xbin][ybin][conebin]= (Float_t)myHoughTransform->GetEntry(xbin, ybin);
			}
		}
	}

	// Now deconvolve and find the peaks
	Int_t nfound;
	const float *** source2 = (const float ***) source;
	std::cout << "Running peak finder" << std::endl;
  nfound = search->SearchHighRes(source2, dest, fHoughX, fHoughY, fConeAngleBins, 3, 35, kTRUE, 2, kFALSE, 5);   
	if(nfound > 0)
	{
		std::cout << "Found " << nfound << " peaks" << std::endl;
	
	  Float_t *PosX = new Float_t[nfound];        
	  Float_t *PosY = new Float_t[nfound];
		Float_t *PosAngle = new Float_t[nfound];
		PosX = search->GetPositionX();
		PosY = search->GetPositionY();
		PosAngle = search->GetPositionZ();
	
		for(int iRings = 0; iRings < nfound; iRings++){	
			std::cout << PosX[iRings] << "   " << PosY[iRings] << std::endl;
	
			Int_t binX = TMath::Nint(PosX[iRings]);
			Int_t binY = TMath::Nint(PosY[iRings]);
			Int_t binAngle = TMath::Nint(PosAngle[iRings]);
			WCSimHoughTransform * ringHoughTransform = (WCSimHoughTransform *)this->GetHoughTransform(binAngle);
			double Angle = this->GetConeAngle(binAngle);
			double Height = ringHoughTransform->GetEntry(binX,binY);

		  Float_t phiradians = ringHoughTransform->GetX(binX);  // This comes from PosX but we need to convert bins -> phi
		  Float_t costheta = ringHoughTransform->GetY(binY);
			Float_t sintheta = sqrt(1.0-costheta*costheta);
			
			std::cout << "Phiradians = " << phiradians << "   sintheta = " << costheta << std::endl;
		  double DirX = sintheta * cos(phiradians);
			double DirY = sintheta * sin(phiradians);
			double DirZ = costheta;
	
			std::cout << DirX << "  " << DirY << "  " << DirZ << "   " << std::endl;
			
			houghDirX.push_back(DirX);
	 		houghDirY.push_back(DirY);
	  	houghDirZ.push_back(DirZ);
			houghAngle.push_back(Angle);
			houghHeight.push_back(Height);
		}
	}

	return;
}





// Use a TSpectrum2 to find multiple peaks in the output of AngleMaxima
//======================================================================
void WCSimHoughTransformArray::FitTSpectrum2( std::vector<Double_t> &houghDirX, std::vector<Double_t> &houghDirY,
																							std::vector<Double_t> &houghDirZ, std::vector<Double_t> &houghAngle, 
																							std::vector<Double_t> &houghHeight )
{

  std::cout << " *** WCSimHoughTransformArray::FitTSpectrum2() *** " << std::endl;

	std::vector<std::vector<Float_t> > fHoughArrayMaxAngle;	
	this->AngleMaxima(fHoughArrayMaxAngle);
 
	// Convert fHoughArrayMaxAngle into the right form to pass to TSpectrum2
	std::cout << "Converting to double pointer" << std::endl;
	float ** source = new float *[fHoughX];
	float ** dest = new float *[fHoughX];     

	for(int xbin=0; xbin < fHoughX; xbin++ ){
	   source[xbin]=new float [fHoughY];
	   dest[xbin]=new float [fHoughY];
	}          
	std::cout << "Initializing TSpectrum2" << std::endl;	
	TSpectrum2 *search = new TSpectrum2();
	
	std::cout << "Populating source" << std::endl;
	for (int xbin = 0; xbin < fHoughX; xbin++){
	  for (int ybin = 0; ybin < fHoughY; ybin++){
	       source[xbin][ybin]= (Float_t)fHoughArrayMaxAngle[ xbin * fHoughY + ybin ][0];
		}
	}

	// Now deconvolve and find the peaks
	Int_t nfound;
	std::cout << "Running peak finder" << std::endl;
  nfound = search->SearchHighRes(source, dest, fHoughX, fHoughY, 4.0, 20, kTRUE, 5, kFALSE, 5);   
	if(nfound > 0)
	{
		std::cout << "Found " << nfound << " peaks" << std::endl;
	
	  Float_t *PosX = new Float_t[nfound];        
	  Float_t *PosY = new Float_t[nfound];
		PosX = search->GetPositionX();
		PosY = search->GetPositionY();
	
		for(int iRings = 0; iRings < nfound; iRings++){	
			std::cout << PosX[iRings] << "   " << PosY[iRings] << std::endl;
	
			Int_t binX = TMath::Nint(PosX[iRings]);
			Int_t binY = TMath::Nint(PosY[iRings]);
			Float_t Angle = fHoughArrayMaxAngle[binX*fHoughY + binY][1];
			Int_t binAngle = this->FindBin(Angle);

			WCSimHoughTransform * ringHoughTransform = (WCSimHoughTransform *)this->GetHoughTransform(binAngle);
			double Height = ringHoughTransform->GetEntry(binX,binY);

		  Float_t phiradians = TMath::Pi() / 180.0 * ringHoughTransform->GetX(binX);  // This comes from PosX but we need to convert bins -> phi
		  Float_t costheta = ringHoughTransform->GetY(binY);
			Float_t sintheta = sqrt(1.0-costheta*costheta);
			
			std::cout << "Phiradians = " << phiradians << "   sintheta = " << costheta << std::endl;
		  double DirX = sintheta * cos(phiradians);
			double DirY = sintheta * sin(phiradians);
			double DirZ = costheta;
	
			std::cout << DirX << "  " << DirY << "  " << DirZ << "   " << std::endl;
			
			houghDirX.push_back(DirX);
	 		houghDirY.push_back(DirY);
	  	houghDirZ.push_back(DirZ);
			houghAngle.push_back(Angle);
			houghHeight.push_back(Height);
		}
/*
		for(int iRings = 0; iRings < nfound; iRings++){	
			std::cout << PosX[iRings] << "   " << PosY[iRings] << std::endl;
	
		  Int_t binX = (fHoughArray.at(0))->GetBinX(PosX[iRings]);
		  Int_t binY = (fHoughArray.at(0))->GetBinY(PosY[iRings]); 
			double Angle = fHoughArrayMaxAngle[binX*fHoughY + binY][1];
			double Height = fHoughArrayMaxAngle[binX*fHoughY + binY][0];

		  Float_t phiradians = (-TMath::Pi() + (PosX[iRings])*TMath::Pi()*2.0/fHoughX);  // This comes from PosX but we need to convert bins -> phi
		  Float_t costheta = -1.0 + PosY[iRings] * 2.0/fHoughY;
			Float_t sintheta = sqrt(1.0-costheta*costheta);
			
			std::cout << "Phiradians = " << phiradians << "   sintheta = " << costheta << std::endl;
		  double DirX = sintheta * cos(phiradians);
			double DirY = sintheta * sin(phiradians);
			double DirZ = costheta;
	
			std::cout << DirX << "  " << DirY << "  " << DirZ << "   " << std::endl;
			
			houghDirX.push_back(DirX);
	 		houghDirY.push_back(DirY);
	  	houghDirZ.push_back(DirZ);
			houghAngle.push_back(Angle);
			houghHeight.push_back(Height);
		}*/
	}

	return;
}

void WCSimHoughTransformArray::FitMultiPeaksSmooth(std::vector<Double_t> &houghDirX, std::vector<Double_t> &houghDirY,
                                                   std::vector<Double_t> &houghDirZ, Double_t seedDirX, Double_t seedDirY, 
                                                   Double_t seedDirZ, std::vector<Double_t> &houghAngle,
                                                   std::vector<Double_t> &houghHeight){
  // Transform for the expected cone angle.
  WCSimHoughTransform* houghTrans = this->GetHoughTransform(this->FindBin(42));

  // If we are at the edge of the histogram space then rotate it.
  double seedPhi = TMath::ATan2(seedDirY,seedDirX);
  double houghRotationPhi = 0.0;
  if(fabs(seedPhi) > 0.75 * TMath::Pi()){
    houghRotationPhi = 180.0;
  }

  // Get the histogram corresponding to the expected cone angle.
  TH2D* houghSpace = houghTrans->GetRotatedTH2D("houghSpace",houghRotationPhi);

  // Smooth slightly to help prevent noise problems
  houghSpace->Smooth(1);

  // Search for peaks in the histogram.
  TSpectrum2 *peakFinder = new TSpectrum2();
  peakFinder->Search(houghSpace,4,"col");

  // Search for ridges around cosTheta = +/- 1
  // These aren't peaks because phi is arbitrary for upward/downward particles
  TH1D * houghProjection1D = new TH1D("houghProjection1D",";cos#theta;height",
                                      houghSpace->GetYaxis()->GetNbins(),
                                      houghSpace->GetYaxis()->GetXmin(),
                                      houghSpace->GetYaxis()->GetXmax());
  for(int yBin = 1; yBin <= houghSpace->GetNbinsY(); ++yBin)
  {
    double total = 0.0;
    double nonZeroBins = 0;
    for(int xBin = 1; xBin <= houghSpace->GetNbinsX(); ++xBin)
    {
      total += houghSpace->GetBinContent(xBin, yBin);
      if(houghSpace->GetBinContent(xBin, yBin) > 0)
      {
        nonZeroBins++;
      }
    }
    if(nonZeroBins)
    {
      houghProjection1D->SetBinContent(yBin, total/nonZeroBins);
    }
  }
  TSpectrum *peakFinder1D = new TSpectrum();
  peakFinder1D->Search(houghProjection1D, 4, "col");

  std::vector<std::pair<double,TVector2> > peakListToSort;
  float *phiArray = peakFinder->GetPositionX();
  float *thetaArray = peakFinder->GetPositionY();
  for(int i = 0; i < peakFinder->GetNPeaks(); ++i){

    TVector2 tempVec(phiArray[i] - houghRotationPhi,thetaArray[i]);
    std::pair<double,TVector2> tempPair(houghSpace->Interpolate(phiArray[i],thetaArray[i]),tempVec);
    peakListToSort.push_back(tempPair);
  }
  
  float * thetaArray1D = peakFinder1D->GetPositionX();
  for(int i = 0; i < peakFinder1D->GetNPeaks(); ++i){
    double max = 0.0;
    double maxPhi = 1;
    double cosTheta = thetaArray1D[i];
    if(fabs(cosTheta) < 0.95){ continue; }
    for(int xBin = 1; xBin < houghSpace->GetNbinsX(); ++xBin)
    {
      double tmp = houghSpace->GetBinContent(xBin, houghSpace->GetYaxis()->FindBin(thetaArray1D[i]));
      if(tmp > max)
      { 
        max = tmp; 
        maxPhi = houghSpace->GetXaxis()->GetBinCenter(xBin);  
      }
    }
    TVector2 tempVec(maxPhi - houghRotationPhi,thetaArray1D[i]);
    std::pair<double,TVector2> tempPair(houghSpace->Interpolate(maxPhi,thetaArray1D[i]),tempVec);
    peakListToSort.push_back(tempPair);
  }
  // Check to make sure we didn't miss a peak at theta = -1.
  if(houghProjection1D->GetMaximumBin() == 1){
    double thetaVal = houghProjection1D->GetXaxis()->GetBinCenter(1);
    double maxPhi = 1;
    double max = 0.0;
    for(int xBin = 1; xBin < houghSpace->GetNbinsX(); ++xBin)
    {
      double tmp = houghSpace->GetBinContent(xBin, houghSpace->GetYaxis()->FindBin(thetaVal));
      if(tmp > max)
      { 
        max = tmp; 
        maxPhi = houghSpace->GetXaxis()->GetBinCenter(xBin);  
      }
    }
    TVector2 tempVec(maxPhi - houghRotationPhi,thetaVal);
    std::pair<double,TVector2> tempPair(houghSpace->Interpolate(maxPhi,thetaVal),tempVec);
    peakListToSort.push_back(tempPair);
  }

  std::sort(peakListToSort.begin(),peakListToSort.end(),PairSort);

  double maxHeight = -999;
  if(peakListToSort.size() != 0){
    maxHeight = peakListToSort[0].first;
  }

  for(unsigned int v = 0; v < peakListToSort.size(); ++v){

    // Only consider fairly decent sized peaks
    if(peakListToSort[v].first < (maxHeight / 2.0)) continue;

    std::cout << "Peak height = " << peakListToSort[v].first << std::endl;
    Float_t phiradians = TMath::Pi() / 180.0 * (peakListToSort[v].second).X();
    Float_t costheta = (peakListToSort[v].second).Y();
    Float_t sintheta = sqrt(1.0-costheta*costheta);

    double dirX = sintheta * cos(phiradians);
    double dirY = sintheta * sin(phiradians);
    double dirZ = costheta;
    std::cout << "Phiradians = " << phiradians << "   costheta = " << costheta << "    (" << dirX << ", " << dirY << ", " << dirZ << ") " << std::endl;

    houghDirX.push_back(dirX);
    houghDirY.push_back(dirY);
    houghDirZ.push_back(dirZ);
    houghAngle.push_back(42.0);
    houghHeight.push_back(peakListToSort[v].first);
  }

  // Tidy up
  delete houghSpace;
  delete peakFinder;
  delete houghProjection1D;
  houghSpace = 0x0;
  peakFinder = 0x0;

}

bool PairSort(const std::pair<double,TVector2> &a, const std::pair<double,TVector2> &b){
  return a.first > b.first;
}
 
