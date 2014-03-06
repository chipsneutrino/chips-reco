#include "WCSimHoughTransformArray.hh"
#include "WCSimHoughTransform.hh"

#include "TMath.h"
#include "TSpectrum2.h"
#include "TSpectrum3.h"

#include <iostream>
#include <vector>

ClassImp(WCSimHoughTransformArray)

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
						<< "  phi      (bins, min, max) =(" << fHoughX << ",0,360)" << std::endl
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

 
 
