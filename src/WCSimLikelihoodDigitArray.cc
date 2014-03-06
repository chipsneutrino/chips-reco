#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimGeometry.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimInterface.hh"
#include "TCollection.h"
#include <cstdlib>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodDigitArray)
#endif

///////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodDigitArray::WCSimLikelihoodDigitArray( )
{   
	return;
}

WCSimLikelihoodDigitArray::WCSimLikelihoodDigitArray( WCSimRootEvent * myEvent )
{
  // Check the number of PMTs and create a TClonesArray with one entry per PMT
  Int_t numPMT = ((WCSimGeometry::Instance())->GetNumPMTs());	 
  WCSimGeometry::PrintGeometry();
  std::cout << "There are " << numPMT << " PMTs" << std::endl;
  std::cout << "here 1" << std::endl;
	fLikelihoodDigitArray = new TClonesArray("WCSimLikelihoodDigit",numPMT);
  std::cout << "here 2" << std::endl;
  TClonesArray &digitArray = *fLikelihoodDigitArray;
  std::cout << "here 3" << std::endl;

  // Get the geometry type and maximum x, y and z coordinates of the detector volume
  fGeomType = -1;
  if((WCSimGeometry::Instance())->IsCylinder()){ fGeomType = 0; }
  else if((WCSimGeometry::Instance())->IsMailBox()){ fGeomType = 1;}
  
  fExtent[0] = 0.;
  fExtent[1] = 0.;
  fExtent[2] = 0.;
  if( this->IsCylinder() ) 
  {
    fExtent[0] = WCSimGeometry::Instance()->GetCylRadius();
    fExtent[1] = fExtent[0];
    fExtent[2] = 0.5 * WCSimGeometry::Instance()->GetCylLength();
  }
  else if(this->IsMailBox())
  {
    fExtent[0] = WCSimGeometry::Instance()->GetMailBoxX();
    fExtent[1] = WCSimGeometry::Instance()->GetMailBoxY();
    fExtent[2] = WCSimGeometry::Instance()->GetMailBoxZ();
  }
  else
  {
    std::cerr << "Error: WCSimLikelihoodDigitArray() - Couldn't determine geometry type" << std::endl;
    exit(EXIT_FAILURE);
  }

		
  // Initialize an array saying that none of the PMTs were hit    
	std::vector<bool> hitPMT;
	for(int iPMT = 0; iPMT < numPMT; ++iPMT)
	{
		hitPMT.push_back(kFALSE);
	}

  // For now, find the trigger with the largest number of hits
  WCSimRootTrigger * myTrigger = WCSimInterface::FilterTrigger( myEvent );

	// Use this trigger to flag the tubes that were hit and fill the ClonesArray with
  // the corresponding WCSimLikelihoodDigit
  TClonesArray * fCherenkovDigiHits = myTrigger->GetCherenkovDigiHits();
  std::cout << "here 4" << std::endl;
  TIter digiHitItr(fCherenkovDigiHits);
	while( WCSimRootCherenkovDigiHit * myChDigiHit = (WCSimRootCherenkovDigiHit *) digiHitItr.Next())
  {	
	  Int_t tubeID = myChDigiHit->GetTubeId();
    Int_t arrayIndex = tubeID - 1;
    
    if(tubeID < 10) std::cout << tubeID << std::endl;
	  hitPMT[arrayIndex] = kTRUE; 
	  new(digitArray[arrayIndex]) WCSimLikelihoodDigit(myChDigiHit);
	}
  std::cout << "here 55" << std::endl;
	
  // Now fill all the unhit PMTs with an empty (q = 0) hit
	for(int arrayIndex = 0; arrayIndex < numPMT; ++arrayIndex)
	{
    if(arrayIndex < 10) std::cout << "arrayIndex = " << arrayIndex << std::endl;
		if(hitPMT[arrayIndex] == kFALSE)
		{
      Int_t tubeID = arrayIndex + 1;
      WCSimRootCherenkovDigiHit * emptyHit = new WCSimRootCherenkovDigiHit(0.,0.,tubeID);
			new (digitArray[arrayIndex]) WCSimLikelihoodDigit(emptyHit);
			delete emptyHit;
    }
	}
  
  // Also save the array size 
  fNLikelihoodDigits = numPMT;

}

// Use hits before digitizing them to help tune the charge likelihood
WCSimLikelihoodDigitArray::WCSimLikelihoodDigitArray( WCSimRootEvent * myRootEvent, Bool_t useUndigitized)
{
	// Check the number of PMTs and create a TClonesArray with one entry per PMT
    Int_t numPMT = ((WCSimGeometry::Instance())->GetNumPMTs());  // WCSimGeometry adds 1 to the true number to prevent overflow errors
  	WCSimGeometry::PrintGeometry();
  	std::cout << "There are " << numPMT << " PMTs" << std::endl;
	fLikelihoodDigitArray = new TClonesArray("WCSimLikelihoodDigit",numPMT);
  std::cout << "here 5" << std::endl;
  	TClonesArray &digitArray = *fLikelihoodDigitArray;
  std::cout << "here 6" << std::endl;

  	// Get the geometry type and maximum x, y and z coordinates of the detector volume
  	fGeomType = -1;
  	if((WCSimGeometry::Instance())->IsCylinder()){ fGeomType = 0; }
  	else if((WCSimGeometry::Instance())->IsMailBox()){ fGeomType = 1;}
  
  	fExtent[0] = 0.;
  	fExtent[1] = 0.;
  	fExtent[2] = 0.;
  	if( this->IsCylinder() ) 
  	{
  	  fExtent[0] = WCSimGeometry::Instance()->GetCylRadius();
    	fExtent[1] = fExtent[0];
    	fExtent[2] = 0.5 * WCSimGeometry::Instance()->GetCylLength();
  	}
  	else if(this->IsMailBox())
  	{
    	fExtent[0] = WCSimGeometry::Instance()->GetMailBoxX();
    	fExtent[1] = WCSimGeometry::Instance()->GetMailBoxY();
    	fExtent[2] = WCSimGeometry::Instance()->GetMailBoxZ();
  	}
  	else
  	{
    	std::cerr << "Error: WCSimLikelihoodDigitArray() - Couldn't determine geometry type" << std::endl;
    	exit(EXIT_FAILURE);
  	}

		
  	// Initialize an array saying that none of the PMTs were hit    
	std::vector<bool> hitPMT;
	for(int iPMT = 0; iPMT < numPMT; ++iPMT)
	{
		hitPMT.push_back(kFALSE);
	}

  	// For now, find the trigger with the largest number of hits
  	WCSimRootTrigger * myTrigger = WCSimInterface::FilterTrigger( myRootEvent );
    std::cout << "Got here" << std::endl;

	// Use this trigger to flag the tubes that were hit and fill the ClonesArray with
  	// the corresponding WCSimLikelihoodDigit
  	TClonesArray * fCherenkovDigiHits = myTrigger->GetCherenkovHits();
 	TClonesArray * fCherenkovHitTimes = myTrigger->GetCherenkovHitTimes();
 	
  	TIter undigiHitItr(fCherenkovDigiHits);
	while( WCSimRootCherenkovHit * myChHit = (WCSimRootCherenkovHit *) undigiHitItr.Next())
  	{	
	  	Int_t tubeID = myChHit->GetTubeID();
      Int_t arrayIndex = tubeID-1;
    	Double_t Q = myChHit->GetTotalPe(1);
    	if(Q <= 0) continue;
 
    	Float_t t = ((WCSimRootCherenkovHitTime *)fCherenkovHitTimes->At(myChHit->GetTotalPe(0)))->GetTruetime();	
		
		  WCSimRootGeom * myGeom = (WCSimRootGeom*)(WCSimGeometry::Instance())->GetWCSimGeometry();
    	WCSimRootPMT myPMT = myGeom->GetPMTFromTubeID(tubeID);
    	Double_t posX = myPMT.GetPosition(0);
    	Double_t posY = myPMT.GetPosition(1);
    	Double_t posZ = myPMT.GetPosition(2);
    	Double_t faceX = myPMT.GetOrientation(0);
    	Double_t faceY = myPMT.GetOrientation(1);
    	Double_t faceZ = myPMT.GetOrientation(2);
	  	hitPMT[arrayIndex] = kTRUE; 
	  	new(digitArray[arrayIndex]) WCSimLikelihoodDigit( posX, posY, posZ, t, Q, tubeID, faceX, faceY, faceZ  );
	}
	

     // Now fill all the unhit PMTs with an empty (q = 0) hit
	for(int arrayIndex = 0; arrayIndex < numPMT; ++arrayIndex)
	{
		if(hitPMT[arrayIndex] == kFALSE)
		{
      Int_t tubeID = arrayIndex + 1;
			WCSimRootCherenkovDigiHit * emptyHit = new WCSimRootCherenkovDigiHit(0,0,tubeID); // These index by tubeID so count from 1
			new (digitArray[arrayIndex]) WCSimLikelihoodDigit(emptyHit); // But the array indices go from zero
			delete emptyHit;
   	}
	}

    // Also save the array size 
    fNLikelihoodDigits = numPMT;
}


///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodDigitArray::~WCSimLikelihoodDigitArray()
{
}


///////////////////////////////////////////////////////////////////////////
// Getters
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodDigit * WCSimLikelihoodDigitArray::GetDigit( int digit )
{
  if(digit < fNLikelihoodDigits) return (WCSimLikelihoodDigit*)fLikelihoodDigitArray->At(digit);
  else
  {
    std::cout << "Error getting digit from fLikelihoodDigitArray: index out of bounds" << std::endl;
    return 0;
  }
}

Int_t WCSimLikelihoodDigitArray::GetNDigits(){ return fNLikelihoodDigits;}

Bool_t WCSimLikelihoodDigitArray::IsCylinder(){ return fGeomType == 0;}
Bool_t WCSimLikelihoodDigitArray::IsMailBox(){ return fGeomType == 1;}

Double_t WCSimLikelihoodDigitArray::GetExtent(Int_t i)
{
  if( 0 <= i && i < 3) return fExtent[i];
  else return -1.0;
}

