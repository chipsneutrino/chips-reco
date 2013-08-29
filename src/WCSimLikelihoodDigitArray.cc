#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimGeometry.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimInterface.hh"
#include "TCollection.h"

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
  Int_t numPMT = ((WCSimGeometry::Instance())->GetNumPMTs()) - 1;	// WCSimGeometry adds 1 to the true number to prevent overflow errors 
  WCSimGeometry::PrintGeometry();
  std::cout << "There are " << numPMT << " PMTs" << std::endl;
	fLikelihoodDigitArray = new TClonesArray("WCSimLikelihoodDigit",numPMT);
  TClonesArray &digitArray = *fLikelihoodDigitArray;

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
  TIter digiHitItr(fCherenkovDigiHits);
	while( WCSimRootCherenkovDigiHit * myChDigiHit = (WCSimRootCherenkovDigiHit *) digiHitItr.Next())
  {	
	  Int_t tubeID = myChDigiHit->GetTubeId();
	  hitPMT[tubeID] = kTRUE; 
	  new(digitArray[tubeID]) WCSimLikelihoodDigit(myChDigiHit);
	}
	

     // Now fill all the unhit PMTs with an empty (q = 0) hit
	for(int jPMT = 0; jPMT < numPMT; ++jPMT)
	{
		if(hitPMT[jPMT] == kFALSE)
		{
			WCSimRootCherenkovDigiHit * emptyHit = new WCSimRootCherenkovDigiHit(0,0,jPMT);
			new (digitArray[jPMT]) WCSimLikelihoodDigit(emptyHit);
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
    Int_t numPMT = ((WCSimGeometry::Instance())->GetNumPMTs()) - 1;  // WCSimGeometry adds 1 to the true number to prevent overflow errors
  	WCSimGeometry::PrintGeometry();
  	std::cout << "There are " << numPMT << " PMTs" << std::endl;
	fLikelihoodDigitArray = new TClonesArray("WCSimLikelihoodDigit",numPMT);
  	TClonesArray &digitArray = *fLikelihoodDigitArray;

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

	// Use this trigger to flag the tubes that were hit and fill the ClonesArray with
  	// the corresponding WCSimLikelihoodDigit
  	TClonesArray * fCherenkovDigiHits = myTrigger->GetCherenkovHits();
 	TClonesArray * fCherenkovHitTimes = myTrigger->GetCherenkovHitTimes();
 	
  	TIter undigiHitItr(fCherenkovDigiHits);
	while( WCSimRootCherenkovHit * myChHit = (WCSimRootCherenkovHit *) undigiHitItr.Next())
  	{	
	  	Int_t tubeID = myChHit->GetTubeID();
    	Double_t Q = myChHit->GetTotalPe(1);
    	if(Q <= 0) continue;
 
    	Float_t t = ((WCSimRootCherenkovHitTime *)fCherenkovHitTimes->At(myChHit->GetTotalPe(0)))->GetTruetime();	
		
		WCSimRootGeom * myGeom = (WCSimRootGeom*)(WCSimGeometry::Instance())->GetWCSimGeometry();
    	WCSimRootPMT myPMT = myGeom->GetPMT(tubeID);
    	Double_t posX = myPMT.GetPosition(0);
    	Double_t posY = myPMT.GetPosition(1);
    	Double_t posZ = myPMT.GetPosition(2);
    	Double_t faceX = myPMT.GetOrientation(0);
    	Double_t faceY = myPMT.GetOrientation(1);
    	Double_t faceZ = myPMT.GetOrientation(2);
	  	hitPMT[tubeID] = kTRUE; 
	  	new(digitArray[tubeID]) WCSimLikelihoodDigit( posX, posY, posZ, t, Q, tubeID, faceX, faceY, faceZ  );
	}
	

     // Now fill all the unhit PMTs with an empty (q = 0) hit
	for(int jPMT = 0; jPMT < numPMT; ++jPMT)
	{
		if(hitPMT[jPMT] == kFALSE)
		{
			WCSimRootCherenkovDigiHit * emptyHit = new WCSimRootCherenkovDigiHit(0,0,jPMT);
			new (digitArray[jPMT]) WCSimLikelihoodDigit(emptyHit);
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

