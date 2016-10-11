#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimGeometry.hh"
#include "WCSimReco.hh"
#include "WCSimRecoFactory.hh"
#include "WCSimRecoDigit.hh"
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
  fNumHitPMTs = 0;
	return;
}

WCSimLikelihoodDigitArray::WCSimLikelihoodDigitArray( WCSimRootEvent * myEvent )
{
  fNumHitPMTs = 0;
  // Check the number of inner detector PMTs and create a TClonesArray with one entry per PMT
//  Int_t numPMT = ((WCSimGeometry::Instance())->GetNumPMTs());	 
  Int_t numPMT = ((WCSimGeometry::Instance())->GetNumInnerPMTs());	 
  WCSimGeometry::PrintGeometry();
	fLikelihoodDigitArray = new TClonesArray("WCSimLikelihoodDigit",numPMT);
  TClonesArray &digitArray = *fLikelihoodDigitArray;

  // Get the geometry type and maximum x, y and z coordinates of the detector volume
  fGeomType = kUnknown;
  if((WCSimGeometry::Instance())->IsCylinder()){ fGeomType = kCylinder; }
  else if((WCSimGeometry::Instance())->IsMailBox()){ fGeomType = kMailBox;}
  
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
  WCSimReco * myReco = WCSimRecoFactory::Instance()->MakeReco(); // This calls new() - delete when done
  WCSimRecoEvent* myRecoEvent = WCSimInterface::RecoEvent();
  myReco->RunFilter(myRecoEvent);

	// Use this trigger to flag the tubes that were hit and fill the ClonesArray with
  // the corresponding WCSimLikelihoodDigit
  std::vector<WCSimRecoDigit*>* filteredDigits = myRecoEvent->GetFilterDigitList();
  std::vector<WCSimRecoDigit*>::iterator filterIter = filteredDigits->begin();
  fNumHitPMTs = 0;
  fFirstTime = -999;
  fLastTime = -999;
	while(filterIter != filteredDigits->end())
  {	
	  Int_t tubeID = (*filterIter)->GetTubeID();
    Int_t arrayIndex = tubeID - 1;

	  hitPMT[arrayIndex] = kTRUE; 

    // Things for likelihood digit constructor
    double q = (*filterIter)->GetRawQPEs();
    double t = (*filterIter)->GetRawTime();
    double pos[3] = {0,0,0};
    WCSimRootGeom * myGeom = (WCSimRootGeom*) (WCSimGeometry::Instance())->GetWCSimGeometry();
    WCSimRootPMT myPMT = myGeom->GetPMTFromTubeID(tubeID);
    pos[0]  = myPMT.GetPosition(0);
    pos[1]  = myPMT.GetPosition(1);
    pos[2]  = myPMT.GetPosition(2);
    double face[3] = {0,0,0};
    face[0] = myPMT.GetOrientation(0);
    face[1] = myPMT.GetOrientation(1);
    face[2] = myPMT.GetOrientation(2);
    TString name = myPMT.GetPMTName();

	  new(digitArray[arrayIndex]) WCSimLikelihoodDigit(pos[0], pos[1], pos[2], t, q, tubeID, face[0], face[1], face[2], name, WCSimDetectorParameters::WavelengthAveragedQE(name.Data()), WCSimDetectorParameters::QEAveragedRefIndex(name.Data()), WCSimDetectorParameters::PMTExposeHeight(name.Data()));
    fNumHitPMTs++;
    filterIter++;

      if(t < fFirstTime || fFirstTime == -999){ fFirstTime = t; }
      if(t > fLastTime || fLastTime == -999){ fLastTime = t; }
	}
	
  // Now fill all the unhit PMTs with an empty (q = 0) hit
	for(int arrayIndex = 0; arrayIndex < numPMT; ++arrayIndex)
	{
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
  std::cout << "There are " << fNLikelihoodDigits << " PMTs and " << fNumHitPMTs << " were hit" << std::endl;

}

/**
 * Use hits before digitizing them to help tune the charge likelihood
 * @TODO Handle this properly using one constructor and some kind of Initialise() function
 */
WCSimLikelihoodDigitArray::WCSimLikelihoodDigitArray( WCSimRootEvent * myRootEvent, Bool_t useUndigitized)
{
    if( ! useUndigitized )
    {
        std::cerr << "Error: bool useUndigitized has to be true" << std::endl;
        exit(EXIT_FAILURE);
    }
    fNumHitPMTs = 0;
	// Check the number of PMTs and create a TClonesArray with one entry per PMT
    Int_t numPMT = ((WCSimGeometry::Instance())->GetNumInnerPMTs());  // WCSimGeometry adds 1 to the true number to prevent overflow errors
  	WCSimGeometry::PrintGeometry();
  	std::cout << "There are " << numPMT << " PMTs" << std::endl;
	fLikelihoodDigitArray = new TClonesArray("WCSimLikelihoodDigit",numPMT);
  std::cout << "here 5" << std::endl;
  	TClonesArray &digitArray = *fLikelihoodDigitArray;
  std::cout << "here 6" << std::endl;

  	// Get the geometry type and maximum x, y and z coordinates of the detector volume
  	fGeomType = kUnknown;
  	if((WCSimGeometry::Instance())->IsCylinder()){ fGeomType = kCylinder; }
  	else if((WCSimGeometry::Instance())->IsMailBox()){ fGeomType = kMailBox;}
  
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
    fFirstTime = -999;
    fLastTime = -999;
	while( WCSimRootCherenkovHit * myChHit = (WCSimRootCherenkovHit *) undigiHitItr.Next())
  	{	
	  	Int_t tubeID = myChHit->GetTubeID();
      Int_t arrayIndex = tubeID-1;
    	Double_t Q = myChHit->GetTotalPe(1);
    	if(Q <= 0) continue;
 
    	Float_t t = ((WCSimRootCherenkovHitTime *)fCherenkovHitTimes->At(myChHit->GetTotalPe(0)))->GetTruetime();	
        if(t < fFirstTime || fFirstTime == -999){ fFirstTime = t; }
        if(t > fLastTime || fLastTime == -999){ fLastTime = t; }
		
    	WCSimRootGeom * myGeom = (WCSimRootGeom*)(WCSimGeometry::Instance())->GetWCSimGeometry();
    	WCSimRootPMT myPMT = myGeom->GetPMTFromTubeID(tubeID);
    	Double_t posX = myPMT.GetPosition(0);
    	Double_t posY = myPMT.GetPosition(1);
    	Double_t posZ = myPMT.GetPosition(2);
    	Double_t faceX = myPMT.GetOrientation(0);
    	Double_t faceY = myPMT.GetOrientation(1);
    	Double_t faceZ = myPMT.GetOrientation(2);
      TString name = myPMT.GetPMTName();
	  	hitPMT[arrayIndex] = kTRUE; 
	  	new(digitArray[arrayIndex]) WCSimLikelihoodDigit( posX, posY, posZ, t, Q, tubeID, faceX, faceY, faceZ, name, WCSimDetectorParameters::WavelengthAveragedQE(name.Data()), WCSimDetectorParameters::QEAveragedRefIndex(name.Data()), WCSimDetectorParameters::PMTExposeHeight(name.Data()) );
      fNumHitPMTs++;
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
Int_t WCSimLikelihoodDigitArray::GetNHits(){ return fNumHitPMTs;}

Bool_t WCSimLikelihoodDigitArray::IsCylinder(){ return fGeomType == kCylinder;}
Bool_t WCSimLikelihoodDigitArray::IsMailBox(){ return fGeomType == kMailBox;}
WCSimLikelihoodDigitArray::GeomType_t WCSimLikelihoodDigitArray::GetGeomType(){ return fGeomType; }

Double_t WCSimLikelihoodDigitArray::GetExtent(Int_t i)
{
  if( 0 <= i && i < 3) return fExtent[i];
  else return -1.0;
}

Double_t WCSimLikelihoodDigitArray::GetDuration() const
{
    return fLastTime - fFirstTime;
}

