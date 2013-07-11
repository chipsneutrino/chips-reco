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
    Int_t numPMT = (WCSimGeometry::Instance())->GetNumPMTs();
    std::cout << "There are " << numPMT << " PMTs" << std::endl;
	fLikelihoodDigitArray = new TClonesArray("WCSimLikelihoodDigit",numPMT);
    TClonesArray &digitArray = *fLikelihoodDigitArray;
		
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
    fNLikelihoodDigits = fLikelihoodDigitArray->GetEntriesFast();

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
