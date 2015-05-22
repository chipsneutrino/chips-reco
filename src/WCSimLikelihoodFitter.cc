#include "WCSimAnalysisConfig.hh"
#include "WCSimChargePredictor.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimFitterTree.hh"
#include "WCSimFitterInterface.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimReco.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoFactory.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimTrackParameterEnums.hh"

#include "TClonesArray.h"
#include "TCollection.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TStopwatch.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodFitter)
#endif

/**
 * @todo Remove hardcoding of track type
 */
WCSimLikelihoodFitter::WCSimLikelihoodFitter()
{
  fFitterPlots = NULL;
  fFitterTree = 0x0;
  fTotalLikelihood = NULL;
  fRootEvent = NULL;
  fLikelihoodDigitArray = NULL;
  fTrueLikelihoodTracks = NULL;

  fStartVals = NULL;
  fMinVals = NULL;
  fMaxVals = NULL;
  fStepSizes = NULL;
  fFixed = NULL;
  fIsEnergy = NULL;
  fNames = NULL;
  fUseHoughFitterForSeed = WCSimAnalysisConfig::Instance()->GetUseHoughFitterForSeed();
  ResetEvent();
}

WCSimLikelihoodFitter::~WCSimLikelihoodFitter()
{
}

UInt_t WCSimLikelihoodFitter::GetNPars()
{
  // Do we know how to fit this number of tracks?
  unsigned int nPars = WCSimFitterConfig::Instance()->GetNumIndependentParameters();
  return nPars;
}



void WCSimLikelihoodFitter::Minimize2LnL()
{
   // Set up the minimizer
   ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

   // Alternatively: use a different algorithm to check the minimizer works
   // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "GSLSimAn");

   if( fTrackParMap.size() == 0)
   {
     SetTrackParameterMap();
   }

   min->SetMaxFunctionCalls(100);
   min->SetMaxIterations(9);
   min->SetPrintLevel(3);
   min->SetTolerance(200);
   min->SetStrategy(2);
   std::cout << " Tolerance = " << min->Tolerance() << std::endl;

  // Convert nTracks into number of parameters
  const unsigned int nPars = this->GetNPars();
  
  // Tell the minimizer the function to minimize
  // We have to wrap it in this functor because it's a member function of this class
  ROOT::Math::Functor func(this,&WCSimLikelihoodFitter::WrapFunc, nPars);
  
  // Print the parameters we're using
  for(UInt_t j = 0; j < nPars; ++j)
  {
      std::cout << j << "   " << fNames[j] << "   " << fStartVals[j] << "   " << fMinVals[j] << "   " << fMaxVals[j] << std::endl;
  }

  // Tell the minimizer the functor and variables to consider
  min->SetFunction(func);

  // Set the fixed and free variables
  Bool_t isAnythingFree = false;

  // Use the old track fitter to seed better start values:
  // The seed algorithm will update fStartVals
  if( fUseHoughFitterForSeed )
  {
	  SeedParams();
  }


  for(UInt_t i = 0; i < nPars; ++i)
  {
    // Sometimes the fast seed goes outside the allowed range:
    if(fStartVals[i] < fMinVals[i])
    {
      fStartVals[i] = fMinVals[i] + 0.01 * ( fMaxVals[i] - fMinVals[i] );
    }
    if(fStartVals[i] > fMaxVals[i])
    {
      fStartVals[i] = fMaxVals[i] - 0.01 * (fMaxVals[i] - fMinVals[i]);
    }


	  if( fFixed[i] || fIsEnergy[i])
	  {
	  	min->SetFixedVariable(i, fNames[i], fStartVals[i]);
	  }
	  else
	  {
		  isAnythingFree = true;
	  	  min->SetLimitedVariable(i, fNames[i], fStartVals[i], fStepSizes[i], fMinVals[i], fMaxVals[i]);
	  }
  }

  // Perform the minimization
  fCalls = 0;
  if( isAnythingFree )
  {
    min->Minimize();
    fMinimum = min->MinValue();
    fStatus = min->Status();
  }
  
  // Now we go again but free up the energy:
  Bool_t doSecondStep = false;

  // This is what we ought to do:
  const Double_t * outPar = min->X();
  for(UInt_t i = 0; i < nPars; ++i)
  {
    if(fIsEnergy[i] && !fFixed[i])
    {
      doSecondStep = true;
      min->SetLimitedVariable(i, fNames[i], fStartVals[i], fStepSizes[i], fMinVals[i], fMaxVals[i]);
    }
    else
    {
      min->SetFixedVariable(i, fNames[i], outPar[i]);
    }
  }
  
  //  // Either we want to minimise over the energy, or we haven't freed anything up at all
  //  // and to perform the calculation at least once
  //  if( doSecondStep || !isAnythingFree )
  //  {
  //    min->Minimize();
  //  	fMinimum = min->MinValue();
  //    outPar = min->X();
  //  }
  
  // Minuit gives us a const double array, but (for now) we want to do the grid search and then modify it
  std::cout << "Here's outPar..." << std::endl;
  std::cout << outPar[0] << "  " << outPar[1] << std::endl;
  std::vector<Double_t> outParVec;
  for(unsigned 	int iPar = 0; iPar < min->NDim(); ++iPar)
  {
    outParVec.push_back(outPar[iPar]);
  }

  // But actually we'll just do a grid search:
  if( doSecondStep || !isAnythingFree )
  {
  	Double_t bestLnL = 0.0;
  	std::vector<Double_t> bestEnergies;
  	PerformEnergyGridSearch(bestLnL, bestEnergies);
    std::cout << "Size of bestEnergies = " << bestEnergies.size() << std::endl;
    for(UInt_t iTrack = 0; iTrack < WCSimFitterConfig::Instance()->GetNumTracks(); ++iTrack)
    {
      std::pair<UInt_t, FitterParameterType::Type> trackEnergy = std::make_pair(iTrack, FitterParameterType::kEnergy);
      std::cout << "Making output parameter vector" << std::endl;
      std::cout << "Size of outParVec = " << outParVec.size() << " and size of fTrackParMap = " << fTrackParMap.size() << std::endl;
      std::cout << "Trying to fill " << fTrackParMap[trackEnergy] << std::endl;
      outParVec.at(fTrackParMap[trackEnergy]) = bestEnergies.at(iTrack);
      std::cout << "Made output vector" << std::endl;
    }
  }


  // Get and print the fit results
  // std::cout << "Best fit track: " << std::endl;
  // RescaleParams(outPar2[0],  outPar2[1],  outPar2[2],  outPar2[3],  outPar2[4],  outPar2[5],  outPar2[6], fType).Print();

  // Now save the best-fit tracks
  fBestFit.clear();
  for(UInt_t iTrack = 0; iTrack < WCSimFitterConfig::Instance()->GetNumTracks(); ++iTrack)
  {
	  std::pair<UInt_t, FitterParameterType::Type> trackParX = std::make_pair(iTrack, FitterParameterType::kVtxX);
	  std::pair<UInt_t, FitterParameterType::Type> trackParY = std::make_pair(iTrack, FitterParameterType::kVtxY);
	  std::pair<UInt_t, FitterParameterType::Type> trackParZ = std::make_pair(iTrack, FitterParameterType::kVtxZ);
	  std::pair<UInt_t, FitterParameterType::Type> trackParT = std::make_pair(iTrack, FitterParameterType::kVtxT);
	  std::pair<UInt_t, FitterParameterType::Type> trackParTh = std::make_pair(iTrack, FitterParameterType::kDirTh);
	  std::pair<UInt_t, FitterParameterType::Type> trackParPhi = std::make_pair(iTrack, FitterParameterType::kDirPhi);
	  std::pair<UInt_t, FitterParameterType::Type> trackParE = std::make_pair(iTrack, FitterParameterType::kEnergy);
	std::cout << "Track " << iTrack << ", parameter " << FitterParameterType::AsString(trackParX.second) << "   outParVec.at(" << fTrackParMap[trackParX] << ") = ";  std::cout << outParVec.at(fTrackParMap[trackParX]) << std::endl;
	std::cout << "Track " << iTrack << ", parameter " << FitterParameterType::AsString(trackParY.second) << "   outParVec.at(" << fTrackParMap[trackParY] << ") = ";  std::cout << outParVec.at(fTrackParMap[trackParY]) << std::endl;
	std::cout << "Track " << iTrack << ", parameter " << FitterParameterType::AsString(trackParZ.second) << "   outParVec.at(" << fTrackParMap[trackParZ] << ") = ";  std::cout << outParVec.at(fTrackParMap[trackParZ]) << std::endl;
	std::cout << "Track " << iTrack << ", parameter " << FitterParameterType::AsString(trackParT.second) << "   outParVec.at(" << fTrackParMap[trackParT] << ") = ";  std::cout << outParVec.at(fTrackParMap[trackParT]) << std::endl;
	std::cout << "Track " << iTrack << ", parameter " << FitterParameterType::AsString(trackParTh.second) << "   outParVec.at(" << fTrackParMap[trackParTh] << ") = ";  std::cout << outParVec.at(fTrackParMap[trackParTh]) << std::endl;
	std::cout << "Track " << iTrack << ", parameter " << FitterParameterType::AsString(trackParPhi.second) << "   outParVec.at(" << fTrackParMap[trackParPhi] << ") = ";  std::cout << outParVec.at(fTrackParMap[trackParPhi]) << std::endl;
	std::cout << "Track " << iTrack << ", parameter " << FitterParameterType::AsString(trackParE.second) << "   outParVec.at(" << fTrackParMap[trackParE] << ") = ";  std::cout << outParVec.at(fTrackParMap[trackParE]) << std::endl;
	  WCSimLikelihoodTrack track(
								  outParVec.at(fTrackParMap[trackParX]),
								  outParVec.at(fTrackParMap[trackParY]),
								  outParVec.at(fTrackParMap[trackParZ]),
								  outParVec.at(fTrackParMap[trackParT]),
								  outParVec.at(fTrackParMap[trackParTh]),
								  outParVec.at(fTrackParMap[trackParPhi]),
								  outParVec.at(fTrackParMap[trackParE]),
								  fTypes.at(iTrack)
								  );
	  fBestFit.push_back(track);
  }

  return;
}

/**
 * Wrapper: constructs the correct number of track objects
 * then works out the total likelihood
 */
Double_t WCSimLikelihoodFitter::WrapFunc(const Double_t * x)
{
  std::cout << "Function evaluated " << ++fCalls << " times for event " << fEvent << std::endl;
  UInt_t nTracks = WCSimFitterConfig::Instance()->GetNumTracks();
  std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes();

  // If we've fixed some track parameters together, then our array of fit parameters won't just be of size
  // n tracks * m parameters per track, and we need to work out which entry corresponds to which parameter
  // of which track(s) - this is probably expensive, so we'll do it once and look it up in a map thereafter

  std::vector<WCSimLikelihoodTrack> tracksToFit;
  for(UInt_t iTrack = 0 ; iTrack < nTracks; ++iTrack)
  {
	  std::pair<UInt_t, FitterParameterType::Type> trackParX = std::make_pair(iTrack, FitterParameterType::kVtxX);
	  std::pair<UInt_t, FitterParameterType::Type> trackParY = std::make_pair(iTrack, FitterParameterType::kVtxY);
	  std::pair<UInt_t, FitterParameterType::Type> trackParZ = std::make_pair(iTrack, FitterParameterType::kVtxZ);
	  std::pair<UInt_t, FitterParameterType::Type> trackParT = std::make_pair(iTrack, FitterParameterType::kVtxT);
	  std::pair<UInt_t, FitterParameterType::Type> trackParTh = std::make_pair(iTrack, FitterParameterType::kDirTh);
	  std::pair<UInt_t, FitterParameterType::Type> trackParPhi = std::make_pair(iTrack, FitterParameterType::kDirPhi);
	  std::pair<UInt_t, FitterParameterType::Type> trackParE = std::make_pair(iTrack, FitterParameterType::kEnergy);

	  WCSimLikelihoodTrack track(
								  x[fTrackParMap[trackParX]],
								  x[fTrackParMap[trackParY]],
								  x[fTrackParMap[trackParZ]],
								  x[fTrackParMap[trackParT]],
								  x[fTrackParMap[trackParTh]],
								  x[fTrackParMap[trackParPhi]],
								  x[fTrackParMap[trackParE]],
								  fTypes.at(iTrack)
								  );
	  tracksToFit.push_back(track);
    std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
  }



  // if(fIsFirstCall ) // TEMP: Print it out all the time for now
  // {
  //   std::cout << "Tracks used for first call:" << std::endl;
  //   std::vector<WCSimLikelihoodTrack>::iterator itr = tracksToFit.begin();
  //   for( ; itr < tracksToFit.end(); ++itr)
  //   {
  //     (*itr).Print();
  //   }
  // }

  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(tracksToFit);
  minus2LnL = fTotalLikelihood->Calc2LnL();

  return minus2LnL;
}

void WCSimLikelihoodFitter::SetTrackParameterMap()
{
  fTrackParMap.clear();
  UInt_t arrayCounter = 0;
  UInt_t nTracks = WCSimFitterConfig::Instance()->GetNumTracks();
  std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes();

  for(UInt_t iTrack = 0; iTrack < nTracks; ++iTrack)
  {
	  for(UInt_t iParam = 0; iParam < allTypes.size(); ++iParam)
	  {
		  // Is it joined?
		  std::pair<UInt_t, FitterParameterType::Type> trackPar = std::make_pair(iTrack, allTypes.at(iParam));
		  UInt_t isJoinedWith = WCSimFitterConfig::Instance()->GetTrackIsJoinedWith(iTrack, FitterParameterType::AsString(trackPar.second).c_str());
		  if( isJoinedWith == iTrack ) 	// Itself, i.e. not joined with any other track
		  {
			  // Then set it to the value from the first track
			  assert(arrayCounter < WCSimFitterConfig::Instance()->GetNumIndependentParameters());
			  fTrackParMap[trackPar] = arrayCounter;
			  arrayCounter++;
		  }
		  else
		  {
			  // Otherwise don't; set it to the next index value instead
			  std::pair<UInt_t, FitterParameterType::Type> joinPar = std::make_pair(isJoinedWith, allTypes.at(iParam));
			  fTrackParMap[trackPar] = fTrackParMap[joinPar];
		  }
	  }
  }
}

// Use the previous Hough transform-based reconstruction algorithm to seed the multi-track one
void WCSimLikelihoodFitter::SeedParams()
{
	std::cout << " *** WCSimLikelihoodFitter::SeedParams() *** " << std::endl;
	WCSimReco * myReco = WCSimRecoFactory::Instance()->MakeReco(); // This calls new() - delete when done
	WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
	// myReco->RunFilter(recoEvent); // Might need to do this first, I'm not sure
	myReco->RunRecoVertex(recoEvent);

	fSeedMap.clear();
	fSeedMap.insert(std::make_pair(FitterParameterType::kVtxX, recoEvent->GetVtxX()));
	fSeedMap.insert(std::make_pair(FitterParameterType::kVtxY, recoEvent->GetVtxY()));
	fSeedMap.insert(std::make_pair(FitterParameterType::kVtxZ, recoEvent->GetVtxZ()));

	Double_t recoDirX = recoEvent->GetDirX();
	Double_t recoDirY = recoEvent->GetDirY();
	Double_t recoDirZ = recoEvent->GetDirZ();
	Double_t seedTheta   = TMath::ACos(recoDirZ);
	Double_t seedPhi     = -999.9;
	if(recoDirY != 0.0){ seedPhi = TMath::ATan2(recoDirY,recoDirX); }// Ensure range is -pi to pi
	else{ seedPhi = (recoDirX < 0.0)? 0.5*TMath::Pi() : -0.5*TMath::Pi(); }


	fSeedMap.insert(std::make_pair(FitterParameterType::kDirTh, seedTheta));
	fSeedMap.insert(std::make_pair(FitterParameterType::kDirPhi, seedPhi));

	//	std::cout << "Seed values:" << std::endl;
	//	std::map<FitterParameterType::Type, Double_t>::const_iterator printer = fSeedMap.begin();
	//	for(; printer != fSeedMap.end(); ++printer)
	//	{
	//		std::cout << "Seed for " << FitterParameterType::AsString((*printer).first) << " = " << (*printer).second << std::endl;
	//	}
	//   fSeedTheta = 0.5*TMath::Pi();
	//   fSeedPhi = 0.25*TMath::Pi();

	// For every fit parameter we've set, find the array index for that parameter in fStartVals
	// And set that array entry to the corresponding seed value we just worked out
	//                 Track             Type            Array index
	std::map<std::pair<UInt_t, FitterParameterType::Type>, UInt_t >::iterator trackParItr;

	// Loop over all parameters for all tracks
	for(trackParItr = fTrackParMap.begin(); trackParItr != fTrackParMap.end(); ++trackParItr)
	{
	   // Find the seed corresponding this parameter type
	   FitterParameterType::Type searchType = ((*trackParItr).first).second;
	   std::map<FitterParameterType::Type, Double_t>::const_iterator seedSearch = fSeedMap.find(searchType);
	   if(seedSearch != fSeedMap.end())
	   {
		   // Then set the start value for this parameter to the seed value
		   UInt_t arrIndex = (*trackParItr).second;
		   if(!fFixed[arrIndex])
		   {
			   std::cout << "Setting seed for " << arrIndex << " to " << (*seedSearch).second << std::endl;
			   fStartVals[arrIndex] = (*seedSearch).second;
		   }
	   }
	}
	delete myReco;
	// std::cout << "Seed dir = (vx, vy, vz) = ( " << recoDirX << "," << recoDirY << "," << recoDirZ << ")" << std::endl
	// std::cout << "-> theta = " << fSeedTheta << "    phi = " << fSeedPhi << std::endl
	//           << "Seed position = (x,y,z) = ( " << fSeedVtxX << "," << fSeedVtxY << "," <<fSeedVtxZ << ")" << std::endl;
}


std::vector<WCSimLikelihoodTrack> WCSimLikelihoodFitter::GetBestFit()
{
  return fBestFit;
}

Double_t WCSimLikelihoodFitter::GetMinimum()
{
  return fMinimum;
}

void WCSimLikelihoodFitter::FitEventNumber(Int_t iEvent) {
	fIsFirstCall = true;
	fEvent = iEvent;
	fRootEvent = WCSimInterface::Instance()->GetWCSimEvent(iEvent);
	fLikelihoodDigitArray = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(iEvent);
  std::cout << "There are " << fLikelihoodDigitArray->GetNDigits() << " digits" << std::endl;
	if(fTotalLikelihood != NULL)
	{
		delete fTotalLikelihood;
    fTotalLikelihood = NULL;
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}
	else
	{
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}


	if(WCSimFitterInterface::Instance()->GetMakeFits() ) {
		WCSimInterface::Instance()->BuildEvent(iEvent);
    std::cout << "Fitting event " << iEvent << std::endl;
		Minimize2LnL();
		fTrueLikelihoodTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
		FillPlots();
		FillTree();
		FillHitComparison();
	}

	return;


}

void WCSimLikelihoodFitter::CreateParameterArrays() {
	Int_t numParams = WCSimFitterConfig::Instance()->GetNumIndependentParameters();
	if(fStartVals != NULL){ delete [] fStartVals; }
	if(fMinVals != NULL){   delete [] fMinVals; }
	if(fMaxVals != NULL){   delete [] fMaxVals; }
	if(fStepSizes != NULL){ delete [] fStepSizes; }
	if(fFixed != NULL){     delete [] fFixed; }
	if(fIsEnergy!= NULL){   delete [] fIsEnergy; }
	if(fNames != NULL){     delete [] fNames; }
	fStartVals = new Double_t[numParams];
	fMinVals   = new Double_t[numParams];
	fMaxVals   = new Double_t[numParams];
	fStepSizes = new Double_t[numParams];
	fFixed	   = new Bool_t[numParams];
	fIsEnergy  = new Bool_t[numParams];
	fNames 	   = new std::string[numParams];
}

void WCSimLikelihoodFitter::SetParameterArrays() {
	std::cout << "Setting parameter arrays" << std::endl;
	Int_t iParam = 0;
	UInt_t numTracks = WCSimFitterConfig::Instance()->GetNumTracks();
	fTypes.clear();
	for(UInt_t jTrack = 0; jTrack < numTracks; ++jTrack)
	{
		fTypes.push_back(WCSimFitterConfig::Instance()->GetTrackType(jTrack));
		std::vector<FitterParameterType::Type> paramTypes = FitterParameterType::GetAllAllowedTypes();
		std::vector<FitterParameterType::Type>::const_iterator typeItr = paramTypes.begin();
		for( ; typeItr != paramTypes.end(); ++typeItr)
		{
			const char * name  = FitterParameterType::AsString(*typeItr).c_str();

			if(WCSimFitterConfig::Instance()->GetIsParameterJoined(jTrack, name)){
				continue;
			}
			else{
				fStartVals[iParam] = WCSimFitterConfig::Instance()->GetParStart(jTrack,name);
				fMinVals[iParam]   = WCSimFitterConfig::Instance()->GetParMin(jTrack,name);
				fMaxVals[iParam]   = WCSimFitterConfig::Instance()->GetParMax(jTrack,name);
				fStepSizes[iParam] = (WCSimFitterConfig::Instance()->GetParMax(jTrack,name) - WCSimFitterConfig::Instance()->GetParMin(jTrack, name)) / 10.0;
				fFixed[iParam]     = WCSimFitterConfig::Instance()->GetIsFixedParameter(jTrack,name);
				TString parName;
				parName.Form("%s_track%d", name, jTrack);
				fNames[iParam]     = std::string(parName.Data());
        if( (*typeItr) == FitterParameterType::kEnergy )
        {
          fIsEnergy[iParam] = kTRUE;
        }
        else
        {
          fIsEnergy[iParam] = kFALSE;
        }

        std::cout << "Setting parameter: " << "name = " << fNames[iParam] << " min = " << fMinVals[iParam] << "   max = " << fMaxVals[iParam] << "  start = " << fStartVals[iParam] << "  fix = " << fFixed[iParam] << std::endl;
				iParam++;
			}
		}
	}
	return;
}

void WCSimLikelihoodFitter::DeleteParameterArrays() {
	if( fStartVals != NULL ) { delete [] fStartVals; }
	if( fMinVals != NULL )   { delete [] fMinVals; }
	if( fMaxVals != NULL )   { delete [] fMaxVals; }
	if( fStepSizes != NULL ) { delete [] fStepSizes; }
	if( fFixed != NULL ) 	 { delete [] fFixed; }
	if( fNames != NULL ) 	 { delete [] fNames; }
}

void WCSimLikelihoodFitter::SetFitterPlots(WCSimFitterPlots* fitterPlots) {
  std::cout << "Setting fFitterPlots" << std::endl;
	fFitterPlots = fitterPlots;
  std::cout << "Number of surface bins" << fFitterPlots->GetNumSurfaceBins() << std::endl;
}

void WCSimLikelihoodFitter::Make1DSurface(
	std::pair<UInt_t, FitterParameterType::Type> trackPar) {
  std::cout << " *** WCSimLikelihoodFitter::Make1DSurface() *** " << std::endl;
 
  SetTrackParameterMap(); 
  std::map<std::pair<UInt_t, FitterParameterType::Type>, UInt_t >::iterator itr = fTrackParMap.find(trackPar);
	if( itr != fTrackParMap.end() )
	{
		unsigned int arrIndex = (*itr).second;
    std::cout << "arrIndex = " << arrIndex << std::endl;
		double min = fMinVals[arrIndex];
		double max = fMaxVals[arrIndex];

		Double_t * x = new Double_t[WCSimFitterConfig::Instance()->GetNumIndependentParameters()];
		for(unsigned int iBin = 0; iBin < WCSimFitterConfig::Instance()->GetNumIndependentParameters(); ++iBin)
		{
			x[iBin] = fStartVals[iBin];
		}

		for(unsigned int iBin = 1; iBin <= fFitterPlots->GetNumSurfaceBins(); ++iBin)
		{
			double stepVal = min + (max - min) * iBin / (double)fFitterPlots->GetNumSurfaceBins();
			x[arrIndex] = stepVal;
			fFitterPlots->Fill1DProfile(trackPar, iBin, WrapFunc(x));
		}
		delete [] x;
	}


}

void WCSimLikelihoodFitter::Make2DSurface(
		std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > trackPar) 
{
  std::cout << " *** WCSimLikelihoodFitter::Make1DSurface() *** " << std::endl;
  
  SetTrackParameterMap(); 
  std::map<std::pair<UInt_t, FitterParameterType::Type>, UInt_t >::iterator itr1 = fTrackParMap.find(trackPar.first);
  std::map<std::pair<UInt_t, FitterParameterType::Type>, UInt_t >::iterator itr2 = fTrackParMap.find(trackPar.second);

  if( itr1 != fTrackParMap.end() && itr2 != fTrackParMap.end())
  {
	  unsigned int arrIndexX = (*itr1).second;
	  double minX = fMinVals[arrIndexX];
	  double maxX = fMaxVals[arrIndexX];

	  unsigned int arrIndexY = (*itr2).second;
	  double minY = fMinVals[arrIndexY];
	  double maxY = fMaxVals[arrIndexY];


	  Double_t * x = new Double_t[WCSimFitterConfig::Instance()->GetNumIndependentParameters()];
	  // Set everything to the default value
	  for(unsigned int iXBin = 0; iXBin < WCSimFitterConfig::Instance()->GetNumIndependentParameters(); ++iXBin)
	  {
	  	x[iXBin] = fStartVals[iXBin];
	  }

	  // Now step over each variable
	  for(unsigned int iXBin = 0; iXBin < fFitterPlots->GetNumSurfaceBins(); ++iXBin)
	  {
	  	double stepValX = minX + (maxX - minX) * iXBin / (double)fFitterPlots->GetNumSurfaceBins();
	  	x[arrIndexX] = stepValX;


	  	for(unsigned int iYBin = 0; iYBin < fFitterPlots->GetNumSurfaceBins(); ++iYBin)
	  	{
	  		double stepValY = minY + (maxY - minY) * iYBin / (double)fFitterPlots->GetNumSurfaceBins();
	  		x[arrIndexY] = stepValY;
	  		fFitterPlots->Fill2DProfile(trackPar, iXBin, iYBin, WrapFunc(x));
	  	}
	  }
  }

  return;
}

void WCSimLikelihoodFitter::PerformEnergyGridSearch(Double_t& best2LnL,
		std::vector<Double_t>& bestEnergies)
{
  SetTrackParameterMap();

	std::pair<unsigned int, FitterParameterType::Type> trackAndEnergy;
	if( WCSimFitterConfig::Instance()->GetNumTracks() > 2 )
	{
		std::cerr << "Error: Performing a grid search with more than two tracks will take an extremely long time" << std::endl;
		std::cerr << "       I'm sorry Dave, I can't let you do that" << std::endl;
		assert(WCSimFitterConfig::Instance()->GetNumTracks() <= 2);
	}


	// Get energies allowed for the first track
	std::pair<UInt_t, FitterParameterType::Type> trackZeroEnergy = std::make_pair(0, FitterParameterType::kEnergy);
	std::map<std::pair<UInt_t, FitterParameterType::Type>, UInt_t >::iterator trackZeroEnergyItr = fTrackParMap.find(trackZeroEnergy);
	UInt_t trackZeroIndex = (*trackZeroEnergyItr).second;
  std::cout << "Could I find it? " << (trackZeroEnergyItr != fTrackParMap.end()) << std::endl;
  std::cout << "Array index for track 0 energy is " << trackZeroIndex << std::endl;
	WCSimLikelihoodTrack::TrackType trackZeroType = WCSimFitterConfig::Instance()->GetTrackType(0);
  
  std::cout << "Track zero type should be " << WCSimLikelihoodTrack::TrackTypeToString(trackZeroType) << std::endl;
	WCSimLikelihoodTrack * tempTrack = new WCSimLikelihoodTrack();
	tempTrack->SetType(trackZeroType);
  tempTrack->SetE(WCSimFitterConfig::Instance()->GetParStart(0, "kEnergy"));

  std::cout << "Making empty emission profile" << std::endl;
	WCSimEmissionProfiles * tempProfileZero = new WCSimEmissionProfiles();
  
  std::cout << "Setting track 0" << std::endl;
  tempProfileZero->SetTrack(tempTrack);
  
  std::cout << "Getting profile energies" << std::endl;
	std::vector<Double_t> energyBinsZero = tempProfileZero->GetProfileEnergies();

  // Emission profiles can only be used for energies between (not equal to) the first and last bins
  assert(energyBinsZero.size() > 2);
  for(unsigned int i = 0; i < energyBinsZero.size(); ++i)
  {
    std::cout << "Energy bin " << i << " = " << energyBinsZero.at(i) << std::endl;
  }
  energyBinsZero.pop_back();
  energyBinsZero.erase(energyBinsZero.begin(), energyBinsZero.begin() + 1);
  
  std::cout << "after fiddling around..." << std::endl;
  for(unsigned int i = 0; i < energyBinsZero.size(); ++i)
  {
    std::cout << "Energy bin " << i << " = " << energyBinsZero.at(i) << std::endl;
  }

	delete tempTrack;
	delete tempProfileZero;

	// Now for the second track
	WCSimLikelihoodTrack::TrackType trackOneType = WCSimLikelihoodTrack::Unknown;
	std::pair<UInt_t, FitterParameterType::Type> trackOneEnergy = std::make_pair(1,FitterParameterType::kEnergy);
  std::map<std::pair<UInt_t, FitterParameterType::Type>, UInt_t >::iterator trackOneEnergyItr = fTrackParMap.find(trackOneEnergy);
	int trackOneIndex = WCSimFitterConfig::Instance()->GetNumTracks() > 1 ? (*trackOneEnergyItr).second : 0;
  std::vector<Double_t> energyBinsOne;
  if( WCSimFitterConfig::Instance()->GetNumTracks() > 1)
	{
		trackOneType = WCSimFitterConfig::Instance()->GetTrackType(1);
		WCSimLikelihoodTrack * tempTrack = new WCSimLikelihoodTrack();
		tempTrack->SetType(trackOneType);
    tempTrack->SetE(WCSimFitterConfig::Instance()->GetParStart(1, "kEnergy"));
  	WCSimEmissionProfiles * tempProfileOne = new WCSimEmissionProfiles();
    tempProfileOne->SetTrack(tempTrack);
		energyBinsOne = tempProfileOne->GetProfileEnergies();
    assert(energyBinsOne.size() > 2);
    energyBinsOne.pop_back();
    energyBinsOne.erase(energyBinsOne.begin(), energyBinsOne.begin() + 1);
		delete tempTrack;
		delete tempProfileOne;

	}

	// This loop assumes we have at most two tracks
	best2LnL = -999.9;

	Double_t * x = new Double_t[WCSimFitterConfig::Instance()->GetNumIndependentParameters()];
  // Set everything to the default value
  for(unsigned int iXBin = 0; iXBin < WCSimFitterConfig::Instance()->GetNumIndependentParameters(); ++iXBin)
  {
  	x[iXBin] = fStartVals[iXBin];
  }


	bool isFirstLoop = true;
  for(unsigned int iBin = 0; iBin < energyBinsZero.size(); ++iBin)
	{
    std::cout << "TrackZeroIndex = " << trackZeroIndex << std::endl;
    std::cout << "Energy bin " << iBin << " is " << energyBinsZero.at(iBin) << std::endl;
  	x[trackZeroIndex] = energyBinsZero.at(iBin);

    std::cout << "about to do the loop" << std::endl;


    unsigned int jBin = 0;
    do
  	{
      std::cout << "jBin = " << jBin << std::endl;
  		if( WCSimFitterConfig::Instance()->GetNumTracks() > 1)
  		{
  			x[trackOneIndex] = energyBinsOne.at(jBin);
  		}
  		double tempMin = WrapFunc(x);
			if( tempMin < best2LnL || isFirstLoop )
			{
				best2LnL = tempMin;
				isFirstLoop = false;
				bestEnergies.clear();
        std::cout << "Pushing back for iBin = " << iBin << std::endl;
				bestEnergies.push_back(energyBinsZero.at(iBin));
				if( WCSimFitterConfig::Instance()->GetNumTracks() == 2 )
				{
          std::cout << "And for jBin = " << std::endl;
					bestEnergies.push_back(energyBinsOne.at(jBin));
				}
        jBin++;
			}
  	} while(jBin < energyBinsOne.size());


	}
  fMinimum = best2LnL;
  std::cout << "Grid search finished" << std::endl;

}

void WCSimLikelihoodFitter::FillHitComparison() {
	std::vector<double> predictedCharges = fTotalLikelihood->GetPredictedChargeVector();
	std::vector<double> measuredCharges = fTotalLikelihood->GetMeasuredChargeVector();
	std::vector<double> best2LnLs = fTotalLikelihood->GetTotal2LnLVector();

	std::vector<WCSimLikelihoodTrack*> *correctTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
	std::vector<WCSimLikelihoodTrack> correctTracksNoPtr;
	for(int i = 0; i < correctTracks->size(); ++i)
	{
		correctTracksNoPtr.push_back(*(correctTracks->at(i)));
	  correctTracksNoPtr.at(i).Print();
  }
  
	fTotalLikelihood->SetTracks(correctTracksNoPtr);
	fTotalLikelihood->Calc2LnL();
	std::vector<double> correctPredictedCharges = fTotalLikelihood->GetPredictedChargeVector();
	std::vector<double> correct2LnLs = fTotalLikelihood->GetTotal2LnLVector();
	fFitterTree->FillHitComparison(fLikelihoodDigitArray, predictedCharges, correctPredictedCharges, measuredCharges, best2LnLs, correct2LnLs);


}

WCSimLikelihoodTrack WCSimLikelihoodFitter::RescaleParams(Double_t x, Double_t y, Double_t z, Double_t t,
                                                          Double_t th, Double_t phi, Double_t E, 
                                                          WCSimLikelihoodTrack::TrackType type)
{
  // Currently do nothing...
  return WCSimLikelihoodTrack(x,y,z,t,th,phi,E,type);

  // But could do this - scale variables in the range (0,1) up to proper track parameers:
  double pi = TMath::Pi();
  double x2,y2,z2,t2,th2,phi2,E2;
  // x2   = (1900.+x)/3800.;
  // y2   = (1900.+y)/3800.;
  // z2   = (940. +z)/1880.;
  // t2   = t;
  // th2  = th / pi;
  // phi2 = (phi+pi)/(2*pi);
  // E2   = E / 1500;

  // Dirty hack: detector dimensions hardcoded for quick testing
  x2 = 3800.0 * x - 1900.;
  y2 = 3800.0 * y - 1900.;
  z2 = 1880.0 * z - 940.;
  th2 = pi * th;
  std::cout << "theta2 = " << th2 << std::endl;
  phi2 = 2*pi*phi - pi;
  E2 = 1500*E;
  WCSimLikelihoodTrack myTrack(x2,y2,z2,t2,th2,phi2,E2,type);
  return myTrack;
}

Int_t WCSimLikelihoodFitter::GetStatus()
{
  return fStatus;
}

void WCSimLikelihoodFitter::RunFits() {
	CreateParameterArrays();

	UInt_t firstEvent = WCSimFitterConfig::Instance()->GetFirstEventToFit();
	UInt_t numEventsToFit = WCSimFitterConfig::Instance()->GetNumEventsToFit();


	for(UInt_t iEvent = firstEvent; iEvent < firstEvent + WCSimFitterConfig::Instance()->GetNumEventsToFit() ; ++iEvent)
	{
		SetParameterArrays();
		FitEventNumber(iEvent);
    WCSimFitterInterface::Instance()->SaveResults();
	}

	DeleteParameterArrays();
}

void WCSimLikelihoodFitter::RunSurfaces() {
	
  CreateParameterArrays();
  SetParameterArrays();
  fRootEvent = WCSimInterface::Instance()->GetWCSimEvent(WCSimFitterConfig::Instance()->GetFirstEventToFit());
  fLikelihoodDigitArray = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(WCSimFitterConfig::Instance()->GetFirstEventToFit());
  std::cout << "There are " << fLikelihoodDigitArray->GetNDigits() << " digits" << std::endl;
	if(fTotalLikelihood != NULL)
	{
    delete fTotalLikelihood;
    fTotalLikelihood = NULL;
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}
	else
	{
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}

	
  std::vector<std::pair<unsigned int, FitterParameterType::Type> > all1DProfile = fFitterPlots->GetAll1DSurfaceKeys();
  std::cout << "Size of 1D surface vector = " << all1DProfile.size() << std::endl;
	for(std::vector<std::pair<unsigned int, FitterParameterType::Type> >::iterator itr1D = all1DProfile.begin() ;
		itr1D != all1DProfile.end(); ++itr1D)
	{
		Make1DSurface(*itr1D);
	}


	std::vector<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > > all2DProfile = fFitterPlots->GetAll2DSurfaceKeys();
	for(std::vector<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > >::iterator itr2D = all2DProfile.begin() ;
			itr2D != all2DProfile.end(); ++itr2D)
	{
			Make2DSurface(*itr2D);
	}

  DeleteParameterArrays();
	return;
}

void WCSimLikelihoodFitter::ResetEvent() {
    fStatus = -999;
    fRootEvent = NULL;
    fLikelihoodDigitArray = NULL; // WCSimInterface cleans up after itself so no need to delete these here
    if(fTotalLikelihood != NULL)
    {
      delete fTotalLikelihood;
      fTotalLikelihood = NULL;
    }

    fStatus = -999;
    fTypes.clear();

    fParMap[1] = 8; // The number of fit parameters for n tracks, plus 1 for nTracks itself (fixed)
    fParMap[2] = 11;
    fMinimum  = 0;
    fSeedMap.clear();
    fIsFirstCall = true;
}

void WCSimLikelihoodFitter::FillPlots() {
	if(fFitterPlots != NULL)
	{
		fFitterPlots->FillPlots(fBestFit);
		fFitterPlots->FillRecoMinusTrue(fBestFit, fTrueLikelihoodTracks);
	}
	return;
}

void WCSimLikelihoodFitter::FillTree() {
	if(fFitterTree != NULL)
	{
		std::vector<Bool_t> bestFitEscapes;
		std::vector<Bool_t> trueTrackEscapes;
		for(unsigned int i = 0; i < fBestFit.size(); ++i)
		{
			bestFitEscapes.push_back(this->GetFitTrackEscapes(i));
		}
		for(unsigned int j = 0; j < fTrueLikelihoodTracks->size(); ++j)
		{
			trueTrackEscapes.push_back(this->GetTrueTrackEscapes(j));
		}
		fFitterTree->Fill(fEvent, fBestFit, bestFitEscapes, *fTrueLikelihoodTracks, trueTrackEscapes, fMinimum);
	}
	return;
}

void WCSimLikelihoodFitter::SetFitterTree(WCSimFitterTree* fitterTree) {
  std::cout << "Setting fFitterTree" << std::endl;
  fFitterTree = fitterTree;
}

Bool_t WCSimLikelihoodFitter::GetTrueTrackEscapes(unsigned int iTrack) const{
	WCSimLikelihoodTrack * track = fTrueLikelihoodTracks->at(iTrack);
	return GetTrackEscapes(track);
}

Bool_t WCSimLikelihoodFitter::GetFitTrackEscapes( unsigned int iTrack) const{
	WCSimLikelihoodTrack track = (fBestFit.at(iTrack));
	return GetTrackEscapes(&track);
}

Bool_t WCSimLikelihoodFitter::GetTrackEscapes(WCSimLikelihoodTrack * track) const{
	Double_t distanceToEdge = WCSimGeometry::Instance()->ForwardProjectionToEdge(track->GetX(), track->GetY(), track->GetZ(),
																	 	 	 	 track->GetDirX(), track->GetDirY(), track->GetDirZ());
  Double_t stoppingDistance = 0.0;
  if( track->GetType() != WCSimLikelihoodTrack::Unknown )
  {
	  WCSimEmissionProfiles ep(track);
	  stoppingDistance = ep.GetStoppingDistance(track);
  }
  else
  {
    std::cerr << "Warning, track was of unknown type" << std::endl;
    std::cerr << "Will assume it was contained" << std::endl; 
  }
	return (stoppingDistance > distanceToEdge);
}
