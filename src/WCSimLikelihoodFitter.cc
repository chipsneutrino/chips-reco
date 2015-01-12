#include "WCSimChargePredictor.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimFitterInterface.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimReco.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"

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
  for(UInt_t i = 0; i < nPars; ++i)
  {
	  if( fFixed[i] || fIsEnergy[i])
	  {
	  	min->SetFixedVariable(i, fNames[i], fStartVals[i]);
	  }
	  else
	  {
		  min->SetLimitedVariable(i, fNames[i], fStartVals[i], fStepSizes[i], fMinVals[i], fMaxVals[i]);
    }
  }

  // Perform the minimization
  min->Minimize();
  fMinimum = min->MinValue();
  fStatus = min->Status();

  
  // Get and print the fit results
  const Double_t * outPar = min->X();

  // Now we go again but free up the energy:
  for(UInt_t i = 0; i < nPars; ++i)
  {
    if(fIsEnergy[i] && !fFixed[i])
    {
      min->SetLimitedVariable(i, fNames[i], fStartVals[i], fStepSizes[i], fMinVals[i], fMaxVals[i]);
    }
    else
    {
      min->SetFixedVariable(i, fNames[i], outPar[i]);
    }
  }

  //  // Fix the vertex
  //  min->SetFixedVariable(1, parName[1], outPar[1]);
  //  min->SetFixedVariable(2, parName[2], outPar[2]);
  //  min->SetFixedVariable(3, parName[3], outPar[3]);
  //  min->SetLimitedVariable(5, parName[5], par[5], stepSize[5], maxVal[5], minVal[5]);
  //  min->SetLimitedVariable(6, parName[6], par[6], stepSize[6], maxVal[6], minVal[6]);
  //  for(UInt_t j = 0; j < nPars; ++j){std::cout << j << "   " << parName[j] << "   " << par[j] << "   " << minVal[j] << "   " << maxVal[j] << std::endl;}
  min->Minimize();


  // Get and print the fit results
  const Double_t * outPar2 = min->X();
  std::cout << "Best fit track: " << std::endl;
  RescaleParams(outPar2[0],  outPar2[1],  outPar2[2],  outPar2[3],  outPar2[4],  outPar2[5],  outPar2[6], fType).Print();

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
	  WCSimLikelihoodTrack track(
								  outPar2[fTrackParMap[trackParX]],
								  outPar2[fTrackParMap[trackParY]],
								  outPar2[fTrackParMap[trackParZ]],
								  outPar2[fTrackParMap[trackParT]],
								  outPar2[fTrackParMap[trackParTh]],
								  outPar2[fTrackParMap[trackParPhi]],
								  outPar2[fTrackParMap[trackParE]],
								  fType
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
  UInt_t nTracks = WCSimFitterConfig::Instance()->GetNumTracks();
  std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes();

  // If we've fixed some track parameters together, then our array of fit parameters won't just be of size
  // n tracks * m parameters per track, and we need to work out which entry corresponds to which parameter
  // of which track(s) - this is probably expensive, so we'll do it once and look it up in a map thereafter
  if( fIsFirstCall )
  { 
    SetTrackParameterMap();
  }

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
								  fType
								  );
	  tracksToFit.push_back(track);
  }



  if(fIsFirstCall || 1) // TEMP: Print it out all the time for now
  {
    std::cout << "Tracks used for first call:" << std::endl;
    std::vector<WCSimLikelihoodTrack>::iterator itr = tracksToFit.begin();
    for( ; itr < tracksToFit.end(); ++itr)
    {
      (*itr).Print();
    }
  }

  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(tracksToFit);
  minus2LnL = fTotalLikelihood->Calc2LnL();

  fIsFirstCall = false;
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
void WCSimLikelihoodFitter::SeedParams( WCSimReco * myReco )
{

   WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
   myReco->Run(recoEvent);
   fSeedVtxX = recoEvent->GetVtxX();
   fSeedVtxY = recoEvent->GetVtxY();
   fSeedVtxZ = recoEvent->GetVtxZ();

   Double_t recoDirX = recoEvent->GetDirX();
   Double_t recoDirY = recoEvent->GetDirY();
   Double_t recoDirZ = recoEvent->GetDirZ();
  

   fSeedTheta   = TMath::ACos(recoDirZ);
   if(recoDirY != 0) fSeedPhi = TMath::ATan2(recoDirY,recoDirX); // Ensure range is -pi to pi
   else fSeedPhi = (recoDirX < 0.0)? 0.5*TMath::Pi() : -0.5*TMath::Pi();

   //   fSeedTheta = 0.5*TMath::Pi();
   //   fSeedPhi = 0.25*TMath::Pi();

  // std::cout << "Seed dir = (vx, vy, vz) = ( " << recoDirX << "," << recoDirY << "," << recoDirZ << ")" << std::endl
  std::cout << "-> theta = " << fSeedTheta << "    phi = " << fSeedPhi << std::endl
            << "Seed position = (x,y,z) = ( " << fSeedVtxX << "," << fSeedVtxY << "," <<fSeedVtxZ << ")" << std::endl;
}

WCSimLikelihoodTrack WCSimLikelihoodFitter::GetSeedParams()
{
  WCSimLikelihoodTrack seedTrack(fSeedVtxX, fSeedVtxY, fSeedVtxZ, fSeedTime, fSeedTheta, fSeedPhi, fSeedEnergy, WCSimLikelihoodTrack::MuonLike);
  return seedTrack;
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

	Minimize2LnL();

	if(WCSimFitterInterface::Instance()->GetMakeFits() ) {
		WCSimInterface::Instance()->BuildEvent(iEvent);
    std::cout << "Getting true likelihood tracks 2" << std::endl;
		fTrueLikelihoodTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
		std::cout << "Fill plots now 3" << std::endl;
    FillPlots();
    std::cout << "Plots filled" << std::endl;
	}

	return;


}

void WCSimLikelihoodFitter::CreateParameterArrays() {
	Int_t numParams = WCSimFitterConfig::Instance()->GetNumIndependentParameters();

	fStartVals = new Double_t[numParams];
	fMinVals   = new Double_t[numParams];
	fMaxVals   = new Double_t[numParams];
	fStepSizes = new Double_t[numParams];
	fFixed	   = new Bool_t[numParams];
	fIsEnergy  = new Bool_t[numParams];
	fNames 	   = new std::string[numParams];
}

void WCSimLikelihoodFitter::SetParameterArrays() {
	Int_t iParam = 0;
	UInt_t numTracks = WCSimFitterConfig::Instance()->GetNumTracks();
	for(UInt_t jTrack = 0; jTrack < numTracks; ++jTrack)
	{
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
				fStepSizes[iParam] = (WCSimFitterConfig::Instance()->GetParMax(jTrack,name) - WCSimFitterConfig::Instance()->GetParMin(jTrack, name)) / 50.0;
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

		for(unsigned int iBin = 0; iBin < fFitterPlots->GetNumSurfaceBins(); ++iBin)
		{
			double stepVal = min + (max - min) * iBin / (double)fFitterPlots->GetNumSurfaceBins();
			x[arrIndex] = stepVal;
			fFitterPlots->Fill1DProfile(trackPar, stepVal, WrapFunc(x));
		}
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
	  		fFitterPlots->Fill2DProfile(trackPar, stepValX, stepValY, WrapFunc(x));
	  	}
	  }
  }

  return;
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

	for(UInt_t iEvent = 0; iEvent < WCSimFitterConfig::Instance()->GetNumEventsToFit() ; ++iEvent)
	{
		SetParameterArrays();
		FitEventNumber(iEvent);
	}

	DeleteParameterArrays();
}

void WCSimLikelihoodFitter::RunSurfaces() {
	
  CreateParameterArrays();
  SetParameterArrays();
  fRootEvent = WCSimInterface::Instance()->GetWCSimEvent(0);
	fLikelihoodDigitArray = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(0);
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
    fType = WCSimLikelihoodTrack::MuonLike;
    fParMap[1] = 8; // The number of fit parameters for n tracks, plus 1 for nTracks itself (fixed)
    fParMap[2] = 11;
    fMinimum  = 0;
    fSeedVtxX = 0.6;
    fSeedVtxY = 0.6;
    fSeedVtxZ = 0.6;
    fSeedTheta = 0.1;
    fSeedPhi = 0.5;
    fSeedTime = -40;
    fSeedEnergy = 1500;
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
