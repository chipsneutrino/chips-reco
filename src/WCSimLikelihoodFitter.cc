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
#include "WCSimChargeLikelihood.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"
#include "TClonesArray.h"
#include "TCollection.h"
#include "TMath.h"
#include "TMinuit.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

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
WCSimLikelihoodFitter::WCSimLikelihoodFitter(WCSimRootEvent * myRootEvent)
{
    fStatus = -999;
    fRootEvent = myRootEvent;
    fLikelihoodDigitArray = new WCSimLikelihoodDigitArray(fRootEvent);
    fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
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



WCSimLikelihoodFitter::~WCSimLikelihoodFitter()
{
  delete fLikelihoodDigitArray;
  delete fTotalLikelihood;
}



Int_t WCSimLikelihoodFitter::GetNPars(Int_t nTracks)
{
  // Do we know how to fit this number of tracks?
  Int_t nPars = -1;
  std::map<Int_t, Int_t>::const_iterator mapItr = fParMap.find(nTracks);
  if(mapItr == fParMap.end()) std::cerr << "WCSimLikelihoodFitter::GetNPars *** Error: could not look up the number of parameters for " << nTracks << " tracks" << std::endl;
  else nPars = fParMap[nTracks];

  return nPars;
}



void WCSimLikelihoodFitter::Minimize2LnL(Int_t nTracks)
{
   // Set up the minimizer
   ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

   // Alternatively: use a different algorithm to check the minimizer works
   // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "GSLSimAn");

   min->SetMaxFunctionCalls(9999999);
   min->SetMaxIterations(9);
   min->SetPrintLevel(3);
   min->SetTolerance(200);
   min->SetStrategy(2);
   std::cout << " Tolerance = " << min->Tolerance() << std::endl;

  // Convert nTracks into number of parameters
  const Int_t nPars = this->GetNPars(nTracks);
  
  // Tell the minimizer the function to minimize
  // We have to wrap it in this functor because it's a member function of this class
  ROOT::Math::Functor func(this,&WCSimLikelihoodFitter::WrapFunc, nPars);
  
  Double_t * par        = new Double_t[nPars];            // the start values
  Double_t * stepSize   = new Double_t[nPars];            // step sizes 
  Double_t * minVal     = new Double_t[nPars];            // minimum bound on parameter 
  Double_t * maxVal     = new Double_t[nPars];            // maximum bound on parameter
  std::cout << "nPars = " << nPars << std::endl;
  std::string * parName = new std::string[nPars];

  // Set parameter start values to the seed output
  if(nTracks == 1 || nTracks ==2)
  {
    
    // par[0] = 80.;    // x
    // par[1] = -200.;    // y
    // par[2] = 320.;    // z
    // par[3] = 0.;    // t
    // par[4] = 0.2*TMath::Pi();    // theta
    // par[5] = 0.05*TMath::Pi();    // phi
    // par[6] = 1500;  // energy

    par[0] = nTracks;    // x
    par[1] = fSeedVtxX;    // x
    par[2] = fSeedVtxY;    // y
    par[3] = fSeedVtxZ;    // z
    par[4] = fSeedTime;    // t
    par[5] = fSeedTheta;    // theta
    par[6] = fSeedPhi;    // phi
    par[7] = fSeedEnergy;  // energy

    minVal[0] = nTracks;
    minVal[1] = -1900.;
    minVal[2] = -1900.;
    minVal[3] = -940.;
    minVal[4] = -250.;
    minVal[5] = 0. * TMath::Pi();
    minVal[6] = -1.0 * TMath::Pi();
    minVal[7] = 1500;

    maxVal[0] = nTracks;
    maxVal[1] = 1900.;
    maxVal[2] = 1900.;
    maxVal[3] = 940.;
    maxVal[4] = 100.;
    maxVal[5] = TMath::Pi();
    maxVal[6] = 1.0 * TMath::Pi();;
    maxVal[7] = 1500;

    stepSize[0] = 0.;
    stepSize[1] = 50.;
    stepSize[2] = 50.;
    stepSize[3] = 50.;
    stepSize[4] = 10.;
    stepSize[5] = 0.02;
    stepSize[6] = 0.02;
    stepSize[7] = 1500;

    parName[0] = "nTracks";
    parName[1] = "vtxX";
    parName[2] = "vtxY";
    parName[3] = "vtxZ";
    parName[4] = "vtxT";
    parName[5] = "theta";
    parName[6] = "phi";
    parName[7] = "energy";
  }
  if(nTracks == 2)
  {
    par[8]  = 0.;    // second track theta
    par[9]  = 0.;    // second track phi
    par[10] = 1500;  // second track energy

    minVal[8]  = 0.;
    minVal[9]  = 0.;
    minVal[10] = 1500.;

    maxVal[8]  = TMath::Pi();
    maxVal[9]  = 2.0 * TMath::Pi();;
    maxVal[10] = 1500;

    stepSize[8]  = 0.05;
    stepSize[9]  = 0.05;
    stepSize[10] = 1500;

    parName[8]  = "theta_2";
    parName[9]  = "phi_2";
    parName[10] = "energy_2";
  }

 /* 
    // If we test the using the simulated annealing algorithm to do the minimizing
    // it converges better if we vary parameters between 0 and 1 and scale them up later
    par[0] = nTracks;    // x
    par[1] = fSeedVtxX;    // x
    par[2] = fSeedVtxY;    // y
    par[3] = fSeedVtxZ;    // z
    par[4] = fSeedTime;    // t
    par[5] = fSeedTheta;    // theta
    par[6] = fSeedPhi;    // phi
    par[7] = fSeedEnergy;  // energy

    minVal[0] = nTracks;
    minVal[1] = 0.0; // -1900.;
    minVal[2] = 0.0; // -1900.;
    minVal[3] = 0.0; // -940.;
    minVal[4] = 0.0; // 0.;
    minVal[5] = 0.0; // 0. * TMath::Pi();
    minVal[6] = 0.0; // -1.0 * TMath::Pi();
    minVal[7] = 1.0; // 1500;

    maxVal[0] = nTracks;
    maxVal[1] = 1.0; // 1900.;
    maxVal[2] = 1.0; // 1900.;
    maxVal[3] = 1.0; // 940.;
    maxVal[4] = 1.0; // 0.;
    maxVal[5] = 1.0; // TMath::Pi();
    maxVal[6] = 1.0; // 1.0 * TMath::Pi();;
    maxVal[7] = 1.0; // 1500;

    stepSize[0] = 0.25; // 0.;
    stepSize[1] = 0.25; // 50.;
    stepSize[2] = 0.25; // 50.;
    stepSize[3] = 0.25; // 50.;
    stepSize[4] = 0.25; // 1.;
    stepSize[5] = 0.25; // 0.02;
    stepSize[6] = 0.25; // 0.02;
    stepSize[7] = 1; // 1500;

    parName[0] = "nTracks";
    parName[1] = "vtxX";
    parName[2] = "vtxY";
    parName[3] = "vtxZ";
    parName[4] = "vtxT";
    parName[5] = "theta";
    parName[6] = "phi";
    parName[7] = "energy";
  }
  if(nTracks == 2)
  {
    par[8]  = 0.;    // second track theta
    par[9]  = 0.;    // second track phi
    par[10] = 1;  // second track energy

    minVal[8]  = 0.;
    minVal[9]  = 0.;
    minVal[10] = 1.;

    maxVal[8]  = TMath::Pi();
    maxVal[9]  = 2.0 * TMath::Pi();;
    maxVal[10] = 1.;

    stepSize[8]  = 0.05;
    stepSize[9]  = 0.05;
    stepSize[10] = 1;

    parName[8]  = "theta_2";
    parName[9]  = "phi_2";
    parName[10] = "energy_2";
  }
  */

  // Print the parameters we're using
  for(UInt_t j = 0; j < nPars; ++j)
  {
      std::cout << j << "   " << parName[j] << "   " << par[j] << "   " << minVal[j] << "   " << maxVal[j] << std::endl;
  }

  // Tell the minimizer the functor and variables to consider
  min->SetFunction(func);
  // min->SetFixedVariable(0,"nTracks",nTracks);  // Lets the wrapper construct the right likelihood

  // Set the fixed and free variables
  for(Int_t i = 0; i < nPars; ++i)
  {
    min->SetVariable(i,parName[i], par[i], stepSize[i]);
    if((1 == i || 2 == i || 3 == i || i ==5 || i==6 ) && (maxVal[i] != minVal[i])){ min->SetLimitedVariable(i,parName[i], par[i], stepSize[i], minVal[i], maxVal[i]);}
    if(7 == i) { min->SetFixedVariable(i, parName[i], par[i]); }
    std::cout << i << "   " << parName[i] << "   " << par[i] << "   " << minVal[i] << "   " << maxVal[i] << std::endl;
  }

  //  // Fix the direction
  //  min->SetFixedVariable(5, parName[5], par[5]);
  //  min->SetFixedVariable(6, parName[6], par[6]);


  // Perform the minimization
  min->Minimize();
  fMinimum = min->MinValue();
  fStatus = min->Status();

  
  // Get and print the fit results
  const Double_t * outPar = min->X();

  //  // Fix the vertex
  //  min->SetFixedVariable(1, parName[1], outPar[1]);
  //  min->SetFixedVariable(2, parName[2], outPar[2]);
  //  min->SetFixedVariable(3, parName[3], outPar[3]);
  //  min->SetLimitedVariable(5, parName[5], par[5], stepSize[5], maxVal[5], minVal[5]);
  //  min->SetLimitedVariable(6, parName[6], par[6], stepSize[6], maxVal[6], minVal[6]);
  //  for(UInt_t j = 0; j < nPars; ++j){std::cout << j << "   " << parName[j] << "   " << par[j] << "   " << minVal[j] << "   " << maxVal[j] << std::endl;}
  //  min->Minimize();


  // Get and print the fit results
  const Double_t * outPar2 = min->X();
  const Double_t * err    = min->Errors();
  std::cout << "Best fit track: " << std::endl;
  RescaleParams(outPar2[1],  outPar2[2],  outPar2[3],  outPar2[4],  outPar2[5],  outPar2[6],  outPar2[7], fType).Print();

  // Now save the best-fit tracks
  fBestFit.clear();
  switch(nTracks)
  {
    case 2:
    {  
      WCSimLikelihoodTrack tempTrack2 = RescaleParams( outPar[1], outPar[2], outPar[3], outPar[4], outPar[8], outPar[9], outPar[10], fType );
      fBestFit.push_back(tempTrack2);
    }
    case 1:
    {
      WCSimLikelihoodTrack tempTrack1 = RescaleParams( outPar[1], outPar[2], outPar[3], outPar[4], outPar[5], outPar[6], outPar[7], fType );
      fBestFit.push_back(tempTrack1);
      std::reverse(fBestFit.begin(), fBestFit.end());
    }
    default:
      break;
  }
  return;
}

/**
 * Wrapper: constructs the correct number of track objects
 * then works out the total likelihood
 */
Double_t WCSimLikelihoodFitter::WrapFunc(const Double_t * x)
{
  std::vector<WCSimLikelihoodTrack> trackVec;
  WCSimLikelihoodTrack myTrack = RescaleParams(x[1],x[2],x[3],x[4],x[5],x[6],x[7], fType);
  trackVec.push_back(myTrack);

  // Check we can deal with the number of tracks provided
  Int_t nTracks = static_cast<Int_t>(x[0]);
  switch(nTracks)
  {
    case 1:
      break;
    case 2:
      break;
    default:
      std::cerr << "WCSimLikelihoodFitter::WrapFunc *** Error: can only fit for 1 or 2 tracks" << std::endl;
      return 0;
  }
  if(x[0] == 2)
  {
    WCSimLikelihoodTrack myTrack2 = RescaleParams(x[1],x[2],x[3],x[4],x[8],x[9],x[10], fType);
    trackVec.push_back(myTrack2);
  }

  if(fIsFirstCall || 1) // TEMP: Print it out all the time for now
  {
    std::cout << "Tracks used for first call:" << std::endl;
    std::vector<WCSimLikelihoodTrack>::iterator itr = trackVec.begin();
    for( ; itr < trackVec.end(); ++itr)
    {
      (*itr).Print();
    }
    fIsFirstCall = false;
  }

  std::cout << " ...  Done" << std::endl;
  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(trackVec);
  minus2LnL = fTotalLikelihood->Calc2LnL();
  return minus2LnL;
}



/*
 * Integrate the time likelihood here:
double WCSimLikelihoodFitter::Time2LnL( )
{
    WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood( fLikelihoodDigitArray, fTrack );
    WCSimTimeLikelihood * myTimeLikelihood = new WCSimTimeLikelihood( fLikelihoodDigitArray, fTrack, myChargeLikelihood );
    Double_t time2LnL = myTimeLikelihood->Calc2LnL();
    delete myTimeLikelihood;
    delete myChargeLikelihood;

    //std::cout << "Returning -2* time LnL" << std::endl;
    std::cout << "-2 * Time2LnL = " << time2LnL << std::endl;
    return time2LnL;
}
*/

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
