<<<<<<< HEAD
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

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodFitter)
#endif

WCSimLikelihoodFitter::WCSimLikelihoodFitter(WCSimRootEvent * myRootEvent)
{
    fRootEvent = myRootEvent;
    fLikelihoodDigitArray = new WCSimLikelihoodDigitArray(fRootEvent);
    fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
    fType = WCSimLikelihoodTrack::MuonLike;
    fParMap[1] = 8; // The number of fit parameters for n tracks, plus 1 for nTracks itself (fixed)
    fParMap[2] = 11;
    fMinimum  = 0;
    fSeedVtxX = 0.;
    fSeedVtxY = 0.;
    fSeedVtxZ = 0.;
    fSeedTheta = 0.0;
    fSeedPhi = 0.0;
    fIsFirstCall = true;
}

WCSimLikelihoodFitter::~WCSimLikelihoodFitter()
{
  delete fLikelihoodDigitArray;
  delete fTotalLikelihood;
}

Int_t WCSimLikelihoodFitter::GetNPars(Int_t nTracks)
{
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
  min->SetMaxFunctionCalls(9999999);
  min->SetMaxIterations(9);
  min->SetPrintLevel(3);
  min->SetTolerance(100);
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
    par[4] = 0.0;    // t
    par[5] = fSeedTheta;    // theta
    par[6] = fSeedPhi;    // phi
    par[7] = 1500;  // energy

    minVal[0] = nTracks;
    minVal[1] = -1900.;
    minVal[2] = -1900.;
    minVal[3] = -940.;
    minVal[4] = 0.;
    minVal[5] = 0. * TMath::Pi();
    minVal[6] = -1.0 * TMath::Pi();
    minVal[7] = 1500;

    maxVal[0] = nTracks;
    maxVal[1] = 1900.;
    maxVal[2] = 1900.;
    maxVal[3] = 940.;
    maxVal[4] = 0.;
    maxVal[5] = TMath::Pi();
    maxVal[6] = 1.0 * TMath::Pi();;
    maxVal[7] = 1500;

    stepSize[0] = 0.;
    stepSize[1] = 50.;
    stepSize[2] = 50.;
    stepSize[3] = 50.;
    stepSize[4] = 1.;
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
    par[8] = 0.;    // second track theta
    par[9] = 0.;    // second track phi
    par[10] = 1500;  // second track energy

    minVal[8] = 0.;
    minVal[9] = 0.;
    minVal[10] = 1500;

    maxVal[8] = TMath::Pi();
    maxVal[9] = 2.0 * TMath::Pi();;
    maxVal[10] = 1500;

    stepSize[8] = 0.05;
    stepSize[9] = 0.05;
    stepSize[10] = 1500;

    parName[8] = "theta_2";
    parName[9] = "phi_2";
    parName[10] = "energy_2";
  }
  for(UInt_t j = 0; j < nPars; ++j)
  {std::cout << j << "   " << parName[j] << "   " << par[j] << "   " << minVal[j] << "   " << maxVal[j] << std::endl;}

  // Tell the minimizer the functor and variables to consider
  min->SetFunction(func);
//  min->SetFixedVariable(0,"nTracks",nTracks);  // Lets the wrapper construct the right likelihood
  for(Int_t i = 0; i < nPars; ++i)
  {
    min->SetVariable(i,parName[i], par[i], stepSize[i]);
    if((1 == i || 2 == i || 3 == i || 5 == i || 6 == i ) && (maxVal[i] != minVal[i])){ min->SetLimitedVariable(i,parName[i], par[i], stepSize[i], minVal[i], maxVal[i]);}
    if(0 == i || 4 == i || 7 == i) min->SetFixedVariable(i, parName[i], par[i]);
    std::cout << i << "   " << parName[i] << "   " << par[i] << "   " << minVal[i] << "   " << maxVal[i] << std::endl;
  }

//  // Fix the direction
//  min->SetFixedVariable(5, parName[5], par[5]);
//  min->SetFixedVariable(6, parName[6], par[6]);
//
  // Perform the minimization
  min->Minimize();
  fMinimum = min->MinValue();
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
  std::cout << "Best fit track: " << std::endl
            << "(x,y,z,t)   = (" << outPar2[1] << "," << outPar2[2] << "," << outPar2[3] << "," << outPar2[4] << ")" << std::endl
            << "(theta,phi) = (" << outPar2[5] << "," << outPar2[6] << ")" << std::endl
            << "Energy      =  " << outPar2[7] << std::endl;
  std::cout << "Best fit track, errors: " << std::endl
            << "(x,y,z,t)   = (" << err[1] << "," << err[2] << "," << err[3] << "," << err[4] << ")" << std::endl
            << "(theta,phi) = (" << err[5] << "," << err[6] << ")" << std::endl
            << "Energy      =  " << err[7] << std::endl;

  fBestFit.clear();
  switch(nTracks)
  {
    case 2:
    {  
      WCSimLikelihoodTrack tempTrack2( outPar[1], outPar[2], outPar[3], outPar[4], outPar[8], outPar[9], outPar[10], fType );
      fBestFit.push_back(tempTrack2);
    }
    case 1:
    {
      WCSimLikelihoodTrack tempTrack1( outPar[1], outPar[2], outPar[3], outPar[4], outPar[5], outPar[6], outPar[7], fType );
      fBestFit.push_back(tempTrack1);
    }
    default:
      break;
  }
  return;
}

////////////////////////////////////////////////////////////
// Wrapper: constructs the correct number of track objects
// then works out the total likelihood
////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodFitter::WrapFunc(const Double_t * x)
{
  std::vector<WCSimLikelihoodTrack> trackVec;
  WCSimLikelihoodTrack myTrack(x[1],x[2],x[3],x[4],x[5],x[6],x[7], fType);
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
    WCSimLikelihoodTrack myTrack2(x[1],x[2],x[3],x[4],x[8],x[9],x[10], fType);
    trackVec.push_back(myTrack2);
  }

  if(fIsFirstCall || 1)
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
   if(recoDirY != 0) fSeedPhi = TMath::ATan2(recoDirY,recoDirX); // range -pi to pi
   else fSeedPhi = (recoDirX < 0.0)? 0.5*TMath::Pi() : -0.5*TMath::Pi();

  // std::cout << "Seed dir = (vx, vy, vz) = ( " << recoDirX << "," << recoDirY << "," << recoDirZ << ")" << std::endl
  std::cout << "-> theta = " << fSeedTheta << "    phi = " << fSeedPhi << std::endl
            << "Seed position = (x,y,z) = ( " << fSeedVtxX << "," << fSeedVtxY << "," <<fSeedVtxZ << ")" << std::endl;

}

WCSimLikelihoodTrack WCSimLikelihoodFitter::GetSeedParams()
{
  WCSimLikelihoodTrack seedTrack(fSeedVtxX, fSeedVtxY, fSeedVtxZ, 0, fSeedTheta, fSeedPhi, 1500, WCSimLikelihoodTrack::MuonLike);
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
=======
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimChargeLikelihood.hh"
#include "WCSimTimeLikelihood.hh"
#include "TClonesArray.h"
#include "TCollection.h"
#include "TMath.h"

WCSimLikelihoodFitter::WCSimLikelihoodFitter(WCSimRootEvent * myRootEvent)
{
    fRootEvent = myRootEvent;
    fLikelihoodDigitArray = new WCSimLikelihoodDigitArray(fRootEvent);
//    fRecoDigitList = *(myRecoEvent->GetDigitList());
    //std::cout << "Calculated LnL()" << std::endl;
    //ctor
}

WCSimLikelihoodFitter::~WCSimLikelihoodFitter()
{
}

void WCSimLikelihoodFitter::Minimize2LnL()
{
    //TODO: should set up some WCSimLikelihoodTrack!
    this->Calc2LnL();
    // Then minimize it
    return ;
}

double WCSimLikelihoodFitter::Calc2LnL()  
{
  //TODO: check if there is a defined track?  //mp
  //  like: if(!fTrack) ShoutOnUser();
  return (this->Charge2LnL() + this->Time2LnL());
}

double WCSimLikelihoodFitter::Calc2LnL(WCSimLikelihoodTrack * myTrack)
{
  return (this->Charge2LnL( myTrack ) + this->Time2LnL());
}

double WCSimLikelihoodFitter::Charge2LnL( WCSimLikelihoodTrack * myTrack )
{
  fTrack = myTrack;
  return (this->Charge2LnL());
}


double WCSimLikelihoodFitter::Charge2LnL( )
{
    if(fTrack == NULL)
    {
      std::cerr << "WCSimLikelihoodFitter::Charge2LnL() - Error: fTrack = NULL" << std::endl;
      exit(EXIT_FAILURE);
    }
    double Charge2LnL = 1.0;

    //std::cout << "Calculating charge LnL" << std::endl;

    WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood( fLikelihoodDigitArray, fTrack );
    Charge2LnL = myChargeLikelihood->Calc2LnL();
    delete myChargeLikelihood;
  //        double myProb = TMath::Poisson( Q, mu );
  //        ChargeLnL += TMath::Log( myProb );
      //}
   // }
    //std::cout << "Returning  -2 * charge LnL" << std::endl;
    std::cout << "-2 * ChargeLnL = " << Charge2LnL << std::endl;
    return Charge2LnL;
}

double WCSimLikelihoodFitter::CalcMuDirect()
{   
    double muDirect = 0.0;
    return muDirect;
}

double WCSimLikelihoodFitter::CalcMuIndirect()
{
    double muIndirect = 1.0;
    return muIndirect;
}

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
>>>>>>> 1d76eafaf5c2e96d8b78405738da852825d0e829
