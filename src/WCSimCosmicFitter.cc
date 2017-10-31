#include "WCSimAnalysisConfig.hh"
#include "WCSimChargePredictor.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimOutputTree.hh"
#include "WCSimFitterInterface.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimCosmicFitter.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimReco.hh"
#include "WCSimRecoSeed.hh"
#include "WCSimCosmicSeed.hh"
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
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TParticlePDG.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include "TCanvas.h"
#include "TPolyMarker.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimCosmicFitter)
#endif

  //bool WCSimCosmicFitter::RingSort(const std::pair<WCSimRecoRing*,double> &a, const std::pair<WCSimRecoRing*,double> &b){
  //  return (a.first)->GetHeight() > (b.first)->GetHeight();
  //}

  /**
   * @todo Remove hardcoding of track type
   */
WCSimCosmicFitter::WCSimCosmicFitter(WCSimFitterConfig * config) : WCSimLikelihoodFitter(config)
{
//  fCalls = 0;
//  fFitterConfig = config;
//  fFitterPlots = NULL;
//  fFitterTree = 0x0;
//  fTotalLikelihood = NULL;
//  fRootEvent = NULL;
//  fLikelihoodDigitArray = NULL;
//  fTrueLikelihoodTracks = NULL;
  fSingleTrackToFit = 0; 

  ResetEvent();

}

WCSimCosmicFitter::~WCSimCosmicFitter()
{
}

void WCSimCosmicFitter::RunFits() {
  std::cout << "WCSimCosmicFitter::RunFits() " << std::endl;

  UInt_t firstEvent = fFitterConfig->GetFirstEventToFit();
  UInt_t numEventsToFit = fFitterConfig->GetNumEventsToFit();

  std::cout << firstEvent << ", " << numEventsToFit << std::endl;

  for(UInt_t iEvent = firstEvent; iEvent < firstEvent + numEventsToFit ; ++iEvent)
  {
    FitEventNumber(iEvent);
    if( fMinimum > 0 && !fFailed )
    {
      std::cout << "Filling successful event" << std::endl;
      FillPlots();
    }
    else
    {
      std::cout << "Failed for some reason: fMinimum = " << fMinimum << " and fFailed = " << fFailed << std::endl;
    }
    FillTree();
    fOutputTree->SaveTree();
    fFitterPlots->SavePlots();
  }
}

void WCSimCosmicFitter::FitEventNumber(Int_t iEvent) {
  fFailed = false;
  fIsFirstCall = true;
  SetEvent(iEvent);
  bool usingCharge = WCSimAnalysisConfig::Instance()->GetUseCharge();
  bool usingTime = WCSimAnalysisConfig::Instance()->GetUseTime();

  if(fFitterConfig->GetMakeFits() && CanFitEvent()) {
    std::cout << "Fitting event " << iEvent << std::endl;


    fCalls = 0;
    fFitterTrackParMap.Set();
    fFitterTrackParMap.Print(); 
    // Run seed
    SeedEvent();

    std::cout << "FITTER STAGE 0 - Fit the energy" << std::endl;
    FixVertex();
    FixTime();
    FixDirection();
    FreeEnergy();
    FitEnergy();

    std::cout << "FITTER STAGE 1 - Take seed, fit time (time-only)" << std::endl;
    if(usingTime)
    {
      WCSimAnalysisConfig::Instance()->SetUseTimeOnly();
      FreeTime();
      FixVertex();
      FixDirection();
      FixEnergy();
      FitTime();
      FreeVertex();
      FreeEnergy();
      if(usingCharge){ WCSimAnalysisConfig::Instance()->SetUseChargeAndTime(); }
    }
    else
    {
      std::cout << "  - Not using time - skipping this stage  - " << std::endl;
    }


    std::cout << "FITTER STAGE 2 - Take time-corrected seed, fit energy" << std::endl;
    FreeEnergy();
    FixVertex();
    FixDirection();
    FixTime();
    FitEnergy();
    FreeTime();
    FreeDirection();
    FreeVertex();

    if(    WCSimAnalysisConfig::Instance()->GetEqualiseChargeAndTime()
        && usingCharge 
        && usingTime)    
    {
      std::cout << "EQUALISING CHARGE AND TIME" << std::endl;
      std::vector<double> time2LnLVec = fTotalLikelihood->GetTime2LnLVector();
      std::vector<double> charge2LnLVec = fTotalLikelihood->GetCharge2LnLVector();
      double time2LnL = std::accumulate(time2LnLVec.begin(), time2LnLVec.end(), 0.0);
      double charge2LnL = std::accumulate(charge2LnLVec.begin(), charge2LnLVec.end(), 0.0);

      double scaleFactor = 1;
      if( charge2LnL > 0 && time2LnL > 0)
      {
        scaleFactor = charge2LnL / time2LnL;
      }
      fTotalLikelihood->SetTimeScaleFactor(scaleFactor);
    }

    std::cout << "FITTER STAGE 3 - Iteratively correct direction and fit along track" << std::endl;
    for(int times = 0; times < 2; ++times)
    {
      FixDirection();
      FitAlongTrack();
      FreeDirection();

      FixVertex();
      FixEnergy();
      FixTime();
      Fit();
      FreeVertex();
      FreeEnergy();
      FreeTime();
    }



    // Now freely fit the vertex and direction    
    std::cout << "FITTER STAGE 4 - Completely free vertex fit" << std::endl;
    FitVertex();

    // Finesse the final fit result along the track.  The time likelihood
    // tends to help perpendicular to the track direction, so we'll lock
    // this in and then use charge for a fit along the track
    if( usingCharge && usingTime )
    {
      // Fit energy and vertex position along track using charge
      WCSimAnalysisConfig::Instance()->SetUseChargeOnly();
      std::cout << "FITTER STAGE 5 - Adjust along the track using charge only" << std::endl;
    }
    else
    {
      std::cout << "FITTER STAGE 5 - Adjust along the track" << std::endl;
    }

    FixDirection();
    FixTime();
    FitAlongTrack();
    FreeDirection();
    FreeTime();

    // Then fit energy
    if(usingCharge && usingTime){ WCSimAnalysisConfig::Instance()->SetUseChargeAndTime(); }
    std::cout << "FITTER STAGE 6 - Final energy fit" << std::endl;
    FixVertex();
    FixDirection();
    FixTime();
    FitEnergy();
    FreeVertex();
    FreeDirection();
    FreeEnergy();

    // Finally fit the time
    std::cout << "FITTER STAGE 7 - Final timing fit" << std::endl;
    FixVertex();
    FixEnergy();
    FixDirection();
    FreeTime();
    FitTime();
    FreeVertex();
    FreeEnergy();
    FreeDirection();
  }

  fTrueLikelihoodTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
  fTotalLikelihood->SetTracks(fBestFit);
  fMinimum = fTotalLikelihood->Calc2LnL();
  fMinimumTimeComponent = fTotalLikelihood->GetLastTime2LnL();
  fMinimumChargeComponent = fTotalLikelihood->GetLastCharge2LnL();

  std::cout << "Fitted event number " << iEvent << std::endl;
  return;
}

void WCSimCosmicFitter::SeedEvent()
{

  // Run the old Hough transform reco
  WCSimRecoSeed * myReco = dynamic_cast<WCSimRecoSeed*>(WCSimRecoFactory::Instance()->MakeReco("seed")); // This calls new() - delete when done
  myReco->SetNumberOfTracks(fFitterConfig->GetNumTracks());
  // Tell the reco seed it needs to use the cosmic seed.
  myReco->SetCosmicFit(fFitterConfig->GetIsCosmicFit());
  WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
  std::vector<WCSimRecoEvent*> slicedEvents = myReco->RunSeed(recoEvent);

  // Also get the cosmic seed
  //  WCSimCosmicSeed cosmicSeed(recoEvent->GetVetoDigitList());
  //  cosmicSeed.CalcSeedVtxAndDir();
  WCSimCosmicSeed* cosmicSeed = myReco->GetCosmicSeed();

  // Make a vector of all of the available rings we have. 
  std::vector<WCSimRecoRing*> ringVec;
  std::vector<std::pair<WCSimRecoRing*,double> > otherRings;
  std::vector<double> ringTime;

  for(unsigned int e = 0; e < slicedEvents.size(); ++e){
    ringVec.push_back(slicedEvents[e]->GetRing(0));
    ringTime.push_back(slicedEvents[e]->GetVtxTime());
    for(int r = 1; r < slicedEvents[e]->GetNRings(); ++r){
      std::pair<WCSimRecoRing*,double> tempPair;
      tempPair = std::make_pair<WCSimRecoRing*,double>(slicedEvents[e]->GetRing(r),slicedEvents[e]->GetVtxTime());
      otherRings.push_back(tempPair);
    }
  }
  // Sort the otherRings vector by height and add to the main vector
  std::sort(otherRings.begin(),otherRings.end(),RingSort);
  for(unsigned int r = 0; r < otherRings.size(); ++r){
    ringVec.push_back(otherRings[r].first);
    ringTime.push_back(otherRings[r].second);
  }

  bool dedicatedCosmicSeed = cosmicSeed->GetSuccess();
  // If we only have one seed after slicing, use the dedicated cosmic seed if we have one.

  if(dedicatedCosmicSeed) std::cout << "COSMIC SEED SUCCESS with " << ringVec.size() << " rings found." << std::endl;

  // If we don't have two slices, rely on the veto to seed things.
  if(dedicatedCosmicSeed){
    if((fFitterTrackParMap.GetTrackType(0) == TrackType::MuonLike && 
          fFitterTrackParMap.GetTrackType(1) == TrackType::ElectronLike) || 
        (fFitterTrackParMap.GetTrackType(1) == TrackType::MuonLike && 
         fFitterTrackParMap.GetTrackType(0) == TrackType::ElectronLike)){

      unsigned int muonTrack = 0;
      unsigned int otherTrack = 1;
      if(fFitterTrackParMap.GetTrackType(1) == TrackType::MuonLike){
        muonTrack = 1;
        otherTrack = 0;
      }

      TVector3 cosmicVtx = cosmicSeed->GetFittedVtx(); 
      TVector3 cosmicDir = cosmicSeed->GetFittedDir();
      double cosmicT = cosmicSeed->GetFittedVtxT();

      std::cout << " USING COSMIC SEED " << cosmicVtx.X() << ", " << cosmicVtx.Y() << ", " << cosmicVtx.Z() << " :: "
        << cosmicDir.X() << ", " << cosmicDir.Y() << ", " << cosmicDir.Z() << std::endl;
      // If the seed happens to be outside the allowed region of the detector, move 
      // the track back along its direction until it re-enters
      double tempX = cosmicVtx.X();
      double tempY = cosmicVtx.Y();
      double tempZ = cosmicVtx.Z();
      double tempDX = cosmicDir.X();
      double tempDY = cosmicDir.Y();
      double tempDZ = cosmicDir.Z();
      if(IsOutsideAllowedRegion(muonTrack, cosmicVtx.X(), cosmicVtx.Y(), cosmicVtx.Z()) )
      {
        MoveBackInside(muonTrack, tempX, tempY, tempZ, tempDX, tempDY, tempDZ);
        cosmicVtx = TVector3(tempX,tempY,tempZ);
      }

      // Fill the muon track with the muon seed
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxX, cosmicVtx.X());     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxY, cosmicVtx.Y());     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxZ, cosmicVtx.Z());     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxT, cosmicT);     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kDirTh,cosmicDir.Theta());     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kDirPhi,cosmicDir.Phi());     
      // And the other track with what is left 

      TVector3 trkVtx = TVector3(ringVec[0]->GetVtxX(),ringVec[0]->GetVtxY(),ringVec[0]->GetVtxZ());
      TVector3 trkDir = TVector3(ringVec[0]->GetDirX(),ringVec[0]->GetDirY(),ringVec[0]->GetDirZ());
      // If the seed happens to be outside the allowed region of the detector, move 
      // the track back along its direction until it re-enters
      tempX = trkVtx.X();
      tempY = trkVtx.Y();
      tempZ = trkVtx.Z();
      tempDX = trkDir.X();
      tempDY = trkDir.Y();
      tempDZ = trkDir.Z();
      if(IsOutsideAllowedRegion(otherTrack, cosmicVtx.X(), cosmicVtx.Y(), cosmicVtx.Z()) )
      {
        MoveBackInside(otherTrack, tempX, tempY, tempZ, tempDX, tempDY, tempDZ);
        trkVtx = TVector3(tempX,tempY,tempZ);
      }

      // Finding a single seed based on a big overlapping mess isn't a good idea. Let's try just pointing
      // the beam particle in the beam direction.
      //TVector3 trkDir = TVector3(1,0,0);
      //TVector3 fakeVertex = TVector3(0,0,0);

      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxX, trkVtx.X());     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxY, trkVtx.Y());     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxZ, trkVtx.Z());     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxT, ringTime[0]);     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kDirTh, trkDir.Theta());     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kDirPhi, trkDir.Phi());

    }
  }
  else{
    // We only need to worry if we have an electron and a muon.
    bool takeAction = false;
    TrackType track0 = fFitterTrackParMap.GetTrackType(0);
    TrackType track1 = fFitterTrackParMap.GetTrackType(1);
    unsigned int muonTrack = 1;
    unsigned int otherTrack = 0;
    if(fFitterTrackParMap.GetTrackType(0) == TrackType::MuonLike){
      if(fFitterTrackParMap.GetTrackType(1) != TrackType::MuonLike){
        takeAction = true;
        muonTrack = 0;
        otherTrack = 1;
      }
    }
    else{
      if(fFitterTrackParMap.GetTrackType(1) == TrackType::MuonLike) takeAction = true;
    }
    // Make sure we actually have two seeds before doing this.
    if(ringVec.size() < 2){takeAction = false;}

    if(takeAction){
      // This means we have two events.
      std::vector<TVector3> seedVtx;
      seedVtx.push_back(TVector3(ringVec[0]->GetVtxX(),ringVec[0]->GetVtxY(),ringVec[0]->GetVtxZ()));
      seedVtx.push_back(TVector3(ringVec[1]->GetVtxX(),ringVec[1]->GetVtxY(),ringVec[1]->GetVtxZ()));

      std::vector<TVector3> seedDir;
      seedDir.push_back(TVector3(ringVec[0]->GetDirX(),ringVec[0]->GetDirY(),ringVec[0]->GetDirZ()));
      seedDir.push_back(TVector3(ringVec[1]->GetDirX(),ringVec[1]->GetDirY(),ringVec[1]->GetDirZ()));

      // Make an assumption then check if it is wrong. 
      unsigned int muonSeed = 0;
      unsigned int beamSeed = 1;

      bool switchSeeds = false;        

      // If the vertices are the same, see what we can do.
      if(seedVtx[0] == seedVtx[1]){
        if((seedDir[1].Z() < seedDir[0].Z()) && (seedDir[1].X() < seedDir[0].X())){
          switchSeeds = true;
        }
      }
      // Is the second track seed steep?
      else if(seedDir[1].Z() < -0.9){
        switchSeeds = true;
      }
      // Does the first seed look like a beam direction track?
      else if(seedDir[0].X() > 0.5){
        switchSeeds = true;
      }   
      // Is the second seed at large R?
      else if(seedVtx[1].Perp() > 0.8*WCSimGeometry::Instance()->GetCylRadius()){
        switchSeeds = true;
      }
      // Is the second seed at large Z?
      else if(fabs(seedVtx[1].Z()) > 0.8*(WCSimGeometry::Instance()->GetCylLength()/2.0)){
        switchSeeds = true;
      }

      if(switchSeeds){
        muonSeed = 1;
        beamSeed = 0;
      }

      std::cout << "Muon seed = " << muonSeed << ", hence beam seed = " << beamSeed << std::endl;

      std::vector<double> seedTheta;
      std::vector<double> seedPhi;

      for(unsigned int i = 0; i < 2; ++i){  
        seedTheta.push_back(TMath::ACos(seedDir[i].Z()));
        if(seedDir[i].Y() != 0.0){ seedPhi.push_back(TMath::ATan2(seedDir[i].Y(),seedDir[i].X())); }// Ensure range is -pi to pi
        else{ seedPhi.push_back((seedDir[i].X() < 0.0)? 0.5*TMath::Pi() : -0.5*TMath::Pi()); }
      }

      // Fill the muon track with the muon seed
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxX, ringVec[muonSeed]->GetVtxX());     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxY, ringVec[muonSeed]->GetVtxY());     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxZ, ringVec[muonSeed]->GetVtxZ());     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kVtxT, ringTime[muonSeed]);     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kDirTh, seedTheta[muonSeed]);     
      fFitterTrackParMap.SetCurrentValue(muonTrack, FitterParameterType::kDirPhi, seedPhi[muonSeed]);     
      // And the other track with what is left 
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxX, ringVec[beamSeed]->GetVtxX());     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxY, ringVec[beamSeed]->GetVtxY());     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxZ, ringVec[beamSeed]->GetVtxZ());     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kVtxT, ringTime[beamSeed]);     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kDirTh, seedTheta[beamSeed]);     
      fFitterTrackParMap.SetCurrentValue(otherTrack, FitterParameterType::kDirPhi, seedPhi[beamSeed]);     
    } // End take action if statement
    else{
      // If no special treatment is required, just set the tracks as in the default fitter.
      // No action occurs if fitting two muons, or if we don't have the required number of seeds.
      for(unsigned int iTrack = 0; iTrack < fFitterConfig->GetNumTracks(); ++iTrack)
      {
        double seedX, seedY, seedZ, seedT;
        double dirX, dirY, dirZ;
        double seedTheta, seedPhi;
        std::cout << "Track number " <<  iTrack << " compared to " << slicedEvents.size() << std::endl;
        if(iTrack < ringVec.size()){
          seedX = ringVec[iTrack]->GetVtxX();
          seedY = ringVec[iTrack]->GetVtxY();
          seedZ = ringVec[iTrack]->GetVtxZ();
          seedT = ringTime[iTrack];
          dirX = ringVec[iTrack]->GetDirX();
          dirY = ringVec[iTrack]->GetDirY();
          dirZ = ringVec[iTrack]->GetDirZ();
        }
        else{
          seedX = ringVec[0]->GetVtxX();
          seedY = ringVec[0]->GetVtxY();
          seedZ = ringVec[0]->GetVtxZ();
          seedT = ringTime[0];
          dirX = ringVec[0]->GetDirX();
          dirY = ringVec[0]->GetDirY();
          dirZ = ringVec[0]->GetDirZ();
        }
        // If the seed happens to be outside the allowed region of the detector, move 
        // the track back along its direction until it re-enters
        if(IsOutsideAllowedRegion(iTrack, seedX, seedY, seedZ))
        {
          MoveBackInside(iTrack, seedX, seedY, seedZ, dirX, dirY, dirZ);
        }

        seedTheta = TMath::ACos(dirZ);
        if(dirY != 0.0){ seedPhi = TMath::ATan2(dirY,dirX); }// Ensure range is -pi to pi
        else{ seedPhi = (dirX < 0.0)? 0.5*TMath::Pi() : -0.5*TMath::Pi(); }

        std::cout << "Track seed " << iTrack << ": " << seedX << ", " << seedY << ", " << seedZ << " :: "
          << dirX << ", " << dirY << ", " << dirZ << ", " << seedTheta << ", " << seedPhi << std::endl; 

        // Only see the parameters if they are not requested to be fixed:
        // Also check if any of the parameters are joined with any other tracks: we should
        // let the lower-numbered track win in this case, because it has the higher Hough peak
        WCSimFitterConfig * fitterConfig = fFitterConfig;
        if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxX)){
          if(fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxX) >= iTrack){
            fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxX, seedX);
          }
        }
        if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxY)){
          if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxY) >= iTrack){
            fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxY, seedY);
          }
        }
        if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxZ)){
          if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxZ) >= iTrack){
            fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxZ, seedZ);
          }
        }
        if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kVtxT)){
          if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kVtxT) >= iTrack){
            fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kVtxT, seedT);
          }
        }
        if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kDirTh)){
          if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kDirTh) >= iTrack){
            fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kDirTh, seedTheta);
          }
        }
        if(!fFitterTrackParMap.GetIsFixed(iTrack,FitterParameterType::kDirPhi)){
          if( fitterConfig->GetTrackIsJoinedWith(iTrack, FitterParameterType::kDirPhi) >= iTrack){
            fFitterTrackParMap.SetCurrentValue(iTrack, FitterParameterType::kDirPhi, seedPhi);
          }
        }
      }
    }
  } // End else statement

  // Need to delete the elements of slicedEvents as they are not destroyed by WCSimRecoSlicer
  for(unsigned int v = 0; v < slicedEvents.size(); ++v){
    delete slicedEvents[v];
  }

  delete myReco;
}

void WCSimCosmicFitter::Fit(const char * minAlgorithm)
{
  double likelihood = FitAndGetLikelihood(minAlgorithm);
  std::cout << "Best-fit -2Ln(likelihood) was " << likelihood << std::endl;
  return;
}

double WCSimCosmicFitter::FitAndGetLikelihood(const char * minAlgorithm)
{
  // Set up the minimizer
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", minAlgorithm);

  // Alternatively: use a different algorithm to check the minimizer works
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "GSLSimAn");


  min->SetMaxFunctionCalls(500);
  min->SetMaxIterations(1);
  min->SetPrintLevel(3);
  //min->SetTolerance(0.);
  min->SetErrorDef(1.0);
  min->SetStrategy(2);
  std::cout << " Tolerance = " << min->Tolerance() << std::endl;

  // Convert nTracks into number of parameters
  const unsigned int nPars = this->GetNPars();

  // Tell the minimizer the function to minimize
  // We have to wrap it in this functor because it's a member function of this class
  ROOT::Math::Functor func(this,&WCSimCosmicFitter::WrapFunc, nPars);


  // Tell the minimizer the functor and variables to consider
  min->SetFunction(func);

  // Set the fixed and free variables
  Bool_t isAnythingFree = false;

  std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
  std::vector<double> minVals = fFitterTrackParMap.GetMinValues();
  std::vector<double> maxVals = fFitterTrackParMap.GetMaxValues();
  std::vector<double> steps = fFitterTrackParMap.GetSteps();
  std::vector<bool> currentlyFixed = fFitterTrackParMap.GetCurrentlyFixed();
  std::vector<std::string> names = fFitterTrackParMap.GetNames();
  for(UInt_t i = 0; i < nPars; ++i)
  {
    // Sometimes the fast seed goes outside the allowed range:
    if( currentlyFixed[i] ){
      std::cout << "Fixed: " << names[i] << "  Start: " << startVals[i] << std::endl;
      min->SetFixedVariable(i, names[i], startVals[i]);
    }
    else{
      isAnythingFree = true;
      std::cout << "Free: " << names[i] << "  Start: " << startVals[i] << "   Min: " << minVals[i] << "   Max: " << maxVals[i] << std::endl;
      min->SetLimitedVariable(i, names[i], startVals[i], steps[i], minVals[i], maxVals[i]);
    }

    if(GetUsePiZeroMassConstraint())
    {
      if(i == fFitterTrackParMap.GetIndex(1, FitterParameterType::kEnergy))
      {
        min->SetFixedVariable(i, names[i], startVals[i]);
      }
    }
  }

  // Print the parameters we're using
  for(UInt_t j = 0; j < nPars; ++j)
  {
    std::cout << j << "   " << names[j] << "   " << startVals[j] << "   " << minVals[j] << "   " << maxVals[j] << std::endl;
  }




  // Perform the minimization
  if( isAnythingFree )
  {
    try
    {
      min->Minimize();
      fMinimum = min->MinValue();
      fStatus = min->Status();
    }
    catch(FitterArgIsNaN &e)
    {
      std::cerr << "Error: Fitter encountered an exception when minimising due to NaN argument" << std::endl;
      std::cerr << "	   Will flag this event and skip this minimiser stage" << std::endl;
      fFailed = true;
    }
    catch(...){
      std::cerr << "Caught some unknown?! error..." << std::endl;
      fFailed = true;
    }
  }
  else{
    std::cout << "Nothing in this world is free" << std::endl;
  }
  // This is what we ought to do:
  const Double_t * outPar = min->X();
  // Future tracks need to start from this best fit
  for( unsigned int i = 0; i < nPars; ++i)
  {
    fFitterTrackParMap.SetCurrentValue(i, outPar[i]);

    if(GetUsePiZeroMassConstraint())
    {
      if(i == fFitterTrackParMap.GetIndex(1, FitterParameterType::kEnergy))
      {
        fFitterTrackParMap.SetCurrentValue(i, GetPiZeroSecondTrackEnergy(outPar));
      }
    }
  }

  // Minuit gives us a const double array, but (for now) we want to do the grid search and then modify it
  // std::cout << "Here's outPar..." << std::endl;
  // std::cout << outPar[0] << "  " << outPar[1] << std::endl;

  // Get and print the fit results
  // std::cout << "Best fit track: " << std::endl;
  // RescaleParams(outPar2[0],  outPar2[1],  outPar2[2],  outPar2[3],  outPar2[4],  outPar2[5],  outPar2[6], fType).Print();
  UpdateBestFits();
  return min->MinValue();
}


void WCSimCosmicFitter::FitAlongTrack()
{
  std::cout << "Fitting along track" << std::endl;
  const unsigned int nTracks = fFitterConfig->GetNumTracks();

  // Make use of the WCSimGeometry class
  WCSimGeometry* wcGeom = WCSimGeometry::Instance();

  // We are only looking for two tracks.
  if(nTracks != 2){
    return;
  }

  for(unsigned int i = 0; i < nTracks; ++i){
    bool canMove = (   !fFitterTrackParMap.GetIsFixed(i, FitterParameterType::kVtxX)
        && !fFitterTrackParMap.GetIsFixed(i, FitterParameterType::kVtxY)
        && !fFitterTrackParMap.GetIsFixed(i, FitterParameterType::kVtxZ)
        );

    TVector3 direction;
    direction.SetMagThetaPhi(1.0, 
        fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirTh), 
        fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirPhi)
        );
    double minE = fFitterTrackParMap.GetMinValue(i, FitterParameterType::kEnergy);
    double maxE = fFitterTrackParMap.GetMaxValue(i, FitterParameterType::kEnergy);
    double currE = fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kEnergy);
    double stepE = fFitterTrackParMap.GetStep(i, FitterParameterType::kEnergy);
    bool fixE = fFitterTrackParMap.GetIsFixed(i, FitterParameterType::kEnergy);

    //double minT = fFitterTrackParMap.GetMinValue(i, FitterParameterType::kVtxT);
    //double maxT = fFitterTrackParMap.GetMaxValue(i, FitterParameterType::kVtxT);
    //double currT = fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxT);
    //double stepT = fFitterTrackParMap.GetStep(i, FitterParameterType::kVtxT);
    bool fixT = (fFitterTrackParMap.GetIsFixed(i, FitterParameterType::kVtxT) || !WCSimAnalysisConfig::Instance()->GetUseTime());

    TVector3 startVertex(fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxX),
        fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxY),
        fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxZ)
        );

    TVector3 stepSizes(fFitterTrackParMap.GetStep(i, FitterParameterType::kVtxX),
        fFitterTrackParMap.GetStep(i, FitterParameterType::kVtxY),
        fFitterTrackParMap.GetStep(i, FitterParameterType::kVtxZ)
        );

    TVector3 mins(fFitterTrackParMap.GetMinValue(i, FitterParameterType::kVtxX),
        fFitterTrackParMap.GetMinValue(i, FitterParameterType::kVtxY),
        fFitterTrackParMap.GetMinValue(i, FitterParameterType::kVtxZ));

    TVector3 maxes(fFitterTrackParMap.GetMaxValue(i, FitterParameterType::kVtxX),
        fFitterTrackParMap.GetMaxValue(i, FitterParameterType::kVtxY),
        fFitterTrackParMap.GetMaxValue(i, FitterParameterType::kVtxZ));  

    std::vector<double> escapeDists;
    for(int e = 0; e < 3; ++e)
    {
      if(fabs(direction[e]) > 1e-6)
      {
        escapeDists.push_back( (maxes[e] - startVertex[e])/direction[e] );
        escapeDists.push_back( (mins[e] - startVertex[e])/direction[e] );
      }
    }
    assert(escapeDists.size() > 1);
    std::sort(escapeDists.begin(), escapeDists.end());
    double maxVal = static_cast<unsigned int>(-1);
    double minVal = -1.0 * maxVal;
    std::cout << "minVal = " << minVal << " and start maxVal = " << maxVal << std::endl;
    /*     for(unsigned int e = 0; e < escapeDists.size(); ++e)
           {
           if(escapeDists[e] < 0 && escapeDists[e] > minVal){ minVal = escapeDists[e]; }
           if(escapeDists[e] > 0 && escapeDists[e] < maxVal){ maxVal = escapeDists[e]; }
           }
           */     
    std::cout << startVertex.X() << "," << startVertex.Y() << "," << startVertex.Z() << " :: "
      << direction.X() << "," << direction.Y() << "," << direction.Z() << std::endl;
    wcGeom->ProjectMinMaxPosition(startVertex, direction, minVal, maxVal);
    std::cout << "MinVal = " << minVal << " and maxVal = " << maxVal << std::endl;

    double stepSize = sqrt(pow((stepSizes.X() * direction.X()), 2) + pow((stepSizes.Y() * direction.Y()), 2) + pow((stepSizes.Z() * direction.Z()), 2));
    stepSize = 250;
    std::cout << "Step size = " << stepSize << std::endl;

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

    min->SetMaxFunctionCalls(2500);
    min->SetMaxIterations(10);
    min->SetPrintLevel(3);
    //min->SetTolerance(0.);
    min->SetErrorDef(1.0);
    min->SetStrategy(2);
    std::cout << " Tolerance = " << min->Tolerance() << std::endl;

    const unsigned int nPars = 3; // Distance along track, energy and time.
    fSingleTrackToFit = i; // Make sure the wrap function knows which track to vary.

    // Tell the minimizer the function to minimize
    // We have to wrap it in this functor because it's a member function of this class
    ROOT::Math::Functor func(this,&WCSimCosmicFitter::WrapFuncAlongSingleTrack, nPars);

    // Tell the minimizer the functor and variables to consider
    min->SetFunction(func);
    if( canMove )
    {
      min->SetLimitedVariable(0, "Distance from seed", 0., stepSize, minVal, maxVal);
    }
    else
    {
      min->SetFixedVariable(0, "Distance from seed", 0);
    }
    // Set up the energy variable
    if(fixE == false)
    {
      min->SetLimitedVariable(1, "Energy", currE, stepE, minE, maxE);
    }
    else
    {
      min->SetFixedVariable(1, "Energy", currE);
    }
    // Set up the time variable    
    if(fixT == false && WCSimAnalysisConfig::Instance()->GetUseTime())
    {
      //        min->SetLimitedVariable(2, "Time shift", 0, 1., -5, 5);
      min->SetLimitedVariable(2, "Time shift", 0., 1., minVal/30.0, maxVal/30.);
    }
    else
    {
      min->SetFixedVariable(2, "Time shift",0);
    }
//    min->Minimize();

    double finalDistance = 0;
    double finalDeltaE = 0;
    try
    {
      min->Minimize();
      finalDistance = min->X()[0];
      finalDeltaE = min->X()[1];
    }
    catch(FitterArgIsNaN &e)
    {
      std::cerr << "Error: Fitter encountered an exception when minimising along track, due to NaN argument" << std::endl;
      std::cerr << "	   Will flag this event and skip this minimiser stage" << std::endl;
      fFailed = true;
    }
    catch(...){
      std::cerr << "Caught some unknown?! error..." << std::endl;
      fFailed = true;
    }

    // Get the final fitted parameters from the minimiser
    TVector3 newVertex = startVertex + finalDistance * direction;
    std::cout << "Started at " << startVertex.X() << ", " << startVertex.Y() << ", " << startVertex.Z() << std::endl;
    std::cout << "Shift track vertex by " << finalDistance << " to" << std::endl;
    const double speed = WCSimLikelihoodTrackBase::GetPropagationSpeedFrac(fFitterTrackParMap.GetTrackType(i)) * TMath::C() / 1e7;
    newVertex.Print();
    WCSimFitterTrackParMap& ftpm = fFitterTrackParMap; // Reference this so our lines aren't giant unreadable garbage
    ftpm.SetCurrentValue(ftpm.GetIndex(i, FitterParameterType::kVtxX), newVertex.X());
    ftpm.SetCurrentValue(ftpm.GetIndex(i, FitterParameterType::kVtxY), newVertex.Y());
    ftpm.SetCurrentValue(ftpm.GetIndex(i, FitterParameterType::kVtxZ), newVertex.Z());
    ftpm.SetCurrentValue(ftpm.GetIndex(i, FitterParameterType::kVtxT), 
        ftpm.GetCurrentValue(i, FitterParameterType::kVtxT) 
        + finalDistance/speed
        + min->X()[2]);
    ftpm.SetCurrentValue(ftpm.GetIndex(i, FitterParameterType::kEnergy), min->X()[1]);
  }
  UpdateBestFits();

  return;
}

double WCSimCosmicFitter::WrapFuncAlongSingleTrack(const Double_t * x)
{
  UInt_t nTracks = fFitterConfig->GetNumTracks();

  // Make sure nothing has gone wrong.
  CheckTrackParametersForNaN(x,3); // Distance along the track, the energy and the time.

  // If we've fixed some track parameters together, then our array of fit parameters won't just be of size
  // n tracks * m parameters per track, and we need to work out which entry corresponds to which parameter
  // of which track(s) - this is probably expensive, so we'll do it once and look it up in a map thereafter

  std::vector<WCSimLikelihoodTrackBase*> tracksToFit;
  TrackType::Type trackType = fFitterTrackParMap.GetTrackType(fSingleTrackToFit);
  std::vector<FitterParameterType::Type> allTypes = FitterParameterType::GetAllAllowedTypes(trackType);
  TrackAndType trackParX = std::make_pair(fSingleTrackToFit, FitterParameterType::kVtxX);
  TrackAndType trackParY = std::make_pair(fSingleTrackToFit, FitterParameterType::kVtxY);
  TrackAndType trackParZ = std::make_pair(fSingleTrackToFit, FitterParameterType::kVtxZ);
  TrackAndType trackParT = std::make_pair(fSingleTrackToFit, FitterParameterType::kVtxT);
  TrackAndType trackParTh = std::make_pair(fSingleTrackToFit, FitterParameterType::kDirTh);
  TrackAndType trackParPhi = std::make_pair(fSingleTrackToFit, FitterParameterType::kDirPhi);
  TrackAndType trackParE = std::make_pair(fSingleTrackToFit, FitterParameterType::kEnergy);

  TVector3 vertex(fFitterTrackParMap.GetCurrentValue(trackParX),
      fFitterTrackParMap.GetCurrentValue(trackParY),
      fFitterTrackParMap.GetCurrentValue(trackParZ));

  TVector3 direction(0,0,0);
  direction.SetMagThetaPhi(1.0,fFitterTrackParMap.GetCurrentValue(trackParTh),fFitterTrackParMap.GetCurrentValue(trackParPhi));

  double time = fFitterTrackParMap.GetCurrentValue(trackParT) + x[0]/30.0 + x[2];

  TVector3 newVertex = vertex + direction * x[0];
  std::cout << "x = " << x[0] << " So vertex goes from (" << vertex.X() << ", " << vertex.Y() << ", " << vertex.Z() << ") to ("
    << newVertex.X() << ", " << newVertex.Y() << ", " << newVertex.Z() << ")" << std::endl;
  std::cout << "Time is now " << time << " and energy is " << x[1] << std::endl;

  for(unsigned int i  = 0 ; i < nTracks; ++i)
  {
    WCSimLikelihoodTrackBase * track;
    if(i == fSingleTrackToFit){
      track = WCSimLikelihoodTrackFactory::MakeTrack(
          fFitterTrackParMap.GetTrackType(i),
          newVertex.X(),
          newVertex.Y(),
          newVertex.Z(),
          time,
          fFitterTrackParMap.GetCurrentValue(trackParTh),
          fFitterTrackParMap.GetCurrentValue(trackParPhi),
          x[1]
          );
    }
    else{
      // Keep the other tracks as they were
      track = WCSimLikelihoodTrackFactory::MakeTrack(
          fFitterTrackParMap.GetTrackType(i),
          fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxX),
          fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxY),
          fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxZ),
          fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kVtxT),
          fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirTh),
          fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kDirPhi),
          fFitterTrackParMap.GetCurrentValue(i, FitterParameterType::kEnergy)
          );
    }
    tracksToFit.push_back(track);
  }

  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(tracksToFit);
  minus2LnL = fTotalLikelihood->Calc2LnL();

  // std::cout << "deleting tracks" << std::endl;
  for(size_t iTrack = 0; iTrack < tracksToFit.size(); ++iTrack)
  {
    // std::cout << iTrack << "/" << tracksToFit.size() << std::endl;
    delete (tracksToFit.at(iTrack));
  }
  // std::cout << "Clearing" << std::endl;
  tracksToFit.clear();
  // std::cout << "Done" << std::endl;
  std::cout << "Function evaluated " << ++fCalls << " times for event " << fEvent << " ... " << " -2Ln(L) = " << minus2LnL << std::endl << std::endl;
  return minus2LnL;
}


