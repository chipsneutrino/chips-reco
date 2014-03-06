<<<<<<< HEAD
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "TArrayD.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TArrayD.h"
#include "TH2D.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"

#include "WCSimChargeLikelihood.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimRootGeom.hh"
#include "WCSimGeometry.hh"

/////////////////////////////////////////////////////
// Constructor - if the Tuner is created without
// specifying the size of the detector then it 
// assumes the PMTs are floating in space and will 
// only prevent a PMT from seeing light if the track
// is behind it, not if it's outside the tank on the 
// opposite side.               
////////////////////////////////////////////////////
WCSimLikelihoodTuner::WCSimLikelihoodTuner( Bool_t calculateIntegrals = kFALSE)
{
  	fConstrainExtent = kFALSE;
  	fExtent[0] = -1.0;
  	fExtent[1] = -1.0;
  	fExtent[2] = -1.0;
  
  
  	fCalculateIntegrals = calculateIntegrals;
  	this->Initialize();

}

/////////////////////////////////////////////////////
// Otherwise if the extent of the detector in x, y 
// and z is given, the tuner will kill tracks as soon
// as they go outside of this region
/////////////////////////////////////////////////////
WCSimLikelihoodTuner::WCSimLikelihoodTuner(Double_t xMax, Double_t yMax, Double_t zMax, Bool_t calculateIntegrals = kFALSE)
{
  	fConstrainExtent = kTRUE;
  	fExtent[0] = xMax;
  	fExtent[1] = yMax;
  	fExtent[2] = zMax;
    std::cout << fExtent[0] << " " << fExtent[1] << " " << fExtent[2] << std::endl;
  
  	fCalculateIntegrals = true;
  	this->Initialize();
//    std::cout << "Done!!" << std::endl;
}

///////////////////////////////////////////////////////////
// These values should all be the same regardless of
// the constructor called.  Later we should teach it to
// distinguish between track types and open the appropriate
// emission profile file
///////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::Initialize()
{
  fIsOpen = WCSimLikelihoodTrack::Unknown;
  fProfileLocation = new TString("./config/emissionProfilesElectron.root");
  fProfiles = new TFile(fProfileLocation->Data());
  fHistArray = 0;
  fAngHistArray= 0;
  fFluxArray= 0;
  fIsOpen = WCSimLikelihoodTrack::ElectronLike;
  fHistArray = (TObjArray *) fProfiles->Get("histArray");
  fHistArray->SetOwner(kTRUE);
  fAngHistArray = (TObjArray *) fProfiles->Get("angHistArray");
  fAngHistArray->SetOwner(kTRUE);
  fFluxArray = (TObjArray *) fProfiles->Get("fluxArray");
  fFluxArray->SetOwner(kTRUE);
  fWhichHisto = (TH1D *) fProfiles->Get("hWhichHisto");
  fAverageQE = 1.0;

  // Pointer to the last track for which we calculated the cutoff, to prevent repetition
  fLastCutoff = 0;
  
  //  std::cout << "Setting up files " << std::endl;
  fRhoGIntegralFile = new TFile();
	fRhoIntegralFile = new TFile();
	fRhoGIntegralTree = 0;
	fRhoIntegralTree = 0;
  fRhoGIntegrals = 0;
  fRhoIntegrals = 0;
    
  // The binning scheme for the direct integral tables
	fNR0Bins = 200;
	fNCosTheta0Bins = 1000;
	fNSBins = 500;
	fNEBins = 1;
  fNBinsRhoG = fNR0Bins * fNCosTheta0Bins * fNSBins * fNEBins;
  fR0Min = 0.0;
	fR0Max = 4000.0;
  fCosTheta0Min = 0.0;
  fCosTheta0Max = 1.0;
	fEMin = 1500.0;
	fEMax = 3000.0;
  fSMin = 0.0;
	fSMax = 2500.0; 
  fIntegralEnergyBin = -1;
  fIntegralSMax = -1;
  fIntegralParticleType = WCSimLikelihoodTrack::Unknown;
  
  // The binning scheme for the indirect integral tables
  fNSBinsRho = fNSBins;
  fNEBinsRho = fNEBins;
  fNBinsRho = fNSBinsRho * fNEBinsRho;
  fSMinRho = fSMin;
  fSMaxRho = fSMax;
  fEMinRho = fEMin;
  fEMaxRho = fEMax;
  
//  std::cout << "Done!" << std::endl;
  return;




}
/////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////
WCSimLikelihoodTuner::~WCSimLikelihoodTuner()
{
    //std::cout << " *** WCSimLikelihoodTuner::~WCSimLikelihoodTuner() *** Cleaning up" << std::endl;
    fProfiles->Close();
    fIsOpen = WCSimLikelihoodTrack::Unknown;
    delete fProfileLocation;
    delete fProfiles;
    if(fHistArray) delete fHistArray;
    if(fAngHistArray) delete fAngHistArray;
    if(fFluxArray) delete fFluxArray;
//    if(fWhichHisto) delete fWhichHisto;
    fProfileLocation = 0;
    fProfiles = 0;
    fHistArray = 0;
    fAngHistArray= 0;
    fFluxArray= 0;
    
    if(fRhoGIntegralTree) delete fRhoGIntegralTree;
    if(fRhoIntegralTree) delete fRhoIntegralTree;
    if(fRhoGIntegralFile) delete fRhoGIntegralFile;
	  if(fRhoIntegralFile) delete fRhoIntegralFile;
    if(fRhoGIntegrals) delete fRhoGIntegrals;
    if(fRhoIntegrals) delete fRhoIntegrals;

}

void WCSimLikelihoodTuner::UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray)
{
  fConstrainExtent = true; 
  fExtent[0] = myDigitArray->GetExtent(0);
  fExtent[1] = myDigitArray->GetExtent(1);
  fExtent[2] = myDigitArray->GetExtent(2);
  return;
}


////////////////////////////////////////////////////////////////////////
// Load the appropriate emission profiles for this track's particle type
////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadEmissionProfiles( WCSimLikelihoodTrack * myTrack )
{
  this->LoadEmissionProfiles(myTrack->GetType());
  return;
}

///////////////////////////////////////////////////////////////
// Load the appropriate emission profile for a given track type
///////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadEmissionProfiles( WCSimLikelihoodTrack::TrackType myType )
{
  //std::cout << " *** WCSimLikelihoodTuner::LoadEmissionProfiles - Loading profile" << std::endl;

  if( myType == fIsOpen )
  {
    //std::cout << "Is already open" << std::endl;
    return;
  }

  if( fProfileLocation != NULL)
  {
    delete fProfileLocation;
    fProfileLocation = 0;
  }

  switch( myType )
  {
    case WCSimLikelihoodTrack::ElectronLike:
      fProfileLocation = new TString("./config/emissionProfilesElectron.root");
      std::cerr << "Track type = " << myType << std::endl;
      break;
    case WCSimLikelihoodTrack::MuonLike:
      fProfileLocation = new TString("./config/emissionProfilesMuon.root");
      std::cerr << "Track type = " << myType << std::endl;
      break;
    default:
      std::cerr << "Error: unkown track type in WCSimLikelihoodTuner::LoadEmissionProfiles" << std::endl;
      std::cerr << "Track type = " << myType << std::endl;
      exit(EXIT_FAILURE);
  }
 
  if( fProfiles != NULL )
  {
    //std::cout << fProfiles << std::endl;
    delete fProfiles;
    delete fHistArray;
    delete fAngHistArray;
    delete fFluxArray;
    //std::cout << " Was null" << fProfiles << std::endl;
  }

  fProfiles = new TFile(fProfileLocation->Data(),"READ");
  fIsOpen = myType;
  fHistArray = (TObjArray *) fProfiles->Get("histArray");
  fHistArray->SetOwner(kTRUE);
  fAngHistArray = (TObjArray *) fProfiles->Get("angHistArray");
  fAngHistArray->SetOwner(kTRUE);
  fFluxArray = (TObjArray *) fProfiles->Get("fluxArray");
  fFluxArray->SetOwner(kTRUE);
  fWhichHisto = (TH1D*) fProfiles->Get("hWhichHisto");
  return;
  

}

/////////////////////////////////////////////////////////////////////////////////////////
// Work out the probability of light surviving to the PMT as it travels through the water
/////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::TransmissionFunction(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{

	
	if( s== 0) return 1;
    // First we need the distance from the photon emission to the PMT 
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    Double_t r = (pmtPos - emissionPos).Mag();


    // We'll use a triple exponential to start with
    Double_t nu[3] = {-1.137e-5,-5.212e-4, -4.359e-3}; // nu = 1/Decay length in mm
    Double_t f[3]      = {0.8827, 0.08162, 0.03515};
    Double_t trans=0.0;
    for(int i = 0; i < 3; ++i){ trans+= f[i]*exp(1.0 * 10 * r * nu[i]);}  //Convert to cm -> factor 10
    return trans;
}

///////////////////////////////////////////////////////////////////////////////////////
// The efficiency of the PMT as a function of the incident angle of the photon
// This uses the same formula as WCSim.  Quantum efficiency (QE) is accounted for here.
// Efficiency of digitization is handled elsewhere
///////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::Efficiency(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
	return 1;
    // We need the angle of incidence at the PMT
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    TVector3 pmtToEm = emissionPos - pmtPos;
    TVector3 pmtFace(myDigit->GetFaceX(), myDigit->GetFaceY(), myDigit->GetFaceZ());

    Double_t cosTheta = pmtFace.Dot(pmtToEm) / (pmtFace.Mag() * pmtToEm.Mag()); 

	// The MiniBooNE method:
    Double_t theta = TMath::ACos(cosTheta) * 180. / TMath::Pi();
    if( theta > 90.0 )
    {
      theta = 180.0 - theta;
    }
	
	Double_t efficiency =  (1 + (-1.182e-4) * pow(theta, 2) + 4.959e-9 * pow(theta, 4) - 7.371e-14 * pow(theta, 6));
	return efficiency;

    // Function for the PMT acceptance's dependence on the angle: 
    // WCSim defines arrays of efficiency at 10 degree intervals and linearly interpolates

 /*   if( cosTheta < 0.0 || cosTheta > 1.0 )
    {

//      std::cout << "Behind the PMT, cosTheta = " << cosTheta << std::endl;
//      pmtPos.Print() ;
//      emissionPos.Print();
      return 0.0;
    }

    Double_t theta = TMath::ACos(cosTheta) * 180.0 / TMath::Pi();
    Double_t collection_angle[10]={0.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
    Double_t collection_eff[10]={100.,100.,99.,95.,90.,85.,80.,69.,35.,13.}; 
    Int_t num_elements = sizeof( collection_angle ) / sizeof( collection_angle[0] );

    Double_t efficiency = 0.0;
    for(int iEntry = 0; iEntry < num_elements-1; ++iEntry)
    {
      if( theta >= collection_angle[iEntry] && theta < collection_angle[iEntry+1])
      { 
        efficiency = collection_eff[iEntry] 
                     + (theta - collection_angle[iEntry]) / (collection_angle[iEntry+1] - collection_angle[iEntry])
                     * (collection_eff[iEntry+1] - collection_eff[iEntry]);
      }
    }
      return efficiency/100.;
*/






}

/*
//////////////////////////////////////////////////////////////////////////////////////////
// Scale the expected photons by the quantum efficiency, which for now we treat as a 
// constant by averaging the QE over the wavelength spectrum for all emitted photons
// and scaling the predicted charge by this number
//////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::QuantumEfficiency(WCSimLikelihoodTrack * myTrack)
{
	WCSimLikelihoodTrack::TrackType type = myTrack->GetType();
	if(type == WCSimLikelihoodTrack::MuonLike) return 0.05;
	else if(type == WCSimLikelihoodTrack::ElectronLike) return 1.0;
	else return 1.0;
}
*/

//////////////////////////////////////////////////////////////////////////////////////////
// Work out the effect of solid angle on the probability of photon incidenc; pure geometry
//////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::SolidAngle(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
    // Now some physical parameters of the phototube
	Double_t WCSimPMTRadius = 12.7;
//    Double_t WCSimPMTExposeHeight = WCSimPMTRadius - 1.0;

	
    // We need the distance from emission to the PMT
    // These are from WCSim so they're in cm
    TVector3 pmtPos(myDigit->GetX() + myDigit->GetFaceX() * WCSimPMTRadius,
    				myDigit->GetY() + myDigit->GetFaceY() * WCSimPMTRadius, 
    				myDigit->GetZ() + myDigit->GetFaceZ() * WCSimPMTRadius);
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    Double_t r = (pmtPos - emissionPos).Mag();

    // Purely geometry: we need the solid angle of a cone whose bottom is a circle of the same radius as the circle of PMT poking
    // through the blacksheet.  This is 2pi( 1 - cos(coneAngle) ) where coneAngle can be deduced from trig and the PMT geometry
	// So the fraction of total solid angle is this/4pi = 0.5(1 - cos(coneAngle))
	Double_t solidAngle = 2.0*TMath::Pi()*(1.0 - (r)/sqrt( (r*r + WCSimPMTRadius * WCSimPMTRadius)));
	
	// I multiply by cosTheta to account for the angle from track to PMT diminishing its effective area
	// This isn't strictly correct (solid angle of a sphere cut off by a plane turns out to be quite complicated
	// but it's a decent approximation to a small effect
//	solidAngle *= TMath::Cos((emissionPos - pmtPos).Angle(TVector3(myDigit->GetFaceX(), myDigit->GetFaceY(), myDigit->GetFaceZ())));

	// See whether the PMT covers more than 1 bin in cosTheta - if it does we should have picked multiple theta bins from g(cosTheta)
	// so scale the solid angle accordingly.
//    this->LoadEmissionProfiles(myTrack);
////	std::cout << fWhichHisto->FindBin(myTrack->GetE());
//	TH2D * hG = (TH2D*)fAngHistArray->At(fWhichHisto->FindBin(myTrack->GetE()) - 1);
////	std::cout << fWhichHisto->FindBin(myTrack->GetE() - 1) << "   " << hG << std::endl;
//	Double_t cosThetaBins = hG->GetNbinsX();
//	Double_t cosThetaRange = (hG->GetXaxis())->GetXmax() - (hG->GetXaxis())->GetXmin();
//	Double_t dCosTheta = cosThetaRange / cosThetaBins;
//	
//	Double_t angleToCenter = pmtPos.Angle(emissionPos);
//	Double_t rangeThetaPMT = TMath::ACos(1.0 - 0.5*solidAngle / TMath::Pi());
////	std::cout << "angle to center " << angleToCenter << "     rangeThetaPMT " << rangeThetaPMT << std::endl;
//	Double_t rangeCosThetaPMT = 2 * TMath::Sin(angleToCenter) * TMath::Sin(rangeThetaPMT);
//
//
//	std::cout << "pmt range = " << rangeCosThetaPMT << " dCosTheta = " << dCosTheta << "  range / dCosTheta = " << rangeCosThetaPMT / dCosTheta << std::endl;
//	solidAngle *= rangeCosThetaPMT / dCosTheta;
	/*if(pmtPos.Z() > 1000)
	{
		std::cout << "cosTheta = " << TMath::Cos(TMath::ATan2(WCSimPMTRadius, r)) << std::endl;
    	std::cout << "solid angle = " << solidAngle << std::endl;
    	std::cout << "r = " << r << "    WCSimPMTRadius = " << WCSimPMTRadius << "   solid angle  = " << solidAngle << "    hit Q = " << myDigit->GetQ() << std::endl;;
    	std::cout << "Pmt pos = " ;
    	pmtPos.Print();
    	std::cout << "Emission pos = " ;
    	emissionPos.Print();
    	std::cout << std::endl;
	}*/
    return solidAngle;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////
// The scattering table assumes that scattered light can be modelled by scaling the
// number of photons you would see from an isotropic source of Cherenkov light by
// a factor dependent on geometry, energy etc.
// It's small and computationally intensive to tabulate, so we're taking a constant initially
/////////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::ScatteringTable(Double_t s)
{
    // 20% scattered light for now
    return 0.001;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//  We define J as the product of the transmission, efficiency and solid angle fuctions (and
//  also the scattering table for indirect light).  However we need it to be zero once the
//  particle has exited the detector, so we check for that before calling the other functions
/////////////////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::CalculateJ( Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit ) 
{
    std::vector<Double_t> J;

    // Check make sure the particle is still inside the detector
	if( fConstrainExtent )
	{
	  TVector3 pos = myTrack->GetVtx() + s*myTrack->GetDir();
	  if( fabs(pos.X()) > fExtent[0] || fabs(pos.Y()) > fExtent[1] || fabs(pos.Z()) > fExtent[2])
      {
        J.push_back(0.0);
        J.push_back(0.0);
        return J;
      }
    }


    // Work out the direct and indirect contributions to J
    // J[0] = J_dir, J[1] = J_ind
    if( s == 0.0)
    {
    	std::cout << "Transmission = " << this->TransmissionFunction(s, myTrack, myDigit) << std::endl
    			  << "Efficiency = "   << this->Efficiency(s, myTrack, myDigit) << std::endl
    			  << "SolidAngle = " << this->SolidAngle(s, myTrack, myDigit) << std::endl;
    }
    J.push_back(   this->TransmissionFunction(s, myTrack, myDigit) 
                 * this->Efficiency(s, myTrack, myDigit)
            //     * this->QuantumEfficiency(myTrack)
                 * this->SolidAngle(s, myTrack, myDigit));
    J.push_back(J.at(0) * this->ScatteringTable(s));
    return J;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The direct and indirect expectations come from a large integral over several factors.  We separate out the
// transmission function, efficiency and solid angle and fit these to a quadratic in s (distance along track)
// This just leaves integrals over emission profiles that can be tabulated in advance.  This function gives
// the coefficients of these quadratics
///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void WCSimLikelihoodTuner::CalculateCoefficients(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit )
{
    
    // std::cout << "*** WCSimLikelihoodTuner::CalculateCoefficients() *** Calculating the coefficients for the integrals from tuning info" << std::endl; 
    for(int i = 0; i < 3; ++i)
    { 
        fDirCoeffs[i] = 0.0;    // Coefficients for direct (Cherenkov) light
        fIndCoeffs[i] = 0.0;    // And for indirect light
     }
    
    // Load the file containing the emission profiles for different energies
    // whichHisto should be a TH1I with the same energy binning as we want for the reconstruction
    // There should also be a TObjArray with the emission profiles for the center value of each energy bin
//    std::cout << "Getting the histogram" << std::endl;
    this->LoadEmissionProfiles(myTrack);
 
    Int_t whichBin = fWhichHisto->GetXaxis()->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    // std::cout << "WhichBin = " << whichBin << "   energy = " << myTrack->GetE() << std::endl;
    if(whichBin < 0 || whichBin > fWhichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateCoefficients() looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }

    // std::cout << "We want number " << whichBin << std::endl;
    // std::cout << "fHistArray= " << fHistArray << std::endl;
    TH1D * hProfile = (TH1D*)(fHistArray->At(whichBin));
    // Calculate the 3 s values we care about: s=0, s where integral(rho(s') ds')_0^{s} = 0.75, and double that
    Double_t runningTotal=0.;
    Double_t s[3];
    
    for( int iBin = 1; iBin <= hProfile->GetNbinsX(); ++iBin )
    {
        runningTotal += hProfile->GetBinContent(iBin) * hProfile->GetBinWidth(iBin);
        if(runningTotal >= 0.75)
        {
            //std::cout << "Setting s[i]" << std::endl;
            s[0] = 10.0;
            s[1] = hProfile->GetBinCenter(iBin);
            s[2] = 2 * s[1];
            break;
        }
    }
   
    // Check these s values are inside the detector:
    this->CalculateCutoff(myTrack);
    if( fCutoffIntegral > 0.0 && s[2] > fCutoffIntegral )
    {
      s[2] = 0.8*fCutoffIntegral;
      s[1] = 0.5 * s[2];
    }

    // Evaluate J at each point
    Double_t JDir[3] = {0.,0.,0.};
    Double_t JInd[3] = {0.,0.,0.};
        
        // std::cout << " Second loop: kk = " << kk << "    s = " << s[kk] << "     JDir = " << JDir[kk] << std::endl;
    for(int k = 0; k < 3; ++k)
    {
        std::vector<Double_t> J = this->CalculateJ( s[k], myTrack, myDigit);
        JDir[k] = J.at(0); 
        JInd[k] = J.at(1);
    }
    // Now do a quadratic fit to the Cherenkov and indirect parts, sampling them at these 3 points
    // We can actually work this out analytically:    
    
    fDirCoeffs[0] = JDir[0];
    fDirCoeffs[1] = (4.*JDir[1] - JDir[2] -3.*JDir[0]) / (2*s[1]);
    fDirCoeffs[2] = (JDir[0] + JDir[2] - 2.*JDir[1]) / (2 * s[1] * s[1]);
//    std::cout << "jDir:        " << JDir[0]       << "," << JDir[1]       << "," << JDir[2]       << std::endl;
//    std::cout << "fDirCoeffs:  " << fDirCoeffs[0] << "," << fDirCoeffs[1] << "," << fDirCoeffs[2] << std::endl; 
//    

    fIndCoeffs[0] = JInd[0];
    fIndCoeffs[1] = (4.*JInd[1] - JInd[2] -3.*JInd[0]) / (2*s[1]);
    fIndCoeffs[2] = (JInd[0] + JInd[2] - 2.*JInd[1]) / (2 * s[1] * s[1]);

//    if((myDigit->GetTubeId() % 10000) == 0)
//    {
//      TGraph * gr = new TGraph();
//      for(int i = 0; i < hProfile->GetNbinsX(); ++i)
//      {
//        double tempS = hProfile->GetXaxis()->GetBinCenter(i);
//        std::vector<Double_t> tempJ = this->CalculateJ(tempS, myTrack, myDigit);
//        gr->SetPoint(gr->GetN(), tempS, tempJ.at(1));
//      }
//      TF1 * jInd = new TF1("jInd","[0] + [1] * x + [2] * x * x", 0, hProfile->GetXaxis()->GetXmax());
//      jInd->SetParameters(fIndCoeffs[0], fIndCoeffs[1], fIndCoeffs[2]);
//  
//      TCanvas * c1 = new TCanvas("c1","c1", 800, 600);
//      gr->SetLineColor(kBlue);
//      gr->SetLineWidth(2);
//      gr->Draw("ALP");
//      jInd->SetLineColor(kRed);
//      jInd->SetLineWidth(2);
//      jInd->Draw("same");
//      TString str;
//      str.Form("jInd%d.png",myDigit->GetTubeId());
//      c1->SaveAs(str.Data());
//      delete c1;
//      delete jInd;
//      delete gr;
//    }  
  
    // std::cout << "Done here" << std::endl;  
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Same as the function above, but returns a vector of the coefficients as well as updating the class member variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::CalculateCoefficientsVector(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit )
{
    this->CalculateCoefficients(myTrack, myDigit);
    
    // Return these in a vector
    std::vector<Double_t> coeffs;
    for(int j = 0; j < 3; ++j){ coeffs.push_back(fDirCoeffs[j]); }
    for(int j = 0; j < 3; ++j){ coeffs.push_back(fIndCoeffs[j]); }
    return coeffs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// We shouldn't integrate over any more of the track once the particle has exited the detector.
// Here, we input the track and work out how far it has to travel before it reaches the edge of the
// detector volume, and return that value
// ////////////////////////////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::CalculateCutoff( WCSimLikelihoodTrack * myTrack )
{
//  std::cout << "Calculating the cutoff" << std::endl;
  if(myTrack == fLastCutoff) return;
  
  Double_t cutoff = fSMax;

  if( fConstrainExtent )
  {
    TVector3 vtx = myTrack->GetVtx();
    TVector3 dir = myTrack->GetDir();
    Double_t sMax[3] = {0., 0.,0.};
    for(int i = 0; i < 3; ++i)
    {
      if(dir[i] > 1e-6)
      {
        sMax[i] = -1. * (vtx(i) / dir(i)) + TMath::Sqrt( fExtent[i] * fExtent[i] / (dir(i)*dir(i)));
        std::cout << sMax[i] << std::endl;
      }
      else sMax[i] = cutoff;
    }
  
    // The world's laziest sorting algorithm:
    if( sMax[0] < cutoff) cutoff = sMax[0];
    if( sMax[1] < cutoff ) cutoff = sMax[1];
    if( sMax[2] < cutoff ) cutoff = sMax[2];
  }
  std::cout << "Cutoff = " << cutoff << ".... returning" << std::endl;
  fCutoffIntegral = cutoff;
  fLastCutoff = myTrack;
  return;
}

///////////////////////////////////////////////////////////////////////////
//	Work out what array entry we need to query to look up the integrals
//////////////////////////////////////////////////////////////////////////
Int_t WCSimLikelihoodTuner::GetEBin(Double_t energy)
{
  if(energy == fEMax) return fNEBins-1;
  return static_cast<Int_t>( (energy - fEMin) / ((fEMax - fEMin)/(Double_t)fNEBins) );
}

Int_t WCSimLikelihoodTuner::GetSBin(Double_t sMax)
{  
  if(sMax >= fSMax) return 3*(fNSBins-1);
  return 3*static_cast<Int_t>( (sMax - 0) / (fSMax/(Double_t)fNSBins ));
}


Int_t WCSimLikelihoodTuner::GetESBin(Double_t energy, Double_t sMax)
{  
  Int_t eBin = (energy >= fEMaxRho)? fNEBinsRho-1 : static_cast<Int_t>( (energy - fEMinRho) / ((fEMaxRho-fEMinRho)/(Double_t)fNEBinsRho ));
  Int_t sBin = (sMax >= fSMaxRho)? fNSBins-1 : static_cast<Int_t>( (sMax - fSMinRho) / ((fSMaxRho-fSMinRho)/(Double_t)fNSBinsRho ));
//  std::cout << "energy = " << energy << "   so eBin = " << eBin << std::endl
//            << "sMax   = " << fSMax  << "   so sBin = " << sBin << std::endl
//            << "ESBin  = " << sBin + fNSBins * eBin << std::endl;
  return 3*(sBin + fNSBinsRho * eBin);
}

Int_t WCSimLikelihoodTuner::GetIntegralBin(Double_t E, Double_t sMax, Double_t R0, Double_t cosTheta0)
{
	  Int_t eBin = (E >= fEMax) ? fNEBins-1 : static_cast<Int_t>( (E - fEMin) / ((fEMax - fEMin)/(Double_t)fNEBins) );
    Int_t sBin = (sMax >= fSMax) ? fNSBins-1 : static_cast<Int_t>( (sMax - fSMin) / ((fSMax-fSMin)/(Double_t)fNSBins ));
    Int_t R0Bin = (R0 >= fR0Max) ? fNR0Bins-1 : static_cast<Int_t>( (R0 - fR0Min) / ((fR0Max-fR0Min)/(Double_t)fNR0Bins ));
    Int_t cosTheta0Bin = ( cosTheta0 >= fCosTheta0Max)? fNCosTheta0Bins-1 : static_cast<Int_t>( (cosTheta0 - fCosTheta0Min) / ((fCosTheta0Max - fCosTheta0Min)/(Double_t)fNCosTheta0Bins ));
//    std::cout << "Values are " << "Energy = " << E << "   sMax = " << sMax << "   R0 = " << R0 << "   cosTheta0" << cosTheta0 << std::endl;
//    std::cout << "Bins are " << eBin << "/" << fNEBins << "    " << sBin << "/" << fNSBins << "    " << R0Bin << "/" << fNR0Bins << "     " << cosTheta0Bin << "/" << fNCosTheta0Bins <<  std::endl;
    
   	return  3*(sBin + fNSBins*(cosTheta0Bin + fNCosTheta0Bins*(R0Bin + fNR0Bins *( eBin ))));
}

///////////////////////////////////////////////////////////////////////////
// We decouple the direct part into a source factor (flux x profile) and 
// an acceptance factor (solid angle x transmittance x PMT acceptance)
// We then integrate over the source factor, but it's much quicker to 
// tabulate this in advance.  This function opens the file containing those
// tables
///////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadTabulatedIntegrals( WCSimLikelihoodTrack * myTrack )
{
  // Check we have the same track type, energy and s cutoff as before
	WCSimLikelihoodTrack::TrackType myType = myTrack->GetType();
  this->CalculateCutoff(myTrack);

  // In which case we don't need to do anything!
	if( myType == fIntegralParticleType )
	{
//			std::cout << "Returning!" << std::endl;
		return;
	}

		
  std::cout << "Looking up integrals" << std::endl;
  // Otherwise open up the right table file
  std::ifstream integralFileRhoG;
  std::ifstream integralFileRho;
	switch(myType)
	{
		case WCSimLikelihoodTrack::MuonLike:
			std::cout << "It's muon-like" << std::endl;
      integralFileRhoG.open("config/integralsMuonRhoG.dat",std::ios::in|std::ios::binary);
      integralFileRho.open("config/integralsMuonRho.dat",std::ios::in|std::ios::binary);
      break;
	  	case WCSimLikelihoodTrack::ElectronLike:
	  		std::cout << "Is it tripping the electronlike flag too?" << WCSimLikelihoodTrack::ElectronLike << std::endl;
      integralFileRhoG.open("config/integralsElectronRhoG.dat",std::ios::in|std::ios::binary);
      integralFileRho.open("config/integralsElectronRho.dat",std::ios::in|std::ios::binary);
	  	break;
		case WCSimLikelihoodTrack::Unknown:
  		std::cout << "Or is it?" << "   " << myType << std::endl;
	  		std::cerr << "Error: could not identify particle type, exiting" << std::endl;	 			
        exit(EXIT_FAILURE);	
			break;
	}

  if(integralFileRhoG.is_open() == false)
  {
    std::cerr << "Could not open tabulated direct integral file" << std::endl;
    exit(EXIT_FAILURE);
  }
  else
  {
     // The first 12 entries tell us how many bins we have in each variable, and what the minimum and maximum is
    
    integralFileRhoG.read((char *) &fNEBins, sizeof(fNEBins));
    integralFileRhoG.read((char *) &fNR0Bins, sizeof(fNR0Bins));
    integralFileRhoG.read((char *) &fNCosTheta0Bins, sizeof(fNCosTheta0Bins));
    integralFileRhoG.read((char *) &fNSBins, sizeof(fNSBins));
    integralFileRhoG.read((char *) &fEMin, sizeof(fEMin));
    integralFileRhoG.read((char *) &fEMax, sizeof(fEMax));
    integralFileRhoG.read((char *) &fR0Min, sizeof(fEMin));
    integralFileRhoG.read((char *) &fR0Max, sizeof(fEMax));
    integralFileRhoG.read((char *) &fCosTheta0Min, sizeof(fCosTheta0Min));
    integralFileRhoG.read((char *) &fCosTheta0Max, sizeof(fCosTheta0Max));
    integralFileRhoG.read((char *) &fSMin, sizeof(fSMin));
    integralFileRhoG.read((char *) &fSMax, sizeof(fSMax));
       std::cout << fNCosTheta0Bins << "  " << fSMax << std::endl;

    fNBinsRhoG = 3 * fNEBins * fNR0Bins * fNCosTheta0Bins * fNSBins;
    std::cout << fNBinsRhoG;
    if(fRhoGIntegrals) delete fRhoGIntegrals;
    fRhoGIntegrals = new Double_t[fNBinsRhoG];
    integralFileRhoG.read((char *) &fRhoGIntegrals[0], fNBinsRhoG * sizeof(fRhoGIntegrals[0]));
//    ofstream outFile("whatWeReadRhoG.txt");
//    for(int i = 0; i < fNBinsRhoG; i+=3)
//    {
//      outFile << fRhoGIntegrals[i] << " ";
//      outFile << fRhoGIntegrals[i+1] << " ";
//      outFile << fRhoGIntegrals[i+2] << " ";
//      if((i % 60) == 0) outFile << std::endl;
//    }
//    outFile.close();

  } 
    
  if(integralFileRhoG.is_open() == false)
  {
    std::cerr << "Could not open tabulated indirect integral file" << std::endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    integralFileRho.read((char *) &fNEBinsRho, sizeof(fNEBinsRho));
    integralFileRho.read((char *) &fNSBinsRho, sizeof(fNSBinsRho));
    integralFileRho.read((char *) &fEMinRho, sizeof(fEMinRho));
    integralFileRho.read((char *) &fEMaxRho, sizeof(fEMaxRho));
    integralFileRho.read((char *) &fSMinRho, sizeof(fSMinRho));
    integralFileRho.read((char *) &fSMaxRho, sizeof(fSMaxRho));
   
    std::cout << fNEBinsRho << "  " << fNSBinsRho << "  " << fEMinRho << "  " << fEMaxRho <<"   " << fSMinRho <<"   " << fSMaxRho << std::endl;
   
    if(fRhoIntegrals) delete fRhoIntegrals;
    fNBinsRho = 3 * fNEBinsRho * fNSBinsRho;
    fRhoIntegrals = new Double_t[fNBinsRho];
    integralFileRho.read((char *) &fRhoIntegrals[0], fNBinsRho * sizeof(fRhoIntegrals[0]));
//    ofstream outFile;
//    outFile.open("whatWeReadRho.txt");
//    for(int i = 0; i < fNBinsRho; i+=3)
//    { 
//      outFile << fRhoIntegrals[i] << " ";
//      outFile << fRhoIntegrals[i+1] << " ";
//      outFile << fRhoIntegrals[i+2] << " ";
//      if((i%60)==0 && i != 0) outFile << std::endl;
//    }
//    outFile.close();
  }

  fIntegralParticleType = myType;

}

///////////////////////////////////////////////////////////////////////////
//	Look up the tabulated integrals - just return 1 power of s
///////////////////////////////////////////////////////////////////////////
double WCSimLikelihoodTuner::LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, Int_t sPower)
{
    if( !(sPower < 3 && sPower >= 0) )
    {
        std::cerr << "WCSimLikelihoodTuner::LookupChIntegrals() ERROR: power of s must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }	
    
    std::vector<Double_t> integralsVec = this->LookupChIntegrals(myTrack, myDigit);
    
  if( (UInt_t)sPower > integralsVec.size() ) std::cerr << "There's a problem with integralsVec!" << std::endl;  
	return integralsVec.at(sPower);
}

///////////////////////////////////////////////////////////////////////////
// Look up the tabulated integrals: return a vector of all 3 powers of s
// (this is the fastest way to work out the likelihood)
///////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
//  	std::cout << "*** WCSimLikelihoodTuner::LookupChIntegrals() *** Looking up the tabulated integrals for direct Cherenkov light" << std::endl; 
	this->LoadTabulatedIntegrals( myTrack );

 	
 	//std::cout << "Doing the vector stuff ..." ;
 	TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
 	Double_t R0 = (pmtPos - (myTrack->GetVtx())).Mag();
 	Double_t cosTheta0 = myTrack->GetDir().Dot(pmtPos - myTrack->GetVtx()) / (R0);
 	//std::cout << "Done! ... ";
 	
 	//std::cout << "Checking which bin...";
  Int_t whichBin = this->GetIntegralBin(myTrack->GetE(), fCutoffIntegral, R0, cosTheta0);
  //std::cout << "  E      = " << myTrack->GetE() << "  sMax      = " << fCutoffIntegral << " R0     = " << R0 << "    cosTheta0 = " << cosTheta0 << std::endl;
//  std::cout << "Int bin = " << this->GetIntegralBin(myTrack->GetE(), fCutoffIntegral, R0, cosTheta0) << std::endl;
//  std::cout << "Integrals are      i[0]    = " << fRhoGIntegrals[whichBin] << "         i[1] = " << fRhoGIntegrals[1+whichBin] << "       i[2] = " << fRhoGIntegrals[whichBin+2] << std::endl;

// 	std::cout << "Bin = " << whichBin << "   Done! ... ";
 	
	//std::cout << "Making the vector...";
	std::vector<Double_t> integralsVec;
  if( whichBin+2 > fNBinsRhoG ) std::cerr << "There's a problem with fRhoGIntegrals!" << std::endl
                                                      << fNBinsRhoG << "   " << whichBin << std::endl;
	integralsVec.push_back(fRhoGIntegrals[whichBin]);
	integralsVec.push_back(fRhoGIntegrals[whichBin+1]);
	integralsVec.push_back(fRhoGIntegrals[whichBin+2]);
	
 	//std::cout << "Done! ... ";
 	//std::cout << std::endl;
  
//  std::cout << "R0 = " << R0 << "    cosTheta0 = " << cosTheta0 << "     E = " << myTrack->GetE() << "    sMax = " << fCutoffIntegral << "     s term = " << integralsVec.at(1) << "    integral bin = " << whichBin << std::endl;
  
 	return integralsVec;


}

///////////////////////////////////////////////////////////////////////////
//	Look up the tabulated integrals
///////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::LookupIndIntegrals(WCSimLikelihoodTrack * myTrack, Int_t sPower)
{
	if( !(sPower < 3 && sPower >= 0) )
    {
        std::cerr << "WCSimLikelihoodTuner::LookupChIntegrals() ERROR: power of s must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::vector<Double_t> integralsVec = this->LookupIndIntegrals(myTrack);
    if( sPower > integralsVec.size() ) std::cerr << "There's a problem with integralsVec!" << std::endl;
	return integralsVec.at(sPower);
}

std::vector<Double_t> WCSimLikelihoodTuner::LookupIndIntegrals(WCSimLikelihoodTrack * myTrack)
{
//  	std::cout << "*** WCSimLikelihoodTuner::LookupIndIntegrals() *** Looking up the tabulated integrals for indirect light" << std::endl; 
 	  this->LoadTabulatedIntegrals( myTrack );
    
    Int_t ESBin = this->GetESBin(myTrack->GetE(), fCutoffIntegral); 
    std::vector<Double_t> integralsVec;
    if( ESBin > fNBinsRho ) std::cerr << "There's a problem with fRhoIntegrals!" << std::endl;
    integralsVec.push_back(fRhoIntegrals[ESBin]);
    integralsVec.push_back(fRhoIntegrals[ESBin+1]);
    integralsVec.push_back(fRhoIntegrals[ESBin+2]);
    return integralsVec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Work out the integrals over the emission profiles numerically.  This is slow and only good for testing - 
// eventually they'll be tabulated and looked-up
///////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, int i)
{
	std::vector<Double_t> integrals = this->CalculateChIntegrals(myTrack, myDigit);
	if( i < 0 || i > 2) 
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateChIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    return integrals.at(i);
}

std::vector<Double_t> WCSimLikelihoodTuner::CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{


//    std::cout << "*** WCSimLikelihoodTuner::CalculateChIntegrals() ***" << std::endl;
    this->LoadEmissionProfiles(myTrack);
    Int_t whichBin = fWhichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin > fWhichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateChIntegrals() is looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }


    //std::cout << "We want number " << whichBin << std::endl;
    TH1D * hProfile = (TH1D*)fHistArray->At(whichBin);
    TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(whichBin);
 
    // Check to see if we need to curtail the integrals:
    std::vector<Double_t> integrals;
    for(int i = 0; i < 3; ++i)
    {
    	integrals.push_back(0.0);
    }
    
    CalculateCutoff( myTrack );
    int cutoffBin = hProfile->GetNbinsX();
    if( fCutoffIntegral > 0.0 && fCutoffIntegral < hProfile->GetXaxis()->GetXmax() && fCutoffIntegral > hProfile->GetXaxis()->GetXmin() )
    {
      cutoffBin = hProfile->FindBin(fCutoffIntegral);
    }

    // Now do the integrals. 
    // Need these to calculate theta; can declare them outside of the loop
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 vtxPos(myTrack->GetVtx());
    TVector3 vtxDir(myTrack->GetDir());
    for(int iBin = 1; iBin <= cutoffBin; ++iBin)
    {
        Double_t rho = 0.;
        Double_t g = 0.; 
        Double_t s = hProfile->GetBinCenter(iBin);
//        Double_t sLow = hProfile->GetBinLowEdge(iBin);
//        Double_t sHigh = sLow + hProfile->GetBinWidth(iBin);
        if( s > fCutoffIntegral && fCutoffIntegral >= 0. ) 
        {
          break;
        }

        Double_t cosTheta = TMath::Cos(vtxDir.Angle( pmtPos - (vtxPos + s * vtxDir)));
//        Double_t cosThetaLow = TMath::Cos(vtxDir.Angle( pmtPos - (vtxPos + sLow * vtxDir) ));
//        Double_t cosThetaHigh = TMath::Cos(vtxDir.Angle( pmtPos - (vtxPos + sHigh * vtxDir) ));

                // Make sure the histograms in s and s, cosTheta(s) have the same binning on the s axis
        rho = hProfile->GetBinContent(iBin);
        g = hAngularProfile->GetBinContent( hAngularProfile->GetXaxis()->FindBin(cosTheta), iBin );

//        Int_t binCosThetaLow = hAngularProfile->GetXaxis()->FindBin(cosThetaLow);
//        Double_t fractionLow = (cosThetaLow - hAngularProfile->GetXaxis()->GetBinLowEdge(binCosThetaLow))/
//                               (hAngularProfile->GetXaxis()->GetBinWidth(binCosThetaLow));
//        Int_t binCosThetaHigh = hAngularProfile->GetXaxis()->FindBin(cosThetaHigh);
//        Double_t fractionHigh = (cosThetaHigh - hAngularProfile->GetXaxis()->GetBinLowEdge(binCosThetaHigh))/
//                               (hAngularProfile->GetXaxis()->GetBinWidth(binCosThetaHigh));
//
//        Int_t binS = hAngularProfile->GetYaxis()->FindBin(sLow);
//        // Int_t gBin = hAngularProfile->GetBin(binCosTheta,binS,0); //see doc of TH1::GetBin
//        if(binCosThetaHigh >= binCosThetaLow) // The high and low refer to s not cosTheta - which way round the theta                                 
//        {                                     // bins are depends on the particular PMT/track geometry being considered
//          for(int sumBins = binCosThetaLow; sumBins <= binCosThetaHigh; ++sumBins)
//          {
//            g += hAngularProfile->GetBinContent(sumBins, binS);
//     //       if(g) std::cout << g << std::endl;
//          }
//    //      if( binCosThetaLow != binCosThetaHigh) std::cout << "It matters!" << binCosThetaHigh - binCosThetaLow << std::endl;
//        }
//        else
//        {
//          for(int sumBins = binCosThetaHigh; sumBins <= binCosThetaLow; ++sumBins)
//          {
//            g += hAngularProfile->GetBinContent(sumBins, binS);
//  //          if(g) std::cout << g << std::endl;
//          }
//    //      if( binCosThetaLow != binCosThetaHigh) std::cout << "It matters!" << binCosThetaLow - binCosThetaHigh << std::endl;
//        }
//          g += hAngularProfile->GetBinContent(binCosThetaLow, binS) * (fractionLow - 1);
//          g += hAngularProfile->GetBinContent(binCosThetaHigh, binS) * (fractionHigh - 1);
//
//


//        if(0.7 < cosTheta && 00.8 > cosTheta && 125 > s)
//        {
//          printf("Getting information... \n s = %f \n cosTheta = %f \n cosThetaBin = %d \n sBin = %d \n rho = %f \n G = %f \n PMT pos = ", s, cosTheta, binCosTheta, binS, rho, g);
//          pmtPos.Print();
//          printf("\n vtxPos = ");
//          vtxPos.Print();
//          printf("\n vtxDir = ");
//          vtxDir.Print();
//        }
//        if(g) std::cout << "g = " << g << "   rho = " << rho << "   binWidth = " << hProfile->GetBinWidth(iBin) << std::endl;



        integrals[0] += rho * g * hProfile->GetBinWidth(iBin);
        integrals[1] += rho * g * s * hProfile->GetBinWidth(iBin);
        integrals[2] += rho * g * s * s * hProfile->GetBinWidth(iBin);
      }
//   std::cout << "s = " << fCutoffIntegral << "   " << integrals[0] << "   " << integrals[1] << "   " << integrals[2] << std::endl;
  
    return integrals;
}

///////////////////////////////////////////////////////////////////////////////////////////
// As for the previous function, but now we're calculating the integrals for indirect light
///////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack, int i)
{
  
	std::vector<Double_t> integrals = this->CalculateIndIntegrals(myTrack);
	if( i < 0 || i > 2) 
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateIndIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    return integrals.at(i);
}

std::vector<Double_t> WCSimLikelihoodTuner::CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack)
{ 
    std::vector<Double_t> integrals;
    for(int i = 0; i < 3; ++i)
    {
    	integrals.push_back(0.0);
    }

    this->LoadEmissionProfiles(myTrack);
   Int_t whichBin = fWhichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin >fWhichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateChIntegrals() is looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }

    //std::cout << "We want number " << whichBin << std::endl;
    TH1D * hProfile = (TH1D*)fHistArray->At(whichBin);
    //std::cout << "We got it" << hProfile << std::endl;
    
    // Check to see if we need to curtail the integrals: 
    CalculateCutoff( myTrack );
    int cutoffBin = hProfile->GetNbinsX();
    if( fCutoffIntegral > 0.0 && fCutoffIntegral < hProfile->GetXaxis()->GetXmax() && fCutoffIntegral > hProfile->GetXaxis()->GetXmin() )
    {
      cutoffBin = hProfile->FindBin(fCutoffIntegral);
    }

    // Now do the integrals
    for(int iBin = 1; iBin <= cutoffBin; ++iBin)
    {
        Double_t rho, s;
        s = hProfile->GetBinCenter(iBin);
        if( s > fCutoffIntegral && fCutoffIntegral >= 0.0)
        {
          break;
        }
        rho = hProfile->GetBinContent(iBin);
//        if( rho ) std::cout << "Indirect: s = " << s << " and rho = " << rho << std::endl;
		integrals[0] += rho * hProfile->GetBinWidth(iBin);
        integrals[1] += rho * s * hProfile->GetBinWidth(iBin);
        integrals[2] += rho * s * s * hProfile->GetBinWidth(iBin);
    }
//	std::cout << "Is it 1? " << integrals[0] << std::endl; 
	return integrals;

}

///////////////////////////////////////////////////////////////////////////////////////
// Tabulate the integrals over the emission profiles.  Indirect light does not need the 
// angular emission profile so this one is only indexed by energy.  This should only
// need running once - from then on we can just look the integrals up in the table.
///////////////////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::TabulateIndirectIntegrals( WCSimLikelihoodTrack::TrackType myType, TString filename)
{

    // This might be quite big, so a vector puts it on the heap
    std::vector<Double_t> integrals; 


    //  The structure of how the data is stored is as follows:
    //      Data saved in a TTree
    //      The tree contains one branch, which is a vector of doubles
    //      We write a new entry to the tree for each track energy
    //      The vector itself contains 3 * number of s bins entries
    //      These are groups of three numbers in increasing order of s
    //      ... which correspond to the integrals of rho, s*rho, s*s*rho
    
    std::ofstream outFile;
    outFile.open("integralsRho.dat", std::ios::out|std::ios::binary|std::ios::app);
    outFile.write(reinterpret_cast<char*>(&fNEBinsRho),sizeof(fNEBinsRho));
    outFile.write(reinterpret_cast<char*>(&fNSBinsRho),sizeof(fNSBinsRho));
    outFile.write(reinterpret_cast<char*>(&fEMinRho),sizeof(fEMinRho));
    outFile.write(reinterpret_cast<char*>(&fEMaxRho),sizeof(fEMaxRho));
    outFile.write(reinterpret_cast<char*>(&fSMinRho),sizeof(fSMinRho));
    outFile.write(reinterpret_cast<char*>(&fSMaxRho),sizeof(fSMaxRho));
    
    std::ofstream checkWriting;
    checkWriting.open("whatWeWroteRho.txt");    
  	for( int iEBin = 0; iEBin < fNEBinsRho; ++iEBin)
    {
      integrals.clear();
      std::cout << iEBin << std::endl;
		  Double_t E = fEMin + iEBin*((fEMax - fEMin)/fNEBins);
      TH1D * hProfile = (TH1D*)fHistArray->At(fWhichHisto->FindBin(E) - 1);
//      hProfile->Print();
      
      // An array to keep track of the intermediate values of the integrals
      Double_t rhoS[3] = {0.0, 0.0, 0.0};
		  for( Int_t iSBin = 0; iSBin < fNSBinsRho; ++iSBin)
		  {   
        std::cout << "iSBin = " << iSBin << "/" << fNSBinsRho << std::endl;
			  Double_t s = hProfile->GetBinCenter(iSBin+1);
			  Double_t sWidth = hProfile->GetBinWidth(iSBin+1);
        Double_t rho = hProfile->GetBinContent(iSBin+1);
        rhoS[0] += rho * sWidth;
        rhoS[1] += rho * s * sWidth;
        rhoS[2] += rho * s * s * sWidth;
        // Save the integrals at every step
        outFile.write(reinterpret_cast<char*>(&rhoS),sizeof(rhoS));
        checkWriting << rhoS[0] << " " << rhoS[1] << " " << rhoS[2] << " ";
        if((iSBin % 20) == 0 && iSBin != 0) checkWriting << std::endl;
      }  	  	
    }
    checkWriting.close();
	  outFile.close();
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tabulate the integrals over the emission profiles.  Direct (Cherenkov) light requires the angular emission
// profile as well so we need to be careful with how we index the array.  The evolution of theta as a function
// of s is completely determined by the distance from the vertex to PMT and the angle of the trajectory to the
// vector from vertex to PMT, so we will bin in particle energy, and these two variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::TabulateDirectIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename)
{

	  Double_t rhoGS[3] = {0.,0.,0.};

    // We need this temporary array to let us loop over s at fixed R0, cosTheta0
    // But tabulate looping over R0 and cosTheta0 at fixed s
    Double_t ** tempArray = new Double_t *[fNSBins*fNEBins];
    for( int i = 0; i < fNSBins*fNEBins; ++i )
    {
      tempArray[i] = new Double_t[3 * fNR0Bins * fNCosTheta0Bins]; 
    }
    
    this->LoadEmissionProfiles(myType);    
    
    std::ofstream outFile;
    outFile.open("integralsRhoG.dat", std::ios::out|std::ios::binary|std::ios::app);
    outFile.write(reinterpret_cast<char*>(&fNEBins),sizeof(fNEBins));
    outFile.write(reinterpret_cast<char*>(&fNR0Bins),sizeof(fNR0Bins));
    outFile.write(reinterpret_cast<char*>(&fNCosTheta0Bins),sizeof(fNCosTheta0Bins));
    outFile.write(reinterpret_cast<char*>(&fNSBins),sizeof(fNSBins));
    outFile.write(reinterpret_cast<char*>(&fEMin),sizeof(fEMin));
    outFile.write(reinterpret_cast<char*>(&fEMax),sizeof(fEMax));
    outFile.write(reinterpret_cast<char*>(&fR0Min),sizeof(fR0Min));
    outFile.write(reinterpret_cast<char*>(&fR0Max),sizeof(fR0Max));
    outFile.write(reinterpret_cast<char*>(&fCosTheta0Min),sizeof(fCosTheta0Min));
    outFile.write(reinterpret_cast<char*>(&fCosTheta0Max),sizeof(fCosTheta0Max));
    outFile.write(reinterpret_cast<char*>(&fSMin),sizeof(fSMin));
    outFile.write(reinterpret_cast<char*>(&fSMax),sizeof(fSMax));
           
    std::ofstream checkWriting;
    checkWriting.open("whatWeWroteRhoG.txt");       
    Int_t theCount = 0;                         
    for( int iEBin = 0; iEBin < fNEBins; ++iEBin)
    {
		  Double_t E = fEMin + iEBin*((fEMax - fEMin)/fNEBins);
      TH1D * hProfile = (TH1D*)fHistArray->At(fWhichHisto->FindBin(E)-1); // -1 to convert from bins to array positions
      TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(fWhichHisto->FindBin(E)-1); // as arrays number from 0 but bins from 1

            
		  for( Int_t iR0Bin = 0; iR0Bin < fNR0Bins; ++iR0Bin)
		  {
		    std::cout << "iEBin = " << iEBin << "    iR0Bin = " << iR0Bin << std::endl;
        std::cout << "fR0Max = " << fR0Max << "   fR0Min = " << fR0Min << std::endl; 
			  Double_t R0 = fR0Min + iR0Bin * ((fR0Max - fR0Min)/(Double_t)fNR0Bins) ;
			  for( Int_t iCosTheta0Bin = 0; iCosTheta0Bin < fNCosTheta0Bins; ++iCosTheta0Bin)
			  {
				  Double_t cosTheta0 = fCosTheta0Min + iCosTheta0Bin*((fCosTheta0Max - fCosTheta0Min)/fNCosTheta0Bins);
				  rhoGS[0] = 0.0;
				  rhoGS[1] = 0.0;
				  rhoGS[2] = 0.0;
				  for( Int_t iSBin = 0; iSBin < fNSBins; ++iSBin)
				  {

            // Work out the co-ordinates needed to look up the emission profiles at this step
					  Double_t s = hProfile->GetBinCenter(iSBin+1);
            Double_t sWidth = hProfile->GetBinWidth(iSBin+1); // +1 is array numbering vs. bin numbering
            Double_t cosTheta = (R0*cosTheta0 - s)/TMath::Sqrt(R0*R0 + s*s - 2*R0*s*cosTheta0);

            if( iEBin == 0 && iR0Bin == 14 && iSBin == 250 && iCosTheta0Bin == 46 )
            {
              std::cout << "It's that bin... " << "   E = " << iEBin << "   sMax = " << s << "    cosTheta0 = " << cosTheta0 << "    R0=  " << R0 << std::endl;
            }
            // Now get the values of the emission profiles
            Double_t rho = hProfile->GetBinContent(iSBin+1);
            Int_t binCosTheta = hAngularProfile->GetXaxis()->FindBin(cosTheta);
            Int_t gBin = hAngularProfile->GetBin(binCosTheta,iSBin+1,0); //see doc of TH1::GetBin
            Double_t g = hAngularProfile->GetBinContent(gBin);

            // Do this step of the integral
            rhoGS[0] += rho * g * sWidth;
            rhoGS[1] += rho * g * s * sWidth;
            rhoGS[2] += rho * g * s * s * sWidth;
            outFile.write(reinterpret_cast<char*>(&rhoGS),sizeof(rhoGS));
            checkWriting << rhoGS[0] << " " << rhoGS[1] << " " << rhoGS[2] << " ";
            if(((iSBin + fNSBins*(iCosTheta0Bin + fNCosTheta0Bins*(iR0Bin + fNR0Bins*(iEBin)))) % 20) == 0) checkWriting << std::endl;
				    theCount += 3;
          }	
			  }
		  }
    }
    checkWriting.close();
   outFile.close();
   std::cout << "Total writes = " << theCount << std::endl;
   return;
}


void WCSimLikelihoodTuner::TabulateIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename)
{
  if(filename == "")
  {
    std::cerr << "WCSimLikelihoodTuner::TabulateIntegrals - Error: A filename is needed" << std::endl;
    return;
  }


  TabulateIndirectIntegrals(myType, filename);
  TabulateDirectIntegrals(myType, filename);
  return; 
}



=======
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "TArrayD.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TArrayD.h"
#include "TH2D.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"

#include "WCSimChargeLikelihood.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimRootGeom.hh"
#include "WCSimGeometry.hh"

/////////////////////////////////////////////////////
// Constructor - if the Tuner is created without
// specifying the size of the detector then it 
// assumes the PMTs are floating in space and will 
// only prevent a PMT from seeing light if the track
// is behind it, not if it's outside the tank on the 
// opposite side.               
////////////////////////////////////////////////////
WCSimLikelihoodTuner::WCSimLikelihoodTuner( Bool_t calculateIntegrals = kFALSE)
{
  	fConstrainExtent = kFALSE;
  	fExtent[0] = -1.0;
  	fExtent[1] = -1.0;
  	fExtent[2] = -1.0;
  
  
  	fCalculateIntegrals = calculateIntegrals;
  	this->Initialize();

}

/////////////////////////////////////////////////////
// Otherwise if the extent of the detector in x, y 
// and z is given, the tuner will kill tracks as soon
// as they go outside of this region
/////////////////////////////////////////////////////
WCSimLikelihoodTuner::WCSimLikelihoodTuner(Double_t xMax, Double_t yMax, Double_t zMax, Bool_t calculateIntegrals = kFALSE)
{
  	fConstrainExtent = kTRUE;
  	fExtent[0] = xMax;
  	fExtent[1] = yMax;
  	fExtent[2] = zMax;
    std::cout << fExtent[0] << " " << fExtent[1] << " " << fExtent[2] << std::endl;
  
  	fCalculateIntegrals = true;
  	this->Initialize();
//    std::cout << "Done!!" << std::endl;
}

///////////////////////////////////////////////////////////
// These values should all be the same regardless of
// the constructor called.  Later we should teach it to
// distinguish between track types and open the appropriate
// emission profile file
///////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::Initialize()
{
  fIsOpen = WCSimLikelihoodTrack::Unknown;
  fProfileLocation = new TString("./config/emissionProfilesMuon.root");
  fProfiles = new TFile(fProfileLocation->Data());
  fHistArray = 0;
  fAngHistArray= 0;
  fFluxArray= 0;
  fIsOpen = WCSimLikelihoodTrack::MuonLike;
  fHistArray = (TObjArray *) fProfiles->Get("histArray");
  fHistArray->SetOwner(kTRUE);
  fAngHistArray = (TObjArray *) fProfiles->Get("angHistArray");
  fAngHistArray->SetOwner(kTRUE);
  fFluxArray = (TObjArray *) fProfiles->Get("fluxArray");
  fFluxArray->SetOwner(kTRUE);
  fWhichHisto = (TH1D *) fProfiles->Get("hWhichHisto");
  fAverageQE = 1.0;

  // Pointer to the last track for which we calculated the cutoff, to prevent repetition
  fLastCutoff = 0;
  
  //  std::cout << "Setting up files " << std::endl;
  fRhoGIntegralFile = new TFile();
	fRhoIntegralFile = new TFile();
	fRhoGIntegralTree = 0;
	fRhoIntegralTree = 0;
  fRhoGIntegrals = 0;
  fRhoIntegrals = 0;
    
  // The binning scheme for the direct integral tables
	fNR0Bins = 400;
	fNCosTheta0Bins = 500;
	fNSBins = 400;
	fNEBins = 1;
  fNBinsRhoG = fNR0Bins * fNCosTheta0Bins * fNSBins * fNEBins;
  fR0Min = 0.0;
	fR0Max = 4000.0;
  fCosTheta0Min = 0.0;
  fCosTheta0Max = 1.0;
	fEMin = 250.0;
	fEMax = 3250.0;
  fSMin = 0.0;
	fSMax = 2500.0; 
  fIntegralEnergyBin = -1;
  fIntegralSMax = -1;
  fIntegralParticleType = WCSimLikelihoodTrack::Unknown;
  
  // The binning scheme for the indirect integral tables
  fNSBinsRho = fNSBins;
  fNEBinsRho = fNEBins;
  fNBinsRho = fNSBinsRho * fNEBinsRho;
  fSMinRho = fSMin;
  fSMaxRho = fSMax;
  fEMinRho = fEMin;
  fEMaxRho = fEMax;
  
//  std::cout << "Done!" << std::endl;
  return;




}
/////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////
WCSimLikelihoodTuner::~WCSimLikelihoodTuner()
{
    //std::cout << " *** WCSimLikelihoodTuner::~WCSimLikelihoodTuner() *** Cleaning up" << std::endl;
    fProfiles->Close();
    fIsOpen = WCSimLikelihoodTrack::Unknown;
    delete fProfileLocation;
    delete fProfiles;
    if(fHistArray) delete fHistArray;
    if(fAngHistArray) delete fAngHistArray;
    if(fFluxArray) delete fFluxArray;
//    if(fWhichHisto) delete fWhichHisto;
    fProfileLocation = 0;
    fProfiles = 0;
    fHistArray = 0;
    fAngHistArray= 0;
    fFluxArray= 0;
    
    if(fRhoGIntegralTree) delete fRhoGIntegralTree;
    if(fRhoIntegralTree) delete fRhoIntegralTree;
    if(fRhoGIntegralFile) delete fRhoGIntegralFile;
	  if(fRhoIntegralFile) delete fRhoIntegralFile;
    if(fRhoGIntegrals) delete fRhoGIntegrals;
    if(fRhoIntegrals) delete fRhoIntegrals;

}

void WCSimLikelihoodTuner::UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray)
{
  fConstrainExtent = true; 
  fExtent[0] = myDigitArray->GetExtent(0);
  fExtent[1] = myDigitArray->GetExtent(1);
  fExtent[2] = myDigitArray->GetExtent(2);
  return;
}


////////////////////////////////////////////////////////////////////////
// Load the appropriate emission profiles for this track's particle type
////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadEmissionProfiles( WCSimLikelihoodTrack * myTrack )
{
  this->LoadEmissionProfiles(myTrack->GetType());
  return;
}

///////////////////////////////////////////////////////////////
// Load the appropriate emission profile for a given track type
///////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadEmissionProfiles( WCSimLikelihoodTrack::TrackType myType )
{
  //std::cout << " *** WCSimLikelihoodTuner::LoadEmissionProfiles - Loading profile" << std::endl;

  if( myType == fIsOpen )
  {
    //std::cout << "Is already open" << std::endl;
    return;
  }

  if( fProfileLocation != NULL)
  {
    delete fProfileLocation;
    fProfileLocation = 0;
  }

  switch( myType )
  {
    case WCSimLikelihoodTrack::ElectronLike:
      fProfileLocation = new TString("./config/emissionProfilesElectron.root");
      break;
    case WCSimLikelihoodTrack::MuonLike:
      fProfileLocation = new TString("./config/emissionProfilesMuon.root");
      break;
    default:
      std::cerr << "Error: unkown track type in WCSimLikelihoodTuner::LoadEmissionProfiles" << std::endl;
      exit(EXIT_FAILURE);
  }
 
  if( fProfiles != NULL )
  {
    //std::cout << fProfiles << std::endl;
    delete fProfiles;
    delete fHistArray;
    delete fAngHistArray;
    delete fFluxArray;
    //std::cout << " Was null" << fProfiles << std::endl;
  }

  fProfiles = new TFile(fProfileLocation->Data(),"READ");
  fIsOpen = myType;
  fHistArray = (TObjArray *) fProfiles->Get("histArray");
  fHistArray->SetOwner(kTRUE);
  fAngHistArray = (TObjArray *) fProfiles->Get("angHistArray");
  fAngHistArray->SetOwner(kTRUE);
  fFluxArray = (TObjArray *) fProfiles->Get("fluxArray");
  fFluxArray->SetOwner(kTRUE);
  return;
  

}

/////////////////////////////////////////////////////////////////////////////////////////
// Work out the probability of light surviving to the PMT as it travels through the water
/////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::TransmissionFunction(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{

	
	if( s== 0) return 1;
    // First we need the distance from the photon emission to the PMT 
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    Double_t r = (pmtPos - emissionPos).Mag();


    // We'll use a triple exponential to start with
    Double_t nu[3] = {-1.137e-5,-5.212e-4, -4.359e-3}; // nu = 1/Decay length in mm
    Double_t f[3]      = {0.8827, 0.08162, 0.03515};
    Double_t trans=0.0;
    for(int i = 0; i < 3; ++i){ trans+= f[i]*exp(1.0 * 10 * r * nu[i]);}  //Convert to cm -> factor 10
    return trans;
}

///////////////////////////////////////////////////////////////////////////////////////
// The efficiency of the PMT as a function of the incident angle of the photon
// This uses the same formula as WCSim.  Quantum efficiency (QE) is accounted for here.
// Efficiency of digitization is handled elsewhere
///////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::Efficiency(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
	return 1;
    // We need the angle of incidence at the PMT
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    TVector3 pmtToEm = emissionPos - pmtPos;
    TVector3 pmtFace(myDigit->GetFaceX(), myDigit->GetFaceY(), myDigit->GetFaceZ());

    Double_t cosTheta = pmtFace.Dot(pmtToEm) / (pmtFace.Mag() * pmtToEm.Mag()); 

	// The MiniBooNE method:
    Double_t theta = TMath::ACos(cosTheta) * 180. / TMath::Pi();
    if( theta > 90.0 )
    {
      theta = 180.0 - theta;
    }
	
	Double_t efficiency =  (1 + (-1.182e-4) * pow(theta, 2) + 4.959e-9 * pow(theta, 4) - 7.371e-14 * pow(theta, 6));
	return efficiency;

    // Function for the PMT acceptance's dependence on the angle: 
    // WCSim defines arrays of efficiency at 10 degree intervals and linearly interpolates

 /*   if( cosTheta < 0.0 || cosTheta > 1.0 )
    {

//      std::cout << "Behind the PMT, cosTheta = " << cosTheta << std::endl;
//      pmtPos.Print() ;
//      emissionPos.Print();
      return 0.0;
    }

    Double_t theta = TMath::ACos(cosTheta) * 180.0 / TMath::Pi();
    Double_t collection_angle[10]={0.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
    Double_t collection_eff[10]={100.,100.,99.,95.,90.,85.,80.,69.,35.,13.}; 
    Int_t num_elements = sizeof( collection_angle ) / sizeof( collection_angle[0] );

    Double_t efficiency = 0.0;
    for(int iEntry = 0; iEntry < num_elements-1; ++iEntry)
    {
      if( theta >= collection_angle[iEntry] && theta < collection_angle[iEntry+1])
      { 
        efficiency = collection_eff[iEntry] 
                     + (theta - collection_angle[iEntry]) / (collection_angle[iEntry+1] - collection_angle[iEntry])
                     * (collection_eff[iEntry+1] - collection_eff[iEntry]);
      }
    }
      return efficiency/100.;
*/






}

/*
//////////////////////////////////////////////////////////////////////////////////////////
// Scale the expected photons by the quantum efficiency, which for now we treat as a 
// constant by averaging the QE over the wavelength spectrum for all emitted photons
// and scaling the predicted charge by this number
//////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::QuantumEfficiency(WCSimLikelihoodTrack * myTrack)
{
	WCSimLikelihoodTrack::TrackType type = myTrack->GetType();
	if(type == WCSimLikelihoodTrack::MuonLike) return 0.05;
	else if(type == WCSimLikelihoodTrack::ElectronLike) return 1.0;
	else return 1.0;
}
*/

//////////////////////////////////////////////////////////////////////////////////////////
// Work out the effect of solid angle on the probability of photon incidenc; pure geometry
//////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::SolidAngle(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
    // Now some physical parameters of the phototube
	Double_t WCSimPMTRadius = 12.7;
//    Double_t WCSimPMTExposeHeight = WCSimPMTRadius - 1.0;

	
    // We need the distance from emission to the PMT
    // These are from WCSim so they're in cm
    TVector3 pmtPos(myDigit->GetX() + myDigit->GetFaceX() * WCSimPMTRadius,
    				myDigit->GetY() + myDigit->GetFaceY() * WCSimPMTRadius, 
    				myDigit->GetZ() + myDigit->GetFaceZ() * WCSimPMTRadius);
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    Double_t r = (pmtPos - emissionPos).Mag();

    // Purely geometry: we need the solid angle of a cone whose bottom is a circle of the same radius as the circle of PMT poking
    // through the blacksheet.  This is 2pi( 1 - cos(coneAngle) ) where coneAngle can be deduced from trig and the PMT geometry
	// So the fraction of total solid angle is this/4pi = 0.5(1 - cos(coneAngle))
	Double_t solidAngle = 2.0*TMath::Pi()*(1.0 - (r)/sqrt( (r*r + WCSimPMTRadius * WCSimPMTRadius)));
	
	// I multiply by cosTheta to account for the angle from track to PMT diminishing its effective area
	// This isn't strictly correct (solid angle of a sphere cut off by a plane turns out to be quite complicated
	// but it's a decent approximation to a small effect
//	solidAngle *= TMath::Cos((emissionPos - pmtPos).Angle(TVector3(myDigit->GetFaceX(), myDigit->GetFaceY(), myDigit->GetFaceZ())));

	// See whether the PMT covers more than 1 bin in cosTheta - if it does we should have picked multiple theta bins from g(cosTheta)
	// so scale the solid angle accordingly.
//    this->LoadEmissionProfiles(myTrack);
////	std::cout << fWhichHisto->FindBin(myTrack->GetE());
//	TH2D * hG = (TH2D*)fAngHistArray->At(fWhichHisto->FindBin(myTrack->GetE()) - 1);
////	std::cout << fWhichHisto->FindBin(myTrack->GetE() - 1) << "   " << hG << std::endl;
//	Double_t cosThetaBins = hG->GetNbinsX();
//	Double_t cosThetaRange = (hG->GetXaxis())->GetXmax() - (hG->GetXaxis())->GetXmin();
//	Double_t dCosTheta = cosThetaRange / cosThetaBins;
//	
//	Double_t angleToCenter = pmtPos.Angle(emissionPos);
//	Double_t rangeThetaPMT = TMath::ACos(1.0 - 0.5*solidAngle / TMath::Pi());
////	std::cout << "angle to center " << angleToCenter << "     rangeThetaPMT " << rangeThetaPMT << std::endl;
//	Double_t rangeCosThetaPMT = 2 * TMath::Sin(angleToCenter) * TMath::Sin(rangeThetaPMT);
//
//
//	std::cout << "pmt range = " << rangeCosThetaPMT << " dCosTheta = " << dCosTheta << "  range / dCosTheta = " << rangeCosThetaPMT / dCosTheta << std::endl;
//	solidAngle *= rangeCosThetaPMT / dCosTheta;
	/*if(pmtPos.Z() > 1000)
	{
		std::cout << "cosTheta = " << TMath::Cos(TMath::ATan2(WCSimPMTRadius, r)) << std::endl;
    	std::cout << "solid angle = " << solidAngle << std::endl;
    	std::cout << "r = " << r << "    WCSimPMTRadius = " << WCSimPMTRadius << "   solid angle  = " << solidAngle << "    hit Q = " << myDigit->GetQ() << std::endl;;
    	std::cout << "Pmt pos = " ;
    	pmtPos.Print();
    	std::cout << "Emission pos = " ;
    	emissionPos.Print();
    	std::cout << std::endl;
	}*/
    return solidAngle;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////
// The scattering table assumes that scattered light can be modelled by scaling the
// number of photons you would see from an isotropic source of Cherenkov light by
// a factor dependent on geometry, energy etc.
// It's small and computationally intensive to tabulate, so we're taking a constant initially
/////////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::ScatteringTable(Double_t s)
{
    // 20% scattered light for now
    return 0.001;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//  We define J as the product of the transmission, efficiency and solid angle fuctions (and
//  also the scattering table for indirect light).  However we need it to be zero once the
//  particle has exited the detector, so we check for that before calling the other functions
/////////////////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::CalculateJ( Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit ) 
{
    std::vector<Double_t> J;

    // Check make sure the particle is still inside the detector
	if( fConstrainExtent )
	{
	  TVector3 pos = myTrack->GetVtx() + s*myTrack->GetDir();
	  if( fabs(pos.X()) > fExtent[0] || fabs(pos.Y()) > fExtent[1] || fabs(pos.Z()) > fExtent[2])
      {
        J.push_back(0.0);
        J.push_back(0.0);
        return J;
      }
    }


    // Work out the direct and indirect contributions to J
    // J[0] = J_dir, J[1] = J_ind
    if( s == 0.0)
    {
    	std::cout << "Transmission = " << this->TransmissionFunction(s, myTrack, myDigit) << std::endl
    			  << "Efficiency = "   << this->Efficiency(s, myTrack, myDigit) << std::endl
    			  << "SolidAngle = " << this->SolidAngle(s, myTrack, myDigit) << std::endl;
    }
    J.push_back(   this->TransmissionFunction(s, myTrack, myDigit) 
                 * this->Efficiency(s, myTrack, myDigit)
            //     * this->QuantumEfficiency(myTrack)
                 * this->SolidAngle(s, myTrack, myDigit));
    J.push_back(J.at(0) * this->ScatteringTable(s));
    return J;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The direct and indirect expectations come from a large integral over several factors.  We separate out the
// transmission function, efficiency and solid angle and fit these to a quadratic in s (distance along track)
// This just leaves integrals over emission profiles that can be tabulated in advance.  This function gives
// the coefficients of these quadratics
///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void WCSimLikelihoodTuner::CalculateCoefficients(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit )
{
    
//    std::cout << "*** WCSimLikelihoodTuner::CalculateCoefficients() *** Calculating the coefficients for the integrals from tuning info" << std::endl; 
    for(int i = 0; i < 3; ++i)
    { 
        fDirCoeffs[i] = 0.0;    // Coefficients for direct (Cherenkov) light
        fIndCoeffs[i] = 0.0;    // And for indirect light
     }
    
    // Load the file containing the emission profiles for different energies
    // whichHisto should be a TH1I with the same energy binning as we want for the reconstruction
    // There should also be a TObjArray with the emission profiles for the center value of each energy bin
//    std::cout << "Getting the histogram" << std::endl;
    this->LoadEmissionProfiles(myTrack);
 
    Int_t whichBin = fWhichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin > fWhichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateCoefficients() looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }

//    std::cout << "We want number " << whichBin << std::endl;
//    std::cout << "fHistArray= " << fHistArray << std::endl;
    TH1D * hProfile = (TH1D*)(fHistArray->At(whichBin));
    // Calculate the 3 s values we care about: s=0, s where integral(rho(s') ds')_0^{s} = 0.75, and double that
    Double_t runningTotal=0.;
    Double_t s[3];
    
    for( int iBin = 1; iBin <= hProfile->GetNbinsX(); ++iBin )
    {
        runningTotal += hProfile->GetBinContent(iBin) * hProfile->GetBinWidth(iBin);
        if(runningTotal >= 0.75)
        {
            //std::cout << "Setting s[i]" << std::endl;
            s[0] = 10.0;
            s[1] = hProfile->GetBinCenter(iBin);
            s[2] = 2 * s[1];
            break;
        }
    }
   
    // Check these s values are inside the detector:
    this->CalculateCutoff(myTrack);
    if( fCutoffIntegral > 0.0 && s[2] > fCutoffIntegral )
    {
      s[2] = 0.8*fCutoffIntegral;
      s[1] = 0.5 * s[2];
    }

    // Evaluate J at each point
    Double_t JDir[3] = {0.,0.,0.};
    Double_t JInd[3] = {0.,0.,0.};
        
        // std::cout << " Second loop: kk = " << kk << "    s = " << s[kk] << "     JDir = " << JDir[kk] << std::endl;
    for(int k = 0; k < 3; ++k)
    {
        std::vector<Double_t> J = this->CalculateJ( s[k], myTrack, myDigit);
        JDir[k] = J.at(0); 
        JInd[k] = J.at(1);
    }
    // Now do a quadratic fit to the Cherenkov and indirect parts, sampling them at these 3 points
    // We can actually work this out analytically:    
    
    fDirCoeffs[0] = JDir[0];
    fDirCoeffs[1] = (4.*JDir[1] - JDir[2] -3.*JDir[0]) / (2*s[1]);
    fDirCoeffs[2] = (JDir[0] + JDir[2] - 2.*JDir[1]) / (2 * s[1] * s[1]);
//    std::cout << "jDir:        " << JDir[0]       << "," << JDir[1]       << "," << JDir[2]       << std::endl;
//    std::cout << "fDirCoeffs:  " << fDirCoeffs[0] << "," << fDirCoeffs[1] << "," << fDirCoeffs[2] << std::endl; 
//    

    fIndCoeffs[0] = JInd[0];
    fIndCoeffs[1] = (4.*JInd[1] - JInd[2] -3.*JInd[0]) / (2*s[1]);
    fIndCoeffs[2] = (JInd[0] + JInd[2] - 2.*JInd[1]) / (2 * s[1] * s[1]);

//    if((myDigit->GetTubeId() % 10000) == 0)
//    {
//      TGraph * gr = new TGraph();
//      for(int i = 0; i < hProfile->GetNbinsX(); ++i)
//      {
//        double tempS = hProfile->GetXaxis()->GetBinCenter(i);
//        std::vector<Double_t> tempJ = this->CalculateJ(tempS, myTrack, myDigit);
//        gr->SetPoint(gr->GetN(), tempS, tempJ.at(1));
//      }
//      TF1 * jInd = new TF1("jInd","[0] + [1] * x + [2] * x * x", 0, hProfile->GetXaxis()->GetXmax());
//      jInd->SetParameters(fIndCoeffs[0], fIndCoeffs[1], fIndCoeffs[2]);
//  
//      TCanvas * c1 = new TCanvas("c1","c1", 800, 600);
//      gr->SetLineColor(kBlue);
//      gr->SetLineWidth(2);
//      gr->Draw("ALP");
//      jInd->SetLineColor(kRed);
//      jInd->SetLineWidth(2);
//      jInd->Draw("same");
//      TString str;
//      str.Form("jInd%d.png",myDigit->GetTubeId());
//      c1->SaveAs(str.Data());
//      delete c1;
//      delete jInd;
//      delete gr;
//    }  
  
  
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Same as the function above, but returns a vector of the coefficients as well as updating the class member variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::CalculateCoefficientsVector(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit )
{
    this->CalculateCoefficients(myTrack, myDigit);
    
    // Return these in a vector
    std::vector<Double_t> coeffs;
    for(int j = 0; j < 3; ++j){ coeffs.push_back(fDirCoeffs[j]); }
    for(int j = 0; j < 3; ++j){ coeffs.push_back(fIndCoeffs[j]); }
    return coeffs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// We shouldn't integrate over any more of the track once the particle has exited the detector.
// Here, we input the track and work out how far it has to travel before it reaches the edge of the
// detector volume, and return that value
// ////////////////////////////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::CalculateCutoff( WCSimLikelihoodTrack * myTrack )
{
//  std::cout << "Calculating the cutoff" << std::endl;
  if(myTrack == fLastCutoff) return;
  
  Double_t cutoff = fSMax;

  if( fConstrainExtent )
  {
    TVector3 vtx = myTrack->GetVtx();
    TVector3 dir = myTrack->GetDir();
    Double_t sMax[3] = {0., 0.,0.};
    for(int i = 0; i < 3; ++i)
    {
      if(dir[i] > 1e-6)
      {
        sMax[i] = -1. * (vtx(i) / dir(i)) + TMath::Sqrt( fExtent[i] * fExtent[i] / (dir(i)*dir(i)));
        std::cout << sMax[i] << std::endl;
      }
      else sMax[i] = cutoff;
    }
  
    // The world's laziest sorting algorithm:
    if( sMax[0] < cutoff) cutoff = sMax[0];
    if( sMax[1] < cutoff ) cutoff = sMax[1];
    if( sMax[2] < cutoff ) cutoff = sMax[2];
  }
  std::cout << "Cutoff = " << cutoff << ".... returning" << std::endl;
  fCutoffIntegral = cutoff;
  fLastCutoff = myTrack;
  return;
}

///////////////////////////////////////////////////////////////////////////
//	Work out what array entry we need to query to look up the integrals
//////////////////////////////////////////////////////////////////////////
Int_t WCSimLikelihoodTuner::GetEBin(Double_t energy)
{
  if(energy == fEMax) return fNEBins-1;
  return static_cast<Int_t>( (energy - fEMin) / ((fEMax - fEMin)/(Double_t)fNEBins) );
}

Int_t WCSimLikelihoodTuner::GetSBin(Double_t sMax)
{  
  if(sMax >= fSMax) return 3*(fNSBins-1);
  return 3*static_cast<Int_t>( (sMax - 0) / (fSMax/(Double_t)fNSBins ));
}


Int_t WCSimLikelihoodTuner::GetESBin(Double_t energy, Double_t sMax)
{  
  Int_t eBin = (energy >= fEMaxRho)? fNEBinsRho-1 : static_cast<Int_t>( (energy - fEMinRho) / ((fEMaxRho-fEMinRho)/(Double_t)fNEBinsRho ));
  Int_t sBin = (sMax >= fSMaxRho)? fNSBins-1 : static_cast<Int_t>( (sMax - fSMinRho) / ((fSMaxRho-fSMinRho)/(Double_t)fNSBinsRho ));
//  std::cout << "energy = " << energy << "   so eBin = " << eBin << std::endl
//            << "sMax   = " << fSMax  << "   so sBin = " << sBin << std::endl
//            << "ESBin  = " << sBin + fNSBins * eBin << std::endl;
  return 3*(sBin + fNSBinsRho * eBin);
}

Int_t WCSimLikelihoodTuner::GetIntegralBin(Double_t E, Double_t sMax, Double_t R0, Double_t cosTheta0)
{
	  Int_t eBin = (E >= fEMax) ? fNEBins-1 : static_cast<Int_t>( (E - fEMin) / ((fEMax - fEMin)/(Double_t)fNEBins) );
    Int_t sBin = (sMax >= fSMax) ? fNSBins-1 : static_cast<Int_t>( (sMax - fSMin) / ((fSMax-fSMin)/(Double_t)fNSBins ));
    Int_t R0Bin = (R0 >= fR0Max) ? fNR0Bins-1 : static_cast<Int_t>( (R0 - fR0Min) / ((fR0Max-fR0Min)/(Double_t)fNR0Bins ));
    Int_t cosTheta0Bin = ( cosTheta0 >= fCosTheta0Max)? fNCosTheta0Bins-1 : static_cast<Int_t>( (cosTheta0 - fCosTheta0Min) / ((fCosTheta0Max - fCosTheta0Min)/(Double_t)fNCosTheta0Bins ));
//    std::cout << "Values are " << "Energy = " << E << "   sMax = " << sMax << "   R0 = " << R0 << "   cosTheta0" << cosTheta0 << std::endl;
//    std::cout << "Bins are " << eBin << "/" << fNEBins << "    " << sBin << "/" << fNSBins << "    " << R0Bin << "/" << fNR0Bins << "     " << cosTheta0Bin << "/" << fNCosTheta0Bins <<  std::endl;
    
   	return  3*(sBin + fNSBins*(cosTheta0Bin + fNCosTheta0Bins*(R0Bin + fNR0Bins *( eBin ))));
}

///////////////////////////////////////////////////////////////////////////
// We decouple the direct part into a source factor (flux x profile) and 
// an acceptance factor (solid angle x transmittance x PMT acceptance)
// We then integrate over the source factor, but it's much quicker to 
// tabulate this in advance.  This function opens the file containing those
// tables
///////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadTabulatedIntegrals( WCSimLikelihoodTrack * myTrack )
{
  // Check we have the same track type, energy and s cutoff as before
	WCSimLikelihoodTrack::TrackType myType = myTrack->GetType();
  this->CalculateCutoff(myTrack);

  // In which case we don't need to do anything!
	if( myType == fIntegralParticleType )
	{
//			std::cout << "Returning!" << std::endl;
		return;
	}

		
  std::cout << "Looking up integrals" << std::endl;
  // Otherwise open up the right table file
  ifstream integralFileRhoG;
  ifstream integralFileRho;
	switch(myType)
	{
		case WCSimLikelihoodTrack::MuonLike:
	//		std::cout << "It's muon-like" << std::endl;
      integralFileRhoG.open("config/integralsMuonRhoG.dat",std::ios::in|std::ios::binary);
      integralFileRho.open("config/integralsMuonRho.dat",std::ios::in|std::ios::binary);
      break;
	  	case WCSimLikelihoodTrack::ElectronLike:
	//  		std::cout << "Is it tripping the electronlike flag too?" << WCSimLikelihoodTrack::ElectronLike << std::endl;
      integralFileRhoG.open("config/integralsElectronRhoG.dat",std::ios::in|std::ios::binary);
      integralFileRho.open("config/integralsElectronRho.dat",std::ios::in|std::ios::binary);
	  	break;
		case WCSimLikelihoodTrack::Unknown:
	//		std::cout << "Or is it?" << "   " << myType << std::endl;
	  		std::cerr << "Error: could not identify particle type, exiting" << std::endl;	 			
        exit(EXIT_FAILURE);	
			break;
	}

  if(integralFileRhoG.is_open() == false)
  {
    std::cerr << "Could not open tabulated direct integral file" << std::endl;
    exit(EXIT_FAILURE);
  }
  else
  {
     // The first 12 entries tell us how many bins we have in each variable, and what the minimum and maximum is
    
    integralFileRhoG.read((char *) &fNEBins, sizeof(fNEBins));
    integralFileRhoG.read((char *) &fNR0Bins, sizeof(fNR0Bins));
    integralFileRhoG.read((char *) &fNCosTheta0Bins, sizeof(fNCosTheta0Bins));
    integralFileRhoG.read((char *) &fNSBins, sizeof(fNSBins));
    integralFileRhoG.read((char *) &fEMin, sizeof(fEMin));
    integralFileRhoG.read((char *) &fEMax, sizeof(fEMax));
    integralFileRhoG.read((char *) &fR0Min, sizeof(fEMin));
    integralFileRhoG.read((char *) &fR0Max, sizeof(fEMax));
    integralFileRhoG.read((char *) &fCosTheta0Min, sizeof(fCosTheta0Min));
    integralFileRhoG.read((char *) &fCosTheta0Max, sizeof(fCosTheta0Max));
    integralFileRhoG.read((char *) &fSMin, sizeof(fSMin));
    integralFileRhoG.read((char *) &fSMax, sizeof(fSMax));
       std::cout << fNCosTheta0Bins << "  " << fSMax << std::endl;

    fNBinsRhoG = 3 * fNEBins * fNR0Bins * fNCosTheta0Bins * fNSBins;
    std::cout << fNBinsRhoG;
    if(fRhoGIntegrals) delete fRhoGIntegrals;
    fRhoGIntegrals = new Double_t[fNBinsRhoG];
    integralFileRhoG.read((char *) &fRhoGIntegrals[0], fNBinsRhoG * sizeof(fRhoGIntegrals[0]));
//    ofstream outFile("whatWeReadRhoG.txt");
//    for(int i = 0; i < fNBinsRhoG; i+=3)
//    {
//      outFile << fRhoGIntegrals[i] << " ";
//      outFile << fRhoGIntegrals[i+1] << " ";
//      outFile << fRhoGIntegrals[i+2] << " ";
//      if((i % 60) == 0) outFile << std::endl;
//    }
//    outFile.close();

  } 
    
  if(integralFileRhoG.is_open() == false)
  {
    std::cerr << "Could not open tabulated indirect integral file" << std::endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    integralFileRho.read((char *) &fNEBinsRho, sizeof(fNEBinsRho));
    integralFileRho.read((char *) &fNSBinsRho, sizeof(fNSBinsRho));
    integralFileRho.read((char *) &fEMinRho, sizeof(fEMinRho));
    integralFileRho.read((char *) &fEMaxRho, sizeof(fEMaxRho));
    integralFileRho.read((char *) &fSMinRho, sizeof(fSMinRho));
    integralFileRho.read((char *) &fSMaxRho, sizeof(fSMaxRho));
   
    std::cout << fNEBinsRho << "  " << fNSBinsRho << "  " << fEMinRho << "  " << fEMaxRho <<"   " << fSMinRho <<"   " << fSMaxRho << std::endl;
   
    if(fRhoIntegrals) delete fRhoIntegrals;
    fNBinsRho = 3 * fNEBinsRho * fNSBinsRho;
    fRhoIntegrals = new Double_t[fNBinsRho];
    integralFileRho.read((char *) &fRhoIntegrals[0], fNBinsRho * sizeof(fRhoIntegrals[0]));
//    ofstream outFile;
//    outFile.open("whatWeReadRho.txt");
//    for(int i = 0; i < fNBinsRho; i+=3)
//    { 
//      outFile << fRhoIntegrals[i] << " ";
//      outFile << fRhoIntegrals[i+1] << " ";
//      outFile << fRhoIntegrals[i+2] << " ";
//      if((i%60)==0 && i != 0) outFile << std::endl;
//    }
//    outFile.close();
  }

  fIntegralParticleType = myType;

}

///////////////////////////////////////////////////////////////////////////
//	Look up the tabulated integrals - just return 1 power of s
///////////////////////////////////////////////////////////////////////////
double WCSimLikelihoodTuner::LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, Int_t sPower)
{
    if( !(sPower < 3 && sPower >= 0) )
    {
        std::cerr << "WCSimLikelihoodTuner::LookupChIntegrals() ERROR: power of s must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }	
    
    std::vector<Double_t> integralsVec = this->LookupChIntegrals(myTrack, myDigit);
    
  if( (UInt_t)sPower > integralsVec.size() ) std::cerr << "There's a problem with integralsVec!" << std::endl;  
	return integralsVec.at(sPower);
}

///////////////////////////////////////////////////////////////////////////
// Look up the tabulated integrals: return a vector of all 3 powers of s
// (this is the fastest way to work out the likelihood)
///////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
//  	std::cout << "*** WCSimLikelihoodTuner::LookupChIntegrals() *** Looking up the tabulated integrals for direct Cherenkov light" << std::endl; 
	this->LoadTabulatedIntegrals( myTrack );

 	
 	//std::cout << "Doing the vector stuff ..." ;
 	TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
 	Double_t R0 = (pmtPos - (myTrack->GetVtx())).Mag();
 	Double_t cosTheta0 = myTrack->GetDir().Dot(pmtPos - myTrack->GetVtx()) / (R0);
 	//std::cout << "Done! ... ";
 	
 	//std::cout << "Checking which bin...";
  Int_t whichBin = this->GetIntegralBin(myTrack->GetE(), fCutoffIntegral, R0, cosTheta0);
  std::cout << "  E      = " << myTrack->GetE() << "  sMax      = " << fCutoffIntegral << " R0     = " << R0 << "    cosTheta0 = " << cosTheta0 << std::endl;
//  std::cout << "Int bin = " << this->GetIntegralBin(myTrack->GetE(), fCutoffIntegral, R0, cosTheta0) << std::endl;
//  std::cout << "Integrals are      i[0]    = " << fRhoGIntegrals[whichBin] << "         i[1] = " << fRhoGIntegrals[1+whichBin] << "       i[2] = " << fRhoGIntegrals[whichBin+2] << std::endl;

// 	std::cout << "Bin = " << whichBin << "   Done! ... ";
 	
	//std::cout << "Making the vector...";
	std::vector<Double_t> integralsVec;
  if( whichBin+2 > fNBinsRhoG ) std::cerr << "There's a problem with fRhoGIntegrals!" << std::endl
                                                      << fNBinsRhoG << "   " << whichBin << std::endl;
	integralsVec.push_back(fRhoGIntegrals[whichBin]);
	integralsVec.push_back(fRhoGIntegrals[whichBin+1]);
	integralsVec.push_back(fRhoGIntegrals[whichBin+2]);
	
 	//std::cout << "Done! ... ";
 	//std::cout << std::endl;
  
//  std::cout << "R0 = " << R0 << "    cosTheta0 = " << cosTheta0 << "     E = " << myTrack->GetE() << "    sMax = " << fCutoffIntegral << "     s term = " << integralsVec.at(1) << "    integral bin = " << whichBin << std::endl;
  
 	return integralsVec;


}

///////////////////////////////////////////////////////////////////////////
//	Look up the tabulated integrals
///////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::LookupIndIntegrals(WCSimLikelihoodTrack * myTrack, Int_t sPower)
{
	if( !(sPower < 3 && sPower >= 0) )
    {
        std::cerr << "WCSimLikelihoodTuner::LookupChIntegrals() ERROR: power of s must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::vector<Double_t> integralsVec = this->LookupIndIntegrals(myTrack);
    if( sPower > integralsVec.size() ) std::cerr << "There's a problem with integralsVec!" << std::endl;
	return integralsVec.at(sPower);
}

std::vector<Double_t> WCSimLikelihoodTuner::LookupIndIntegrals(WCSimLikelihoodTrack * myTrack)
{
//  	std::cout << "*** WCSimLikelihoodTuner::LookupIndIntegrals() *** Looking up the tabulated integrals for indirect light" << std::endl; 
 	  this->LoadTabulatedIntegrals( myTrack );
    
    Int_t ESBin = this->GetESBin(myTrack->GetE(), fCutoffIntegral); 
    std::vector<Double_t> integralsVec;
    if( ESBin > fNBinsRho ) std::cerr << "There's a problem with fRhoIntegrals!" << std::endl;
    integralsVec.push_back(fRhoIntegrals[ESBin]);
    integralsVec.push_back(fRhoIntegrals[ESBin+1]);
    integralsVec.push_back(fRhoIntegrals[ESBin+2]);
    return integralsVec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Work out the integrals over the emission profiles numerically.  This is slow and only good for testing - 
// eventually they'll be tabulated and looked-up
///////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, int i)
{
	std::vector<Double_t> integrals = this->CalculateChIntegrals(myTrack, myDigit);
	if( i < 0 || i > 2) 
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateChIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    return integrals.at(i);
}

std::vector<Double_t> WCSimLikelihoodTuner::CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{


    std::cout << "*** WCSimLikelihoodTuner::CalculateChIntegrals() ***" << std::endl;
    this->LoadEmissionProfiles(myTrack);
    Int_t whichBin = fWhichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin > fWhichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateChIntegrals() is looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }


    //std::cout << "We want number " << whichBin << std::endl;
    TH1D * hProfile = (TH1D*)fHistArray->At(whichBin);
    TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(whichBin);
 
    // Check to see if we need to curtail the integrals:
    std::vector<Double_t> integrals;
    for(int i = 0; i < 3; ++i)
    {
    	integrals.push_back(0.0);
    }
    
    CalculateCutoff( myTrack );
    int cutoffBin = hProfile->GetNbinsX();
    if( fCutoffIntegral > 0.0 && fCutoffIntegral < hProfile->GetXaxis()->GetXmax() && fCutoffIntegral > hProfile->GetXaxis()->GetXmin() )
    {
      cutoffBin = hProfile->FindBin(fCutoffIntegral);
    }

    // Now do the integrals. 
    // Need these to calculate theta; can declare them outside of the loop
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 vtxPos(myTrack->GetVtx());
    TVector3 vtxDir(myTrack->GetDir());
    for(int iBin = 1; iBin <= cutoffBin; ++iBin)
    {
        Double_t rho = 0.;
        Double_t g = 0.; 
        Double_t s = hProfile->GetBinCenter(iBin);
//        Double_t sLow = hProfile->GetBinLowEdge(iBin);
//        Double_t sHigh = sLow + hProfile->GetBinWidth(iBin);
        if( s > fCutoffIntegral && fCutoffIntegral >= 0. ) 
        {
          break;
        }

        Double_t cosTheta = TMath::Cos(vtxDir.Angle( pmtPos - (vtxPos + s * vtxDir)));
//        Double_t cosThetaLow = TMath::Cos(vtxDir.Angle( pmtPos - (vtxPos + sLow * vtxDir) ));
//        Double_t cosThetaHigh = TMath::Cos(vtxDir.Angle( pmtPos - (vtxPos + sHigh * vtxDir) ));

                // Make sure the histograms in s and s, cosTheta(s) have the same binning on the s axis
        rho = hProfile->GetBinContent(iBin);
        g = hAngularProfile->GetBinContent( hAngularProfile->GetXaxis()->FindBin(cosTheta), iBin );

//        Int_t binCosThetaLow = hAngularProfile->GetXaxis()->FindBin(cosThetaLow);
//        Double_t fractionLow = (cosThetaLow - hAngularProfile->GetXaxis()->GetBinLowEdge(binCosThetaLow))/
//                               (hAngularProfile->GetXaxis()->GetBinWidth(binCosThetaLow));
//        Int_t binCosThetaHigh = hAngularProfile->GetXaxis()->FindBin(cosThetaHigh);
//        Double_t fractionHigh = (cosThetaHigh - hAngularProfile->GetXaxis()->GetBinLowEdge(binCosThetaHigh))/
//                               (hAngularProfile->GetXaxis()->GetBinWidth(binCosThetaHigh));
//
//        Int_t binS = hAngularProfile->GetYaxis()->FindBin(sLow);
//        // Int_t gBin = hAngularProfile->GetBin(binCosTheta,binS,0); //see doc of TH1::GetBin
//        if(binCosThetaHigh >= binCosThetaLow) // The high and low refer to s not cosTheta - which way round the theta                                 
//        {                                     // bins are depends on the particular PMT/track geometry being considered
//          for(int sumBins = binCosThetaLow; sumBins <= binCosThetaHigh; ++sumBins)
//          {
//            g += hAngularProfile->GetBinContent(sumBins, binS);
//     //       if(g) std::cout << g << std::endl;
//          }
//    //      if( binCosThetaLow != binCosThetaHigh) std::cout << "It matters!" << binCosThetaHigh - binCosThetaLow << std::endl;
//        }
//        else
//        {
//          for(int sumBins = binCosThetaHigh; sumBins <= binCosThetaLow; ++sumBins)
//          {
//            g += hAngularProfile->GetBinContent(sumBins, binS);
//  //          if(g) std::cout << g << std::endl;
//          }
//    //      if( binCosThetaLow != binCosThetaHigh) std::cout << "It matters!" << binCosThetaLow - binCosThetaHigh << std::endl;
//        }
//          g += hAngularProfile->GetBinContent(binCosThetaLow, binS) * (fractionLow - 1);
//          g += hAngularProfile->GetBinContent(binCosThetaHigh, binS) * (fractionHigh - 1);
//
//


//        if(0.7 < cosTheta && 00.8 > cosTheta && 125 > s)
//        {
//          printf("Getting information... \n s = %f \n cosTheta = %f \n cosThetaBin = %d \n sBin = %d \n rho = %f \n G = %f \n PMT pos = ", s, cosTheta, binCosTheta, binS, rho, g);
//          pmtPos.Print();
//          printf("\n vtxPos = ");
//          vtxPos.Print();
//          printf("\n vtxDir = ");
//          vtxDir.Print();
//        }
//        if(g) std::cout << "g = " << g << "   rho = " << rho << "   binWidth = " << hProfile->GetBinWidth(iBin) << std::endl;



        integrals[0] += rho * g * hProfile->GetBinWidth(iBin);
        integrals[1] += rho * g * s * hProfile->GetBinWidth(iBin);
        integrals[2] += rho * g * s * s * hProfile->GetBinWidth(iBin);
      }
   std::cout << "s = " << fCutoffIntegral << "   " << integrals[0] << "   " << integrals[1] << "   " << integrals[2] << std::endl;
  
    return integrals;
}

///////////////////////////////////////////////////////////////////////////////////////////
// As for the previous function, but now we're calculating the integrals for indirect light
///////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack, int i)
{
  
	std::vector<Double_t> integrals = this->CalculateIndIntegrals(myTrack);
	if( i < 0 || i > 2) 
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateIndIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    return integrals.at(i);
}

std::vector<Double_t> WCSimLikelihoodTuner::CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack)
{ 
    std::vector<Double_t> integrals;
    for(int i = 0; i < 3; ++i)
    {
    	integrals.push_back(0.0);
    }

    this->LoadEmissionProfiles(myTrack);
   Int_t whichBin = fWhichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin >fWhichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateChIntegrals() is looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }

    //std::cout << "We want number " << whichBin << std::endl;
    TH1D * hProfile = (TH1D*)fHistArray->At(whichBin);
    //std::cout << "We got it" << hProfile << std::endl;
    
    // Check to see if we need to curtail the integrals: 
    CalculateCutoff( myTrack );
    int cutoffBin = hProfile->GetNbinsX();
    if( fCutoffIntegral > 0.0 && fCutoffIntegral < hProfile->GetXaxis()->GetXmax() && fCutoffIntegral > hProfile->GetXaxis()->GetXmin() )
    {
      cutoffBin = hProfile->FindBin(fCutoffIntegral);
    }

    // Now do the integrals
    for(int iBin = 1; iBin <= cutoffBin; ++iBin)
    {
        Double_t rho, s;
        s = hProfile->GetBinCenter(iBin);
        if( s > fCutoffIntegral && fCutoffIntegral >= 0.0)
        {
          break;
        }
        rho = hProfile->GetBinContent(iBin);
//        if( rho ) std::cout << "Indirect: s = " << s << " and rho = " << rho << std::endl;
		integrals[0] += rho * hProfile->GetBinWidth(iBin);
        integrals[1] += rho * s * hProfile->GetBinWidth(iBin);
        integrals[2] += rho * s * s * hProfile->GetBinWidth(iBin);
    }
//	std::cout << "Is it 1? " << integrals[0] << std::endl; 
	return integrals;

}

///////////////////////////////////////////////////////////////////////////////////////
// Tabulate the integrals over the emission profiles.  Indirect light does not need the 
// angular emission profile so this one is only indexed by energy.  This should only
// need running once - from then on we can just look the integrals up in the table.
///////////////////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::TabulateIndirectIntegrals( WCSimLikelihoodTrack::TrackType myType, TString filename)
{

    // This might be quite big, so a vector puts it on the heap
    std::vector<Double_t> integrals; 


    //  The structure of how the data is stored is as follows:
    //      Data saved in a TTree
    //      The tree contains one branch, which is a vector of doubles
    //      We write a new entry to the tree for each track energy
    //      The vector itself contains 3 * number of s bins entries
    //      These are groups of three numbers in increasing order of s
    //      ... which correspond to the integrals of rho, s*rho, s*s*rho
    
    std::ofstream outFile;
    outFile.open("integralsRho.dat", std::ios::out|std::ios::binary|std::ios::app);
    outFile.write(reinterpret_cast<char*>(&fNEBinsRho),sizeof(fNEBinsRho));
    outFile.write(reinterpret_cast<char*>(&fNSBinsRho),sizeof(fNSBinsRho));
    outFile.write(reinterpret_cast<char*>(&fEMinRho),sizeof(fEMinRho));
    outFile.write(reinterpret_cast<char*>(&fEMaxRho),sizeof(fEMaxRho));
    outFile.write(reinterpret_cast<char*>(&fSMinRho),sizeof(fSMinRho));
    outFile.write(reinterpret_cast<char*>(&fSMaxRho),sizeof(fSMaxRho));
    
    ofstream checkWriting;
    checkWriting.open("whatWeWroteRho.txt");    
  	for( int iEBin = 0; iEBin < fNEBinsRho; ++iEBin)
    {
      integrals.clear();
      std::cout << iEBin << std::endl;
      Double_t E = 1500;
      TH1D * hProfile = (TH1D*)fHistArray->At(fWhichHisto->FindBin(E) - 1);
//      hProfile->Print();
      
      // An array to keep track of the intermediate values of the integrals
      Double_t rhoS[3] = {0.0, 0.0, 0.0};
		  for( Int_t iSBin = 0; iSBin < fNSBinsRho; ++iSBin)
		  {   
        std::cout << "iSBin = " << iSBin << "/" << fNSBinsRho << std::endl;
			  Double_t s = hProfile->GetBinCenter(iSBin+1);
			  Double_t sWidth = hProfile->GetBinWidth(iSBin+1);
        Double_t rho = hProfile->GetBinContent(iSBin+1);
        rhoS[0] += rho * sWidth;
        rhoS[1] += rho * s * sWidth;
        rhoS[2] += rho * s * s * sWidth;
        // Save the integrals at every step
        outFile.write(reinterpret_cast<char*>(&rhoS),sizeof(rhoS));
        checkWriting << rhoS[0] << " " << rhoS[1] << " " << rhoS[2] << " ";
        if((iSBin % 20) == 0 && iSBin != 0) checkWriting << std::endl;
      }  	  	
    }
    checkWriting.close();
	  outFile.close();
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tabulate the integrals over the emission profiles.  Direct (Cherenkov) light requires the angular emission
// profile as well so we need to be careful with how we index the array.  The evolution of theta as a function
// of s is completely determined by the distance from the vertex to PMT and the angle of the trajectory to the
// vector from vertex to PMT, so we will bin in particle energy, and these two variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::TabulateDirectIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename)
{

	  Double_t rhoGS[3] = {0.,0.,0.};

    // We need this temporary array to let us loop over s at fixed R0, cosTheta0
    // But tabulate looping over R0 and cosTheta0 at fixed s
    Double_t ** tempArray = new Double_t *[fNSBins*fNEBins];
    for( int i = 0; i < fNSBins*fNEBins; ++i )
    {
      tempArray[i] = new Double_t[3 * fNR0Bins * fNCosTheta0Bins]; 
    }
    
    this->LoadEmissionProfiles(myType);    
    
    std::ofstream outFile;
    outFile.open("integralsRhoG.dat", std::ios::out|std::ios::binary|std::ios::app);
    outFile.write(reinterpret_cast<char*>(&fNEBins),sizeof(fNEBins));
    outFile.write(reinterpret_cast<char*>(&fNR0Bins),sizeof(fNR0Bins));
    outFile.write(reinterpret_cast<char*>(&fNCosTheta0Bins),sizeof(fNCosTheta0Bins));
    outFile.write(reinterpret_cast<char*>(&fNSBins),sizeof(fNSBins));
    outFile.write(reinterpret_cast<char*>(&fEMin),sizeof(fEMin));
    outFile.write(reinterpret_cast<char*>(&fEMax),sizeof(fEMax));
    outFile.write(reinterpret_cast<char*>(&fR0Min),sizeof(fR0Min));
    outFile.write(reinterpret_cast<char*>(&fR0Max),sizeof(fR0Max));
    outFile.write(reinterpret_cast<char*>(&fCosTheta0Min),sizeof(fCosTheta0Min));
    outFile.write(reinterpret_cast<char*>(&fCosTheta0Max),sizeof(fCosTheta0Max));
    outFile.write(reinterpret_cast<char*>(&fSMin),sizeof(fSMin));
    outFile.write(reinterpret_cast<char*>(&fSMax),sizeof(fSMax));
           
    ofstream checkWriting;
    checkWriting.open("whatWeWroteRhoG.txt");       
    Int_t theCount = 0;                         
    for( int iEBin = 0; iEBin < fNEBins; ++iEBin)
    {
      Double_t E = 1500;      
      TH1D * hProfile = (TH1D*)fHistArray->At(fWhichHisto->FindBin(E)-1); // -1 to convert from bins to array positions
      TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(fWhichHisto->FindBin(E)-1); // as arrays number from 0 but bins from 1

            
		  for( Int_t iR0Bin = 0; iR0Bin < fNR0Bins; ++iR0Bin)
		  {
		    std::cout << "iEBin = " << iEBin << "    iR0Bin = " << iR0Bin << std::endl;
        std::cout << "fR0Max = " << fR0Max << "   fR0Min = " << fR0Min << std::endl; 
			  Double_t R0 = fR0Min + iR0Bin * ((fR0Max - fR0Min)/(Double_t)fNR0Bins) ;
			  for( Int_t iCosTheta0Bin = 0; iCosTheta0Bin < fNCosTheta0Bins; ++iCosTheta0Bin)
			  {
				  Double_t cosTheta0 = fCosTheta0Min + iCosTheta0Bin*((fCosTheta0Max - fCosTheta0Min)/fNCosTheta0Bins);
				  rhoGS[0] = 0.0;
				  rhoGS[1] = 0.0;
				  rhoGS[2] = 0.0;
				  for( Int_t iSBin = 0; iSBin < fNSBins; ++iSBin)
				  {

            // Work out the co-ordinates needed to look up the emission profiles at this step
					  Double_t s = hProfile->GetBinCenter(iSBin+1);
            Double_t sWidth = hProfile->GetBinWidth(iSBin+1); // +1 is array numbering vs. bin numbering
            Double_t cosTheta = (R0*cosTheta0 - s)/TMath::Sqrt(R0*R0 + s*s - 2*R0*s*cosTheta0);

            if( iEBin == 0 && iR0Bin == 14 && iSBin == 250 && iCosTheta0Bin == 46 )
            {
              std::cout << "It's that bin... " << "   E = " << iEBin << "   sMax = " << s << "    cosTheta0 = " << cosTheta0 << "    R0=  " << R0 << std::endl;
            }
            // Now get the values of the emission profiles
            Double_t rho = hProfile->GetBinContent(iSBin+1);
            Int_t binCosTheta = hAngularProfile->GetXaxis()->FindBin(cosTheta);
            Int_t gBin = hAngularProfile->GetBin(binCosTheta,iSBin+1,0); //see doc of TH1::GetBin
            Double_t g = hAngularProfile->GetBinContent(gBin);

            // Do this step of the integral
            rhoGS[0] += rho * g * sWidth;
            rhoGS[1] += rho * g * s * sWidth;
            rhoGS[2] += rho * g * s * s * sWidth;
            outFile.write(reinterpret_cast<char*>(&rhoGS),sizeof(rhoGS));
            checkWriting << rhoGS[0] << " " << rhoGS[1] << " " << rhoGS[2] << " ";
            if(((iSBin + fNSBins*(iCosTheta0Bin + fNCosTheta0Bins*(iR0Bin + fNR0Bins*(iEBin)))) % 20) == 0) checkWriting << std::endl;
				    theCount += 3;
          }	
			  }
		  }
    }
    checkWriting.close();
   outFile.close();
   std::cout << "Total writes = " << theCount << std::endl;
   return;
}


void WCSimLikelihoodTuner::TabulateIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename)
{
  if(filename == "")
  {
    std::cerr << "WCSimLikelihoodTuner::TabulateIntegrals - Error: A filename is needed" << std::endl;
    return;
  }


  TabulateIndirectIntegrals(myType, filename);
  TabulateDirectIntegrals(myType, filename);
  return; 
}



>>>>>>> 1d76eafaf5c2e96d8b78405738da852825d0e829
