#include <math.h>
#include <sstream>

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
WCSimLikelihoodTuner::WCSimLikelihoodTuner()
{
  fConstrainExtent = kFALSE;
  fExtent[0] = -1.0;
  fExtent[1] = -1.0;
  fExtent[2] = -1.0;
  this->Initialize();

}

/////////////////////////////////////////////////////
// Otherwise if the extent of the detector in x, y 
// and z is given, the tuner will kill tracks as soon
// as they go outside of this region
/////////////////////////////////////////////////////
WCSimLikelihoodTuner::WCSimLikelihoodTuner(Double_t xMax, Double_t yMax, Double_t zMax)
{
  fConstrainExtent = kTRUE;
  fExtent[0] = xMax;
  fExtent[1] = yMax;
  fExtent[2] = zMax;
  this->Initialize();
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
  fWhichHisto->Draw();
  fAverageQE = 1.0;
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
    delete fHistArray;
    delete fAngHistArray;
    delete fFluxArray;
    delete fWhichHisto;
    fProfileLocation = 0;
    fProfiles = 0;
    fHistArray = 0;
    fAngHistArray= 0;
    fFluxArray= 0;
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
    Double_t theta = TMath::ACos(cosTheta);
    if( theta > 90.0 )
    {
      theta = 180.0 - theta;
    }
	
	return (1 + (-1.182e-4) * pow(theta, 2) + 4.959e-9 * pow(theta, 4) - 7.371e-14 * pow(theta, 6));


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

    if((myDigit->GetTubeId() % 10000) == 0)
    {
      TGraph * gr = new TGraph();
      for(int i = 0; i < hProfile->GetNbinsX(); ++i)
      {
        double tempS = hProfile->GetXaxis()->GetBinCenter(i);
        std::vector<Double_t> tempJ = this->CalculateJ(tempS, myTrack, myDigit);
        gr->SetPoint(gr->GetN(), tempS, tempJ.at(1));
      }
      TF1 * jInd = new TF1("jInd","[0] + [1] * x + [2] * x * x", 0, hProfile->GetXaxis()->GetXmax());
      jInd->SetParameters(fIndCoeffs[0], fIndCoeffs[1], fIndCoeffs[2]);
  
      TCanvas * c1 = new TCanvas("c1","c1", 800, 600);
      gr->SetLineColor(kBlue);
      gr->SetLineWidth(2);
      gr->Draw("ALP");
      jInd->SetLineColor(kRed);
      jInd->SetLineWidth(2);
      jInd->Draw("same");
      TString str;
      str.Form("jInd%d.png",myDigit->GetTubeId());
      c1->SaveAs(str.Data());
      delete c1;
      delete jInd;
      delete gr;
    }  
  
  
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
  Double_t cutoff = -1.0;
  if( fConstrainExtent )
  {
    TVector3 vtx = myTrack->GetVtx();
    TVector3 dir = myTrack->GetDir();
    Double_t sMax[3];
    for(int i = 0; i < 3; ++i)
    {
      if(dir[i])
      {
        sMax[i] = -1. * (vtx(i) / dir(i)) + TMath::Sqrt( fExtent[i] * fExtent[i] / (dir(i)*dir(i)));
      }
      else sMax[i] = 0.0;
    }
  
    // The world's laziest sorting algorithm:
    cutoff = sMax[0];
    if( sMax[1] > cutoff ) cutoff = sMax[1];
    if( sMax[2] > cutoff ) cutoff = sMax[2];
  }
  fCutoffIntegral = cutoff;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Work out the integrals over the emission profiles numerically.  This is slow and only good for testing - 
// eventually they'll be tabulated and looked-up
///////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateChIntegrals(int i, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{


//    std::cout << "*** WCSimLikelihoodTuner::CalculateChIntegrals() ***" << std::endl;
    if( !(i < 3 && i >= 0) )
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateChIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }

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
    Double_t integrals[3]={0.,0.,0.};
    CalculateCutoff( myTrack );
    int cutoffBin = hProfile->GetNbinsX();
    if( fCutoffIntegral > 0.0 && fCutoffIntegral < hProfile->GetXaxis()->GetXmax() && fCutoffIntegral > hProfile->GetXaxis()->GetXmin() )
    {
      cutoffBin = hProfile->FindBin(fCutoffIntegral);
    }

    // Now do the integrals
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

        // Need to calculate theta
        TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
        TVector3 vtxPos(myTrack->GetVtx());
        TVector3 vtxDir(myTrack->GetDir());
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
        //std::cout << "s = " << s << std::endl;
    }

    Double_t ret = 0;
    if(i==0){ ret = integrals[0];}
    if(i==1){ ret = integrals[1];}
    if(i==2){ret = integrals[2];}
    
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////
// As for the previous function, but now we're calculating the integrals for indirect light
///////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateIndIntegrals(int i, WCSimLikelihoodTrack * myTrack)
{
    if( !(i < 3 && i >= 0) )
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateIndIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    Double_t integrals[3]={0.0,0.0,0.0};

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
        integrals[1] += rho * s * hProfile->GetBinWidth(iBin);
        integrals[2] += rho * s * s * hProfile->GetBinWidth(iBin);
    }
    Double_t ret = 0;
    if(i==0){ ret = 1;}
    if(i==1){ ret = integrals[1];}
    if(i==2){ ret = integrals[2];}
//    std::cout << " Indirect: i = " << i << " and ret = " << ret << std::endl;

    return ret;

}

///////////////////////////////////////////////////////////////////////////////////////
// Tabulate the integrals over the emission profiles.  Indirect light does not need the 
// angular emission profile so this one is only indexed by energy.  This should only
// need running once - from then on we can just look the integrals up in the table.
///////////////////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::TabulateIndirectIntegrals( WCSimLikelihoodTrack::TrackType myType, TString filename)
{
	  // Make a file for saving the integrals.  Require a name to avoid overwriting old ones inadvertantly
    /*if(filename == "")
    {
      std::cerr << "WCSimLikelihoodTuner::TabulateDirectIntegrals - Error: A filename is needed" << std::endl;
      return;
    }
    WCSimLikelihoodTrack * dummyTrack = new WCSimLikelihoodTrack();
    filename += dummyTrack->TrackTypeToString( myType );
    filename += TString(".root");
    std::cout << filename.Data() << std::endl;
    delete dummyTrack;*/
    TFile * f = new TFile("integralsRho.root","RECREATE");
//    TFile * f = new TFile(filename.Data());
   
////////////////////////////////////////////////////////////////////////////////////////////////

for( int iEBin = 0; iEBin < fWhichHisto->GetNbinsX(); ++iEBin)
    {
    	std::cout << "Filling tree for energy bin " << iEBin << std::endl;
        TH1D * hProfile = (TH1D*)fHistArray->At(iEBin);
		std::stringstream ss;
		ss << "energyBin" << iEBin;
		std::string treeName = ss.str();	
		f->cd();
		TTree * t = new TTree(treeName.c_str(), treeName.c_str());
		TArrayD integrals(3);
		Double_t rhoS[3] = {0.0, 0.0, 0.0};
		t->Branch("integrals",&(rhoS[0]), "integrals[3]/D");
		
		
		Int_t sBins = hProfile->GetNbinsX();


		for( Int_t iSBin = 0; iSBin < sBins; ++iSBin)
		{
			Double_t s = hProfile->GetBinCenter(iSBin+1);
			Double_t sWidth = hProfile->GetBinWidth(iSBin+1);
            Double_t rho = hProfile->GetBinContent(iSBin+1);
            rhoS[0] += rho * sWidth;
            rhoS[1] += rho * s * sWidth;
            rhoS[2] += rho * s * s * sWidth;
					
			std::cout << "s = " << s << " Integral[0] = " << rhoS[0] << std::endl;
			t->Fill();
		}	
	

		t->Write();
		gDirectory->Purge();
		delete t;
    }
////////////////////////////////////////////////////////////////////////////////////////////////   

    fWhichHisto->SetName("hWhichEnergyBin");
    fWhichHisto->Write();
	f->Purge();
    f->Close();


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
    // Make a file for saving the integrals.  Require a name to avoid overwriting old ones inadvertantly
/*    if(filename == "")
    {
      std::cerr << "WCSimLikelihoodTuner::TabulateDirectIntegrals - Error: A filename is needed" << std::endl;
      return;
    }
    WCSimLikelihoodTrack * dummyTrack = new WCSimLikelihoodTrack();
    filename += dummyTrack->TrackTypeToString( myType );
    filename += TString(".root");
    std::cout << filename.Data() << std::endl;
    delete dummyTrack;
    */
    TFile * f = new TFile("integralsRhoG.root","RECREATE");

    WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack();
    WCSimLikelihoodDigitArray * myLDArray = new WCSimLikelihoodDigitArray();
    WCSimChargeLikelihood * myCL = new WCSimChargeLikelihood( myLDArray, myTrack ) ;

    Double_t theR0Interval = myCL->GetR0Interval();
    Int_t nR0Bins = (Int_t)(ceil((40000-0)/theR0Interval));
    TH1F * hR0Bins = new TH1F("hR0Bins","hR0Bins", nR0Bins, 0.,40000.);

    Double_t theCosTheta0Interval = myCL->GetCosTheta0Interval();
    Int_t nCosTheta0Bins = (Int_t)(ceil(2./theCosTheta0Interval)); 
    TH1F * hCosTheta0Bins = new TH1F("hCosTheta0Bins","hCosTheta0Bins", nCosTheta0Bins, -1.,1.);

    this->LoadEmissionProfiles(myType);    
    for( int iEBin = 0; iEBin < fWhichHisto->GetNbinsX(); ++iEBin)
    {
    	std::cout << "Filling tree for energy bin " << iEBin << std::endl;
        TH1D * hProfile = (TH1D*)fHistArray->At(iEBin);
        TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(iEBin);
		std::stringstream ss;
		ss << "energyBin" << iEBin;
		std::string treeName = ss.str();	
		f->cd();
		TTree * t = new TTree(treeName.c_str(), treeName.c_str());
		Double_t rhoGS[3] = {0.,0.,0.};
		t->Branch("integrals",&rhoGS,"integrals[3]/D");
		
		
		Int_t sBins = hProfile->GetNbinsX();

		for( Int_t iR0Bin = 0; iR0Bin < nR0Bins; ++iR0Bin)
		{
			std::cout << "R0Bin = " << iR0Bin << std::endl;
			Double_t R0 = hR0Bins->GetBinCenter( iR0Bin + 1) ;
			for( Int_t iCosTheta0Bin = 0; iCosTheta0Bin < nCosTheta0Bins; ++iCosTheta0Bin)
			{
				Double_t cosTheta0 = hCosTheta0Bins->GetBinCenter( iCosTheta0Bin + 1);
				TVector3 dir(-1.0 * cosTheta0, TMath::Sin(TMath::ACos(cosTheta0)),0);
				rhoGS[0] = 0.0;
				rhoGS[1] = 0.0;
				rhoGS[2] = 0.0;
				for( Int_t iSBin = 0; iSBin < sBins; ++iSBin)
				{
					Double_t s = hProfile->GetBinCenter(iSBin+1);
					Double_t sWidth = hProfile->GetBinWidth(iSBin+1);
            		TVector3 toPMT(R0 - s*cosTheta0, R0 + s*TMath::Sin(TMath::ACos(cosTheta0)),0);
            		Double_t cosTheta = TMath::Cos(dir.Angle(toPMT));
            		Double_t rho = hProfile->GetBinContent(iSBin+1);
            		Int_t binCosTheta = hAngularProfile->GetXaxis()->FindBin(cosTheta);
            		Int_t gBin = hAngularProfile->GetBin(binCosTheta,iSBin+1,0); //see doc of TH1::GetBin
            		Double_t g = hAngularProfile->GetBinContent(gBin);
            		rhoGS[0] += rho * g * sWidth;
            		rhoGS[1] += rho * g * s * sWidth;
            		rhoGS[2] += rho * g * s * s * sWidth;
					
					t->Fill();
				}	
			}
		}
			
		t->Write();
		delete t;
    }

////////////////////////////////////////////////////////////////////////////////////////////

 /*   this->LoadEmissionProfiles(myType);
        
    for( int iEBin = 0; iEBin < fWhichHisto->GetNbinsX(); ++iEBin) // loop over energy bins
    {
    	std::cout << "Tabulating energy bin " << iEBin << std::endl;
    	// Get the profiles
    	//std::cout << "We want number " << whichBin << std::endl;
        TH1D * hProfile = (TH1D*)fHistArray->At(iEBin);
        TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(iEBin);
    	
    	R0Array->Clear();

		std::cout << "Creating the arrays" << std::endl;
    	for(int i = 0; i < nR0Bins; ++i)
    	{
    		if( i % 100 == 0) std::cout << "R0 bin = " << i << std::endl;
     		TObjArray * CosTheta0Array = new TObjArray(nCosTheta0Bins);
      		for(int j = 0; j < nCosTheta0Bins;  ++j)
      		{
      	 		TObjArray * sArray = new TObjArray(hProfile->GetNbinsX());
         		for(int k = 0; k < hProfile->GetNbinsX(); ++k)
         		{
	         		TArrayD * integrals = new TArrayD(3);
    		 		sArray->AddAt( (TObject*)integrals, k );
         		}
        		CosTheta0Array->AddAt( sArray, j );
      		}
      		R0Array->AddAt( CosTheta0Array, i );
    	}
    	
      for(int iR0Bin = 0; iR0Bin < nR0Bins; ++ iR0Bin)
      {
        Double_t R0 = hR0Bins->GetBinCenter( iR0Bin + 1) ;
        std::cout << iR0Bin << std::endl;

        for(int iCosTheta0Bin = 0; iCosTheta0Bin < nCosTheta0Bins; ++iCosTheta0Bin)
        {
          Double_t gRhoIntegral[3];
          gRhoIntegral[0] = 0.0;
          gRhoIntegral[1] = 0.0;
          gRhoIntegral[2] = 0.0;
     

          Double_t cosTheta0 = hCosTheta0Bins->GetBinCenter( iCosTheta0Bin + 1); // bin 0 is the underflow bin in the histogram
                                                                                 // but the first bin in the array...
          TVector3 dir(-1.0 * cosTheta0, TMath::Sin(TMath::ACos(cosTheta0)),0);
 
          for(int iBin = 1; iBin <= hProfile->GetNbinsX(); ++iBin)
          {
              Double_t s = hProfile->GetBinCenter(iBin);
              TVector3 toPMT(R0 - s*cosTheta0, R0 + s*TMath::Sin(TMath::ACos(cosTheta0)),0);
              Double_t cosTheta = TMath::Cos(dir.Angle(toPMT));
              Double_t rho = hProfile->GetBinContent(iBin);
              Int_t binCosTheta = hAngularProfile->GetXaxis()->FindBin(cosTheta);
              Int_t gBin = hAngularProfile->GetBin(binCosTheta,iBin,0); //see doc of TH1::GetBin
              Double_t g = hAngularProfile->GetBinContent(gBin);
              gRhoIntegral[0] += rho * g;
              gRhoIntegral[1] += rho * g * s;
              gRhoIntegral[2] += rho * g * s * s;
          
              for(int iIntegral = 0; iIntegral < 3; ++iIntegral)
          	  {
              	((TArrayD*)((TObjArray*)((TObjArray*)R0Array->At(iR0Bin))->At(iCosTheta0Bin))->At(iBin-1))->AddAt(gRhoIntegral[iIntegral], iIntegral); 
            	// This line is horrible. We're getting the TArrayD for a given R0 and CosTheta0 bin and setting each of its 3 entries
           	 	// to the integrals we just calculate.  But there's lots of unpleasant casting because we have arrays of arrays
          	  }
          
          }
          

        }
      }
      std::cout << "Filling the tree" << std::endl;
      t->Fill();
      std::cout << "Filled the tree" << std::endl;
    }
	delete R0Array;  
*/
    fWhichHisto->SetName("hWhichEnergyBin");
//    t->Write();
    fWhichHisto->Write();
    hR0Bins->Write();
    hCosTheta0Bins->Write();
    f->Purge();
    f->Close();


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



