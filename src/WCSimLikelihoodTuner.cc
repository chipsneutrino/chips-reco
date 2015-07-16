#include <cassert>
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
#include "THnSparse.h"
#include "TArrayD.h"
#include "TH2D.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"

#include "WCSimChargePredictor.hh"
#include "WCSimIntegralLookupReader.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimRootGeom.hh"
#include "WCSimGeometry.hh"
#include "WCSimAnalysisConfig.hh"
#include "WCSimPMTManager.hh"

/**
* Constructor - if the Tuner is created without
* specifying the size of the detector then it
* assumes the PMTs are floating in space and will
* only prevent a PMT from seeing light if the track
* is behind it, not if it's outside the tank on the
* opposite side.
*/
WCSimLikelihoodTuner::WCSimLikelihoodTuner()
{
  	fConstrainExtent = false;
  	fExtent[0] = -1.0;
  	fExtent[1] = -1.0;
  	fExtent[2] = -1.0;
    fGeomType  = WCSimLikelihoodDigitArray::kUnknown;
    fLastCutoff = 0x0;
    std::cout << "Make PMT manager" << std::endl;
    fPMTManager = new WCSimPMTManager();
    std::cout << "Made PMT manager" << std::endl;
  
    fAverageQE  = 1.0;
    fIntegralParticleType = TrackType::Unknown;
    fCutoffIntegral = 0.0;
  	fCalculateIntegrals = WCSimAnalysisConfig::Instance()->GetCalculateIntegrals();
  	this->Initialize();

}

/**
 * Otherwise if the extent of the detector in x, y
 * and z is given, the tuner will kill tracks as soon
 * as they go outside of this region
 * @param myDigitArray Array of PMT responses for the event to process
 */
WCSimLikelihoodTuner::WCSimLikelihoodTuner(WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfiles * myEmissionProfiles)
{
    fLastCutoff = 0x0;
	  fCalculateIntegrals = WCSimAnalysisConfig::Instance()->GetCalculateIntegrals();
    this->UpdateDigitArray(myDigitArray);
  	this->Initialize();
  	fPMTManager = new WCSimPMTManager();
    fEmissionProfiles = myEmissionProfiles;
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
  std::cout << "Re-initialising" << std::endl;
  fAverageQE       = 1.0;

  // Pointer to the last track for which we calculated the cutoff, to prevent repetition

  // The binning scheme for the direct integral tables
  fIntegralParticleType = TrackType::Unknown;

  return;
}


/////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////
WCSimLikelihoodTuner::~WCSimLikelihoodTuner()
{
    //std::cout << " *** WCSimLikelihoodTuner::~WCSimLikelihoodTuner() *** Cleaning up" << std::endl;
	delete fPMTManager;
}



void WCSimLikelihoodTuner::UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray)
{

  fConstrainExtent = WCSimAnalysisConfig::Instance()->GetConstrainExtent();; 
  fExtent[0] = myDigitArray->GetExtent(0);
  fExtent[1] = myDigitArray->GetExtent(1);
  fExtent[2] = myDigitArray->GetExtent(2);
  fGeomType  = myDigitArray->GetGeomType();
  // std::cout << "fConstrainExtent = " << fConstrainExtent << " and extent = (" 
            //<< fExtent[0] << "," << fExtent[1] << "," << fExtent[2] 
            //<< ")" << std::endl;
  return;
}



////////////////////////////////////////////////////////////////////////
// Load the appropriate emission profiles for this track's particle type
////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadEmissionProfiles( WCSimLikelihoodTrackBase * myTrack )
{
  fEmissionProfiles->SetTrack(myTrack);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Work out the probability of light surviving to the PMT as it travels through the water
/////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::TransmissionFunction(Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
	if( s== 0 ) { return 1.0; }

  Double_t trans = 1.0;

  if( WCSimAnalysisConfig::Instance()->GetUseTransmission() ) 
  {
    // First we need the distance from the photon emission to the PMT
    TVector3 pmtPos      = myDigit->GetPos();
    TVector3 emissionPos = myTrack->GetPropagatedPos(s);
    Double_t r           = (pmtPos - emissionPos).Mag();

    // We'll use a triple exponential to parameterise the transmission probability
//    Double_t nu[3]     = {-1.137e-5,-5.212e-4, -4.359e-3}; // nu = 1/Decay length in mm
//    Double_t f[3]      = {0.8827, 0.08162, 0.03515};
    
    
    // Assumes scattering length has not been edited down (abwff = 0.625) 
    Double_t nu[3] = {-1.01526e-05, -6.67164e-06, -1.41964e-04  };
    Double_t f[3] = {2.10578e-01, 7.40266e-01, 4.92259e-02 };



    trans = 0.0;
    for(int i = 0; i < 3; ++i){ trans+= f[i]*exp(1.0 * 10 * r * nu[i]);}  //Convert to cm -> factor 10
  }
  //std::cout << "trans = " << trans << std::endl;
  return trans;
}

///////////////////////////////////////////////////////////////////////////////////////
// The efficiency of the PMT as a function of the incident angle of the photon
// This uses the same formula as WCSim.  Quantum efficiency (QE) is accounted for here.
// Efficiency of digitization is handled elsewhere
///////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::Efficiency(Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
  Double_t glassReflect = 0.0;
  if(WCSimAnalysisConfig::Instance()->GetUseGlassCathodeReflection()) { glassReflect = 0.24; }

  Double_t efficiency = 1.0;
  if( WCSimAnalysisConfig::Instance()->GetUseAngularEfficiency() )
  { 

    // We need the angle of incidence at the PMT
    TVector3 pmtPos      = myDigit->GetPos();
    TVector3 emissionPos = myTrack->GetPropagatedPos(s);
    TVector3 pmtToEm     = emissionPos - pmtPos;
    TVector3 pmtFace     = myDigit->GetFace();

    Double_t cosTheta = pmtFace.Dot(pmtToEm) / (pmtFace.Mag() * pmtToEm.Mag()); 
	  // The MiniBooNE method:
    Double_t theta = TMath::ACos(cosTheta) * 180. / TMath::Pi();
    if( theta > 90.0 )
    {
      theta = 180.0 - theta;
    }
	  efficiency =  (1 + (-1.182e-4) * pow(theta, 2) + 4.959e-9 * pow(theta, 4) - 7.371e-14 * pow(theta, 6));
	}
  //std::cout << "Efficiency = " << efficiency * (1-glassReflect) << std::endl;
  return efficiency * (1-glassReflect);

    // Function for the PMT acceptance's dependence on the angle: 
    // WCSim defines arrays of efficiency at 10 degree intervals and linearly interpolates

	/*
    if( cosTheta < 0.0 || cosTheta > 1.0 )
    {
    //  std::cout << "Behind the PMT, cosTheta = " << cosTheta << std::endl;
    //  pmtPos.Print() ;
    //  emissionPos.Print();
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
Double_t WCSimLikelihoodTuner::QuantumEfficiency(WCSimLikelihoodTrackBase * myTrack)
{
	WCSimTrackParameterEnums::TrackType type = myTrack->GetType();
	if(type == TrackType::MuonLike) return 0.05;
	else if(type == TrackType::ElectronLike) return 1.0;
	else return 1.0;
}
*/

//////////////////////////////////////////////////////////////////////////////////////////
// Work out the effect of solid angle on the probability of photon incidenc; pure geometry
//////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::SolidAngle(Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
  // Now some physical parameters of the phototube
  WCSimRootPMT myPMT     = ((WCSimRootGeom*) (WCSimGeometry::Instance())->GetWCSimGeometry())->GetPMTFromTubeID(myDigit->GetTubeId());
  // 
  Double_t mm_to_cm = 0.1;
	Double_t WCSimPMTRadius = myPMT.GetRadius() * mm_to_cm;

	// // Double_t WCSimPMTExposeHeight = WCSimPMTRadius - 1.0;
	

  // We need the distance from emission to the PMT
  // These are from WCSim so they're in cm
  TVector3 pmtDomePos  = myDigit->GetPos() + WCSimPMTRadius * myDigit->GetFace();
  TVector3 emissionPos = myTrack->GetPropagatedPos(s);
  Double_t r = (pmtDomePos - emissionPos).Mag();

  // Purely geometry: we need the solid angle of a cone whose bottom is a circle of the same radius as the circle of PMT poking
  // through the blacksheet.  This is 2pi( 1 - cos(coneAngle) ) where coneAngle can be deduced from trig and the PMT geometry
	// So the fraction of total solid angle is this/4pi = 0.5(1 - cos(coneAngle))
	Double_t solidAngle = 2.0*TMath::Pi()*(1.0 - (r)/sqrt( (r*r + WCSimPMTRadius * WCSimPMTRadius)));
  // std::cout << "Solid angle = " << solidAngle << " = " << solidAngle / (0.04*(M_PI)) << std::endl;
	
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
    // A fixed percentage of scattered light for now
    return 0.001;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//  We define J as the product of the transmission, efficiency and solid angle fuctions (and
//  also the scattering table for indirect light).  However we need it to be zero once the
//  particle has exited the detector, so we check for that before calling the other functions
/////////////////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::CalculateJ( Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit ) 
{
  std::vector<Double_t> J;

  // Check make sure the particle is still inside the detector
	if( fConstrainExtent )
	{
	  TVector3 pos = myTrack->GetPropagatedPos(s);
	  if( IsOutsideDetector(pos))
	  {
        J.push_back(0.0);
        J.push_back(0.0);
        return J;
      }
    }

    // Work out the direct and indirect contributions to J
    // J[0] = J_dir, J[1] = J_ind
//    if( s == 0.0)
//    {
//    	std::cout << "Transmission = " << this->TransmissionFunction(s, myTrack, myDigit) << std::endl
//    			  		<< "Efficiency =   " << this->Efficiency(s, myTrack, myDigit) << std::endl
//    			  		<< "SolidAngle =   " << this->SolidAngle(s, myTrack, myDigit) << std::endl;
//    }

    J.push_back(   this->TransmissionFunction(s, myTrack, myDigit) 
                 * this->Efficiency(s, myTrack, myDigit)
                 * this->QuantumEfficiency(myTrack)
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
void WCSimLikelihoodTuner::CalculateCoefficients(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit )
{
    
    // std::cout << "*** WCSimLikelihoodTuner::CalculateCoefficients() *** Calculating the coefficients for the integrals from tuning info" << std::endl; 
    for(int i = 0; i < 3; ++i)
    { 
        fDirCoeffs[i] = 0.0;    // Coefficients for direct (Cherenkov) light
        fIndCoeffs[i] = 0.0;    // And for indirect light
     }
    
 
    // std::cout << "We want number " << whichBin << std::endl;
    // std::cout << "fHistArray= " << fHistArray << std::endl;
    // Calculate the 3 s values we care about: s=0, s where integral(rho(s') ds')_0^{s} = 0.75, and double that
    Double_t s[3];
    s[0] = 10.0;
    s[1] = GetTrackLengthForPercentile(myTrack, 0.75);
    s[2] = 2 * s[1];

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
    // std::cout << "jDir:        " << JDir[0]       << "," << JDir[1]       << "," << JDir[2]       << std::endl;
    // std::cout << "fDirCoeffs:  " << fDirCoeffs[0] << "," << fDirCoeffs[1] << "," << fDirCoeffs[2] << std::endl;


    fIndCoeffs[0] = JInd[0];
    fIndCoeffs[1] = (4.*JInd[1] - JInd[2] -3.*JInd[0]) / (2*s[1]);
    fIndCoeffs[2] = (JInd[0] + JInd[2] - 2.*JInd[1]) / (2 * s[1] * s[1]);

    //    if((myDigit->GetZ() > 900))
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
std::vector<Double_t> WCSimLikelihoodTuner::CalculateCoefficientsVector(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit )
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
void WCSimLikelihoodTuner::CalculateCutoff( WCSimLikelihoodTrackBase * myTrack )
{
  //  std::cout << "Calculating the cutoff" << std::endl;
  if(fLastCutoff != 0x0 && myTrack->IsSameTrack(fLastCutoff)) { return; }
  
  Double_t cutoff = fEmissionProfiles->GetStoppingDistance(myTrack);
  // std::cout << "From profile, cutoff = " << cutoff << std::endl;

  if( fConstrainExtent )
  {
    assert(fGeomType != WCSimLikelihoodDigitArray::kUnknown);

    if( fGeomType == WCSimLikelihoodDigitArray::kCylinder )
    {
      cutoff = CalculateCylinderCutoff( myTrack );
    }
    else if( fGeomType == WCSimLikelihoodDigitArray::kMailBox )
    {
      cutoff = CalculateMailBoxCutoff( myTrack );
    }
    else assert(false);

  }
  // std::cout << "Energy = " << myTrack->GetE() << "   Cutoff = " << cutoff << ".... returning" << std::endl;
  fCutoffIntegral = cutoff;
  fLastCutoff     = myTrack;
  return;
}

/// Calculate where the integral should be cut off if the detector is a cylinder
Double_t WCSimLikelihoodTuner::CalculateCylinderCutoff(WCSimLikelihoodTrackBase * myTrack)
{
    TVector3 vtx     = myTrack->GetVtx();
    TVector3 dir     = myTrack->GetDir();

    // std::cout << "Vertex is " << std::endl;
    vtx.Print();

    // std::cout << "Direction is " << std::endl;
    dir.Print();
    Double_t cutoff = fEmissionProfiles->GetStoppingDistance(myTrack);

  // std::cout << "cutoff = fSMax = " << fSMax << std::endl;
  Double_t sMax[3] = {0., 0.,0.};
 
  TVector3 flatDir( dir(0), dir(1), 0);
  TVector3 flatVtx( vtx(0), vtx(1), 0); 


  // Escapes in r
  Double_t xDotD = flatVtx.Dot(flatDir);
  Double_t part1 = -1 * xDotD;
  Double_t part2 = TMath::Sqrt( xDotD * xDotD  - flatDir.Mag2() * (flatVtx.Mag2() - fExtent[0]*fExtent[0])  );

  Double_t s1    = (part1 + part2) / flatDir.Mag2();
  Double_t s2    = (part1 - part2) / flatDir.Mag2();
  if( s1 > 0) {      sMax[0] = s1; }
  else if( s2 > 0) { sMax[0] = s2; }

  sMax[1] = sMax[0];
  //std::cout << "It's a cylinder, sR = " << sMax[0] <<std::endl;

  // The z coordinate
  if( dir(2) > 1e-6 )
  {
      sMax[2] = (  fExtent[2] - vtx(2) ) / ( dir(2) );
      //std::cout << "Direction " << 2 << " cutoff is " << sMax[2] << std::endl;
  }
  else if( dir(2) < -1e-6 )
  {
      sMax[2] = ( -fExtent[2] - vtx(2) ) / ( dir(2) );
      // std::cout << "Direction " << 2 << " (negative) cutoff is " << sMax[2] << std::endl;
  }
  else { sMax[2] = cutoff; }
  // std::cout << sMax[2] << std::endl;
  
  // The world's laziest sorting algorithm:
  if( sMax[0] > 0 && sMax[0] < cutoff) { cutoff = sMax[0]; }
  if( sMax[2] > 0 && sMax[2] < cutoff) { cutoff = sMax[2]; }

  return cutoff;
}

/// Calculate where the integral should be cut off if the detector is a box
Double_t WCSimLikelihoodTuner::CalculateMailBoxCutoff(WCSimLikelihoodTrackBase * myTrack)
{
    TVector3 vtx     = myTrack->GetVtx();
    TVector3 dir     = myTrack->GetDir();

    // std::cout << "Vertex is " << std::endl;
    vtx.Print();

    // std::cout << "Direction is " << std::endl;
    dir.Print();
    Double_t cutoff = fEmissionProfiles->GetStoppingDistance(myTrack);
  Double_t sMax[3] = {0., 0.,0.};
  
  // Escapes in r
  for(Int_t i = 0 ; i < 3; ++i )
  {
    if( dir(i) > 1e-6 )
    {
        sMax[i] = (  fExtent[i] - vtx(i) ) / ( dir(i) );
        // std::cout << "Direction " << i << " cutoff is " << sMax[i] << std::endl;
    }
    else if( dir(i) < -1e-6 )
    {
        sMax[i] = ( -fExtent[i] - vtx(i) ) / ( dir(i) );
        // std::cout << "Direction " << i << " (negative) cutoff is " << sMax[i] << std::endl;
    }
    else { sMax[i] = cutoff; }
    // std::cout << sMax[i] << std::endl;
  }
  // The world's laziest sorting algorithm:
  if( sMax[0] < cutoff) { cutoff = sMax[0]; }
  if( sMax[1] < cutoff) { cutoff = sMax[1]; }
  if( sMax[2] < cutoff) { cutoff = sMax[2]; }

  return cutoff;
}

///////////////////////////////////////////////////////////////////////////
// We decouple the direct part into a source factor (flux x profile) and 
// an acceptance factor (solid angle x transmittance x PMT acceptance)
// We then integrate over the source factor, but it's much quicker to 
// tabulate this in advance.  This function opens the file containing those
// tables
///////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::LoadTabulatedIntegrals( WCSimLikelihoodTrackBase * myTrack )
{
  WCSimIntegralLookupReader::Instance()->LoadIntegrals(myTrack);

  return;
}


///////////////////////////////////////////////////////////////////////////
//  Get the integrals along the track, using the config file to decide if
//  they should be looked-up or calculated numerically
///////////////////////////////////////////////////////////////////////////
// First the contribution from Cherenkov light
double WCSimLikelihoodTuner::GetChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit, Int_t sPower)
{
    // Get the integral with a single power of s under the integral sign
    if(fCalculateIntegrals == true)
    {
      std::cout << "Calculating the Cherenkov integrals" << std::endl;
      return this->CalculateChIntegrals(myTrack, myDigit, sPower);
    }
    else 
    {
      std::cout << "Looking up the Cherenkov integrals" << std::endl;
      return this->LookupChIntegrals(myTrack, myDigit, sPower);
    }
}

std::vector<Double_t> WCSimLikelihoodTuner::GetChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
    // If we want all 3 powers it's quicker to get them all at once as it only involves one loop
    if(fCalculateIntegrals == true)
    {
      // std::cout << "Calculating the Cherenkov integrals" << std::endl;
      return this->CalculateChIntegrals(myTrack, myDigit);
    }
    else
    {
      // std::cout << "Looking up the Cherenkov integrals" << std::endl;
      return this->LookupChIntegrals(myTrack, myDigit);
    }
}

// Now the same for indirect light
double WCSimLikelihoodTuner::GetIndIntegrals(WCSimLikelihoodTrackBase * myTrack, Int_t sPower)
{
    // Get the integral with a single power of s under the integral sign
    if(fCalculateIntegrals == true)
    {
      return this->CalculateIndIntegrals(myTrack, sPower);
    }
    else return this->LookupIndIntegrals(myTrack, sPower);
}

std::vector<Double_t> WCSimLikelihoodTuner::GetIndIntegrals(WCSimLikelihoodTrackBase * myTrack)
{
    // If we want all 3 powers it's quicker to get them all at once as it only involves one loop
    if(fCalculateIntegrals == true)
    {
      return this->CalculateIndIntegrals(myTrack);
    }
    else return this->LookupIndIntegrals(myTrack);
}


///////////////////////////////////////////////////////////////////////////
//	Look up the tabulated integrals - just return 1 power of s
///////////////////////////////////////////////////////////////////////////
double WCSimLikelihoodTuner::LookupChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit, Int_t sPower)
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
std::vector<Double_t> WCSimLikelihoodTuner::LookupChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
   // std::cout << "*** WCSimLikelihoodTuner::LookupChIntegrals() *** Looking up the tabulated integrals for direct Cherenkov light" << std::endl;
  TVector3 pmtPos = myDigit->GetPos();
  TVector3 vtxPos = myTrack->GetVtx();
  TVector3 vtxDir = myTrack->GetDir();

	TVector3 toPMT = pmtPos - vtxPos;
  Double_t R0 = toPMT.Mag();
	Double_t cosTheta0 = TMath::Cos(vtxDir.Angle( toPMT ));
 	
 	WCSimIntegralLookupReader * myReader = WCSimIntegralLookupReader::Instance();
 	double E = myTrack->GetE();
 	TrackType::Type type = myTrack->GetType();

	std::vector<Double_t> integralsVec;
	integralsVec.push_back(myReader->GetRhoGIntegral(type, E, fCutoffIntegral, R0, cosTheta0));
	integralsVec.push_back(myReader->GetRhoGSIntegral(type, E, fCutoffIntegral, R0, cosTheta0));
	integralsVec.push_back(myReader->GetRhoGSSIntegral(type, E, fCutoffIntegral, R0, cosTheta0));

 	// std::cout << "Done! ... ";
 	// std::cout << std::endl;
//  if(integralsVec.at(1) != 0){
//  std::cout << "tubeID = " << myDigit->GetTubeId() << "    R0 = " << R0 << "    cosTheta0 = " << cosTheta0 << "     E = " << myTrack->GetE() << "    sMax = " << fCutoffIntegral << "     s term = " << integralsVec.at(1) << std::endl;
//  }
 	return integralsVec;
}

///////////////////////////////////////////////////////////////////////////
//	Look up the tabulated integrals
///////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::LookupIndIntegrals(WCSimLikelihoodTrackBase * myTrack, Int_t sPower)
{
	if( !(sPower < 3 && sPower >= 0) )
    {
        std::cerr << "WCSimLikelihoodTuner::LookupChIntegrals() ERROR: power of s must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::vector<Double_t> integralsVec = this->LookupIndIntegrals(myTrack);
    if( (UInt_t)sPower > integralsVec.size() ) std::cerr << "There's a problem with integralsVec!" << std::endl;
	return integralsVec.at(sPower);
}

std::vector<Double_t> WCSimLikelihoodTuner::LookupIndIntegrals(WCSimLikelihoodTrackBase * myTrack)
{
    // std::cout << "*** WCSimLikelihoodTuner::LookupChIntegrals() *** Looking up the tabulated integrals for direct Cherenkov light" << std::endl;
 	WCSimIntegralLookupReader * myReader = WCSimIntegralLookupReader::Instance();
 	double E = myTrack->GetE();
 	TrackType::Type type = myTrack->GetType();

	std::vector<Double_t> integralsVec;
	integralsVec.push_back(myReader->GetRhoIntegral(type, E, fCutoffIntegral));
	integralsVec.push_back(myReader->GetRhoSIntegral(type, E, fCutoffIntegral));
	integralsVec.push_back(myReader->GetRhoSSIntegral(type, E, fCutoffIntegral));


    return integralsVec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Work out the integrals over the emission profiles numerically.  This is slow and only good for testing - 
// eventually they'll be tabulated and looked-up
///////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit, int i)
{
	std::vector<Double_t> integrals = this->CalculateChIntegrals(myTrack, myDigit);
	if( i < 0 || i > 2) 
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateChIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    return integrals.at(i);
}

std::vector<Double_t> WCSimLikelihoodTuner::CalculateChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
    // std::cout << "*** WCSimLikelihoodTuner::CalculateChIntegrals() ***" << std::endl;
    this->LoadEmissionProfiles(myTrack);

    // Check to see if we need to curtail the integrals:
    
    CalculateCutoff( myTrack );
    
    std::vector<Int_t> sPowers;
    sPowers.push_back(0);
    sPowers.push_back(1);
    sPowers.push_back(2);
    std::vector<Double_t> integrals = fEmissionProfiles->GetRhoGIntegrals(myTrack, myDigit, sPowers, fCutoffIntegral, kTRUE);

    // std::cout << "WCSimLikelihoodTuner::CalculateChIntegrals()   s = " << fCutoffIntegral << "   " << integrals[0] << "   " << integrals[1] << "   " << integrals[2] << std::endl;
  
    return integrals;
}

///////////////////////////////////////////////////////////////////////////////////////////
// As for the previous function, but now we're calculating the integrals for indirect light
///////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateIndIntegrals(WCSimLikelihoodTrackBase * myTrack, int i)
{
  
	std::vector<Double_t> integrals = this->CalculateIndIntegrals(myTrack);
	if( i < 0 || i > 2) 
    {
        std::cerr << "WCSimLikelihoodTuner::CalculateIndIntegrals() ERROR: i must be 0, 1 or 2" << std::endl;
        exit(EXIT_FAILURE);
    }
    return integrals.at(i);
}

std::vector<Double_t> WCSimLikelihoodTuner::CalculateIndIntegrals(WCSimLikelihoodTrackBase * myTrack)
{ 

    this->LoadEmissionProfiles(myTrack);
    std::vector<Int_t> sPowers;
    sPowers.push_back(0);
    sPowers.push_back(1);
    sPowers.push_back(2);

    std::vector<Double_t> integrals = fEmissionProfiles->GetRhoIntegrals(sPowers, myTrack->GetE(), 0, fCutoffIntegral, kTRUE);
    //	std::cout << "Is it 1? " << integrals[0] << std::endl; 
	  return integrals;

}

// DEBUGGING: override the configuration file and choose to calculate or lookup integrals
void WCSimLikelihoodTuner::SetCalculateIntegrals(Bool_t calc)
{
  fCalculateIntegrals = calc;
  return;
}

Bool_t WCSimLikelihoodTuner::GetCalculateIntegrals() const
{
  return fCalculateIntegrals;
}

Bool_t WCSimLikelihoodTuner::IsOutsideDetector(const TVector3 &pos)
{
    return (fabs(pos.X()) > fExtent[0] || fabs(pos.Y()) > fExtent[1] || fabs(pos.Z()) > fExtent[2]);
}

Double_t WCSimLikelihoodTuner::GetLightFlux(WCSimLikelihoodTrackBase * myTrack)
{
  return fEmissionProfiles->GetLightFlux(myTrack);
}

Double_t WCSimLikelihoodTuner::GetTrackLengthForPercentile(WCSimLikelihoodTrackBase * myTrack, const double &percentile)
{
  double length = -999.9;
  if(fCalculateIntegrals || !(WCSimAnalysisConfig::Instance()->GetTruncateIntegrals())) { 
    this->LoadEmissionProfiles(myTrack);
    return fEmissionProfiles->GetTrackLengthForPercentile(percentile); 
  }
  else {
    WCSimIntegralLookupReader * myLookupReader = WCSimIntegralLookupReader::Instance();
    return myLookupReader->GetTrackLengthForPercentile(myTrack->GetType(), myTrack->GetE(), percentile);
  }
  assert(length != -999.9);
  return length;
}

Double_t WCSimLikelihoodTuner::QuantumEfficiency(WCSimLikelihoodTrackBase* myTrack) {
	double qe=  0.0;
	if(myTrack->GetType() == TrackType::MuonLike || myTrack->GetType() == TrackType::ElectronLike)
	{
		qe = 0.436129;
	}
	else assert(myTrack->GetType() == TrackType::ElectronLike || myTrack->GetType() == TrackType::MuonLike);
	return qe;
}

double WCSimLikelihoodTuner::GetCutoff(WCSimLikelihoodTrackBase* myTrack) {
    double cutoff = 0.0;
	if( fGeomType == WCSimLikelihoodDigitArray::kCylinder )
    {
      cutoff = CalculateCylinderCutoff( myTrack );
    }
    else if( fGeomType == WCSimLikelihoodDigitArray::kMailBox )
    {
      cutoff = CalculateMailBoxCutoff( myTrack );
    }
	return cutoff;
}
