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
#include "TMatrixD.h"
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
#include "WCSimTransmissionFunctionLookup.hh"
#include "WCSimScatteringTableManager.hh"

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
    fPMTManager = new WCSimPMTManager();
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
WCSimLikelihoodTuner::WCSimLikelihoodTuner(WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager)
{
    fLastCutoff = 0x0;
	fCalculateIntegrals = WCSimAnalysisConfig::Instance()->GetCalculateIntegrals();
    this->UpdateDigitArray(myDigitArray);
  	this->Initialize();
  	fPMTManager = new WCSimPMTManager();
    fEmissionProfileManager = myEmissionProfileManager;
}

///////////////////////////////////////////////////////////
// These values should all be the same regardless of
// the constructor called.  Later we should teach it to
// distinguish between track types and open the appropriate
// emission profile file
///////////////////////////////////////////////////////////
void WCSimLikelihoodTuner::Initialize()
{
  fAverageQE       = 1.0;
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
  if(fGeomType == WCSimLikelihoodDigitArray::kCylinder)
  {
	  fMaxDistance = 2.0 * sqrt( fExtent[0] * fExtent[0] + fExtent[2] * fExtent[2]);
  }
  else if(fGeomType == WCSimLikelihoodDigitArray::kMailBox)
  {
	  fMaxDistance = 2 * sqrt( fExtent[0] * fExtent[0] + fExtent[1]*fExtent[1] + fExtent[2] * fExtent[2]);
  }
  // std::cout << "fConstrainExtent = " << fConstrainExtent << " and extent = (" 
            //<< fExtent[0] << "," << fExtent[1] << "," << fExtent[2] 
            //<< ")" << std::endl;
  return;
}



/////////////////////////////////////////////////////////////////////////////////////////
// Work out the probability of light surviving to the PMT as it travels through the water
/////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::TransmissionFunction()
{
	if( fCache.GetS() == 0 ) { return 1.0; }

  Double_t trans = (double)(!WCSimAnalysisConfig::Instance()->GetUseTransmission());
  // Default to 1 if we're not using transmission

  if( WCSimAnalysisConfig::Instance()->GetUseTransmission() ) 
  {

    // First we need the distance from the photon emission to the PMT
    Double_t r           = fCache.GetDistanceToPMT();

    // We'll use a triple exponential to parameterise the transmission probability

    //    Double_t nu[3]     = {-1.137e-5,-5.212e-4, -4.359e-3}; // nu = 1/Decay length in mm
    //    Double_t f[3]      = {0.8827, 0.08162, 0.03515};
        
    // Assumes scattering length has not been edited down (abwff = 0.625) 
    // These are in inverse cm
    Double_t nu[3] = {-1.01526e-04, -6.67164e-05, -1.41964e-03  };
    Double_t f[3] = {2.10578e-01, 7.40266e-01, 4.92259e-02 };

    // Scale these down temporarily because I've scaled down the transmission length to 30m max
    // Should eventually re-fit this properly
    for(int i = 0; i < 3; ++i)
    {
      nu[i] = 1.0 / (1.0/nu[i] * 0.0737);
    }
  
    for(int i = 0; i < 3; ++i)
    { 
    	trans += WCSimTransmissionFunctionLookup::Instance()->GetExp(0, fMaxDistance, f[i], nu[i], r);
    }
  }
  return trans;
}



///////////////////////////////////////////////////////////////////////////////////////
// The efficiency of the PMT as a function of the incident angle of the photon
// This uses the same formula as WCSim.  Quantum efficiency (QE) is accounted for here.
// Efficiency of digitization is handled elsewhere
///////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::Efficiency()
{
  // Whether the reflectivity of the PMT glass has been included in the WCSim simulation
  // It's usually switched on there by default
  Double_t glassReflect = 0.0;
  if(WCSimAnalysisConfig::Instance()->GetUseGlassCathodeReflection()) { glassReflect = 0.24; }

  Double_t efficiency = 1.0;
  if( WCSimAnalysisConfig::Instance()->GetUseAngularEfficiency() )
  { 

	  // We need the angle of incidence at the PMT
	  TVector3 pmtToEm     = -fCache.GetToPMT();	// Vector from PMT to photon emission point
	  TVector3 pmtFace     = fCache.GetDigit()->GetFace();	// Direction PMT is facing (unit vector)
	  Double_t cosTheta    = pmtFace.Dot(pmtToEm) / fCache.GetDistanceToPMT();  // Quicker dot product because |pmtFace| = 1
	  Double_t theta 	   = TMath::ACos(cosTheta) * TMath::RadToDeg();
	  if( theta > 90.0 )
	  {
		  theta = 180.0 - theta;
	  }

	  // The PMT angular efficiency taken directly from WCSim
	  const int num_elements = 10;
	  Double_t collection_angle[num_elements]={0.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
	  Double_t collection_eff[num_elements]={1.,1.,0.99,0.95,0.90,0.85,0.80,0.69,0.35,0.13};

    // Interpolate between the efficiencies for the closes two angles
	  for(int iEntry = 0; iEntry < num_elements-1; ++iEntry)
	  {
		  if( theta >= collection_angle[iEntry] && theta < collection_angle[iEntry+1])
		  {
			  efficiency =   collection_eff[iEntry]
			               + (theta - collection_angle[iEntry]) / (collection_angle[iEntry+1] - collection_angle[iEntry])
			               * (collection_eff[iEntry+1] - collection_eff[iEntry]);
		  }
	  }
  }
  return (efficiency) * (1.0 - glassReflect);
}

//////////////////////////////////////////////////////////////////////////////////////////
// Work out the effect of solid angle on the probability of photon incidenc; pure geometry
//////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::SolidAngleFraction()
{
  // Some physical parameters of the phototube
  WCSimRootPMT myPMT     = ((WCSimRootGeom*) (WCSimGeometry::Instance())->GetWCSimGeometry())->GetPMTFromTubeID(fCache.GetDigit()->GetTubeId());
  const Double_t mm_to_cm = 0.1;
  const Double_t WCSimPMTRadius = myPMT.GetRadius() * mm_to_cm;  // Radius of the PMT sphere
  const Double_t exposeHeight = fCache.GetDigit()->GetExposeHeight() * mm_to_cm; // Height of the PMT sphere poking through the blacksheet


  // We need the distance from emission to the PMT
  // These are from WCSim so they're in cm
  TVector3 pmtDomePos  = fCache.GetDigit()->GetPos() + WCSimPMTRadius * fCache.GetDigit()->GetFace();
  TVector3 emissionPos = fCache.GetEmissionPos();
  const Double_t r = (pmtDomePos - emissionPos).Mag();
  const Double_t rPlusExpHeight = r + exposeHeight;
  // Purely geometry: we need the solid angle of a cone whose bottom is a circle of the same radius as the circle of PMT poking
  // through the blacksheet.  This is 2pi( 1 - cos(coneAngle) ) where coneAngle can be deduced from trig and the PMT geometry

  // So the fraction of total solid angle is this/4pi = 0.5(1 - cos(coneAngle))
  const Double_t solidAngleFrac = 0.5 * (1.0 - ((rPlusExpHeight)/sqrt( (rPlusExpHeight*rPlusExpHeight) + WCSimPMTRadius*WCSimPMTRadius)))  ;

  // std::cout << "Solid angle = " << solidAngle << "   c.f. " << 2.0 * TMath::Pi()*(1.0 - (r/sqrt(r*r + WCSimPMTRadius * WCSimPMTRadius))) << std::endl;
  return solidAngleFrac;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// The scattering table assumes that scattered light can be modelled by scaling the
// number of photons you would see from an isotropic source of Cherenkov light by
// a factor dependent on geometry, energy etc.
// It's small and computationally intensive to tabulate, so we're taking a constant initially
/////////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::ScatteringTable()
{
  // If we are not using the scattering table, just use a value of 1%
  if(WCSimAnalysisConfig::Instance()->GetUseScatteringTable() == false){
    return 0.01;
  }

  // Calculate the values of the six required coordinates.
  std::vector<float> coordinates(6, 0.0);

  // Vertex R and Z components first.
  TVector3 vtx = fCache.GetEmissionPos();
  TVector3 vtx2D(vtx.X(),vtx.Y(),0.);
  coordinates[0] = sqrt(vtx.X()*vtx.X() + vtx.Y()*vtx.Y()); // Vertex R position first
  coordinates[1] = vtx.Z(); // Vertex Z position next.

  // Get the PMT and its position (convert from cm -> mm).
  WCSimRootPMT myPMT = ((WCSimRootGeom*) (WCSimGeometry::Instance())->GetWCSimGeometry())->GetPMTFromTubeID(fCache.GetDigit()->GetTubeId());
  TVector3 pmtPos(myPMT.GetPosition(0)*10., myPMT.GetPosition(1)*10., myPMT.GetPosition(2)*10.);
  TVector3 pmtPos2D(pmtPos.X(),pmtPos.Y(),0.);
  coordinates[2] = vtx2D.Angle(pmtPos2D); // Angle zeta between the 2D vtx position and 2D pmt position

  // The PMT position (this is R for caps, and Z for the barrel).
  float pmtCoord = 0.;
  if(myPMT.GetCylLoc() != 1){
    pmtCoord = sqrt(pmtPos.X()*pmtPos.X() + pmtPos.Y()*pmtPos.Y());
  }
  else{
    pmtCoord = pmtPos.Z();
  }
  coordinates[3] = pmtCoord;   

  // Now for the theta and phi directions
  coordinates[4] = (fCache.GetTrack()->GetTheta()); // Theta direction
  coordinates[5] = (fCache.GetTrack()->GetPhi()); // Phi direction

  WCSimScatteringTableManager* scatteringManager = WCSimScatteringTableManager::Instance();

  Double_t scatteringValue = scatteringManager->GetScatteringValue(coordinates,myPMT.GetCylLoc());
//  std::cout << scatteringValue << std::endl;

  if(scatteringValue < 1e-4) scatteringValue = 1e-4;

  return scatteringValue;

}

/////////////////////////////////////////////////////////////////////////////////////////////
//  We define J as the product of the transmission, efficiency and solid angle fuctions (and
//  also the scattering table for indirect light).  However we need it to be zero once the
//  particle has exited the detector, so we check for that before calling the other functions
/////////////////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::CalculateJ( Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit ) 
{
	std::vector<Double_t> J(2, 0.0);

	// Load in this combination of s, digit and track so we don't keep recalculating distances and angles
	MakeCache(s, myTrack, myDigit);


	// Check make sure the particle is still inside the detector
	if( fConstrainExtent )
	{
		TVector3 pos = fCache.GetEmissionPos();;
		if( IsOutsideDetector(pos))
		{
			return J; // It already contains zeros
		}
	}

  // Work out the direct and indirect contributions to J
  // J[0] = J_dir, J[1] = J_ind
  // if( myDigit->GetQ() > 10 && myDigit->GetTubeId() > 4000 && myDigit->GetTubeId() < 4100)
  // {
  // 	std::cout << "Transmission = " << this->TransmissionFunction(s, myTrack, myDigit) << std::endl
  // 			  		<< "Efficiency = " << this->Efficiency(s, myTrack, myDigit) << std::endl
  // 			  		<< "SolidAngle = " << this->SolidAngleFraction(s, myTrack, myDigit) << std::endl
  // 			  		<< "QE           = " << this->QuantumEfficiency(s, myTrack, myDigit) << std::endl;
  // }

  J[0] = (this->TransmissionFunction() 
          * this->Efficiency()
          * this->QuantumEfficiency()
          * this->SolidAngleFraction());
  J[1] = (J.at(0) * this->ScatteringTable());
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
    // Calculate the 3 s values we care about: s=0, s where integral(rho(s') ds')_0^{s} = 0.75, and double that
    Double_t s[3];
    s[0] = 10.0;
    s[1] = GetTrackLengthForPercentile(myTrack, 0.75);
    s[2] = 2 * s[1];

    // Check these s values are inside the detector and adjust them if not:
    this->CalculateCutoff(myTrack);
    if( fCutoffIntegral > 0.0 && s[2] > fCutoffIntegral )
    {
      s[2] = 0.8*fCutoffIntegral;
      s[1] = 0.5 * s[2];
    }

    // Evaluate J at each point
    Double_t JDir[3] = {0.,0.,0.};
    Double_t JInd[3] = {0.,0.,0.};
    for(int k = 0; k < 3; ++k)
    {
        std::vector<Double_t> J = this->CalculateJ( s[k], myTrack, myDigit);
        JDir[k] = J.at(0); 
        JInd[k] = J.at(1);
    }

    /* Now construct a quadratic that fits the Cherenkov and indirect parts, sampling them at these 3 points
    //
    // We want as^2 + bs + c = J
    // So solve:
    // /J0\   / 1  s0  s0*s0 \ /c\
    // |J1| = | 1  s1  s1*s1 | |b|
    // \J2/   \ 1  s2  s2*s2 / \a/
    // By inverting the 3x3 matrix
    */

    TMatrixD jDistMat(3,3);
    jDistMat[0][0] = 1.0;
    jDistMat[0][1] = s[0];
    jDistMat[0][2] = s[0]*s[0];
    jDistMat[1][0] = 1.0;
    jDistMat[1][1] = s[1];
    jDistMat[1][2] = s[1]*s[1];
    jDistMat[2][0] = 1.0;
    jDistMat[2][1] = s[2];
    jDistMat[2][2] = s[2]*s[2];
    if( s[1] < 1e-5){
    	jDistMat.Print("jDistMat");
    	assert(0);
    }

    jDistMat.Invert();
  
    // Multiple the vector of J[i] by the matrix inverse to get the coefficients
    TMatrixD jVec(3,1);
    jVec[0][0] = JDir[0];
    jVec[1][0] = JDir[1];
    jVec[2][0] = JDir[2];
    TMatrixD coeffsDir = jDistMat * jVec;

    // Set the direct coefficients
    fDirCoeffs[0] = coeffsDir[0][0];
    fDirCoeffs[1] = coeffsDir[1][0];
    fDirCoeffs[2] = coeffsDir[2][0];

    // Now multiple the indirect J values by the same matrix to get the other coefficients
    TMatrixD jIndVec(3,1);
    jIndVec[0][0] = JInd[0];
    jIndVec[1][0] = JInd[1];
    jIndVec[2][0] = JInd[2];
    TMatrixD coeffsInd = jDistMat * jIndVec;

    // Finally set the indirect coeffients
    fIndCoeffs[0] = coeffsInd[0][0];
    fIndCoeffs[1] = coeffsInd[1][0];
    fIndCoeffs[2] = coeffsInd[2][0];

    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Same as the function above, but returns a vector of the coefficients as well as updating the class member variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimLikelihoodTuner::CalculateCoefficientsVector(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit )
{
    this->CalculateCoefficients(myTrack, myDigit);
    
    // Return these in a vector
    std::vector<Double_t> coeffs(6, 0.0);
    for(int j = 0; j < 3; ++j){ coeffs[j] = fDirCoeffs[j]; }
    for(int j = 3; j < 6; ++j){ coeffs[j] = fIndCoeffs[j-3]; }
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
  
  Double_t cutoff = fEmissionProfileManager->GetStoppingDistance(myTrack);
  //  std::cout << "From profile, cutoff = " << cutoff << std::endl;

  if( fConstrainExtent )
  {
    assert(fGeomType != WCSimLikelihoodDigitArray::kUnknown);

    if( fGeomType == WCSimLikelihoodDigitArray::kCylinder )
    {
      cutoff = CalculateCylinderCutoff( myTrack );
      //      std::cout << "Now cylinder cutoff = " << cutoff << std::endl;
    }
    else if( fGeomType == WCSimLikelihoodDigitArray::kMailBox )
    {
      cutoff = CalculateMailBoxCutoff( myTrack );
    }
    else assert(false);

  }
  //  std::cout << "Energy = " << myTrack->GetE() << "   Cutoff = " << cutoff << ".... returning" << std::endl;
  fCutoffIntegral = cutoff;
  fLastCutoff     = myTrack;
  return;
}

/// Calculate where the integral should be cut off if the detector is a cylinder
Double_t WCSimLikelihoodTuner::CalculateCylinderCutoff(WCSimLikelihoodTrackBase * myTrack)
{
  TVector3 vtx     = myTrack->GetVtx();
  TVector3 dir     = myTrack->GetDir();
  Double_t cutoff = fEmissionProfileManager->GetStoppingDistance(myTrack);  // Default to where the track stops

  Double_t sMax[3] = {0., 0.,0.};
 
  // Make some pseudo-2D vectors that only contain the x and y components of the vertex and direction
  TVector3 flatDir( dir.X(), dir.Y(), 0);
//  std::cout << "dir = " << std::endl;
//  dir.Print();
  double flatDirMag2 = flatDir.Mag2(); // This gets used a few times, so remember it

  // Make sure there's a component in the (x,y) plane and find the distance where the particle escapes radially
  if(flatDirMag2 > 1e-6)
  {
	  TVector3 flatVtx( vtx.X(), vtx.Y(), 0);

	  // Find where it escapes in r - this is solving the quadratic described by |flatVtx + s * flatDir| = fExtent[0]
	  Double_t xDotD = flatVtx.Dot(flatDir);
	  Double_t part1 = -1 * xDotD;
	  Double_t part2 = TMath::Sqrt( xDotD * xDotD  - flatDirMag2 * (flatVtx.Mag2() - fExtent[0]*fExtent[0])  );

	  Double_t s1    = (part1 + part2) / flatDirMag2;
	  Double_t s2    = (part1 - part2) / flatDirMag2;
	  // For a vertex inside the detector only one solution can be positive
	  if( s1 > 0) {      sMax[0] = s1; }
	  else if( s2 > 0) { sMax[0] = s2; }

	  // With cylinders we return (r,r,z)
	  sMax[1] = sMax[0];
//	  std::cout << "It's a cylinder, sR = " << sMax[0] <<std::endl;
  }
  else
  {
	  sMax[0] = cutoff;
	  sMax[1] = cutoff;
  }

  // The z coordinate
  if( dir.Z() > 1e-6 )
  {
	  // It escapes upwards
      sMax[2] = (  fExtent[2] - vtx(2) ) / ( dir(2) );
  }
  else if( dir.Z() < -1e-6 )
  {
	  // Or downwards
	  sMax[2] = ( -fExtent[2] - vtx(2) ) / ( dir(2) );
  }
  else { sMax[2] = cutoff; } // Or it doesn't go along Z at all

  
  // Figure out which of the r and z escape distance it smaller and choose that one
  if( sMax[0] > 0 && sMax[0] < cutoff) { cutoff = sMax[0]; }
  if( sMax[2] > 0 && sMax[2] < cutoff) { cutoff = sMax[2]; }

  return cutoff;
}

/// Calculate where the integral should be cut off if the detector is a box
Double_t WCSimLikelihoodTuner::CalculateMailBoxCutoff(WCSimLikelihoodTrackBase * myTrack)
{
    TVector3 vtx     = myTrack->GetVtx();
    TVector3 dir     = myTrack->GetDir();

    Double_t cutoff = fEmissionProfileManager->GetStoppingDistance(myTrack);
    Double_t sMax[3] = {0., 0.,0.};
  
  // Iterate over the x, y, and z directions to figure out when it escapes
  for(Int_t i = 0 ; i < 3; ++i )
  {
	// Particle is going forwards in this direction
	if( dir(i) > 1e-6 )
    {
        sMax[i] = (  fExtent[i] - vtx(i) ) / ( dir(i) );
    }
	// Or it's going backwards
    else if( dir(i) < -1e-6 )
    {
        sMax[i] = ( -fExtent[i] - vtx(i) ) / ( dir(i) );
        // std::cout << "Direction " << i << " (negative) cutoff is " << sMax[i] << std::endl;
    }
	// Or it doesn't have a direction component in this coordinate at all
    else { sMax[i] = cutoff; }
  }

  // The world's laziest sorting algorithm:
  if( sMax[0] < cutoff) { cutoff = sMax[0]; }
  if( sMax[1] < cutoff) { cutoff = sMax[1]; }
  if( sMax[2] < cutoff) { cutoff = sMax[2]; }

  return cutoff;
}

///////////////////////////////////////////////////////////////////////////
//  Get the integrals along the track, using the config file to decide if
//  they should be looked-up or calculated numerically
///////////////////////////////////////////////////////////////////////////

// First the contribution from Cherenkov light: get the integral with a single power of s under the integral sign
double WCSimLikelihoodTuner::GetChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit, Int_t sPower)
{
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

// Contribution from Cherenkov light: if we want all s powers (we usually do) then this is quicker as there's only one loop over digits
std::vector<Double_t> WCSimLikelihoodTuner::GetChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit)
{
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

// Contribution from indirect light, considering only one power of s
double WCSimLikelihoodTuner::GetIndIntegrals(WCSimLikelihoodTrackBase * myTrack, Int_t sPower)
{
    // Get the integral with a single power of s under the integral sign
    if(fCalculateIntegrals == true)
    {
      return this->CalculateIndIntegrals(myTrack, sPower);
    }
    else return this->LookupIndIntegrals(myTrack, sPower);
}

// Contribution from indirect light considering all powers of s (the usual case)
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
	// These are the integrals we'll return
  std::vector<Double_t> integralsVec;

  // Our track might have a conversion distance where it doesn't emit photons until it's travelled some way
  // The integral lookup doesn't know about this because it happens for photons, and photons share their 
  // lookup tables with electrons, which don't have a nonzero conversion distance
  // So we'll take care of it here instead:
  TVector3 emissionStart = myTrack->GetFirstEmissionVtx();
  double emissionCutoff = fCutoffIntegral - myTrack->GetConversionDistance();
  emissionCutoff *= (emissionCutoff > 0); // We don't want it to be negative if the particle doesn't convert until after it's escaped

  if(emissionCutoff == 0){ 
    integralsVec.push_back(0);
    integralsVec.push_back(0);
    integralsVec.push_back(0);
  }
  else{
    TVector3 pmtPos = myDigit->GetPos();
    TVector3 vtxDir = myTrack->GetDir();
	TVector3 toPMT = pmtPos - emissionStart; // The vector from where the track starts emitting Cherenkov light to the PMT
    Double_t R0 = toPMT.Mag();
    Double_t cosTheta0 = (vtxDir.Dot( toPMT ))/R0;
 	  
    double E = myTrack->GetE();
    TrackType::Type type = myTrack->GetType();
    integralsVec.push_back(WCSimIntegralLookupReader::Instance()->GetRhoGIntegral(type, E, emissionCutoff, R0, cosTheta0));
    integralsVec.push_back(WCSimIntegralLookupReader::Instance()->GetRhoGSIntegral(type, E, emissionCutoff, R0, cosTheta0));
    integralsVec.push_back(WCSimIntegralLookupReader::Instance()->GetRhoGSSIntegral(type, E, emissionCutoff, R0, cosTheta0));
    //  std::cout << "tubeID = " << myDigit->GetTubeId() << "    R0 = " << R0 << "    cosTheta0 = " << cosTheta0 << "     E = " << myTrack->GetE() << "    sMax = " << fCutoffIntegral << "  EmissionCutoff = " << emissionCutoff << std::endl;
  }
  //  if(integralsVec.at(1) != 0){
  //    std::cout << "tubeID = " << myDigit->GetTubeId() << "    s^0 term = " << integralsVec.at(0) << "   s^1 term = " << integralsVec.at(1) << "   s^2 term = " << integralsVec.at(2) << std::endl;
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
	std::vector<Double_t> integralsVec(3, 0.0);
  
  // Photons have a conversion distance so don't start emitting until they've travelled a certain length
  // Need to allow for this in where we cut off the integral
  double emissionCutoff = fCutoffIntegral - myTrack->GetConversionDistance();
  emissionCutoff *= (emissionCutoff > 0);

  if(emissionCutoff != 0){ 
 	  double E = myTrack->GetE();
 	  TrackType::Type type = myTrack->GetType();
  	  integralsVec[0] = WCSimIntegralLookupReader::Instance()->GetRhoIntegral(type, E, fCutoffIntegral);
	  integralsVec[1] = WCSimIntegralLookupReader::Instance()->GetRhoSIntegral(type, E, fCutoffIntegral);
	  integralsVec[2] = WCSimIntegralLookupReader::Instance()->GetRhoSSIntegral(type, E, fCutoffIntegral);
  }

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
    CalculateCutoff( myTrack );
    
    std::vector<Int_t> sPowers;
    sPowers.push_back(0);
    sPowers.push_back(1);
    sPowers.push_back(2);

    // The track might have a conversion distance, but because the emission profile manager knows about the TrackType, it can take care of this
    std::vector<Double_t> integrals = fEmissionProfileManager->GetRhoGIntegrals(myTrack, myDigit, sPowers, fCutoffIntegral, kTRUE);

    // std::cout << "WCSimLikelihoodTuner::CalculateChIntegrals()   s = " << fCutoffIntegral << "   " << integrals[0] << "   " << integrals[1] << "   " << integrals[2] << std::endl;
  
    return integrals;
}

///////////////////////////////////////////////////////////////////////////////////////////
// As for the previous function, but now we're calculating the integrals for indirect light
///////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::CalculateIndIntegrals(WCSimLikelihoodTrackBase * myTrack, int i)
{
  
  // The track might have a conversion distance, but because the emission profile manager knows about the TrackType, it can take care of this
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

    std::vector<Int_t> sPowers;
    sPowers.push_back(0);
    sPowers.push_back(1);
    sPowers.push_back(2);

    std::vector<Double_t> integrals = fEmissionProfileManager->GetRhoIntegrals(sPowers, myTrack, 0, fCutoffIntegral, kTRUE);
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
	if( fGeomType == WCSimLikelihoodDigitArray::kCylinder )
	{
		return ( (pos.Perp() > fExtent[0]) || (pos.Z() > fExtent[2]) );
	}
	else
	{
		return (fabs(pos.X()) > fExtent[0] || fabs(pos.Y()) > fExtent[1] || fabs(pos.Z()) > fExtent[2]);
	}
	assert(std::cerr << "Error: I don't think the detector is a cylinder or a box" << std::endl && false );
	return false;
}

Double_t WCSimLikelihoodTuner::GetLightFlux(WCSimLikelihoodTrackBase * myTrack)
{
  return fEmissionProfileManager->GetLightFlux(myTrack);
}

Double_t WCSimLikelihoodTuner::GetTrackLengthForPercentile(WCSimLikelihoodTrackBase * myTrack, const double &percentile)
{
  double length = -999.9;
  return fEmissionProfileManager->GetTrackLengthForPercentile(myTrack, percentile); 
}

Double_t WCSimLikelihoodTuner::QuantumEfficiency() {
  
  double distToPMT = fCache.GetDistanceToPMT();
  return (fCache.GetDigit()->GetAverageQE(distToPMT));
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


void WCSimLikelihoodTuner::MakeCache(const double & s, WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit)
{
	if(!fCache.Contains(s, myDigit, myTrack))
	{
		fCache = WCSimLikelihoodTunerCache(s, myDigit, myTrack);
	}
	return;
}


///////////////////////////////////////////////////////////////////////////////////////
// WCSimLikelihoodTunerCache: A lot of functions need the same vector calculations, e.g.
// take a track, travel some distance and work out how far from a given PMT it is - this
// class works them out once per uniqe combination of track, distance and PMT and then
// looks them up instead of recalculating them
WCSimLikelihoodTunerCache::WCSimLikelihoodTunerCache( double const &s, WCSimLikelihoodDigit const * digit, WCSimLikelihoodTrackBase const * track)
{
	fDigit = digit;
	fTrack = track;
	fS = s;
	fEmissionPos = fTrack->GetPropagatedPos(s);
	fToPMT = digit->GetPos() - fEmissionPos;
	fDistanceToPMT = fToPMT.Mag();
	fCosAngleToPMT = fToPMT.Dot(fTrack->GetDir()) / fDistanceToPMT;
	fAngleToPMT = TMath::ACos(fCosAngleToPMT);
}

WCSimLikelihoodTunerCache::WCSimLikelihoodTunerCache()
{
	fDigit = NULL;
	fTrack = NULL;
	fS = -999.9;
	fEmissionPos = TVector3();
	fToPMT = TVector3();
	fDistanceToPMT = -999.9;
	fCosAngleToPMT = -999.9;
	fAngleToPMT = -999.9;
}

bool WCSimLikelihoodTunerCache::Contains(double const &s, WCSimLikelihoodDigit const * digit, WCSimLikelihoodTrackBase const * track)
{
	return ((s == fS) && (*digit == *fDigit) && (track->GetVtx() == fTrack->GetVtx()) && (track->GetDir() == fTrack->GetDir() ));
}

TVector3 WCSimLikelihoodTunerCache::GetEmissionPos()
{
	return fEmissionPos;
}

TVector3 WCSimLikelihoodTunerCache::GetToPMT()
{
	return fToPMT;
}

double WCSimLikelihoodTunerCache::GetAngleToPMT()
{
	return fAngleToPMT;
}

double WCSimLikelihoodTunerCache::GetCosAngleToPMT()
{
	return fCosAngleToPMT;
}

double WCSimLikelihoodTunerCache::GetDistanceToPMT()
{
	return fDistanceToPMT;
}

WCSimLikelihoodDigit const * WCSimLikelihoodTunerCache::GetDigit()
{
	return fDigit;
}

WCSimLikelihoodTrackBase const * WCSimLikelihoodTunerCache::GetTrack()
{
	return fTrack;
}

double WCSimLikelihoodTunerCache::GetS()
{
	return fS;
}





