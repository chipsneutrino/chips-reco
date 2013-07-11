#include <math.h>

#include "TArrayD.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1I.h"
#include "TH1D.h"
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

////////////////////////////////////////////
// Constructor               
////////////////////////////////////////////
WCSimLikelihoodTuner::WCSimLikelihoodTuner()
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
// This uses a formula from MiniBooNE.  Efficiency of digitization is handled elsewhere
///////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::Efficiency(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
    // We need the angle of incidence at the PMT
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    TVector3 toPMT = pmtPos - emissionPos;
    Double_t theta = TMath::ACos( ( -1.0 * myDigit->GetFaceX() * toPMT[0] +
                                    -1.0 * myDigit->GetFaceY() * toPMT[1] + 
                                    -1.0 * myDigit->GetFaceZ() * toPMT[2] 
                                  ) / toPMT.Mag()) ; // angle of incidence wrt PMT normal
    theta *= 180.0 / TMath::Pi();
//    if( theta > 90.0 )
//    {
//      theta = 180.0 - theta;
//    }

    // Function for the PMT acceptance's dependence on the angle: 
    // WCSim defines arrays of efficiency at 10 degree intervals and linearly interpolates

    if( theta < 0.0 || theta >= 90.0 )
    {
      std::cout << "theta = " << theta << std::endl;
      return 0.0;
    }

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





}

//////////////////////////////////////////////////////////////////////////////////////////
// Work out the effect of solid angle on the probability of photon incidenc; pure geometry
//////////////////////////////////////////////////////////////////////////////////////////
Double_t WCSimLikelihoodTuner::SolidAngle(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit)
{
    // We need the distance from emission to the PMT
    // These are from WCSim so they're in cm
    TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
    TVector3 emissionPos(myTrack->GetX() + s*sin(myTrack->GetTheta())*cos(myTrack->GetPhi()),
                         myTrack->GetY() + s*sin(myTrack->GetTheta())*sin(myTrack->GetPhi()), 
                         myTrack->GetZ() + s*cos(myTrack->GetTheta()));
    Double_t r = (pmtPos - emissionPos).Mag();
    
    // And now some physical parameters of the phototube
    Double_t globeRadius = ((WCSimRootGeom*)(WCSimGeometry::Instance())->GetWCSimGeometry())->GetWCPMTRadius();
    Double_t globeHalfHeight = globeRadius - 0.1;
    // This is hardcoded in the detector construction and there doesn't seem to 
    // be a class that contains WCPMTExposeHeight explicitly to query it from

    // Purely geometry: we need the solid angle of a cone whose bottom is a circle of the same radius as the circle of PMT poking
    // through the blacksheet.  This is 2pi( 1 - cos(coneAngle) ) where coneAngle can be deduced from trig and the PMT geometry
    Double_t solidAngle = (1.0 - (r+globeHalfHeight)/sqrt( (r + globeHalfHeight)*(r + globeHalfHeight) + globeRadius*globeRadius));
//    std::cout << "solid angle = " << solidAngle << std::endl;
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
    return 0.1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The direct and indirect expectations come from a large integral over several factors.  We separate out the
// transmission function, efficiency and solid angle and fit these to a quadratic in s (distance along track)
// This just leaves integrals over emission profiles that can be tabulated in advance.  This function gives
// the coefficients of these quadratics
///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void WCSimLikelihoodTuner::CalculateCoefficients(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit )
{
    //std::cout << "*** WCSimLikelihoodTuner::CalculateCoefficients() *** Calculating the coefficients for the integrals from tuning info" << std::endl; 
    for(int i = 0; i < 3; ++i)
    { 
        fDirCoeffs[i] = 0.0;    // Coefficients for direct (Cherenkov) light
        fIndCoeffs[i] = 0.0;    // And for indirect light
     }
    
    // Load the file containing the emission profiles for different energies
    // whichHisto should be a TH1I with the same energy binning as we want for the reconstruction
    // There should also be a TObjArray with the emission profiles for the center value of each energy bin
    //std::cout << "Getting the histogram" << std::endl;
    this->LoadEmissionProfiles(myTrack);
    TH1D * whichHisto = (TH1D*)fProfiles->Get("hWhichHisto");
    Int_t whichBin = whichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin > whichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateCoefficients() looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }

    //std::cout << "We want number " << whichBin << std::endl;
    //std::cout << "fHistArray= " << fHistArray << std::endl;
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
            s[0] = 0.0;
            s[1] = hProfile->GetBinCenter(iBin);
            s[2] = 2 * s[1];
            break;
        }
    }
    
    // Evaluate J at each point
    Double_t JDir[3] = {0.,0.,0.};
    Double_t JInd[3] = {0.,0.,0.};
        
        // std::cout << " Second loop: kk = " << kk << "    s = " << s[kk] << "     JDir = " << JDir[kk] << std::endl;
    for(int k = 0; k < 3; ++k)
    {
        JDir[k] =  this->TransmissionFunction(s[k], myTrack, myDigit) 
                  * this->Efficiency(s[k], myTrack, myDigit)
                  * this->SolidAngle(s[k], myTrack, myDigit);
        JInd[k] = JDir[k] * this->ScatteringTable(s[k]);
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
    TH1D * whichHisto = (TH1D*)fProfiles->Get("hWhichHisto");
    Int_t whichBin = whichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin > whichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateChIntegrals() is looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }


    //std::cout << "We want number " << whichBin << std::endl;
    TH1D * hProfile = (TH1D*)fHistArray->At(whichBin);
    TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(whichBin);
 
    
    Double_t integrals[3]={0.,0.,0.};
    for(int iBin = 1; iBin <= hProfile->GetNbinsX(); ++iBin)
    {
        Double_t rho, g, s, cosTheta;
        s = hProfile->GetBinCenter(iBin);
        
        // Need to calculate theta
        TVector3 pmtPos(myDigit->GetX(), myDigit->GetY(), myDigit->GetZ());
        TVector3 vtxPos(myTrack->GetVtx());
        TVector3 vtxDir(myTrack->GetDir());
        cosTheta = TMath::Cos(vtxDir.Angle( pmtPos - (vtxPos + s * vtxDir) ));
                // Make sure the histograms in s and s, cosTheta(s) have the same binning on the s axis
        rho = hProfile->GetBinContent(iBin);

        Int_t binCosTheta = hAngularProfile->GetXaxis()->FindBin(cosTheta);
        Int_t binS = hAngularProfile->GetYaxis()->FindBin(s);
        Int_t gBin = hAngularProfile->GetBin(binCosTheta,binS,0); //see doc of TH1::GetBin
        g = hAngularProfile->GetBinContent(gBin);

//        if(0.7 < cosTheta && 00.8 > cosTheta && 125 > s)
//        {
//          printf("Getting information... \n s = %f \n cosTheta = %f \n cosThetaBin = %d \n sBin = %d \n rho = %f \n G = %f \n PMT pos = ", s, cosTheta, binCosTheta, binS, rho, g);
//          pmtPos.Print();
//          printf("\n vtxPos = ");
//          vtxPos.Print();
//          printf("\n vtxDir = ");
//          vtxDir.Print();
//        }

        integrals[0] += rho * g * hProfile->GetBinWidth(iBin);
        integrals[1] += rho * g * s * hProfile->GetBinWidth(iBin);
        integrals[2] += rho * g * s * s * hProfile->GetBinWidth(iBin);
        //std::cout << "s = " << s << std::endl;
    }

    Double_t ret = 0;
    if(i==0){ ret = integrals[0];}
    if(i==1){ ret = integrals[1];}
    if(i==2){ret = integrals[2];}
    
//    std::cout << "i = " << i << " ret = " << ret << std::endl;
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
    TH1D * whichHisto = (TH1D*)fProfiles->Get("hWhichHisto");
    Int_t whichBin = whichHisto->FindBin(myTrack->GetE()) - 1; // Histogram bins count from 1, arrays from 0
    if(whichBin < 0 || whichBin > whichHisto->GetNbinsX())
    {
        std::cerr << "Error: WCSimLikelihoodTuner::CalculateChIntegrals() is looking for a histogram bin that isn't in the array" << std::endl;
        exit(EXIT_FAILURE);
    }

    //std::cout << "We want number " << whichBin << std::endl;
    TH1D * hProfile = (TH1D*)fHistArray->At(whichBin);
    //std::cout << "We got it" << hProfile << std::endl;
 
    for(int iBin = 1; iBin <= hProfile->GetNbinsX(); ++iBin)
    {
        Double_t rho, s;
        s = hProfile->GetBinCenter(iBin);
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
/*    // Make a file for saving the integrals.  Require a name to avoid overwriting old ones inadvertantly
    if(filename == "")
    {
      std::cerr << "WCSimLikelihoodTuner::TabulateIndirectIntegrals - Error: A filename is needed" << std::endl;
      return;
    }
    filename += "";
    WCSimLikelihoodTrack * dummyTrack = new WCSimLikelihoodTrack();
    filename += dummyTrack->TrackTypeToString( myType );
    delete dummyTrack;
    TFile * f = new TFile(filename.Data());
    
    // Save the integrals in a tree indexed over energy bins 
    // so we can retrieve entries one at a time
    Int_t iEBin=0;
    Double_t rhoIntegral[3];
    TTree * t = new TTree("integrals","integrals");
    t->Branch("EnergyBin",&iEBin,"EnergyBin/I");
    t->Branch("rhoIntegral",&rhoIntegral,"rhoIntegral[3]/D");

    this->LoadEmissionProfiles(myType);
    TH1D * whichHisto = (TH1D*)fProfiles->Get("hWhichHisto");
    for( iEBin = 0; iEBin < whichHisto->GetNbinsX(); ++iEBin) // loop over energy bins
    {
      rhoIntegral[0] = 1.0; // by definition
      rhoIntegral[1] = 0.0;
      rhoIntegral[2] = 0.0;
     
      //std::cout << "We want number " << whichBin << std::endl;
      TH1D * hProfile = (TH1D*)fHistArray->At(iEBin);
 
      for(int iBin = 1; iBin <= hProfile->GetNbinsX(); ++iBin)
      {
          Double_t rho, s;
          s = hProfile->GetBinCenter(iBin);
          rho = hProfile->GetBinContent(iBin);
          rhoIntegral[1] += rho * s;
          rhoIntegral[2] += rho * s * s;
      }

      t->Fill();
    }

    whichHisto->SetName("hWhichEnergyBin");
    t->Write();
    whichHisto->Write();

    f->Close();

*/
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
/*    // Make a file for saving the integrals.  Require a name to avoid overwriting old ones inadvertantly
    if(filename == "")
    {
      std::cerr << "WCSimLikelihoodTuner::TabulateDirectIntegrals - Error: A filename is needed" << std::endl;
      return;
    }
    WCSimLikelihoodTrack * dummyTrack = new WCSimLikelihoodTrack();
    filename += dummyTrack->TrackTypeToString( myType );
    delete dummyTrack;
    TFile * f = new TFile(filename.Data());

    WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack();
    WCSimLikelihoodDigitArray * myLDArray = new WCSimLikelihoodDigitArray();
    WCSimChargeLikelihood * myCL = new WCSimChargeLikelihood( myLDArray, myTrack ) ;

    Double_t theR0Interval = myCL->GetR0Interval();
    Int_t nR0Bins = (Int_t)(ceil((40000-0)/theR0Interval));
    TH1F * hR0Bins = new TH1F("hR0Bins","hR0Bins", nR0Bins, 0.,40000.);

    Double_t theCosTheta0Interval = myCL->GetCosTheta0Interval();
    Int_t nCosTheta0Bins = (Int_t)(ceil(2./theCosTheta0Interval)); 
    TH1F * hCosTheta0Bins = new TH1F("hCosTheta0Bins","hCosTheta0Bins", nCosTheta0Bins, -1.,1.);

    // Make an array of TArrayD objects so we can quickly find the one we want
    TObjArray * R0Array = new TObjArray(nR0Bins);
    for(int i = 0; i < nR0Bins; ++i)
    {
      TObjArray * CosTheta0Array = new TObjArray(nCosTheta0Bins);
      for(int j = 0; j < nCosTheta0Bins;  ++j)
      {
         TArrayD * integrals = new TArrayD(3);
         CosTheta0Array->AddAt( (TObject*)integrals, j );
      }
      R0Array->AddAt( CosTheta0Array, i );
    }


    // Save the integrals in a tree indexed over energy bins 
    // so we can retrieve entries one at a time
    TTree * t = new TTree("integrals","integrals");
    t->Branch("R0Array","TObjArray",&R0Array);



    this->LoadEmissionProfiles(myType);
    TH1D * whichHisto = (TH1D*)fProfiles->Get("hWhichHisto");
    Double_t * gRhoIntegral = new Double_t[3];
        
    for( int iEBin = 0; iEBin < whichHisto->GetNbinsX(); ++iEBin) // loop over energy bins
    {
      for(int iR0Bin = 0; iR0Bin < nR0Bins; ++ iR0Bin)
      {
        Double_t R0 = hR0Bins->GetBinCenter( iR0Bin + 1) ;

        for(int iCosTheta0Bin = 0; iCosTheta0Bin < nCosTheta0Bins; ++iCosTheta0Bin)
        {
          Double_t gRhoIntegral[3];
          gRhoIntegral[0] = 0.0;
          gRhoIntegral[1] = 0.0;
          gRhoIntegral[2] = 0.0;
     
          //std::cout << "We want number " << whichBin << std::endl;
          TH1D * hProfile = (TH1D*)fHistArray->At(iEBin);
          TH2D * hAngularProfile = (TH2D*)fAngHistArray->At(iEBin);
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
              Int_t gBin = hAngularProfile->GetBin(binCosTheta,iBin+1,0); //see doc of TH1::GetBin
              Double_t g = hAngularProfile->GetBinContent(gBin);
              gRhoIntegral[0] += rho * g;
              gRhoIntegral[1] += rho * g * s;
              gRhoIntegral[2] += rho * g * s * s;
          }
          
          for(int iIntegral = 0; iIntegral < 3; ++iIntegral)
          {
            ((TArrayD*)((TObjArray*)R0Array->At(iR0Bin))->At(iCosTheta0Bin))->AddAt(gRhoIntegral[iIntegral], iIntegral); 
            // This line is horrible. We're getting the TArrayD for a given R0 and CosTheta0 bin and setting each of its 3 entries
            // to the integrals we just calculate.  But there's lots of unpleasant casting because we have arrays of arrays
          }
        }
      }
      t->Fill();
    }
  

    whichHisto->SetName("hWhichEnergyBin");
    t->Write();
    whichHisto->Write();
    hR0Bins->Write();
    hCosTheta0Bins->Write();
    
    f->Close();

*/
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



