#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <exception>

#include "WCSimChargeLikelihood.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"

#include "TAxis.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTree.h"


///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimChargeLikelihood::WCSimChargeLikelihood(WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack )
{   
    //std::cout << "*** WCSimChargeLikelihood::WCSimChargeLikelihood() *** Created charge likelihood object" << std::endl;
    this->Initialize( myDigitArray, myTrack );
}

///////////////////////////////////////////////////////////////////////////
// Set starting values.  Eventually WCSimLikelihoodTuner will be used
// to work these out, but for now they're hardcoded
///////////////////////////////////////////////////////////////////////////
void WCSimChargeLikelihood::Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack )
{
//    std::cout << "*** WCSimChargeLikelihood::Initialize() *** Initializing charge likelihood with tuned values" << std::endl; 
    fCalculateIntegrals = true;    // Setting true will calculate integrals numerically instead of looking them up

    fCoeffsCh.push_back(1);
    fCoeffsCh.push_back(2);
    fCoeffsCh.push_back(3);

    fCoeffsInd.push_back(4);
    fCoeffsInd.push_back(5);
    fCoeffsInd.push_back(6);
 
    fTrack = myTrack;
    fGotTrackParameters = false;
    fDigitArray = myDigitArray;
    
    

    fEnergy = fTrack->GetE();
    
    fEfficiency = 0.985;


    fTuner = new WCSimLikelihoodTuner(fDigitArray->GetExtent(0),
                                      fDigitArray->GetExtent(1),
                                      fDigitArray->GetExtent(2),
                                      fCalculateIntegrals);

	
    TFile * f = new TFile("./config/pePDFs.root","READ");
    fDigiPDF = (TH2D*)(f->Get("digiPDF"));
    fDigiPDF->SetDirectory(gROOT);
    delete f;
  
    return;
}   

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimChargeLikelihood::~WCSimChargeLikelihood()
{
	delete fTuner;   
  delete fDigiPDF; 
}

///////////////////////////////////////////////////////////////////////////
// Calculate the (log) likelihood of the measured PMT hit
// for a given set of track parameters.  First we have to calculate the
// expected charge given the parameters, then it's just Poisson statistics
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::Calc2LnL()
{

//  std::cout << "fDigiPDF = " << fDigiPDF << "   " << fCoeffsCh.at(2) << std::endl;
//  fDigiPDF->Draw("COLZ");
  //std::cout << "*** WCSimChargeLikelihood::Calc2LnL() *** Calculating the charge log likelihood" << std::endl; 
  Double_t Charge2LnL=0.0;

  Double_t X, Y, Z, Qobs, Qpred;
  std::string str("likelihoodDoubleBins.root");
//  str->Form("likelihood_%f_%f.root", fEnergy, fTrack->GetZ());
  TFile * file = new TFile(str.c_str(),"RECREATE");
  TTree * t = new TTree("likelihood","likelihood");
  t->Branch("X",&X,"X/D");
  t->Branch("Y",&Y,"Y/D");
  t->Branch("Z",&Z,"Z/D");
  t->Branch("Qobs",&Qobs,"Qobs/D");
  t->Branch("Qpred",&Qpred,"Qpred/D");


  for(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
    {	
      fDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(iDigit);
      fGotTrackParameters = false;
      Double_t chargeObserved = fDigit->GetQ();
      //std::cout << "chargeObserved = " << chargeObserved << std::endl;
      // This is where we need to work out the digit to vertex parameters
      this->GetTrackParameters();



  /*    //std::cout << "Source position:   (" << sourcePos[0] << "," << sourcePos[1] << "," << sourcePos[2] << ")" << std::endl
                << "Source direction:  (" << TMath::Sin(fTrack->GetTheta()) * TMath::Cos(fTrack->GetPhi()) << "," << TMath::Sin(fTrack->GetTheta()) * TMath::Sin(fTrack->GetPhi()) << "," << TMath::Cos(fTrack->GetTheta()) << ")" << std::endl
                << "Digit position:    (" << myDigit->GetX() << "," << myDigit->GetY() << "," << myDigit->GetZ() << ")" << std::endl
                << "R0 = " << fR0 << std::endl << "CosTheta0 = " << fCosTheta0 << std::endl;
  */
      
  //    likelihood *= TMath::Poisson( chargeObserved, this->ChargeExpectation());
      Double_t chargeExpected = this->ChargeExpectation();
      
      
//  		Double_t likelihood = TMath::Poisson(chargeExpected, chargeObserved);
      Charge2LnL += this->DigitizerMinus2LnL( chargeExpected, chargeObserved);
    
//      if((iDigit % 1) == 0){ std::cout << "iDigit  = " << iDigit << "/" << fDigitArray->GetNDigits() << "  Observed charge = " << chargeObserved << "   Expected charge = " << chargeExpected << "    ... So -2lnL now = " << Charge2LnL << std::endl;  }
      X = fDigit->GetX();
      Y = fDigit->GetY();
      Z = fDigit->GetZ();
      Qobs = chargeObserved;
      Qpred = chargeExpected;//Sthis->DigitizerExpectation(chargeExpected);
      t->Fill();
	 }
     std::cout << "Writing file to " << str.c_str() << std::endl; 
     file->cd();
     t->Write();
     file->Close();
     if(file) delete file;
  return Charge2LnL;
}


///////////////////////////////////////////////////////////////////////////
// Charge at a given PMT is Poisson distributed.  This calculates the 
// expected mean by summing direct and indirect light contributions.
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::ChargeExpectation()
{

    //std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the PMT" << std::endl; 
    double muDir=0, muIndir=0;
    if(fGotTrackParameters)
    {
      muDir = this->GetMuDirect();
      muIndir = this->GetMuIndirect(); 
 //     if(muDir) std::cout << " Direct charge expectation = " << muDir << " and indirect expectation = " << muIndir << "    total = " << muDir+muIndir << std::endl;
    }
    else{std::cout << "Error: did not get track parameters first" << std::endl; }

    // Sometimes the integral is very slightly negative, we'll return zero in these cases:
    if( muDir < 0 )
    {
      std::cerr << "Warning: predicted direct charge was less than zero (" << muDir <<") - returning zero" << std::endl;
      muDir = 0;
    }
    if( muIndir < 0 )
    {
      std::cerr << "Warning: predicted indirect charge was less than zero (" << muIndir <<") - returning zero" << std::endl;
      muIndir = 0;
    }

  return (muDir+muIndir);
}

///////////////////////////////////////////////////////////////////////////
// Calculate the expected contribution to the Poisson mean from direct
// Cherenkov light
///////////////////////////////////////////////////////////////////////////

double WCSimChargeLikelihood::GetMuDirect( )
{
  //std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl; 
    double muDirect = 0.0 ;
  if(!fGotTrackParameters) this->GetTrackParameters();
  if(fGotTrackParameters)
  {
    double lightFlux = this->GetLightFlux();
//    std::cout << "Energy = " << fEnergy << "  and light flux = " << lightFlux <<std::endl;
    //std::cout << " **** CALCULATING COEFFICIENTS ***** " << std::endl;
    fTuner->CalculateCoefficients(fTrack, fDigit);

    std::vector<Double_t> integrals;
    if(!fCalculateIntegrals)
    {
        integrals = fTuner->LookupChIntegrals(fTrack, fDigit);
    }
    else
    {
        integrals = fTuner->CalculateChIntegrals(fTrack, fDigit);
    }

    //std::cout << "i0 = " << i0 << "  i1 = " << i1 << "   i2 = " << i2 << std::endl
    //          << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl;
    
    muDirect = lightFlux * ( integrals[0] * fCoeffsCh[0] + integrals[1] * fCoeffsCh[1] + integrals[2] * fCoeffsCh[2] );
    if(muDirect < 0)
    {
        std::cout << "muDir is NEGATIVE! " << "  i0 = " << integrals[0] << "   i1 = " << integrals[1] << "   i2 = " << integrals[2] << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
        muDirect = 0;        
    }

  
  }
  else std::cout << "Error: did not get track parameters first, returning 0" << std::endl; 
    //std::cout << "muDirect = " << muDirect << std::endl;
  return muDirect;

}

///////////////////////////////////////////////////////////////////////////
// Calculate the expected contribution to the Poisson mean from indirect
// (eg. scatter, reflected) Cherenkov light
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::GetMuIndirect()
{

  //std::cout << "*** WCSimChargeLikelihood::GetMuIndirect() *** Calculating the indirect light contribution to the expected charge" << std::endl; 
  double muIndirect = 0.0;
  if(!fGotTrackParameters) this->GetTrackParameters(  );
  if(fGotTrackParameters)
  { 
    double lightFlux = this->GetLightFlux();
	std::vector<Double_t> integrals;
    if(!fCalculateIntegrals)
    {
		integrals = fTuner->LookupIndIntegrals(fTrack);
    }
    else
    {   
		integrals = fTuner->CalculateIndIntegrals(fTrack);
    }
    muIndirect = lightFlux * ( integrals[0] * fCoeffsInd[0] + integrals[1] * fCoeffsInd[1] + integrals[2] * fCoeffsInd[2] );
    if(muIndirect < 0)
    {
        std::cout << "IT'S NEGATIVE! " << "  i0 = " << integrals[0] << "   i1 = " << integrals[1] << "   i2 = " << integrals[2] << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
        muIndirect = 0;        
    }
//    std::cout << "Mu indirect = " << muIndirect << std::endl;
  }
  else std::cout << "Error: did not get track parameters first, returning 0" << std::endl;
  return muIndirect;
}

///////////////////////////////////////////////////////////////////////////
// Cherenkov flux is a function of energy and is prefactor to the expected
// number of PE - this calculates it
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::GetLightFlux()
{
  //std::cout << "*** WCSimChargeLikelihood::GetLightFlux() *** Getting the light flux as a function of the track KE" << std::endl; 
  
  
   Double_t factorFromG = 1/(4.0*TMath::Pi());
  // We normalize g to 1, but it's angular we've swept away a phi term, so need 1/(2pi)
  
   Double_t factorFromSolidAngle = 1; //2*TMath::Pi(); 
  // There's an overall normalization of 2pi neglected because the solid angle subtended by a cone is 2pi(1 - cosTheta)
  // I've absorbed this into SolidAngle now

   // Double_t nPhotons = 341.726*fEnergy;  
//   Double_t nPhotons = 36.9427*fEnergy - 3356.71;
  Double_t nPhotons = 17.91 + 39.61*fEnergy;
  if(nPhotons < 0) return 0;  
  // From a linear fit to nPhotons as a function of energy (it's very linear)
  // 
  // Calculating these is a waste of time because they're so simple, but this is where the factors come from
  
  return factorFromG * factorFromSolidAngle * nPhotons;
}




///////////////////////////////////////////////////////////////////////////
// Instead of working out another set of integrals for scattered light, we
// note that it must be proportional to the amount of direct light, and
// depend on the position in the tank of its emission, and on its source
///////////////////////////////////////////////////////////////////////////
Double_t WCSimChargeLikelihood::ScatteringTable()
{
  // THIS ISN'T USED YET - THERE'S AN EQUIVALENT IN WCSIMLIKELIHOODTUNER
  // We'll turn this off for the time being...
  return 0.05;
}


void WCSimChargeLikelihood::GetTrackParameters( )
{
  // Get the coefficients to multiply the integrals by
  std::vector<Double_t> coeffs = fTuner->CalculateCoefficientsVector(fTrack, fDigit);
  
  fCoeffsCh[0] = coeffs[0];
  fCoeffsCh[1] = coeffs[1];
  fCoeffsCh[2] = coeffs[2];

  fCoeffsInd[0] = coeffs[3];
  fCoeffsInd[1] = coeffs[4];
  fCoeffsInd[2] = coeffs[5];
//  std::cout << "*** WCSimShcargeLikelihood::GetTrackParameters() *** " << std::endl
//            << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl
//            << "fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;

  fGotTrackParameters = true;
  return;

}

//////////////////////////////////////////////////////////////////////////////
// Calculate the likelihood for a given PMT hit Q and a known track, using the 
// values of the individual functions rather than tabulated integrals.  
// We won't use this in the fit, but it's useful for validation
//////////////////////////////////////////////////////////////////////////////
Double_t WCSimChargeLikelihood::CalculateExactLikelihood()
{
    if(!fGotTrackParameters){ this->GetTrackParameters(); }
    fCalculateIntegrals = true;
    this->Calc2LnL();
    // Charge expectation = integral (flux * rho * solid angle * transmission * acceptance * angular emission * ds)
    return 0; 

}


// Depending on the value of the undigitized charge, we will mimic the WCSim digitizer in 3 separate ways
// For low (<10pe) charges we repeatedly sample the built-in 1pe distribution from WCSim.  This is essentially
// Poisson statistics with a pedestal from the PMT
// For mid-range (10-200pe) charges, this becomes time consuming.  The probability distribution has not yet washed
// out to a Gaussian very well, so we use a Gauss(x)Expo function, with fitted coefficients
// And above 200pe a pure Gaussian works well enough
Double_t WCSimChargeLikelihood::DigitizerMinus2LnL( const Double_t &myUndigiHit, const Double_t &myDigiHit )
{
  Double_t likelihood = this->DigitizerLikelihood(myUndigiHit, myDigiHit);
  Double_t lnL = 0.0;

  // Some of our PDFs are read out of a histogram, so it can only handle probabilities above 
  // roughly 1 / nEntries, or else it will return 0.  If that happens, we'll assume a Poisson-like
  // tail and take the log of that, as we can use a functional form for its log:
  if(likelihood == 0)
  {
//    std::cout << "Oh no! Likelihood = 0!" << myDigiHit << "   " << myUndigiHit << std::endl;
    lnL = (myDigiHit * TMath::Log(myUndigiHit) - myUndigiHit) - TMath::LnGamma(myDigiHit+1);
  }
  else lnL = TMath::Log( likelihood );
   
  // WCSim applies a threshold function: it picks a random uniform number from 0 to 1 and compares this
  // to a polynomial function of the digitized charge (if it is < 1.1pe).  If the function is lower than
  // the random number it throws away the hit.  Hence the probability of the hit surviving is equal to its value
  // The unpleasantly-nested function comes directly from WCSim so it has to look like this, unfortunately
  lnL += (myDigiHit > 1.1) ? 1.0 : TMath::Log(1.0-(-0.06374+myDigiHit*(3.748+myDigiHit*(-63.23+myDigiHit*(452.0+myDigiHit*(-1449.0+myDigiHit*(2513.0+myDigiHit*(-2529.+myDigiHit*(1472.0+myDigiHit*(-452.2+myDigiHit*(51.34+myDigiHit*2.370)))))))))));
//  if(-2.0*lnL < 0) std::cout << "Oh no!  -2.0 LnL = " << -2.0*lnL << "    myDigiHit = " << myDigiHit << "   myUndigiHit =  " << myUndigiHit << std::endl;
//  if(myUndigiHit > 100) std::cout << myUndigiHit << "   " << myDigiHit << "   " << likelihood << "   " << lnL << std::endl;
  return -2.0 * lnL;  
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper to work out the likelihood associated with digitization using different methods
// depending on the predicted charge.  Note that WCSim applies a threshold function after having
// calculated this number.  That threshold is not applied here, but in DigitizerMinus2LnL, which
// is the function called by Calc2LnL
////////////////////////////////////////////////////////////////////////////////////////////////
 Double_t WCSimChargeLikelihood::DigitizerLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit )
 {
 	Double_t digitizerLikelihood = 0.0;
 		
 	// WCSim puts the undigitized charge through a function that repeatedly
 	// samples a 1pe distribution.  Then it goes through a threshold function that
 	// can throw away some lower hits, then it scaled the value by 0.985 and returns
 	// this as the final result.  To get the likelihood, we therefore need the probability
 	// that the sampler gives (measured value/0.985) * the probability (measured value/0.985)
 	// passes the threshold.
 	Double_t myDigiHitEff = myDigiHit/fEfficiency; 						
 	
// 	if(myUndigiHit == 0.0 && myDigiHit == 0.0) digitizerLikelihood = 1.0;
 	
 	// Get the likelihood of an undigitized hit resulting in a given digitized one
 	if( myUndigiHit < 0) digitizerLikelihood = 0.;
 	else if(myUndigiHit == 0.0 && myDigiHitEff == 0.0) digitizerLikelihood = 1;
 	else if(myUndigiHit == 0.0 && myDigiHitEff != 0.0) std::cerr << "Error: measured a charge when predicted charge was 0 - we can't take a likelihood of this" << std::endl;
 	else if( myUndigiHit < 10.0 ) digitizerLikelihood = this->DigitizePickerLikelihood(myUndigiHit, myDigiHitEff);
 	else if( myUndigiHit < 200. ) digitizerLikelihood = this->DigitizeGausExpoLikelihood(myUndigiHit, myDigiHitEff);
 	else digitizerLikelihood = this->DigitizeGausLikelihood(myUndigiHit, myDigiHitEff);
 	
 	if(digitizerLikelihood < 0.0) digitizerLikelihood = 0.0;
	// This is just a safeguard
	
 	if(digitizerLikelihood > 1.0) digitizerLikelihood = 1.0;
 	// There's a threshold below which the digitized PE is forced to be zero.  In this case, the above
 	// function gives 1.0 + 0.06374 so we have to force the likelihood to 1.  When hitQ gets above
 	// threshold, the function will be below 1 and the threshold probability will take effect properly
 
 
 	// std::cout << digitizerLikelihood << std::endl;
 	return digitizerLikelihood;
 }

///////////////////////////////////////////////////////////////////////////////////////////////
// For charges predicted charges below 10pe we read the likelihood from histograms. These are 
// created by repeatedly sampling the 1pe PDF coded into WCSim.  For expected charges below 1pe 
// we first sample a Poisson distribution with the expected charge as its mean to work out how
// many photons arrive at the PMT.  The WCSim simulation is a bit ropey, so this is subject
// to change (ie. to include effects like saturation and nonlinearity)
///////////////////////////////////////////////////////////////////////////////////////////////
 Double_t WCSimChargeLikelihood::DigitizePickerLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit )
 {	
  // 	std::cout << " *** WCSimChargeLikelihood::DigitizePickerLikelihood *** " << std::endl
  //			      << "Undigitized hit = " << myUndigiHit <<"    Digitized hit = " << myDigiHit << std::endl; 
 	Double_t maxQ = 10.0;
 	if( myUndigiHit >= maxQ) return this->DigitizeGausExpoLikelihood(myUndigiHit, myDigiHit);

 	
//	std::cout << "Likelihood = " << digiPDF->GetBinContent(digiPDF->FindBin(myDigiHit)) << std::endl;	
  Int_t whichXBin = ((TAxis*)fDigiPDF->GetXaxis())->FindBin(myDigiHit);
  Int_t whichYBin = ((TAxis*)fDigiPDF->GetYaxis())->FindBin(myUndigiHit);
  Double_t likelihood = fDigiPDF->GetBinContent(whichXBin, whichYBin);	
 	return likelihood;
 }
 
 Double_t WCSimChargeLikelihood::DigitizeGausExpoLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit)
 {	
 	Double_t minQ = 10.0, maxQ = 200.0;
 	if(myUndigiHit < minQ) return this->DigitizePickerLikelihood(myUndigiHit, myDigiHit);
 	if(myUndigiHit >= maxQ) return this->DigitizeGausLikelihood(myUndigiHit, myDigiHit); 

//	std::cout << " *** WCSimChargeLikelihood::DigitizeGausExpoLikelihood *** " << std::endl
//			  << "Undigitized hit = " << myUndigiHit <<"    Digitized hit = " << myDigiHit << std::endl; 


 	TF1 fGausExpo("fGausExpo","[0]*TMath::Exp((pow(([2]/[3]),2)/2)-((x-[1])/[3]))*TMath::Erfc((1/TMath::Sqrt(2))*(([2]/[3])-((x-[1])/[2])))",0.0, 200.0);
 	Double_t norm = 1.055e-2 + 4.091e-4 * myUndigiHit + 5.269e-6 * myUndigiHit * myUndigiHit - 7.559e-8 * TMath::Power(myUndigiHit, 3) - 3.988e-10 * TMath::Power(myUndigiHit, 4) -7.437e-13 * TMath::Power(myUndigiHit, 5); 
 	Double_t gausMean = -3.087 + 0.9785* myUndigiHit;
 	Double_t gausSigma = 0.6997 + 0.1167 * myUndigiHit -6.255e-4 * myUndigiHit * myUndigiHit + 2.605e-6 * TMath::Power(myUndigiHit, 3) - 5.363e-9 * TMath::Power(myUndigiHit, 4) + 3.447e-12 * TMath::Power(myUndigiHit, 5);
 	Double_t expDec = 0.7293 + 0.2400 * myUndigiHit - 5.969e-3 * myUndigiHit * myUndigiHit + 8.371e-5 * TMath::Power(myUndigiHit, 3) - 6.275e-7 * TMath::Power(myUndigiHit, 4) + 2.370e-9 * TMath::Power(myUndigiHit, 5) - 3.538e-12 * TMath::Power(myUndigiHit,6);
 	fGausExpo.SetParameters(norm, gausMean, gausSigma, expDec); 
//	std::cout << "Likelihood = " << fGausExpo.Eval(myDigiHit) << "   myUndigiHit = " << myUndigiHit << "    myDigiHit = " << myDigiHit << std::endl;	
 	return fGausExpo.Eval(myDigiHit);
 }
 
 Double_t WCSimChargeLikelihood::DigitizeGausLikelihood( const Double_t &myUndigiHit, const Double_t &myDigiHit )
 {	
 	Double_t minQ = 200.0;
 	if( myUndigiHit < minQ ) return this->DigitizeGausExpoLikelihood(myUndigiHit, myDigiHit);
// 	std::cout << " *** WCSimChargeLikelihood::DigitizeGausLikelihood *** " << std::endl
//			  << "Undigitized hit = " << myUndigiHit <<"    Digitized hit = " << myDigiHit << std::endl; 

 	TF1 fGaus("fGaus","gaus(0)",200.0, 1.5*myUndigiHit);
 	Double_t gausMean = -3.087 + 0.9785* myUndigiHit;
 	Double_t gausSigma = sqrt(gausMean);
 	fGaus.SetParameters(gausMean, gausSigma);
//	std::cout << "Likelihood = " << fGaus.Eval(myDigiHit) << std::endl;	
 	return fGaus.Eval(myDigiHit);
 }
 
 Double_t WCSimChargeLikelihood::DigitizerExpectation( const Double_t &myUndigiHit )
 {
 	Double_t predictedCharge = 0.0;
 	if( myUndigiHit < 0.5 ) predictedCharge = 0.0;
 	else if( 0.5 <= myUndigiHit && myUndigiHit < 10.0)
 	{
// 		std::cout << "Using the PDFs:	undigi charge = " << myUndigiHit;
 		TFile f("./config/pePDFs.root","READ");
 		TH1D * digiPDF = (TH1D*)((TObjArray*)f.Get("peArray"))->At((Int_t)round(myUndigiHit)-1);
 		predictedCharge = digiPDF->GetMean();
// 		std::cout << "    predicted charge = " << predictedCharge << std::endl;
 		f.Close();
 	}
	else if( 10 <= myUndigiHit && myUndigiHit < 200)
	{

//		std::cout << "Using the gaussExpo:   undigi charge = " << myUndigiHit;
 		TF1 * fGausExpo = new TF1("fGausExpo","[0]*TMath::Exp((pow(([2]/[3]),2)/2)-((x-[1])/[3]))*TMath::Erfc((1/TMath::Sqrt(2))*(([2]/[3])-((x-[1])/[2])))",0.0, 200.0);
	 	Double_t norm = 1.055e-2 + 4.091e-4 * myUndigiHit + 5.269e-6 * myUndigiHit * myUndigiHit - 7.559e-8 * TMath::Power(myUndigiHit, 3) - 3.988e-10 * TMath::Power(myUndigiHit, 4) -7.437e-13 * TMath::Power(myUndigiHit, 5); 
 		Double_t gausMean = -3.087 + 0.9785* myUndigiHit;
 		Double_t gausSigma = 0.6997 + 0.1167 * myUndigiHit - 6.255e-4 * myUndigiHit * myUndigiHit + 2.605e-6 * TMath::Power(myUndigiHit, 3) - 5.363e-9 * TMath::Power(myUndigiHit, 4) + 3.447e-12 * TMath::Power(myUndigiHit, 5);
 		Double_t expDec = 0.7293 + 0.2400 * myUndigiHit - 5.969e-3 * myUndigiHit * myUndigiHit + 8.371e-5 * TMath::Power(myUndigiHit, 3) - 6.275e-7 * TMath::Power(myUndigiHit, 4) + 2.370e-9 * TMath::Power(myUndigiHit, 5) - 3.538e-12 * TMath::Power(myUndigiHit,6);
 		fGausExpo->SetParameters(norm, gausMean, gausSigma, expDec); 
 		predictedCharge =  fGausExpo->Mean(0.0, 200.0);
    delete fGausExpo;
//		std::cout << "    predicted charge = " << predictedCharge << std::endl;
	}
	else
	{
//		std::cout << "Using the gauss mean:    undigi charge = " << myUndigiHit << "    predicted charge = " << 3.087 + 0.9785* myUndigiHit << std::endl;
		predictedCharge = 3.087 + 0.9785* myUndigiHit;
	}	
	return predictedCharge;
		
 }
