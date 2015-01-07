#include "WCSimDigitizerLikelihood.hh"
#include "WCSimAnalysisConfig.hh"

#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TH2.h"

#include <iostream>
#include <cassert>



#ifndef REFLEX_DICTIONARY
ClassImp(WCSimDigitizerLikelihood)
#endif

// static WCSimDigitizerLikelihood* fgWCSDigitizer = NULL;

// Instance function: if a digitizer object already exists, return a pointer to it
// Otherwise create one and return the pointer to that
// Use this once WCSA can read config files and we don't need to set the digitizer type in the constructor
// WCSimiDigitizer* WCSimDigitizerLikelihood::Instance()
// {
//   if( fgWCSDigitizer == NULL ){ fgWCSDigitizer = new WCSimDigitizerLikelihood(); }
//   if( !fgWCSDigitizer ){ assert(fgWCSDigitizer); }
//   if( fgWCSDigitizer ){ }
// 
//   return fgWCSDigitizer;
// }


/////////////////////////////////////////////////////////////////////
/// Constructor - need to initialize the DigiType_t to string map here
WCSimDigitizerLikelihood::WCSimDigitizerLikelihood() : fType(kUnknown), fPDFs(NULL)
{

  fDigiTypeNames["kSimple"]  = kSimple;
  fDigiTypeNames["kWCSim"]   = kWCSim;
  fDigiTypeNames["kUnknown"] = kUnknown;

  // Set the digitizer type
  fType = this->StringToDigiType(WCSimAnalysisConfig::Instance()->GetDigiType());
  if( fType == WCSimDigitizerLikelihood::kWCSim ) { this->OpenPDFs(); }
  std::cout << "Using the digitizer type " << fType << std::endl;/////////////////////////////////////////////////////////////////////

}



WCSimDigitizerLikelihood::~WCSimDigitizerLikelihood()
{
  if(fPDFs != NULL){ delete fPDFs; }
}



/// Parse a string into a DigiType_t
WCSimDigitizerLikelihood::DigiType_t WCSimDigitizerLikelihood::StringToDigiType( const std::string &str ) const
{
  std::map< std::string, DigiType_t>::const_iterator mapItr = fDigiTypeNames.find(str); 
  if( mapItr == fDigiTypeNames.end() ) return kUnknown;
  else return (*mapItr).second;
}



void WCSimDigitizerLikelihood::SetDigiType( WCSimDigitizerLikelihood::DigiType_t type )
{
  if( fType != type )
  {
    fType = type;
    if(type == WCSimDigitizerLikelihood::kWCSim) this->OpenPDFs();
  }
  return;
}



void WCSimDigitizerLikelihood::OpenPDFs()
{
  if( fPDFs != NULL ) delete fPDFs;
  std::string pdfFile = getenv("WCSIMANAHOME");
  pdfFile.append("/config/pePDFs.root");
  fPDFs = new TFile(pdfFile.c_str() ,"READ");
  fDigiPDF = (TH2D*)(fPDFs->Get("digiPDF"));
  return;
}



Double_t WCSimDigitizerLikelihood::GetMinus2LnL( const Double_t &undigi, const Double_t &digi )
{
  // Bail out if the digitizer type is bad
  switch( fType )
  {
    case WCSimDigitizerLikelihood::kSimple:
    {
      return this->GetSimpleMinus2LnL( undigi, digi );
      break;
    }
    case WCSimDigitizerLikelihood::kWCSim:
    {
      return this->GetWCSimMinus2LnL( undigi, digi );
      break;
    }
    case WCSimDigitizerLikelihood::kUnknown:
    {
      assert(fType != WCSimDigitizerLikelihood::kUnknown);
      break;
    }
    default:
    {
      std::cout << "Digitizer type is " << fType << std::endl;
      assert(fType == WCSimDigitizerLikelihood::kWCSim || fType == WCSimDigitizerLikelihood::kSimple);
      break;
    }
  }

  return -1.0;
}


/// Wrapper to call the appropriate likelihood methods
Double_t WCSimDigitizerLikelihood::GetLikelihood( const Double_t &undigi, const Double_t &digi )
{
  switch( fType )
  {
    case WCSimDigitizerLikelihood::kSimple:
    {
      return this->GetSimpleLikelihood( undigi, digi );
      break;
    }
    case WCSimDigitizerLikelihood::kWCSim:
    {
      return this->GetWCSimLikelihood( undigi, digi );
      break;
    }
    default:
    {
      assert(fType == WCSimDigitizerLikelihood::kWCSim || fType == WCSimDigitizerLikelihood::kSimple);
    }
  }

  return -1.0;
}

/// Wrapper to call the appropriate expectation methods
Double_t WCSimDigitizerLikelihood::GetExpectation( const Double_t &undigi )
{
  switch( fType )
  {
    case WCSimDigitizerLikelihood::kSimple:
    {
      return this->GetSimpleExpectation( undigi );
      break;
    }
    case WCSimDigitizerLikelihood::kWCSim:
    {
      return this->GetWCSimExpectation( undigi );
      break;
    }
    default:
    {
      assert(fType == WCSimDigitizerLikelihood::kWCSim || fType == WCSimDigitizerLikelihood::kSimple);
    }
  }

  return -1.0;
}



/**
 * Depending on the value of the undigitized charge, we will mimic the WCSim digitizer in 3 separate ways
 * For low (<10pe) charges we repeatedly sample the built-in 1pe distribution from WCSim.  This is essentially
 * Poisson statistics with a pedestal from the PMT
 * For mid-range (10-200pe) charges, this becomes time consuming.  The probability distribution has not yet washed
 * out to a Gaussian very well, so we use a Gauss(x)Expo function, with fitted coefficients
 * And above 200pe a pure Gaussian works well enough
 * This selects which option to use
 */
Double_t WCSimDigitizerLikelihood::GetWCSimMinus2LnL( const Double_t &undigi, const Double_t &digi )
{
  Double_t likelihood = this->GetWCSimLikelihood(undigi, digi);
  Double_t lnL = 0.0;

  // Some of our PDFs are read out of a histogram, so it can only handle probabilities above
  // roughly 1 / nEntries, or else it will return 0.  If that happens, we'll assume a Poisson-like
  // tail and take the log of that, as we can use a functional form for its log (Stirling's approx):
  if(likelihood == 0)
  {
    //std::cout << "Oh no! Likelihood = 0!" << digi << "   " << undigi << std::endl;
    lnL = (digi * TMath::Log(undigi) - undigi) - TMath::LnGamma(digi+1);
  }
  else lnL = TMath::Log( likelihood );
  //std::cout << "Taking the log:     likelihood = " << likelihood << "    lnL = " << lnL << std::endl;

  // WCSim applies a threshold function: it picks a random uniform number from 0 to 1 and compares this
  // to a polynomial function of the digitized charge (if it is < 1.1pe).  If the function is lower than
  // the random number it throws away the hit.  Hence the probability of the hit surviving is equal to its value
  // The unpleasantly-nested function comes directly from WCSim so it has to look like this, unfortunately
  Double_t digiTest = digi + 0.1;
  Double_t pKeep    = (digiTest > 1.1) ? 1 : std::min(1.0,(-0.06374+digiTest*(3.748+digiTest*(-63.23+digiTest*(452.0+digiTest*(-1449.0+digiTest*(2513.0+digiTest*(-2529.+digiTest*(1472.0+digiTest*(-452.2+digiTest*(51.34+digiTest*2.370)))))))))));
  if(pKeep < 0.) { pKeep = 0.; }
  if(pKeep > 1.) { pKeep = 1.; }
 
  // If we measured a charge we need ln(pKeep) and if we didn't measure a charge, we need ln(1-pKeep):
  
  //  std::cout << "Taking the log of ";
  //  if(digi > 0.) std::cout << "P(keep) = ln(" ;
  //  std::cout << pKeep  << ")" << std::endl;  
  lnL += (digi>0.0) ? TMath::Log(pKeep) : 0; 

  //  std::cout << "Taking the log of P(seen) = ln(" << std::min(1.0,(-0.06374+digiTest*(3.748+digiTest*(-63.23+digiTest*(452.0+digiTest*(-1449.0+digiTest*(2513.0+digiTest*(-2529.+digiTest*(1472.0+digiTest*(-452.2+digiTest*(51.34+digiTest*2.370))))))))))) << ")" << std::endl;
  //  lnL += (digiTest > 1.1) ? 0 : TMath::Log(std::min(1.0,(-0.06374+digiTest*(3.748+digiTest*(-63.23+digiTest*(452.0+digiTest*(-1449.0+digiTest*(2513.0+digiTest*(-2529.+digiTest*(1472.0+digiTest*(-452.2+digiTest*(51.34+digiTest*2.370))))))))))));
  //  if(-2.0*lnL < 0) std::cout << "Oh no!  -2.0 LnL = " << -2.0*lnL << "    digi = " << digi << "   undigi =  " << undigi << std::endl;
  //  if(undigi > 100) std::cout << undigi << "   " << digi << "   " << likelihood << "   " << lnL << std::endl;
  return -2.0 * lnL;
}


/**
 * Wrapper to work out the likelihood associated with digitization using different methods
 * depending on the predicted charge.  Note that WCSim applies a threshold function after having
 * calculated this number.  That threshold is not applied here, but in DigitizerMinus2LnL, which
 * is the function called by Calc2LnL
 */
Double_t WCSimDigitizerLikelihood::GetWCSimLikelihood( const Double_t &undigi, const Double_t &digi )
{
 	Double_t digitizerLikelihood = 0.0;

 	// WCSim puts the undigitized charge through a function that repeatedly
 	// samples a 1pe distribution.  Then it goes through a threshold function that
 	// can throw away some lower hits, then it scaled the value by 0.985 and returns
 	// this as the final result.  To get the likelihood, we therefore need the probability
 	// that the sampler gives (measured value/0.985) * the probability (measured value/0.985)
 	// passes the threshold.
 	Double_t digiEff = digi/fEfficiency;

    // if(undigi == 0.0 && digi == 0.0) digitizerLikelihood = 1.0;

 	// Get the likelihood of an undigitized hit resulting in a given digitized one
 	if( undigi < 0) { digitizerLikelihood = 0.; }
 	else if(undigi == 0.0 && digiEff == 0.0) { digitizerLikelihood = 1; }
 	else if(undigi == 0.0 && digiEff != 0.0)
    {
      // std::cerr << "Error: measured a charge when predicted charge was 0 - we can't take a likelihood of this" << std::endl;
      digitizerLikelihood = 0.00000001;
    }
    else if( undigi < 10.0 ) { digitizerLikelihood = this->GetWCSimPickerLikelihood(undigi, digiEff); }
    else if( undigi < 200. ) { digitizerLikelihood = this->GetWCSimGausExpoLikelihood(undigi, digiEff); }
    else { digitizerLikelihood = this->GetWCSimGausLikelihood(undigi, digiEff); }

 	if(digitizerLikelihood < 0.0) { digitizerLikelihood = 0.0; }
	// This is just a safeguard

 	if(digitizerLikelihood > 1.0) { digitizerLikelihood = 1.0; }
 	// There's a threshold below which the digitized PE is forced to be zero.  In this case, the above
 	// function gives 1.0 + 0.06374 so we have to force the likelihood to 1.  When hitQ gets above
 	// threshold, the function will be below 1 and the threshold probability will take effect properly

    // std::cout << "Digitizer likelihood  = " << digitizerLikelihood << std::endl;
 	return digitizerLikelihood;
}



/**
 * For charges predicted charges below 10pe we read the likelihood from histograms. These are
 * created by repeatedly sampling the 1pe PDF coded into WCSim.  For expected charges below 1pe
 * we first sample a Poisson distribution with the expected charge as its mean to work out how
 * many photons arrive at the PMT.  The WCSim simulation is a bit ropey, so this is subject
 * to change (ie. to include effects like saturation and nonlinearity)
*/
Double_t WCSimDigitizerLikelihood::GetWCSimPickerLikelihood( const Double_t &undigi, const Double_t &digi )
{
  // std::cout << " *** WCSimDigitizerLikelihood::DigitizePickerLikelihood *** " << std::endl
  //           << "Undigitized hit = " << undigi <<"    Digitized hit = " << digi << std::endl;
  Double_t maxQ = 10.0;
  if( undigi >= maxQ){ return this->GetWCSimGausExpoLikelihood(undigi, digi); }

  // std::cout << "Likelihood = " << digiPDF->GetBinContent(digiPDF->FindBin(digi)) << std::endl;
  Int_t whichXBin     = ((TAxis*)fDigiPDF->GetXaxis())->FindBin(digi);
  Int_t whichYBin     = ((TAxis*)fDigiPDF->GetYaxis())->FindBin(undigi);
  Double_t likelihood = fDigiPDF->GetBinContent(whichXBin, whichYBin);
  // std::cout << "X bin = " << whichXBin << "    Y bin = " << whichYBin << "      likelihood = " << likelihood << std::endl;

  return likelihood;
}

Double_t WCSimDigitizerLikelihood::GetWCSimGausExpoLikelihood( const Double_t &undigi, const Double_t &digi)
{
 	Double_t minQ = 10.0, maxQ = 200.0;
 	if(undigi < minQ) { return this->GetWCSimPickerLikelihood(undigi, digi); }
 	if(undigi >= maxQ) { return this->GetWCSimGausLikelihood(undigi, digi); }

    // std::cout << " *** WCSimDigitizerLikelihood::DigitizeGausExpoLikelihood *** " << std::endl
    //           << "Undigitized hit = " << undigi <<"    Digitized hit = " << digi << std::endl;

 	// Coefficients come from a fit to a repeated sample of the WCSim 1pe distribution
 	TF1 fGausExpo("fGausExpo","[0]*TMath::Exp((pow(([2]/[3]),2)/2)-((x-[1])/[3]))*TMath::Erfc((1/TMath::Sqrt(2))*(([2]/[3])-((x-[1])/[2])))",0.0, 200.0);
 	Double_t norm      = 1.055e-2 + 4.091e-4 * undigi + 5.269e-6 * undigi * undigi - 7.559e-8 * TMath::Power(undigi, 3) - 3.988e-10 * TMath::Power(undigi, 4) -7.437e-13 * TMath::Power(undigi, 5);
 	Double_t gausMean  = -3.087 + 0.9785* undigi;
 	Double_t gausSigma = 0.6997 + 0.1167 * undigi -6.255e-4 * undigi * undigi + 2.605e-6 * TMath::Power(undigi, 3) - 5.363e-9 * TMath::Power(undigi, 4) + 3.447e-12 * TMath::Power(undigi, 5);
 	Double_t expDec    = 0.7293 + 0.2400 * undigi - 5.969e-3 * undigi * undigi + 8.371e-5 * TMath::Power(undigi, 3) - 6.275e-7 * TMath::Power(undigi, 4) + 2.370e-9 * TMath::Power(undigi, 5) - 3.538e-12 * TMath::Power(undigi,6);
 	fGausExpo.SetParameters(norm, gausMean, gausSigma, expDec);
    // std::cout << "Likelihood = " << fGausExpo.Eval(digi) << "   undigi = " << undigi << "    digi = " << digi << std::endl;

 	return fGausExpo.Eval(digi);
 }

Double_t WCSimDigitizerLikelihood::GetWCSimGausLikelihood( const Double_t &undigi, const Double_t &digi )
{
  Double_t minQ = 200.0;
  if( undigi < minQ ) return this->GetWCSimGausExpoLikelihood(undigi, digi);
  // std::cout << " *** WCSimDigitizerLikelihood::DigitizeGausLikelihood *** " << std::endl
  //           << "Undigitized hit = " << undigi <<"    Digitized hit = " << digi << std::endl;

  // Coefficients come from a fit to a repeated sample of the WCSim 1pe distribution
  TF1 fGaus("fGaus","gaus(0)",200.0, 1.5*undigi);
  Double_t gausMean  = -3.087 + 0.9785* undigi;
  Double_t gausSigma = sqrt(gausMean);
  fGaus.SetParameters(gausMean, gausSigma);
  // std::cout << "Likelihood = " << fGaus.Eval(digi) << std::endl;

  return fGaus.Eval(digi);
}



Double_t WCSimDigitizerLikelihood::GetWCSimExpectation( const Double_t &undigi )
{
  Double_t predictedCharge = 0.0;
  if( undigi < 0.5 ) predictedCharge = 0.0;
  else if( 0.5 <= undigi && undigi < 10.0)
  {
    // std::cout << "Using the PDFs:	undigi charge = " << undigi;

    //FIXME: So after a long discussion with Andy, I have still no idea
    //  what this should actually do. The old implementation was surely
    //  wrong (outdated), but the function never gets actually called.
    //  But it also means that this is not a priority to fix...
    //  /mpf 27.11.14

    //--Old code--
    //TFile f("./config/pePDFs.root","READ");
    //TH1D * digiPDF = (TH1D*)((TObjArray*)f.Get("peArray"))->At((Int_t)round(undigi)-1);

    //--New code--
    this->OpenPDFs();
    Int_t whichBin = ((TAxis*)fDigiPDF->GetYaxis())->FindBin(undigi);
    //FIXME: this is quite probably wrong
    TH1D *digiPDF = fDigiPDF->ProjectionX("", whichBin, whichBin);
    predictedCharge = digiPDF->GetMean();
    // std::cout << "    predicted charge = " << predictedCharge << std::endl;
  }
  else if( 10 <= undigi && undigi < 200)
  {
    // std::cout << "Using the gaussExpo:   undigi charge = " << undigi;
    TF1 * fGausExpo    = new TF1("fGausExpo","[0]*TMath::Exp((pow(([2]/[3]),2)/2)-((x-[1])/[3]))*TMath::Erfc((1/TMath::Sqrt(2))*(([2]/[3])-((x-[1])/[2])))",0.0, 200.0);
    Double_t norm      = 1.055e-2 + 4.091e-4 * undigi + 5.269e-6 * undigi * undigi - 7.559e-8 * TMath::Power(undigi, 3) - 3.988e-10 * TMath::Power(undigi, 4) -7.437e-13 * TMath::Power(undigi, 5);
    Double_t gausMean  = -3.087 + 0.9785* undigi;
    Double_t gausSigma = 0.6997 + 0.1167 * undigi - 6.255e-4 * undigi * undigi + 2.605e-6 * TMath::Power(undigi, 3) - 5.363e-9 * TMath::Power(undigi, 4) + 3.447e-12 * TMath::Power(undigi, 5);
    Double_t expDec    = 0.7293 + 0.2400 * undigi - 5.969e-3 * undigi * undigi + 8.371e-5 * TMath::Power(undigi, 3) - 6.275e-7 * TMath::Power(undigi, 4) + 2.370e-9 * TMath::Power(undigi, 5) - 3.538e-12 * TMath::Power(undigi,6);
    fGausExpo->SetParameters(norm, gausMean, gausSigma, expDec);

    predictedCharge    =  fGausExpo->Mean(0.0, 200.0);
    delete fGausExpo;
    // std::cout << "    predicted charge = " << predictedCharge << std::endl;
  }
  else
  {
    // std::cout << "Using the gauss mean:    undigi charge = " << undigi << "    predicted charge = " << 3.087 + 0.9785* undigi << std::endl;
    predictedCharge = 3.087 + 0.9785* undigi;
  }

  return predictedCharge;
}



Double_t WCSimDigitizerLikelihood::GetSimpleLikelihood(const Double_t &undigi, const Double_t &digi)
{
  // std::cout << "Digi = " << digi << "   and undigi = " << undigi << std::endl;
  Double_t prob = TMath::Poisson( digi, undigi );
  assert(prob >= 0.);
  return prob;
}

Double_t WCSimDigitizerLikelihood::GetSimpleMinus2LnL( const Double_t &undigi, const Double_t &digi )
{
  Double_t prob = this->GetSimpleLikelihood( undigi, digi );
  if(prob > 1e-40) return -2.0 * TMath::Log(prob);
  else return -(digi * TMath::Log(undigi) - undigi - TMath::LnGamma(digi+1));
}

Double_t WCSimDigitizerLikelihood::GetSimpleExpectation( const Double_t &undigi )
{
  return undigi;
}
