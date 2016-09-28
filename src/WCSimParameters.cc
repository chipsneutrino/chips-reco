#include "WCSimParameters.hh"

#include "WCSimGeometry.hh"

#include <cmath>
#include <iostream>
#include <cassert>

ClassImp(WCSimParameters)

static WCSimParameters* fgParameters = 0;

WCSimParameters* WCSimParameters::Instance()
{
  if( !fgParameters ){
    fgParameters = new WCSimParameters();
  }

  if( !fgParameters ){
    assert(fgParameters);
  }

  if( fgParameters ){

  }

  return fgParameters;
}
  
WCSimParameters::WCSimParameters()
{
  fUseSimpleTimeResolution = 0;
  fUseSimpleTimeSlew = 0;
  fUseSimpleRefractiveIndex = 0;

  // Leigh: Slicer parameters
  fSlicerClusterDistance = 250.; 
  fSlicerMinSize = 25; 
  fSlicerChargeCut = 2.0;
  fSlicerTimeCut = 30.;
  fIterateSlicing = false;

  // Leigh: Veto slicing parameters
  fVetoClusterDistance = 500.;
  fVetoMinSize = 5;
  fVetoMinChargeCut = 2.0;
  fVetoMaxChargeCut = 20.0;
  fVetoTimeCut = 100.0;
}

WCSimParameters::~WCSimParameters()
{

}
  
void WCSimParameters::UseSimpleParameters()
{
  WCSimParameters::UseSimpleTimeResolution();
  WCSimParameters::UseSimpleTimeSlew();
  WCSimParameters::UseSimpleRefractiveIndex();
}

void WCSimParameters::UseSimpleTimeResolution()
{
  WCSimParameters::Instance()->SetSimpleTimeResolution();
}

void WCSimParameters::UseSimpleTimeSlew()
{
  WCSimParameters::Instance()->SetSimpleTimeSlew();
}

void WCSimParameters::UseSimpleRefractiveIndex()
{
  WCSimParameters::Instance()->SetSimpleRefractiveIndex();
}

Double_t WCSimParameters::TimeResolution(Double_t Q)
{
  if( WCSimParameters::Instance()->SimpleTimeResolution() ){
    return WCSimParameters::Instance()->GetSimpleTimeResolution(Q);
  }
  else {
    return WCSimParameters::Instance()->GetTimeResolution(Q);
  }
}

Double_t WCSimParameters::TimeSlew(Double_t Q)
{
  if( WCSimParameters::Instance()->SimpleTimeSlew() ){
    return WCSimParameters::Instance()->GetSimpleTimeSlew();
  }
  else{
    return WCSimParameters::Instance()->GetTimeSlew(Q);
  }
}

Double_t WCSimParameters::RefractiveIndex(Double_t L)
{
  if( WCSimParameters::Instance()->SimpleRefractiveIndex() ){
    return WCSimParameters::Instance()->GetSimpleRefractiveIndex();
  }
  else{
    return WCSimParameters::Instance()->GetRefractiveIndex(L);
  }
}

void WCSimParameters::PrintParameters()
{
  WCSimParameters::Instance()->RunPrintParameters();
}

void WCSimParameters::RunPrintParameters()
{
  std::cout << " *** WCSimParameters::PrintParameters() *** " << std::endl;

  std::cout << "  Reco Parameters: " << std::endl
            << "   UseSimpleTimeResolution = " << fUseSimpleTimeResolution << std::endl
            << "   UseSimpleTimeSlew = " << fUseSimpleTimeSlew << std::endl
            << "   UseSimpleRefractiveIndex = " << fUseSimpleRefractiveIndex << std::endl;

  return;
}

Double_t WCSimParameters::SpeedOfLight()
{
  return 29.9792458;  // velocity of light [cm/ns]
}

Double_t WCSimParameters::CherenkovAngle()
{
  return 42.0;  // degrees
}

Double_t WCSimParameters::ThetaC()
{
  return 42.0;  // degrees
}

Double_t WCSimParameters::CosThetaC()
{
  return 0.743144825477394244;  // return TMath::Cos(42.0*TMath::Pi()/180.0);
}

Double_t WCSimParameters::GetTimeResolution(Double_t Q)
{  
  /*
   // Old Parameterisation (lifted from WCSim)
   // ========================================
   Double_t qpes = Q;
   if( qpes<0.5 ) qpes = 0.5;
   if( qpes>32.0 ) qpes = 32.0;
   Double_t res = 0.33 + sqrt(2.0/qpes);  
  */

  /* 
   // Sep'2010: parameterisation, including scattered light:
   // ======================================================
   Double_t c0 = +0.271, c1 = +3.037, c2 = +2.543;

   // Aug'2011: re-parameterisation, excluding scattered light:
   // ========================================================
   Double_t c0 = +0.013, c1 = +3.592, c2 = -1.635; // [c2 is -ve, take care]

   // Nov'2011: re-parameterisation, for 200 kton geometry:
   // =====================================================
   Double_t c0 = -0.005, c1 = +3.634, c2 = -1.458; // [c2 is -ve, take care]
  */

  Double_t qpes = Q;
  Double_t qpesLow = 0.0;
  if( qpes<0.25 ) qpes = 0.25;
  if( qpes>40.0 ) qpes = 40.0;

  Double_t c0 = -0.005;
  Double_t c1 = +3.634;
  Double_t c2 = -1.458;

  if( c2<0.0 ){
    qpesLow = (2.0*c2/c1)*(2.0*c2/c1);
    if( qpes<qpesLow ) qpes = qpesLow;
  }

  Double_t res = c0 
               + c1/sqrt(qpes) 
               + c2/qpes;
  
  return res;
}

Double_t WCSimParameters::WCSimTimeResolution(Double_t q, Double_t timeConst){

  // Copied from WCSim
  Double_t Q = (q > 0.5) ? q : 0.5;
  Double_t timingResolution = 0.33 + sqrt(timeConst/Q);
  // looking at SK's jitter function for 20" tubes
  if (timingResolution < 0.58) timingResolution=0.58;

  return timingResolution;
}

Double_t WCSimParameters::GetTimeSlew(Double_t Q)
{   
  /*
   // Sep'2010: parameterisation, including scattered light:
   // ======================================================
   Double_t c0 = +3.406, c1 = -2.423, c2 = +0.335;

   // Aug'2011: re-parameterisation, excluding scattered light:
   // =========================================================
   Double_t c0 = +2.234, c1 = -1.362, c2 = +0.125;  

   // Nov'2011: re-parameterisation, for 200 kton geometry:
   // =====================================================
   Double_t c0 = +2.436, c1 = -1.291, c2 = +0.089;  
  */

  Double_t qpes = Q;
  if( qpes<0.25 ) qpes = 0.25;
  if( qpes>40.0 ) qpes = 40.0;

  Double_t c0 = +2.436;
  Double_t c1 = -1.291;
  Double_t c2 = +0.089;

  Double_t dt = c0 
              + c1*log(qpes) 
              + c2*log(qpes)*log(qpes);

  return dt;
}

Double_t WCSimParameters::GetRefractiveIndex(Double_t r)
{
  Double_t c = 29.98;       
  Double_t n0 = 1.33;        // Old Attempt:
  Double_t L0 = 0.0;         // 40.0     
  Double_t dndx = 0.000123;  // 0.00015  

  Double_t L = r/c;

  Double_t n = n0*(1.0+dndx*(L-L0));
  
  return n;
}

Double_t WCSimParameters::GetSimpleTimeResolution(Double_t Q)
{  
  Double_t qpes = Q;
  if( qpes<0.25 ) qpes = 0.25;
  if( qpes>64.0 ) qpes = 64.0;

  Double_t res = 2.0/sqrt(qpes);

  return res;
}
