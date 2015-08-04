#ifndef WCSIMPARAMETERS_HH
#define WCSIMPARAMETERS_HH

#include "TObject.h"

class WCSimParameters : public TObject {

 public:
  static WCSimParameters* Instance();

  static void UseSimpleParameters();
  static void UseSimpleTimeResolution();
  static void UseSimpleTimeSlew();
  static void UseSimpleRefractiveIndex();

  static Double_t SpeedOfLight();
  static Double_t CherenkovAngle();
  static Double_t TimeResolution(Double_t Q);
  static Double_t WCSimTimeResolution(Double_t Q, Double_t timeConst);
  static Double_t TimeSlew(Double_t Q);
  static Double_t RefractiveIndex(Double_t r);

  static Double_t ThetaC();     // Cherenkov Angle
  static Double_t CosThetaC();  // Cosine of Cherenkov Angle

  static void PrintParameters();
  void RunPrintParameters();

  void SetSimpleTimeResolution() { fUseSimpleTimeResolution = 1; }
  Bool_t SimpleTimeResolution() { return fUseSimpleTimeResolution; }

  void SetSimpleTimeSlew() { fUseSimpleTimeSlew = 1; }
  Bool_t SimpleTimeSlew() { return fUseSimpleTimeSlew; }

  void SetSimpleRefractiveIndex(){ fUseSimpleRefractiveIndex = 1; }
  Bool_t SimpleRefractiveIndex() { return fUseSimpleRefractiveIndex; }

  Double_t GetTimeResolution(Double_t Q);
  Double_t GetTimeSlew(Double_t Q);
  Double_t GetRefractiveIndex(Double_t r);  
  
  Double_t GetSimpleTimeResolution(Double_t Q);
  Double_t GetSimpleTimeSlew() { return 0.0; };
  Double_t GetSimpleRefractiveIndex() { return 1.33; };

  // Slicer parameters
  void SetSlicerClusterDistance(Double_t val) {fSlicerClusterDistance = val;};
  Double_t GetSlicerClusterDistance() {return fSlicerClusterDistance;};

  void SetSlicerMinSize(UInt_t val) {fSlicerMinSize = val;};
  UInt_t GetSlicerMinSize() {return fSlicerMinSize;};

  void SetSlicerChargeCut(Double_t val) {fSlicerChargeCut = val;};
  Double_t GetSlicerChargeCut() {return fSlicerChargeCut;};

 private:
  WCSimParameters();
  ~WCSimParameters();

  Bool_t fUseSimpleTimeResolution;
  Bool_t fUseSimpleTimeSlew;
  Bool_t fUseSimpleRefractiveIndex;

  // Slicer parameters
  Double_t fSlicerClusterDistance; // Max distance in cm between hits in the slices
  UInt_t fSlicerMinSize; // Minimum number of hits for a slice to be formed
  Double_t fSlicerChargeCut; // Only consider digits above the charge threshold

  ClassDef(WCSimParameters,0)

};

#endif







