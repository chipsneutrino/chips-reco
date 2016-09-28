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

  void SetSlicerTimeCut(Double_t val) {fSlicerTimeCut = val;};
  Double_t GetSlicerTimeCut() {return fSlicerTimeCut;};

  void SetIterateSlicing(Bool_t val) {fIterateSlicing = val;};
  Bool_t GetIterateSlicing() {return fIterateSlicing;};

  // Veto slicing parameters
  void SetVetoClusterDistance(Double_t val) {fVetoClusterDistance = val;};
  Double_t GetVetoClusterDistance() {return fVetoClusterDistance;};

  void SetVetoMinSize(UInt_t val) {fVetoMinSize = val;};
  UInt_t GetVetoMinSize() {return fVetoMinSize;};

  void SetVetoMinChargeCut(Double_t val) {fVetoMinChargeCut = val;};
  Double_t GetVetoMinChargeCut() {return fVetoMinChargeCut;};

  void SetVetoMaxChargeCut(Double_t val) {fVetoMaxChargeCut = val;};
  Double_t GetVetoMaxChargeCut() {return fVetoMaxChargeCut;};

  void SetVetoTimeCut(Double_t val) {fVetoTimeCut = val;};
  Double_t GetVetoTimeCut() {return fVetoTimeCut;};

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
  Double_t fSlicerTimeCut; // Maximum gap allowed between hits when ordered by time in order
                           // to be associated with the previous hit.
  Double_t fIterateSlicing; // Iterate the charge cut to make as many slices as we have tracks

  // Veto slicing parameters
  Double_t fVetoClusterDistance; // Max distance in cm between hits in the slices
  UInt_t fVetoMinSize; // Minimum number of hits for a slice to be formed
  Double_t fVetoMinChargeCut; // Only consider digits above the charge threshold, set initially
  Double_t fVetoMaxChargeCut; // at the min value, and iterates up to the max value.
  Double_t fVetoTimeCut; // Maximum gap allowed between hits when ordered by time in order
                           // to be associated with the previous hit.

  ClassDef(WCSimParameters,0)

};

#endif







