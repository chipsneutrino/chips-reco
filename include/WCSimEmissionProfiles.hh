/*
 * WCSimEmissionProfiles.h
 *
 *  Created on: 13 Nov 2014
 *      Author: andy
 */

#ifndef WCSIMEMISSIONPROFILES_H_
#define WCSIMEMISSIONPROFILES_H_

#include "WCSimLikelihoodTrack.hh"
#include <string>

class TH1D;
class TH2D;
class TFile;
class TObjArray;
class WCSimLikelihoodDigit;

struct EmissionProfile_t{
	enum Type{
		kUnknown = 0,
		kRho = 1,
		kRhoTimesG = 2
	};
	Type fType;
	EmissionProfile_t(Type t) : fType(t){}
	EmissionProfile_t(const EmissionProfile_t &ep) : fType(ep.fType) {};
	operator Type () const { return fType; };

	std::string AsString(){
		return AsString(fType);
	}

	static std::string AsString(Type type){
		std::string str("");
		if( type == kUnknown )  { str = "kUnknown"; }
		if( type == kRho )  { str = "kRho"; }
		if( type == kRhoTimesG )  { str = "kRhoTimesG"; }
		return str;
	};
};

class WCSimEmissionProfiles {
public:
	WCSimEmissionProfiles();
	WCSimEmissionProfiles(WCSimLikelihoodTrack * myTrack);
  void SaveProfiles();
	void SetTrack(WCSimLikelihoodTrack * myTrack);
	Double_t GetRho(Double_t s);
	Double_t GetG(Double_t s, Double_t cosTheta);

	Double_t GetIntegral(EmissionProfile_t::Type type, WCSimLikelihoodTrack* myTrack,
			WCSimLikelihoodDigit * myDigit, Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth);
	
	Double_t GetRhoIntegral(Int_t sPower, Double_t energy, Double_t startS, Double_t endS, Bool_t multiplyByWidth = kTRUE);
	std::vector<Double_t> GetRhoIntegrals(std::vector<Int_t> sPowers, Double_t energy, Double_t startS, Double_t endS, Bool_t multiplyByWidth = kTRUE);
	Double_t GetRhoGIntegral(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit,
			                     Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth);
	std::vector<Double_t> GetRhoGIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit,
			                                  std::vector<Int_t> sPowers, Double_t cutoffS, Bool_t multiplyByWidth);
  Double_t GetTrackLengthForPercentile(Double_t percentile);

	virtual ~WCSimEmissionProfiles();
  Double_t GetLightFlux(WCSimLikelihoodTrack * myTrack) const;

  std::vector<Double_t> GetProfileEnergies() const;


private:
	UInt_t GetArrayBin(Double_t energy) const;
	UInt_t GetHistBin(Double_t energy) const;
	Double_t GetFractionThroughBin(Double_t energy) const;

	void LoadFile(WCSimLikelihoodTrack * myTrack);
	void InterpolateProfiles(WCSimLikelihoodTrack * myTrack);
	void InterpolateRho(WCSimLikelihoodTrack * myTrack);
	void InterpolateG(WCSimLikelihoodTrack * myTrack);

	Double_t GetCriticalDistance(WCSimLikelihoodTrack::TrackType type, Double_t energy) const;
	Double_t GetLightFlux(WCSimLikelihoodTrack::TrackType type, Double_t energy) const;

	Double_t fLastEnergy;
	WCSimLikelihoodTrack::TrackType fLastType;

	TString fProfileFileName;
	TFile * fProfileFile;

	TObjArray * fRhoArray;
	TObjArray * fGArray;
	TH1D * fBinningHistogram;

	TH1D * fRhoInterp;
	TH2D * fG;

  Bool_t fDebug;
	TH1D * fRhoInterpLo;
	TH1D * fRhoInterpHi;


  ClassDef(WCSimEmissionProfiles,1)

};

#endif /* WCSIMEMISSIONPROFILES_H_ */
