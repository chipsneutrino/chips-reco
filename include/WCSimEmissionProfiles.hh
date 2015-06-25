/*
 * WCSimEmissionProfiles.h
 *
 *  Created on: 13 Nov 2014
 *      Author: andy
 */

#ifndef WCSIMEMISSIONPROFILES_H_
#define WCSIMEMISSIONPROFILES_H_

#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimTrackParameterEnums.hh"
#include <string>


class TH1F;
class TH2F;
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
	WCSimEmissionProfiles(WCSimLikelihoodTrackBase * myTrack);
  void SaveProfiles();
	void SetTrack(WCSimLikelihoodTrackBase * myTrack);
	Double_t GetRho(const Double_t &s);
  Double_t GetRhoWidth(const Double_t &s);
	Double_t GetG(const Double_t &s, const Double_t &cosTheta);
  Double_t GetGWidth(const Double_t &cosTheta);

	Double_t GetIntegral(EmissionProfile_t::Type type, WCSimLikelihoodTrackBase* myTrack,
			WCSimLikelihoodDigit * myDigit, Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth);
	
	Double_t GetRhoIntegral(Int_t sPower, Double_t energy, Double_t startS, Double_t endS, Bool_t multiplyByWidth = kTRUE);
	std::vector<Double_t> GetRhoIntegrals(std::vector<Int_t> sPowers, Double_t energy, Double_t startS, Double_t endS, Bool_t multiplyByWidth = kTRUE);
	Double_t GetRhoGIntegral(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit,
			                     Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth);
	std::vector<Double_t> GetRhoGIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit,
			                                  std::vector<Int_t> sPowers, Double_t cutoffS, Bool_t multiplyByWidth);
  Double_t GetTrackLengthForPercentile(Double_t percentile);

	virtual ~WCSimEmissionProfiles();
  Double_t GetLightFlux(WCSimLikelihoodTrackBase * myTrack) const;

  std::vector<Double_t> GetProfileEnergies() const;
  Double_t GetStoppingDistance(WCSimLikelihoodTrackBase * track);

  TH1F * GetRho(TrackType::Type particle, double energy);
  std::pair<TH2F*,TH2F*> GetG(TrackType::Type particle, double energy);
  TH2F * GetGCoarse(TrackType::Type particle, double energy);
  TH2F * GetGFine(TrackType::Type particle, double energy);
  TH1F * GetEnergyHist(TrackType::Type particle);

private:
	UInt_t GetArrayBin(Double_t energy) const;
	UInt_t GetHistBin(Double_t energy) const;
	Double_t GetFractionThroughBin(Double_t energy) const;

	void LoadFile(WCSimLikelihoodTrackBase * myTrack);
	void InterpolateProfiles(WCSimLikelihoodTrackBase * myTrack);
	void InterpolateRho(WCSimLikelihoodTrackBase * myTrack);
	void InterpolateG(WCSimLikelihoodTrackBase * myTrack);

	Double_t GetCriticalDistance(TrackType::Type type, Double_t energy) const;
	Double_t GetLightFlux(TrackType::Type type, Double_t energy) const;

	Double_t fLastEnergy;
	TrackType::Type fLastType;

	TString fProfileFileName;
	TFile * fProfileFile;

	TObjArray * fRhoArray;
	TObjArray * fGCoarseArray;
	TObjArray * fGFineArray;
	TH1F * fBinningHistogram;

	TH1F * fRhoInterp;
	TH2F * fGCoarse;
	TH2F * fGFine;

  Bool_t fDebug;
	TH1F * fRhoInterpLo;
	TH1F * fRhoInterpHi;


  ClassDef(WCSimEmissionProfiles,1)

};

#endif /* WCSIMEMISSIONPROFILES_H_ */
