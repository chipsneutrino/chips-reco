/*
 * WCSimEmissionProfiles.h
 *
 *  Created on: 13 Nov 2014
 *      Author: andy
 */

#pragma once

#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimTrackParameterEnums.hh"
#include <string>

class TH1F;
class TH2F;
class TFile;
class TTree;
class WCSimLikelihoodDigit;

struct EmissionProfile_t
{
	enum Type
	{
		kUnknown = 0,
		kRho = 1,
		kRhoTimesG = 2
	};
	Type fType;
	EmissionProfile_t(Type t) : fType(t)
	{
	}
	EmissionProfile_t(const EmissionProfile_t &ep) : fType(ep.fType){};
	operator Type() const
	{
		return fType;
	};

	std::string AsString()
	{
		return AsString(fType);
	}

	static std::string AsString(Type type)
	{
		std::string str("");
		if (type == kUnknown)
		{
			str = "kUnknown";
		}
		if (type == kRho)
		{
			str = "kRho";
		}
		if (type == kRhoTimesG)
		{
			str = "kRhoTimesG";
		}
		return str;
	};
};

class WCSimEmissionProfiles
{
public:
	WCSimEmissionProfiles();
	WCSimEmissionProfiles(const TrackType::Type &type, const double &energy);
	void SetEnergy(const double &energy);
	void SaveProfiles(const TrackType::Type &type, const double &energy);
	Double_t GetRho(const Double_t &s);
	Double_t GetRhoWidth(const Double_t &s);
	Double_t GetG(const Double_t &s, const Double_t &cosTheta);
	Double_t GetGWidth(const Double_t &cosTheta);

	Double_t GetIntegral(WCSimLikelihoodTrackBase *myTrack, EmissionProfile_t::Type type,
						 WCSimLikelihoodDigit *myDigit, Int_t sPower, Double_t cutoffS, Bool_t multiplyByWidth);

	Double_t GetRhoIntegral(WCSimLikelihoodTrackBase *myTrack, Int_t sPower, Double_t startS, Double_t endS,
							Bool_t multiplyByWidth = kTRUE);
	std::vector<Double_t> GetRhoIntegrals(WCSimLikelihoodTrackBase *myTrack, std::vector<Int_t> sPowers,
										  Double_t startS, Double_t endS, Bool_t multiplyByWidth = kTRUE);
	Double_t GetRhoGIntegral(WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit, Int_t sPower,
							 Double_t cutoffS, Bool_t multiplyByWidth);
	std::vector<Double_t> GetRhoGIntegrals(WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit,
										   std::vector<Int_t> sPowers, Double_t cutoffS, Bool_t multiplyByWidth);
	Double_t GetTrackLengthForPercentile(Double_t percentile);

	virtual ~WCSimEmissionProfiles();
	Double_t GetLightFlux(const TrackType::Type &type, const double &energy) const;

	std::vector<Double_t> GetProfileEnergies() const;
	Double_t GetStoppingDistance();

	TH1F *GetRho();
	std::pair<TH2F *, TH2F *> GetG();
	TH2F *GetGCoarse();
	TH2F *GetGFine();
	TH1F *GetEnergyHist();

	TH2F *GetSCosThetaForTime();
	TH1F *GetSForTime();
	double GetTimeCosThetaMin();
	double GetTimeCosThetaMax();

private:
	UInt_t GetTreeEntry(Double_t energy) const;
	UInt_t GetHistBin(Double_t energy) const;
	Double_t GetFractionThroughBin(Double_t energy) const;

	void LoadFile(const TrackType::Type &type, const double &energy);
	void InterpolateProfiles(const double &energy);
	void InterpolateRho(const double &energy);
	void InterpolateG(const double &energy);

	Double_t GetCriticalDistance(TrackType::Type type, Double_t energy) const;

	TString fProfileFileName;
	TFile *fProfileFile;

	TTree *fProfileTree;
	TH1F *fBinningHistogram;

	TH1F *fRho;
	TH2F *fGCoarse;
	TH2F *fGFine;

	TH1F *fSForTime;
	TH2F *fSCosThetaForTime;
	double fTimeCosThetaMin;
	double fTimeCosThetaMax;

	Bool_t fDebug;

	double fLastPercentile;
	double fPercentileTrackLength;

	double fStoppingDistance;

	TrackType::Type fType;
	double fEnergy;

	ClassDef(WCSimEmissionProfiles, 1)
};
