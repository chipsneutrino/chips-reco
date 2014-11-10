/*
 * WCSimFitterParameters.hh
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#ifndef WCSIMFITTERPARAMETERS_HH_
#define WCSIMFITTERPARAMETERS_HH_
#include <string>
#include <cassert>
#include <map>
#include <vector>
#include <iostream>

struct FitterParameterType{
	enum Type{
		kUnknown = 0,
		kVtxX = 1,
		kVtxY = 2,
		kVtxZ = 3,
		kVtxT = 4,
		kDirTh = 5,
		kDirPhi = 6,
		kEnergy = 7
	};
	Type fType;
	FitterParameterType(Type t) : fType(t){}
	FitterParameterType(const FitterParameterType &fpt) : fType(fpt.fType) {};
	operator Type () const { return fType; };

	std::string AsString(){
		return AsString(fType);
	}

	static std::string AsString(Type type){
		std::string str("");
		if( type == kVtxX )  { str = "kVtxX"; }
		if( type == kVtxY )  { str = "kVtxY"; }
		if( type == kVtxZ )  { str = "kVtxZ"; }
		if( type == kVtxT )  { str = "kVtxT"; }
		if( type == kDirTh ) { str = "kDirTh"; }
		if( type == kDirPhi ){ str = "kDirPhi"; }
		if( type == kEnergy ){ str = "kEnergy"; }
		return str;
	};
	static FitterParameterType FromName(const char * name)
	{
		std::string nameStr(name);

		if( nameStr.compare(std::string("kVtxX")) == 0){ return FitterParameterType(kVtxX); }
		else if( nameStr.compare(std::string("kVtxY")) == 0){ return FitterParameterType(kVtxY); }
		else if( nameStr.compare(std::string("kVtxZ")) == 0){ return FitterParameterType(kVtxZ); }
		else if( nameStr.compare(std::string("kVtxT")) == 0){ return FitterParameterType(kVtxT); }
		else if( nameStr.compare(std::string("kDirTh")) == 0){ return FitterParameterType(kDirTh); }
		else if( nameStr.compare(std::string("kDirPhi")) == 0){ return FitterParameterType(kDirPhi); }
		else if( nameStr.compare(std::string("kEnergy")) == 0){ return FitterParameterType(kEnergy); }
		else{
			std::cerr << "Error, " << name << " does not correspond to a known fitter parameter" << std::endl;
			assert(0);
		}
		return FitterParameterType(kUnknown);

	}
	static std::vector<Type> GetAllAllowedTypes()
	{
		std::vector<Type> types;
		types.push_back(kVtxX);
		types.push_back(kVtxY);
		types.push_back(kVtxZ);
		types.push_back(kVtxT);
		types.push_back(kDirTh);
		types.push_back(kDirPhi);
		types.push_back(kEnergy);
		return types;
	}

};

class WCSimFitterParameter{
public:
  WCSimFitterParameter();

	WCSimFitterParameter(FitterParameterType::Type type, bool isFixed, double start,
						 double min, double max);

	virtual ~WCSimFitterParameter();

	bool GetIsFixed() const {
		return fIsFixed;
	}

	void SetIsFixed(bool isFixed) {
    std::cout << "Setting isFixed = " << isFixed << std::endl;
		fIsFixed = isFixed;
	}

	double GetMax() const {
		return fMax;
	}

	void SetMax(double max) {
    std::cout << "Setting max = " << max << std::endl;
		fMax = max;
	}

	double GetMin() const {
		return fMin;
	}

	void SetMin(double min) {
    std::cout << "Setting min = " << min << std::endl;
		fMin = min;
	}

	double GetStart() const {
		return fStart;
	}

	void SetStart(double start) {
    std::cout << "Setting start = " << start << std::endl;
		fStart = start;
	}

	FitterParameterType::Type GetType() const {
		return fType;
	}

	void SetType(FitterParameterType::Type type) {
		fType = type;
	}

  void Print();

private:
	FitterParameterType::Type fType;
	bool fIsFixed;
	double fStart;
	double fMin;
	double fMax;
};


class WCSimFitterSingleTrackParameters{

public:
	WCSimFitterSingleTrackParameters();
	virtual ~WCSimFitterSingleTrackParameters();
	void SetDefaultParameters();
	WCSimFitterParameter GetParameter(FitterParameterType::Type type);
	void SetParameter(FitterParameterType::Type type, bool isFixed, double start,
						 double min, double max);

  void SetParMin(FitterParameterType::Type type, double min);
  void SetParMax(FitterParameterType::Type type, double max);
  void SetParStart(FitterParameterType::Type type, double start);
  void SetParRange(FitterParameterType::Type type, double min, double max);
  void SetParIsFixed(FitterParameterType::Type type, bool fixIt);

	bool GetParIsFixed(FitterParameterType type);
	double GetParMax(FitterParameterType type);
	double GetParMin(FitterParameterType type);
	double GetParStart(FitterParameterType type);

	unsigned int GetNumParameters();

private:
	std::map<FitterParameterType::Type, WCSimFitterParameter> fParameters;
};

class WCSimFitterParameters {
public:
	WCSimFitterParameters();
	virtual ~WCSimFitterParameters();
  void SetNumTracks(unsigned int nTracks);
	unsigned int GetNumTracks() const { return fNumTracks;};
	unsigned int GetNumParameters() const;
	unsigned int GetNumIndependentParameters() const;
	void AddTrack(WCSimFitterSingleTrackParameters trackPars);
	WCSimFitterSingleTrackParameters * GetTrackParameters(unsigned int trackNum);

	void JoinParametersTogether(unsigned int track1, unsigned int track2, FitterParameterType::Type type);
	bool GetJoinParametersTogether(unsigned int track1, unsigned int track2, FitterParameterType::Type type);
	bool GetIsParameterJoined(unsigned int track, FitterParameterType::Type type);
	unsigned int GetTrackIsJoinedWith( unsigned int track, FitterParameterType::Type);

private:
	unsigned int fNumTracks;
	unsigned int fNumParameters;
	std::vector<WCSimFitterSingleTrackParameters> fTrackPars;
	std::map<std::pair<unsigned int, unsigned int>, std::vector<FitterParameterType::Type> > fJoinedParams;
};

#endif /* WCSIMFITTERPARAMETERS_HH_ */
