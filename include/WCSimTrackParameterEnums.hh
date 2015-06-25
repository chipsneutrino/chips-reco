/*
 * WCSimTrackParameterEnums.hh
 *
 *  Created on: 15 Jan 2015
 *      Author: ajperch
 */

#ifndef WCSIMTRACKPARAMETERENUMS_HH_
#define WCSIMTRACKPARAMETERENUMS_HH_
#include <cassert>
#include <string>
#include <iostream>
#include <vector>

struct TrackType{
	enum Type{
		Unknown = 0,
    ElectronLike = 1,
    MuonLike = 2,
    PhotonLike = 3
	};
	Type fType;
	TrackType(Type t) : fType(t){}
	TrackType(const TrackType &fpt) : fType(fpt.fType) {};
	operator Type () const { return fType; };

	std::string AsString(){
		return AsString(fType);
	}

	static std::string AsString(Type type){
		std::string str("");
		if( type == Unknown )  { str = "Unknown"; }
		if( type == ElectronLike )  { str = "ElectronLike"; }
		if( type == MuonLike )  { str = "MuonLike"; }
		if( type == PhotonLike )  { str = "PhotonLike"; }
		return str;
	};
	static TrackType FromName(const char * name)
	{
    std::cout << "Getting type from name: " << name << std::endl;
		std::string nameStr(name);

		if( nameStr.compare(std::string("Unknown")) == 0){ return TrackType(Unknown); }
		else if( nameStr.compare(std::string("ElectronLike")) == 0){ return TrackType(ElectronLike); }
		else if( nameStr.compare(std::string("MuonLike")) == 0){ return TrackType(MuonLike); }
		else if( nameStr.compare(std::string("PhotonLike")) == 0){ return TrackType(PhotonLike); }
		else{
			std::cerr << "Error, " << name << " does not correspond to a known fitter parameter" << std::endl;
			assert(0);
		}
		return TrackType(Unknown);

	}
	static std::vector<Type> GetAllAllowedTypes(TrackType trackType)
	{

		std::vector<Type> types;
		types.push_back(Unknown);
		types.push_back(ElectronLike);
		types.push_back(MuonLike);
		types.push_back(PhotonLike);
		return types;
	}

  static Type GetTypeFromPDG(const int &pdg)
  {
    if( abs(pdg) == 11 ) { return ElectronLike; }
    else if( abs(pdg) == 13 ) { return MuonLike; }
    else if( abs(pdg) == 22 ) { return PhotonLike; }
    else
    {
      std::cerr << "Error: unknown track type for PDG code " << pdg << std::endl;
      assert(0);
    }
    return Unknown;
  }

  static int GetPDGFromType(const Type &type)
  {
    if(type == ElectronLike) { return 11; }
    if(type == MuonLike) { return 13; }
    if(type == PhotonLike) { return 22; }
    else
    {
      std::cerr << "Did not recognise track type " << TrackType(type).AsString() << std::endl;
      assert(0);
    }
    return 0;
  }

};

struct FitterParameterType{
	enum Type{
		kUnknown = 0,
		kVtxX = 1,
		kVtxY = 2,
		kVtxZ = 3,
		kVtxT = 4,
		kDirTh = 5,
		kDirPhi = 6,
		kEnergy = 7,
		kConversionDistance = 8
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
		if( type == kConversionDistance ) { str = "kConversionDistance" ; }
		return str;
	};
	static FitterParameterType FromName(const char * name)
	{
    std::cout << "Getting type from name: " << name << std::endl;
		std::string nameStr(name);

		if( nameStr.compare(std::string("kVtxX")) == 0){ return FitterParameterType(kVtxX); }
		else if( nameStr.compare(std::string("kVtxY")) == 0){ return FitterParameterType(kVtxY); }
		else if( nameStr.compare(std::string("kVtxZ")) == 0){ return FitterParameterType(kVtxZ); }
		else if( nameStr.compare(std::string("kVtxT")) == 0){ return FitterParameterType(kVtxT); }
		else if( nameStr.compare(std::string("kDirTh")) == 0){ return FitterParameterType(kDirTh); }
		else if( nameStr.compare(std::string("kDirPhi")) == 0){ return FitterParameterType(kDirPhi); }
		else if( nameStr.compare(std::string("kEnergy")) == 0){ return FitterParameterType(kEnergy); }
		else if( nameStr.compare(std::string("kConversionDistance")) == 0){ return FitterParameterType(kConversionDistance); }
		else{
			std::cerr << "Error, " << name << " does not correspond to a known fitter parameter" << std::endl;
			assert(0);
		}
		return FitterParameterType(kUnknown);

	}

	static std::vector<Type> GetAllAllowedTypes(TrackType trackType)
	{

		std::vector<Type> types;
		types.push_back(kVtxX);
		types.push_back(kVtxY);
		types.push_back(kVtxZ);
		types.push_back(kVtxT);
		types.push_back(kDirTh);
		types.push_back(kDirPhi);
		types.push_back(kEnergy);
		if(trackType == TrackType::PhotonLike) { types.push_back(kConversionDistance); } 
		return types;
	}


};




#endif /* WCSIMTRACKPARAMETERENUMS_HH_ */
