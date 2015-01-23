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
    std::cout << "Getting type from name: " << name << std::endl;
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




#endif /* WCSIMTRACKPARAMETERENUMS_HH_ */
