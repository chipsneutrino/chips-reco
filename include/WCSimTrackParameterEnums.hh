/*
 * WCSimTrackParameterEnums.hh
 *
 *  Created on: 15 Jan 2015
 *      Author: ajperch
 */

#pragma once

#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include "TObject.h"

class TrackType : public TObject
{
public:
	enum Type
	{
		Unknown = 0,
		ElectronLike = 1,
		MuonLike = 2,
		PhotonLike = 3
	};
	Type fType;
	TrackType() : fType(Unknown)
	{
	}
	TrackType(Type t) : fType(t)
	{
	}
	TrackType(const TrackType &fpt) : fType(fpt.fType){};
	virtual ~TrackType()
	{
	}
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
		if (type == Unknown)
		{
			str = "Unknown";
		}
		if (type == ElectronLike)
		{
			str = "ElectronLike";
		}
		if (type == MuonLike)
		{
			str = "MuonLike";
		}
		if (type == PhotonLike)
		{
			str = "PhotonLike";
		}
		return str;
	};
	static TrackType FromName(const char *name)
	{
		std::string nameStr(name);

		if (nameStr.compare(std::string("Unknown")) == 0)
		{
			return TrackType(Unknown);
		}
		else if (nameStr.compare(std::string("ElectronLike")) == 0)
		{
			return TrackType(ElectronLike);
		}
		else if (nameStr.compare(std::string("MuonLike")) == 0)
		{
			return TrackType(MuonLike);
		}
		else if (nameStr.compare(std::string("PhotonLike")) == 0)
		{
			return TrackType(PhotonLike);
		}
		else
		{
			std::cerr << "Error, " << name << " does not correspond to a known fitter parameter" << std::endl;
			assert(0);
		}
		return TrackType(Unknown);
	}
	static std::vector<TrackType::Type> GetAllAllowedTypes(TrackType::Type trackType)
	{

		std::vector<TrackType::Type> types;
		types.push_back(TrackType::Unknown);
		types.push_back(TrackType::ElectronLike);
		types.push_back(TrackType::MuonLike);
		types.push_back(TrackType::PhotonLike);
		return types;
	}

	static Type GetTypeFromPDG(const int &pdg)
	{
		if (abs(pdg) == 11)
		{
			return ElectronLike;
		}
		else if (abs(pdg) == 13)
		{
			return MuonLike;
		}
		else if (abs(pdg) == 22)
		{
			return PhotonLike;
		}
		else
		{
			//std::cerr << "Error: unknown track type for PDG code " << pdg
			//		<< " (only 11, 13, 22 currently implemented)" << std::endl;
			return Unknown;
			//assert(0);
		}
	}

	static int GetPDGFromType(const Type &type)
	{
		if (type == ElectronLike)
		{
			return 11;
		}
		if (type == MuonLike)
		{
			return 13;
		}
		if (type == PhotonLike)
		{
			return 22;
		}
		else
		{
			//std::cerr << "Error: Did not recognise track type " << TrackType(type).AsString()
			//		<< " (only ElectronLike, MuonLike, PhotonLike currently implemented)" << std::endl;
			assert(0);
		}
		return 0;
	}
	ClassDef(TrackType, 1);

	friend std::ostream &operator<<(std::ostream &os, const TrackType &tt)
	{
		std::string name = TrackType(tt).AsString();
		return os << name;
	}
	friend std::ostream &operator<<(std::ostream &os, const Type &tt)
	{
		return os << TrackType(tt).AsString();
	}
};

class FitterParameterType : public TObject
{
public:
	enum Type
	{
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
	FitterParameterType() : fType(kUnknown){};
	FitterParameterType(Type t) : fType(t){};
	virtual ~FitterParameterType(){};
	FitterParameterType(const FitterParameterType &fpt) : fType(fpt.fType){};
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
		if (type == kVtxX)
		{
			str = "kVtxX";
		}
		if (type == kVtxY)
		{
			str = "kVtxY";
		}
		if (type == kVtxZ)
		{
			str = "kVtxZ";
		}
		if (type == kVtxT)
		{
			str = "kVtxT";
		}
		if (type == kDirTh)
		{
			str = "kDirTh";
		}
		if (type == kDirPhi)
		{
			str = "kDirPhi";
		}
		if (type == kEnergy)
		{
			str = "kEnergy";
		}
		if (type == kConversionDistance)
		{
			str = "kConversionDistance";
		}
		return str;
	};
	static FitterParameterType FromName(const char *name)
	{
		std::string nameStr(name);

		if (nameStr.compare(std::string("kVtxX")) == 0)
		{
			return FitterParameterType(kVtxX);
		}
		else if (nameStr.compare(std::string("kVtxY")) == 0)
		{
			return FitterParameterType(kVtxY);
		}
		else if (nameStr.compare(std::string("kVtxZ")) == 0)
		{
			return FitterParameterType(kVtxZ);
		}
		else if (nameStr.compare(std::string("kVtxT")) == 0)
		{
			return FitterParameterType(kVtxT);
		}
		else if (nameStr.compare(std::string("kDirTh")) == 0)
		{
			return FitterParameterType(kDirTh);
		}
		else if (nameStr.compare(std::string("kDirPhi")) == 0)
		{
			return FitterParameterType(kDirPhi);
		}
		else if (nameStr.compare(std::string("kEnergy")) == 0)
		{
			return FitterParameterType(kEnergy);
		}
		else if (nameStr.compare(std::string("kConversionDistance")) == 0)
		{
			return FitterParameterType(kConversionDistance);
		}
		else
		{
			std::cerr << "Error, " << name << " does not correspond to a known fitter parameter" << std::endl;
			assert(0);
		}
		return FitterParameterType(kUnknown);
	}

	static std::vector<FitterParameterType::Type> GetAllAllowedTypes(TrackType::Type trackType)
	{
		std::vector<FitterParameterType::Type> types;
		types.push_back(FitterParameterType::kVtxX);
		types.push_back(FitterParameterType::kVtxY);
		types.push_back(FitterParameterType::kVtxZ);
		types.push_back(FitterParameterType::kVtxT);
		types.push_back(FitterParameterType::kDirTh);
		types.push_back(FitterParameterType::kDirPhi);
		types.push_back(FitterParameterType::kEnergy);
		if (trackType == TrackType::PhotonLike || 1)
		{
			types.push_back(FitterParameterType::kConversionDistance);
		}
		return types;
	}
	ClassDef(FitterParameterType, 1);
};
