/*
 * WCSimTransmissionFunctionLookup.hh
 *
 *  Created on: 4 Nov 2015
 *      Author: ajperch
 */

#pragma once

#include <vector>
#include <iostream>
//
class WCSimTransmissionFunctionLookupTable
{
public:
	// Stores a lookup table for norm * exp(nu * x) with x in the range xMin <= x < xMax
	WCSimTransmissionFunctionLookupTable(double xMin = 0, double xMax = 0, double norm = 0, double nu = 0)
	{
		fXMin = xMin;
		fXMax = xMax;
		fNorm = norm;
		fNu = nu;
		MakeLUT();
	}
	float GetExp(double x);
	inline double const GetXMin()
	{
		return fXMin;
	}
	inline double const GetXMax()
	{
		return fXMax;
	}
	inline double const GetNorm()
	{
		return fNorm;
	}
	inline double const GetNu()
	{
		return fNu;
	}

private:
	void MakeLUT();
	double fXMin;
	double fXMax;
	double fNorm;
	double fNu;
	float fLUT[2048];
};

class WCSimTransmissionFunctionLookup
{
public:
	virtual ~WCSimTransmissionFunctionLookup();

	static WCSimTransmissionFunctionLookup *Instance();
	static float GetExp(double xMin, double xMax, double norm, double nu, double x);
	float GetExponential(double xMin, double xMax, double norm, double nu, double x);

private:
	WCSimTransmissionFunctionLookup();

	int GetTableNum(double xMin, double xMax, double norm, double nu);
	std::vector<WCSimTransmissionFunctionLookupTable> fTables;
};
