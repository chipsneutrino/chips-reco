/*
 * WCSimTransmissionFunctionLookup.cc
 *
 *  Created on: 4 Nov 2015
 *      Author: ajperch
 */

#include "WCSimTransmissionFunctionLookup.hh"
#include <cassert>
#include <cmath>

WCSimTransmissionFunctionLookup * fgWCSimTransmissionFunctionLookup = 0x0;


////////////////////////////////////////////////////////////////
// Struct to hold the tables
////////////////////////////////////////////////////////////////
void WCSimTransmissionFunctionLookupTable::MakeLUT()
{
	double step = (fXMax - fXMin)/2048.;
	for(int i = 0; i < 2048; ++i)
	{
		double x = fXMin + i*step;
		fLUT[i] = fNorm * exp(fNu * x);
	}
	return;
}

float WCSimTransmissionFunctionLookupTable::GetExp(double x)
{
	float width = (fXMax - fXMin) / 2048.;
	int index = (int)((x - fXMin) / width);

	if( !(0 <= index && index < 2048))
	{
		std::cerr << "Warning: you're looking up the transmission function for a distance "
				  << x << std::endl
				  << "         But the table only goes from " << fXMin << " to " << fXMax << std::endl
				  << "         I'll calculate the slow way" << std::endl;
		return exp(-x);
	}
	if(index == 0) { return fLUT[0];}
	if(index == 2047) { return fLUT[2047]; }

	float low = fXMin + width * (int)(x/width);
	float high = low + width;
	return ( (fLUT[index] * (high - x) + fLUT[index+1] * (x - low)) / width );
}

////////////////////////////////////////////////////////////////
// Class to access the tables
////////////////////////////////////////////////////////////////

float WCSimTransmissionFunctionLookup::GetExponential(double xMin, double xMax, double norm, double nu, double x)
{
	int whichTable = this->GetTableNum(xMin, xMax, norm, nu);
	if( whichTable < 0 )
	{
		fTables.push_back(WCSimTransmissionFunctionLookupTable(xMin, xMax, norm, nu));
    whichTable = fTables.size()-1;
	}
	assert(fTables.size() < 4);
	return fTables[whichTable].GetExp(x);
}

float WCSimTransmissionFunctionLookup::GetExp(double xMin, double xMax, double norm, double nu, double x)
{
	return WCSimTransmissionFunctionLookup::Instance()->GetExponential(xMin, xMax, norm, nu, x);
}

WCSimTransmissionFunctionLookup::WCSimTransmissionFunctionLookup()
{
}

WCSimTransmissionFunctionLookup *WCSimTransmissionFunctionLookup::Instance()
{
	if(fgWCSimTransmissionFunctionLookup == 0x0)
	{
		fgWCSimTransmissionFunctionLookup = new WCSimTransmissionFunctionLookup();
	}
	return fgWCSimTransmissionFunctionLookup;
}

WCSimTransmissionFunctionLookup::~WCSimTransmissionFunctionLookup()
{
	// TODO Auto-generated destructor stub
}

int WCSimTransmissionFunctionLookup::GetTableNum(double xMin, double xMax, double norm, double nu)
{
	int tableNum = -1;
	unsigned int nTables = fTables.size();
	for(unsigned int i = 0; i < nTables; ++i)
	{
		if(fTables[i].GetNorm() == norm && fTables[i].GetNu() == nu && fTables[i].GetXMin() == xMin && fTables[i].GetXMax() == xMax)
		{
			tableNum = i;
			break;
		}
	}
	return tableNum;
}



