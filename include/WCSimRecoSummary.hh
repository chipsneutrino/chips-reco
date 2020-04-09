/*
 * WCSimRecoSummary.hh
 *
 *  Created on: 16 Feb 2015
 *      Author: ajperch
 */

#pragma once

#include "TVector3.h"
#include "TObject.h"
#include <vector>

class WCSimRecoSummary : public TObject
{
public:
	// Constructors and Destructor
	WCSimRecoSummary();
	WCSimRecoSummary(const WCSimRecoSummary &ts);
	~WCSimRecoSummary();

	// Reset the parameters to their dummy values
	void ResetValues();

	// Get and set the vertex information
	void AddVertex(TVector3 vtx, double time);
	void AddVertex(double x, double y, double z, double time);

	TVector3 GetVertex(unsigned int p) const;
	std::vector<TVector3> GetVertices() const;
	double GetVertexX(unsigned int p) const;
	double GetVertexY(unsigned int p) const;
	double GetVertexZ(unsigned int p) const;
	double GetVertexT(unsigned int p) const;
	std::vector<double> GetVerticesT() const;

	// Primary particle functions
	void AddPrimary(int pdg, double en, TVector3 dir);
	void AddPrimary(int pdg, double en, double dx, double dy, double dz);
	int GetPrimaryPDG(unsigned int p) const;	   // Index starts at 0
	double GetPrimaryEnergy(unsigned int p) const; // Index starts at 0
	TVector3 GetPrimaryDir(unsigned int p) const;  // Index starts at 0
	std::vector<int> GetPrimaryPDGs() const;
	std::vector<double> GetPrimaryEnergies() const;
	std::vector<TVector3> GetPrimaryDirs() const;
	unsigned int GetNPrimaries() const;

	bool HasCommonVertex() const;
	// Best-fit log likelihoods
	void SetTimeMinus2LnL(double minus2LnL);
	void SetChargeMinus2LnL(double minus2LnL);
	double GetTimeMinus2LnL() const;
	double GetChargeMinus2LnL() const;
	double GetMinus2LnL() const;

	int GetEventNumber() const;
	void SetEventNumber(int event);

private:
	int fEventNumber; // The event number in the WCSim file being reconstructed

	// Likelihoods
	double fTimeMinus2LnL;
	double fChargeMinus2LnL;

	// Vertex position
	std::vector<TVector3> fVertices;
	std::vector<double> fVerticesT;

	// Also store information about the primary particles (only those escaping the nucleus)
	std::vector<int> fPrimaryPDGs;
	std::vector<double> fPrimaryEnergies;
	std::vector<TVector3> fPrimaryDirs;

	ClassDef(WCSimRecoSummary, 2);
};
