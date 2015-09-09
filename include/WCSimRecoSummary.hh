/*
 * WCSimRecoSummary.hh
 *
 *  Created on: 16 Feb 2015
 *      Author: ajperch
 */

#ifndef WCSIMRECOSUMMARY_HH_
#define WCSIMRECOSUMMARY_HH_

#include "TVector3.h"
#include "TObject.h"
#include <vector>


class WCSimRecoSummary : public TObject {
public:
	  // Constructors and Destructor
	  WCSimRecoSummary();
	  WCSimRecoSummary(const WCSimRecoSummary &ts);
	  ~WCSimRecoSummary();

	  // Reset the parameters to their dummy values
	  void ResetValues();

	  // Get and set the vertex information
	  TVector3 GetVertex() const;
	  void SetVertex(TVector3 vtx);
	  void SetVertex(double x, double y, double z);
	  void SetVertex(double x, double y, double z, double t);
    void SetVertexT(double t);

	  double GetVertexX() const;
	  double GetVertexY() const;
	  double GetVertexZ() const;
    double GetVertexT() const;

	  // Primary particle functions
	  void AddPrimary(int pdg, double en, TVector3 dir);
	  void AddPrimary(int pdg, double en, double dx, double dy, double dz);
	  int GetPrimaryPDG(unsigned int p) const; // Index starts at 0
	  double GetPrimaryEnergy(unsigned int p) const; // Index starts at 0
	  TVector3 GetPrimaryDir(unsigned int p) const; // Index starts at 0
	  std::vector<int> GetPrimaryPDGs() const;
	  std::vector<double> GetPrimaryEnergies() const;
	  std::vector<TVector3> GetPrimaryDirs() const;
	  unsigned int GetNPrimaries() const;


	private:

	  // Vertex position
	  TVector3 fVertex;
    double fVertexT;

	  // Also store information about the primary particles (only those escaping the nucleus)
	  std::vector<int> fPrimaryPDGs;
	  std::vector<double> fPrimaryEnergies;
	  std::vector<TVector3> fPrimaryDirs;

	  ClassDef(WCSimRecoSummary,1);
};

#endif /* WCSIMRECOSUMMARY_HH_ */
