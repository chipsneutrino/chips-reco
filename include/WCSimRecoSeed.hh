#pragma once

#include "TVector3.h"
#include "TObject.h"
#include "WCSimRecoEvent.hh"

#include <vector>

class WCSimRecoDigit;
class WCSimRecoVertex;
class WCSimRecoRing;
class WCSimCosmicSeed;

class WCSimRecoSeed : public TObject
{
public:
	WCSimRecoSeed();
	~WCSimRecoSeed();

	// Reconstruction Methods
	// ======================
	void Run();
	void Run(WCSimRecoEvent *evt);
	std::vector<WCSimRecoEvent *> RunSeed(WCSimRecoEvent *evt);
	void RunFilter(WCSimRecoEvent *evt);
	std::vector<WCSimRecoEvent *> RunSlicer(WCSimRecoEvent *evt);
	void RunRecoVertex(WCSimRecoEvent *evt);
	void RunRecoRings(WCSimRecoEvent *evt);

	void SetNumberOfTracks(unsigned int val)
	{
		fNTracks = val;
	}

	unsigned int GetNumberOfTracks()
	{
		return fNTracks;
	}

	TVector3 GetDirBeforeRings() const
	{
		return fDirBeforeRings;
	}

	void SetCosmicFit(bool val)
	{
		fCosmicFit = val;
	}

	bool GetCosmicFit()
	{
		return fCosmicFit;
	}

	WCSimCosmicSeed *GetCosmicSeed()
	{
		return fCosmicSeed;
	}

private:
	TVector3 fDirBeforeRings;
	unsigned int fNTracks;

	WCSimCosmicSeed *fCosmicSeed;

	bool fCosmicFit;

	ClassDef(WCSimRecoSeed, 0)
};
