#ifndef WCSIMRECOAB_HH
#define WCSIMRECOAB_HH

#include "WCSimReco.hh"

#include <vector>

class WCSimRecoDigit;
class WCSimRecoVertex;
class WCSimRecoRing;

class WCSimRecoAB: public WCSimReco {
	public:
		WCSimRecoAB();
		~WCSimRecoAB();

		// Reconstruction Methods
		// ======================
		void Run();
		void Run(WCSimRecoEvent* evt);
		void RunFilter(WCSimRecoEvent* evt);
		void RunRecoVertex(WCSimRecoEvent* evt);
		void RunRecoRings(WCSimRecoEvent* evt);

	private:
		ClassDef(WCSimRecoAB,0)
};

#endif
