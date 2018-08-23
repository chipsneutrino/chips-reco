#ifndef WCSIMRECO_HH
#define WCSIMRECO_HH

#include "TObject.h"

#include "WCSimRecoEvent.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoVertex.hh"
#include "WCSimRecoRing.hh"

class WCSimReco: public TObject {
	public:
		// Reconstruction Methods
		// ======================
		virtual void Run() = 0;
		virtual void Run(WCSimRecoEvent*) = 0;
		virtual void RunFilter(WCSimRecoEvent*) = 0;
		virtual void RunRecoVertex(WCSimRecoEvent*) = 0;
		virtual void RunRecoRings(WCSimRecoEvent*) = 0;

	private:
		ClassDef(WCSimReco,0)

};

#endif
