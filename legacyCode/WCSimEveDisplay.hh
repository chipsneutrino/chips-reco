#ifndef WCSIMEVEDISPLAY_HH
#define WCSIMEVEDISPLAY_HH

#include "WCSimDisplay.hh"
#include "WCSimTruthSummary.hh"
#include "TEvePointSet.h"

class WCSimEveDisplay: public WCSimDisplay {

		enum {
			fNumBins = 10
		};

	public:
		WCSimEveDisplay();
		~WCSimEveDisplay();

		void DrawDisplay(WCSimRecoEvent* recoevent);
		void DrawCleanDisplay(WCSimRecoEvent* recoevent);

		// nothing else implemented
		void DrawRecoEvent(WCSimRecoEvent*);
		void DrawTrueEvent(WCSimRootEvent*);
		void ResetDisplay();
		void PrintDisplay();

	private:
		int fEventNum;
		void Initialize();
		void BuildGeometry();
		void BuildEvePointSets();
		void ClearEvePointSets();

		void BuildTrueTracks(WCSimTruthSummary * mySummary);

		std::vector<std::vector<TEvePointSet*> > fEvePointSets;

		ClassDef(WCSimEveDisplay,0)

};

#endif

