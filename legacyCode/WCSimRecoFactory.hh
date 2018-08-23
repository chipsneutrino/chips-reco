#ifndef WCSIMRECOFACTORY_HH
#define WCSIMRECOFACTORY_HH

#include "TObject.h"

class WCSimReco;

class WCSimRecoFactory: public TObject {
	public:
		static WCSimRecoFactory* Instance();

		WCSimReco* MakeReco();
		WCSimReco* MakeReco(const char *name);

	private:
		WCSimRecoFactory();
		~WCSimRecoFactory();

		ClassDef(WCSimRecoFactory,0)
};

#endif
