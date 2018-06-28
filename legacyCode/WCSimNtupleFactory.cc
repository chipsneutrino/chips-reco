#include "WCSimNtupleFactory.hh"
#include "WCSimNtuple.hh"

#include "WCSimRecoNtuple.hh"
#include "WCSimVertexNtuple.hh"
#include "WCSimVertexSeedNtuple.hh"

#include <string>
#include <iostream>
#include <cassert>

ClassImp (WCSimNtupleFactory)

static WCSimNtupleFactory* fgNtupleFactory = 0;

WCSimNtupleFactory* WCSimNtupleFactory::Instance() {
	if (!fgNtupleFactory) {
		fgNtupleFactory = new WCSimNtupleFactory();
	}

	// die if finder hasn't actually been created
	if (!fgNtupleFactory) {
		assert (fgNtupleFactory);
	}

	// can do re-setting here
	if (fgNtupleFactory) {

	}

	return fgNtupleFactory;
}

WCSimNtuple* WCSimNtupleFactory::MakeNtuple(const char *name) {
	std::cout << " *** WCSimNtupleFactory::MakeNtuple(...) *** " << std::endl;

	if (strcmp(name, "Reco") == 0) {
		if (fRecoNtuple == 0) {
			std::cout << "   Building Ntuple: " << name << std::endl;
			fRecoNtuple = new WCSimRecoNtuple();
		}
		std::cout << "   Returning Ntuple: " << name << std::endl;
		return fRecoNtuple;
	}

	else if (strcmp(name, "Vertex") == 0) {
		if (fVertexNtuple == 0) {
			std::cout << "   Building Ntuple: " << name << std::endl;
			fVertexNtuple = new WCSimVertexNtuple();
		}
		std::cout << "   Returning Ntuple: " << name << std::endl;
		return fVertexNtuple;
	}

	else if (strcmp(name, "VertexSeed") == 0) {
		if (fVertexSeedNtuple == 0) {
			std::cout << "   Building Ntuple: " << name << std::endl;
			fVertexSeedNtuple = new WCSimVertexSeedNtuple();
		}
		std::cout << "   Returning Ntuple: " << name << std::endl;
		return fVertexSeedNtuple;
	}

	else {
		std::cout << "   <warning> unknown ntuple: " << name << std::endl;
	}

	return 0;
}

WCSimNtupleFactory::WCSimNtupleFactory() {
	// list of possible displays
	fRecoNtuple = 0;
	fVertexNtuple = 0;
	fVertexSeedNtuple = 0;
}

WCSimNtupleFactory::~WCSimNtupleFactory() {
	if (fRecoNtuple)
		delete fRecoNtuple;
	if (fVertexNtuple)
		delete fVertexNtuple;
	if (fVertexSeedNtuple)
		delete fVertexSeedNtuple;
}
