#include "WCSimNtupleWriter.hh"
#include "WCSimNtupleFactory.hh"
#include "WCSimNtuple.hh"

#include "WCSimRecoFactory.hh"
#include "WCSimReco.hh"

#include "WCSimRecoEvent.hh"
#include "WCSimTrueEvent.hh"

#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"

#include <iostream>
#include <string>

ClassImp (WCSimNtupleWriter)

WCSimNtupleWriter::WCSimNtupleWriter() :
		fReco(0), fRecoEvent(0), fTrueEvent(0), fNtuple(0), fNtupleName(""), fNtpFileName("wcsim.ntuple.root") {

}

WCSimNtupleWriter::~WCSimNtupleWriter() {
	if (fReco)
		delete fReco;
}

void WCSimNtupleWriter::BuildNtuple() {
	std::cout << " *** WCSimNtupleWriter::BuildNtuple() *** " << std::endl;

	// Building Ntuple
	if (fNtuple) {
		fNtuple->Reset();
	} else {
		fNtuple = WCSimNtupleFactory::Instance()->MakeNtuple(fNtupleName.Data());
	}

	// Set File Name
	if (fNtuple) {
		fNtuple->SetFileName(fNtpFileName.Data());
	}

	// Building Reconstruction
	if (!fReco) {
		fReco = WCSimRecoFactory::Instance()->MakeReco();
	}

	return;
}

void WCSimNtupleWriter::Run(Int_t entries) {
	this->Run(0, entries);
}

void WCSimNtupleWriter::Run(Int_t i1, Int_t i2) {
	std::cout << " *** WCSimNtupleWriter::Run() *** " << std::endl;

	// check for events
	// ================
	if (WCSimInterface::TouchData() == 0) {
		return;
	}

	if (WCSimGeometry::TouchGeometry() == 0) {
		return;
	}

	// number of entries
	// =================
	Int_t ilast = i2;
	if (ilast < 0 || ilast > WCSimInterface::GetNumEvents()) {
		ilast = WCSimInterface::GetNumEvents();
	}

	Int_t ifirst = i1;
	if (ifirst > ilast) {
		return;
	}

	// build ntuple
	// ============
	this->BuildNtuple();

	// check ntuple
	// ============
	if (fNtuple == 0) {
		std::cout << "  <warning> failed to build ntuple! " << std::endl;
		return;
	}

	// loop over events
	// ================
	std::cout << "  running from " << ifirst << " to " << ilast << " (" << ilast - ifirst << " entries) " << std::endl;

	for (Int_t ievent = ifirst; ievent < ilast; ievent++) {
		std::cout << " ===== Next Event: [" << ievent << "/" << ilast << "] ===== " << std::endl;

		WCSimInterface::LoadEvent(ievent);

		// process this event
		// ==================
		if (WCSimInterface::TouchEvent()) {

			// create new event
			fRecoEvent = WCSimInterface::RecoEvent();
			fTrueEvent = WCSimInterface::TrueEvent();

			// run reconstruction
			fReco->Run(fRecoEvent);

			// write to ntuple
			fNtuple->Fill(fTrueEvent, fRecoEvent);

		}
	}

	// close file
	// ==========
	fNtuple->CloseFile();

	return;
}
