#include <iostream>
#include <vector>
#include <algorithm> 
#include <cassert>
#include <cmath>
#include <sstream>

#include "WCSimRecoSlicer.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimParameters.hh"
#include "WCSimCosmicSeed.hh"
#include "WCSimRecoClusteringUtil.hh"

#include <TFile.h>
#include <TDirectory.h>
#include <TH3D.h>

// A little sort function for events
bool SortTheSlices(WCSimRecoEvent* a, WCSimRecoEvent *b);

WCSimRecoSlicer::WCSimRecoSlicer() {
	fInputEvent = 0x0;
	fCosmicSeed = 0x0;
	this->UpdateParameters();
	fClusteringUtil = 0x0;
}

WCSimRecoSlicer::~WCSimRecoSlicer() {
	// Created events must be deleted by the user at a later point

	if (fClusteringUtil != 0x0) {
		delete fClusteringUtil;
	}

}

void WCSimRecoSlicer::UpdateParameters() {
	WCSimParameters* pars = WCSimParameters::Instance();
	fDistanceCut = pars->GetSlicerClusterDistance();
	fMinSliceSize = pars->GetSlicerMinSize();
	fChargeCut = pars->GetSlicerChargeCut();
	fTimeCut = pars->GetSlicerTimeCut();
	fIterateSlicing = pars->GetIterateSlicing();
}

void WCSimRecoSlicer::ResetSlices() {
	// Only want to reset the slice vectors.
	for (unsigned int i = 0; i < fSlicedDigits.size(); ++i) {
		fSlicedDigits[i].clear();
	}
	fSlicedDigits.clear();
}

void WCSimRecoSlicer::SetInputEvent(WCSimRecoEvent* evt) {
	fInputEvent = evt;

}

void WCSimRecoSlicer::Run() {

	if (fInputEvent == 0x0) {
		std::cerr << "WCSimRecoSlicer: Can't slice an empty event!" << std::endl;
		assert(0);
	}

	if (!fInputEvent->IsFilterDone()) {
		std::cerr << "WCSimRecoSlicer: Need to apply filter to digits before slicing." << std::endl;
		assert(0);
	}

	// Make the slices
	this->SliceTheEvent();

	// Do we need to increase the charge cut and try again?
	if (fIterateSlicing && (fSlicedDigits.size() < fNTracks)) {
		std::cout << "Now iterating the slicing" << std::endl;
		while ((fSlicedDigits.size() < fNTracks) && (fChargeCut < 9.01)) {
			++fChargeCut;
			std::cout << " - Charge cut updated to " << fChargeCut << " pe." << std::endl;
			this->ResetSlices();
			this->SliceTheEvent();
		}
		// Do we have the right number of slices?
		if (fSlicedDigits.size() < fNTracks) {
			// It wasn't possible. Go back to default values.
			std::cout << "No addition slicing possible using iteration. Returing to default." << std::endl;
			this->UpdateParameters();
			this->ResetSlices();
			this->SliceTheEvent();
		}
	}

	// Now we need to turn our slices back into WCSimRecoEvents
	this->BuildEvents();

}

void WCSimRecoSlicer::SliceTheEvent() {

	std::vector<WCSimRecoDigit*> digitVector;
	// Make a charge sorted vector of those digits above the charge threshold
	int shadowHits = 0;
	/*
	 TDirectory* tmpd = gDirectory;
	 TFile *shadowFile = new TFile("shadowFile.root","recreate");
	 shadowFile->cd();
	 TH3F* hShadowHits = new TH3F("hShadow",";x;y;z",100,-13000,13000,100,-13000,13000,100,-11000,11000);
	 TH3F* hNormalHits = new TH3F("hNormal",";x;y;z",100,-13000,13000,100,-13000,13000,100,-11000,11000);
	 */
	for (unsigned int v = 0; v < fInputEvent->GetFilterDigitList()->size(); ++v) {
		WCSimRecoDigit *dig = (*(fInputEvent->GetFilterDigitList()))[v];
		// If we are using the cosmic seed, make sure we haven't shadowed this hit
		if (fCosmicSeed != 0x0) {
			if (fCosmicSeed->IsHitShadowed(dig->GetX(), dig->GetY(), dig->GetZ())) {
//        hShadowHits->Fill(dig->GetX(),dig->GetY(),dig->GetZ());
				++shadowHits;
				continue;
			}
//      hNormalHits->Fill(dig->GetX(),dig->GetY(),dig->GetZ());
		}
		digitVector.push_back(dig);
	}

	if (fCosmicSeed != 0x0) {
		std::cout << "Shadowed hits above threshold = " << shadowHits << std::endl;
	}

	if (fClusteringUtil != 0x0) {
		delete fClusteringUtil;
	}
	fClusteringUtil = new WCSimRecoClusteringUtil(digitVector, fChargeCut, fTimeCut, fDistanceCut, fMinSliceSize);

	fClusteringUtil->PerformClustering();
	fSlicedDigits = fClusteringUtil->GetFullSlicedDigits();
	/*
	 // Tidy up the files
	 hShadowHits->Write();
	 hNormalHits->Write();
	 shadowFile->Close();
	 delete hShadowHits;
	 delete hNormalHits;
	 delete shadowFile;
	 gDirectory = tmpd;
	 */
}

void WCSimRecoSlicer::BuildEvents() {

//  std::vector<std::pair<unsigned int, WCSimRecoEvent*> > sortVec;
	for (unsigned int v = 0; v < fSlicedDigits.size(); ++v) {
		WCSimRecoEvent* newEvt = new WCSimRecoEvent();
		// Initialise the event from the input event
		newEvt->SetHeader(fInputEvent->GetRun(), fInputEvent->GetEvent(), fInputEvent->GetTrigger());
		newEvt->SetVertex(fInputEvent->GetVtxX(), fInputEvent->GetVtxY(), fInputEvent->GetVtxZ(),
				fInputEvent->GetVtxTime());
		newEvt->SetDirection(fInputEvent->GetVtxX(), fInputEvent->GetVtxY(), fInputEvent->GetVtxZ());
		newEvt->SetConeAngle(fInputEvent->GetConeAngle());
		newEvt->SetTrackLength(fInputEvent->GetTrackLength());
		newEvt->SetVtxFOM(fInputEvent->GetVtxFOM(), fInputEvent->GetVtxIterations(), fInputEvent->GetVtxPass());
		newEvt->SetVtxStatus(fInputEvent->GetVtxStatus());

		if (fInputEvent->IsFilterDone()) {
			newEvt->SetFilterDone();
		}
		if (fInputEvent->IsVertexFinderDone()) {
			newEvt->SetVertexFinderDone();
		}
		if (fInputEvent->IsRingFinderDone()) {
			newEvt->SetRingFinderDone();
		}

		// Add the digits
		for (unsigned int d = 0; d < fSlicedDigits[v].size(); ++d) {
			newEvt->AddDigit(fSlicedDigits[v][d]);
			newEvt->AddFilterDigit(fSlicedDigits[v][d]);
		}
		fSlicedEvents.push_back(newEvt);
		std::cout << "Created event with " << newEvt->GetNFilterDigits() << " hits." << std::endl;
	}

	if (fSlicedEvents.size() == 0) {
		std::cout << "Found no slices: taking the whole event instead" << std::endl;
		WCSimRecoEvent* newEvt = new WCSimRecoEvent();
		// Initialise the event from the input event
		newEvt->SetHeader(fInputEvent->GetRun(), fInputEvent->GetEvent(), fInputEvent->GetTrigger());
		newEvt->SetVertex(fInputEvent->GetVtxX(), fInputEvent->GetVtxY(), fInputEvent->GetVtxZ(),
				fInputEvent->GetVtxTime());
		newEvt->SetDirection(fInputEvent->GetVtxX(), fInputEvent->GetVtxY(), fInputEvent->GetVtxZ());
		newEvt->SetConeAngle(fInputEvent->GetConeAngle());
		newEvt->SetTrackLength(fInputEvent->GetTrackLength());
		newEvt->SetVtxFOM(fInputEvent->GetVtxFOM(), fInputEvent->GetVtxIterations(), fInputEvent->GetVtxPass());
		newEvt->SetVtxStatus(fInputEvent->GetVtxStatus());

		if (fInputEvent->IsFilterDone()) {
			newEvt->SetFilterDone();
		}
		if (fInputEvent->IsVertexFinderDone()) {
			newEvt->SetVertexFinderDone();
		}
		if (fInputEvent->IsRingFinderDone()) {
			newEvt->SetRingFinderDone();
		}

		// Add the digits
		for (int d = 0; d < fInputEvent->GetNDigits(); ++d) {
			newEvt->AddDigit(fInputEvent->GetDigit(d));
		}
		for (int fd = 0; fd < fInputEvent->GetNFilterDigits(); ++fd) {
			newEvt->AddDigit(fInputEvent->GetFilterDigit(fd));
		}
		fSlicedEvents.push_back(newEvt);
	}

	// Now we want to sort the vector by the number of digits.
	std::sort(fSlicedEvents.begin(), fSlicedEvents.end(), SortTheSlices);

}

// A quick sort function for the reco events.
bool SortTheSlices(WCSimRecoEvent* a, WCSimRecoEvent *b) {
	unsigned int aSize = a->GetFilterDigitList()->size();
	unsigned int bSize = b->GetFilterDigitList()->size();
	return aSize > bSize;
}

