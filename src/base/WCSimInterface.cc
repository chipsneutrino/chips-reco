#include "WCSimInterface.hh"
#include "WCSimGeometry.hh"
#include "WCSimParameters.hh"

#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackFactory.hh"

#include "WCSimTrueEvent.hh"
#include "WCSimTrueTrack.hh"

#include "WCSimRecoEvent.hh"
#include "WCSimRecoDigit.hh"

#include <TChainElement.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TIterator.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParticlePDG.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <set>

ClassImp(WCSimInterface)

static WCSimInterface *fgInterface = 0;

Bool_t WCSimInterface::TouchData()
{
	if (WCSimInterface::Instance()->GetEntries() > 0)
	{
		return true;
	}
	else
	{
		std::cout << " *** WCSimInterface::TouchData() *** " << std::endl;
		std::cout << "  <error> need to input data... " << std::endl
				  << "    Call: WCSimInterface::LoadData(...) "
				  << std::endl;
		return false;
	}
}

WCSimInterface *WCSimInterface::Instance()
{
	if (!fgInterface)
	{
		fgInterface = new WCSimInterface();
	}

	return fgInterface;
}

void WCSimInterface::LoadData(const char *file)
{
	return WCSimInterface::Instance()->AddFile(file);
}

void WCSimInterface::LoadEvent(Int_t ievent)
{
	return WCSimInterface::Instance()->BuildEvent(ievent);
}

Int_t WCSimInterface::GetNumEvents()
{
	return WCSimInterface::Instance()->GetEntries();
}

Bool_t WCSimInterface::TouchEvent()
{
	return WCSimInterface::Instance()->CheckEvent();
}

Int_t WCSimInterface::GetRunNumber()
{
	return WCSimInterface::RecoEvent()->GetRun();
}

Int_t WCSimInterface::GetEventNumber()
{
	return WCSimInterface::RecoEvent()->GetEvent();
}

Int_t WCSimInterface::GetTriggerNumber()
{
	return WCSimInterface::RecoEvent()->GetTrigger();
}

void WCSimInterface::SetEnergyThreshold(Double_t input_mev)
{
	return WCSimInterface::Instance()->SetTrueEnergyThreshold(input_mev);
}

void WCSimInterface::SetRangeThreshold(Double_t input_cm)
{
	return WCSimInterface::Instance()->SetTrueRangeThreshold(input_cm);
}

void WCSimInterface::Reset()
{
	return WCSimInterface::Instance()->ResetForNewSample();
}

WCSimInterface::WCSimInterface()
{
	fTrueEvent = new WCSimTrueEvent();
	fRecoEvent = new WCSimRecoEvent();

	fDigitList = new std::vector<WCSimRecoDigit *>;
	fVetoDigitList = new std::vector<WCSimRecoDigit *>;
	fTrackList = new std::vector<WCSimTrueTrack *>;
	fTrueLikelihoodTracks = new std::vector<WCSimLikelihoodTrackBase *>;

	fTrigger = 0;
	fEvent = 0;
	fGeometry = 0;

	fEnergyThreshold = 25.0; // 25 MeV
	fRangeThreshold = 15.0;	 // 15 cm

	// create event chain
	fChain = new TChain("wcsimT", "chain");
	fChain->SetBranchAddress("wcsimrootevent", &fEvent);

	// create geometry chain
	fChainGeom = new TChain("wcsimGeoT", "chainGeom");
	fChainGeom->SetBranchAddress("wcsimrootgeom", &fGeometry);

	this->ResetForNewSample();
}

WCSimInterface::~WCSimInterface()
{
	this->ResetTrueEvent();
	this->ResetRecoEvent();

	delete fTrueEvent;
	delete fRecoEvent;

	delete fDigitList;
	delete fTrackList;
}

void WCSimInterface::AddFile(const char *file)
{
	std::cout << " *** WCSimInterface::LoadData(...) *** " << std::endl;
	std::cout << "  adding: " << file << std::endl;

	fChain->Add(file);

	std::cout << "   ... total entries=" << fChain->GetEntries() << std::endl;

	if (!fGeometry && fChain->GetEntries() > 0)
	{
		std::cout << " *** WCSimInterface::LoadGeometry(...) *** " << std::endl;
		fChainGeom->Add(file);
		if (fChainGeom->GetEntries() > 0)
		{
			fChainGeom->GetEntry(0);
			if (fGeometry)
				WCSimGeometry::BuildGeometry(fGeometry);
		}
	}

	return;
}

void WCSimInterface::ResetForNewSample()
{
	if (fChain->GetEntries() > 0)
	{
		std::cout << " *** WCSimInterface::Reset() *** " << std::endl;
	}

	// event chain
	if (fChain)
	{
		fChain->Reset();
		fChain->SetBranchAddress("wcsimrootevent", &fEvent);
		if (fEvent)
			delete fEvent;
		fEvent = 0;
	}

	// geometry chain
	if (fChainGeom)
	{
		fChainGeom->Reset();
		fChainGeom->SetBranchAddress("wcsimrootgeom", &fGeometry);
		if (fGeometry)
			delete fGeometry;
		fGeometry = 0;
	}

	// reset true event
	this->ResetTrueEvent();

	// reset true tracks
	this->ResetTrueLikelihoodTracks();

	//  reset reco event
	this->ResetRecoEvent();

	// reset geometry
	WCSimGeometry::Reset();

	return;
}

void WCSimInterface::ResetTrueEvent()
{
	// Reset
	// =====
	fTrueEvent->Reset();

	// delete list of tracks
	// =====================
	for (UInt_t i = 0; i < fTrackList->size(); i++)
	{
		delete (WCSimTrueTrack *)(fTrackList->at(i));
	}

	// clear list of tracks
	// ====================
	fTrackList->clear();

	return;
}

void WCSimInterface::ResetTrueLikelihoodTracks()
{
	// delete list of tracks
	// =====================
	if (fTrueLikelihoodTracks != 0x0)
	{
		for (UInt_t i = 0; i < fTrueLikelihoodTracks->size(); i++)
		{
			delete (WCSimLikelihoodTrackBase *)(fTrueLikelihoodTracks->at(i));
		}
		fTrueLikelihoodTracks->clear();
		delete fTrueLikelihoodTracks;
		fTrueLikelihoodTracks = 0x0;
	}
	return;
}

void WCSimInterface::BuildTrueLikelihoodTracks()
{
	std::cout << " *** WCSimInterface::BuildTrueLikelihoodTracks() *** " << std::endl;

	TDatabasePDG myDatabase;

	WCSimTruthSummary sum = fEvent->GetTruthSummary();

	fTrueLikelihoodTracks = new std::vector<WCSimLikelihoodTrackBase *>;
	if (sum.IsParticleGunEvent())
	{
		TParticlePDG *beamParticle = myDatabase.GetParticle(sum.GetBeamPDG());

		// If we have a pi-zero gun, make sure we treat it properly.
		if (sum.GetBeamPDG() == 111)
		{
			std::vector<WCSimLikelihoodTrackBase *> tracks = this->GetPi0PhotonTracks(sum, sum.GetBeamEnergy(),
																					  fTrigger->GetTracks());
			for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack)
			{
				fTrueLikelihoodTracks->push_back(tracks.at(iTrack));
			}
		}
		else
		{
			double mm_to_cm = 0.1;
			WCSimLikelihoodTrackBase *track = WCSimLikelihoodTrackFactory::MakeTrack(
				TrackType::GetTypeFromPDG(sum.GetBeamPDG()), sum.GetVertexX() * mm_to_cm,
				sum.GetVertexY() * mm_to_cm, sum.GetVertexZ() * mm_to_cm, sum.GetVertexT(),
				sum.GetBeamDir().Theta(), sum.GetBeamDir().Phi(),
				sum.GetBeamEnergy() - beamParticle->Mass() * 1000); // Mass comes in GeV but we work in MeV
			fTrueLikelihoodTracks->push_back(track);
		}
	}
	else
	{
		double mm_to_cm = 0.1;
		for (unsigned int i = 0; i < sum.GetNPrimaries(); ++i)
		{
			// We had a final state pi0, so get the photons it decays to
			if (sum.GetPrimaryPDG(i) == 111)
			{
				std::vector<WCSimLikelihoodTrackBase *> tracks = this->GetPi0PhotonTracks(sum, sum.GetPrimaryEnergy(i),
																						  fTrigger->GetTracks());
				for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack)
				{
					fTrueLikelihoodTracks->push_back(tracks.at(iTrack));
				}
			}
			else if (TrackType::GetTypeFromPDG(sum.GetPrimaryPDG(i)) != TrackType::Unknown)
			{
				TParticlePDG *primaryParticle = myDatabase.GetParticle(sum.GetPrimaryPDG(i));
				if (sum.GetPrimaryEnergy(i) < 20)
				{
					continue;
				}
				WCSimLikelihoodTrackBase *track = WCSimLikelihoodTrackFactory::MakeTrack(
					TrackType::GetTypeFromPDG(sum.GetPrimaryPDG(i)), sum.GetVertexX() * mm_to_cm,
					sum.GetVertexY() * mm_to_cm, sum.GetVertexZ() * mm_to_cm, sum.GetVertexT(),
					sum.GetPrimaryDir(i).Theta(), sum.GetPrimaryDir(i).Phi(),
					sum.GetPrimaryEnergy(i) - primaryParticle->Mass() * 1000 // Mass comes in GeV but we want MeV
				);
				fTrueLikelihoodTracks->push_back(track);
			}
		}
		// Leigh: If this is an overlay event, add the overlays too.
		if (sum.IsOverlayEvent())
		{
			for (unsigned int i = 0; i < sum.GetNOverlays(); ++i)
			{
				// Make the track
				TParticlePDG *overlayParticle = myDatabase.GetParticle(sum.GetOverlayPDG(i));
				WCSimLikelihoodTrackBase *track = WCSimLikelihoodTrackFactory::MakeTrack(
					TrackType::GetTypeFromPDG(sum.GetOverlayPDG(i)), sum.GetOverlayVertexX() * mm_to_cm,
					sum.GetOverlayVertexY() * mm_to_cm, sum.GetOverlayVertexZ() * mm_to_cm, sum.GetOverlayVertexT(),
					sum.GetOverlayDir(i).Theta(), sum.GetOverlayDir(i).Phi(),
					sum.GetOverlayEnergy(i) - overlayParticle->Mass() * 1000 // Mass comes in GeV but we want MeV
				);
				fTrueLikelihoodTracks->push_back(track);
			}
		}
	}
	return;
}

void WCSimInterface::ResetRecoEvent()
{
	// Reset
	// =====
	fRecoEvent->Reset();

	// delete list of digits
	// =====================
	for (UInt_t i = 0; i < fDigitList->size(); i++)
	{
		delete (WCSimRecoDigit *)(fDigitList->at(i));
	}
	for (UInt_t i = 0; i < fVetoDigitList->size(); i++)
	{
		delete (WCSimRecoDigit *)(fVetoDigitList->at(i));
	}

	// clear list of digits
	// ====================
	fDigitList->clear();
	fVetoDigitList->clear();

	return;
}

Int_t WCSimInterface::GetEntries()
{
	return fChain->GetEntries();
}

Bool_t WCSimInterface::CheckEvent()
{
	if (fTrigger)
		return true;
	else
		return false;
}

void WCSimInterface::BuildEvent(Int_t ievent)
{
	std::cout << " *** WCSimInterface::BuildEvent(" << ievent << ") *** " << std::endl;

	WCSimRootTrigger *myTrigger = (WCSimRootTrigger *)(this->GetWCSimTrigger(ievent));

	this->BuildEvent(myTrigger);

	return;
}

WCSimRootEvent *WCSimInterface::GetWCSimEvent(Int_t ievent)
{
	if (fEvent)
	{
		delete fEvent;
		fEvent = 0;
		fTrigger = 0;
	} // Added the braces - was resetting everything to zero without them - AJP 18/04/13

	if (ievent >= 0 && ievent < fChain->GetEntries())
	{ // This was || not && - don't think it should be - AJP 18/04/13
		fChain->GetEntry(ievent);
	}
	else
	{
		std::cout << "  <warning>: this event doesn't exist! " << std::endl;
	}

	return fEvent;
}

WCSimRootTrigger *WCSimInterface::GetWCSimTrigger(Int_t ievent)
{
	WCSimRootEvent *myEvent = (WCSimRootEvent *)(this->GetWCSimEvent(ievent));

	fTrigger = WCSimInterface::FilterTrigger(myEvent);

	return fTrigger;
}

WCSimRootTrigger *WCSimInterface::FilterTrigger(WCSimRootEvent *myEvent)
{
	// Sanity Check
	// ============
	if (myEvent == 0)
		return 0;

	// make a WCSim Trigger
	// ====================
	WCSimRootTrigger *myTrigger = 0;

	// for now, pick the largest trigger
	// =================================

	// Change this so we require there to be atleast 1 CherenkovDigiHit for there to be a trigger...
	/*
	Int_t triggerDigits = -1;
	Int_t triggerNumber = -1;
	for (Int_t nTrigger = 0; nTrigger < myEvent->GetNumberOfEvents(); nTrigger++) {
		WCSimRootTrigger* tempTrigger = (WCSimRootTrigger*) (myEvent->GetTrigger(nTrigger));
		Int_t nDigits = 1 + tempTrigger->GetCherenkovDigiHits()->GetLast();

		if (nDigits > triggerDigits) {
			triggerDigits = nDigits;
			triggerNumber = nTrigger;
		}
	}
	*/

	Int_t triggerDigits = 0;
	Int_t triggerNumber = -1;
	for (Int_t nTrigger = 0; nTrigger < myEvent->GetNumberOfEvents(); nTrigger++)
	{
		WCSimRootTrigger *tempTrigger = (WCSimRootTrigger *)(myEvent->GetTrigger(nTrigger));
		Int_t nDigits = tempTrigger->GetCherenkovDigiHits()->GetLast();

		if (nDigits > triggerDigits)
		{
			triggerDigits = nDigits;
			triggerNumber = nTrigger;
		}
	}

	if (triggerNumber >= 0)
	{
		myTrigger = (WCSimRootTrigger *)(myEvent->GetTrigger(triggerNumber));
	}

	if (myTrigger)
	{
		std::cout << " *** WCSimInterface::FilterTrigger(...) *** " << std::endl
				  << "  Run = "
				  << myTrigger->GetHeader()->GetRun() << std::endl
				  << "  Event = " << myTrigger->GetHeader()->GetEvtNum()
				  << std::endl
				  << "  Trigger = " << myTrigger->GetHeader()->GetSubEvtNumber()
				  << std::endl
				  << "  Ncherenkovhits = " << myTrigger->GetNcherenkovhits() << std::endl;
	}
	else
	{
		std::cout << " <warning>: failed to find trigger for this event " << std::endl;
	}

	return myTrigger;
}

void WCSimInterface::BuildEvent(WCSimRootTrigger *myTrigger)
{
	// Reset everything
	// =====
	this->ResetTrueEvent();

	this->ResetRecoEvent();

	this->ResetTrueLikelihoodTracks();

	// check for trigger
	// =================
	if (myTrigger == 0)
	{
		return;
	}

	// Build True Event
	// ================
	BuildTrueEvent(myTrigger);

	// Build Reco Event
	// ================
	BuildRecoEvent(myTrigger);

	// Build vector of true tracks
	// ===========================
	BuildTrueLikelihoodTracks();

	return;
}

WCSimTrueEvent *WCSimInterface::TrueEvent()
{
	return WCSimInterface::Instance()->GetTrueEvent();
}

WCSimRecoEvent *WCSimInterface::RecoEvent()
{
	return WCSimInterface::Instance()->GetRecoEvent();
}

WCSimRootEvent *WCSimInterface::WCSimEvent()
{
	return WCSimInterface::Instance()->GetWCSimEvent();
}

WCSimRootTrigger *WCSimInterface::WCSimTrigger()
{
	return WCSimInterface::Instance()->GetWCSimTrigger();
}

void WCSimInterface::BuildTrueEvent(WCSimRootTrigger *myTrigger)
{
	// Build True Event
	// ================
	std::cout << " *** WCSimInterface::BuildTrueEvent(...) *** " << std::endl;

	// Get Array of Tracks
	// ===================
	std::cout << "Interface" << std::endl;
	TClonesArray *fTrackArray = (TClonesArray *)(myTrigger->GetTracks());
	std::cout << "Intefrace done" << std::endl;

	Int_t ipdg = 0;
	Int_t ipdgparent = 0;
	Int_t iflag = 0;
	Double_t trkP = 0.0;
	Double_t trkE = 0.0;
	Double_t trkM = 0.0;
	Double_t trkK = 0.0;

	Double_t g4vx = 0.0;
	Double_t g4vy = 0.0;
	Double_t g4vz = 0.0;
	Double_t g4ex = 0.0;
	Double_t g4ey = 0.0;
	Double_t g4ez = 0.0;
	Double_t g4ds = 0.0;

	Double_t vx = 0.0;
	Double_t vy = 0.0;
	Double_t vz = 0.0;
	Double_t ex = 0.0;
	Double_t ey = 0.0;
	Double_t ez = 0.0;
	Double_t ds = 0.0;

	Double_t px = 0.0;
	Double_t py = 0.0;
	Double_t pz = 0.0;
	Double_t qx = 0.0;
	Double_t qy = 0.0;
	Double_t qz = 0.0;

	Int_t vtxReg = -1;
	Int_t endReg = -1;

	Int_t primeIpdg = 0;
	Double_t primeG4Vx = 0.0;
	Double_t primeG4Vy = 0.0;
	Double_t primeG4Vz = 0.0;
	Double_t primeG4Ex = 0.0;
	Double_t primeG4Ey = 0.0;
	Double_t primeG4Ez = 0.0;
	Double_t primeVx = 0.0;
	Double_t primeVy = 0.0;
	Double_t primeVz = 0.0;
	Double_t primeEx = 0.0;
	Double_t primeEy = 0.0;
	Double_t primeEz = 0.0;
	Double_t primePx = 0.0;
	Double_t primePy = 0.0;
	Double_t primePz = 0.0;
	Double_t primeTrkE = 0.0;
	Double_t primeTrkP = 0.0;

	Double_t fom = 0.0;
	Double_t maxfom = 0.0;
	Double_t beta = 0.0;

	Bool_t possible_track = 0;
	Bool_t new_track = 0;

	std::cout << "  looping over tracks: " << std::endl;
	std::cout << " TRACKS -> " << fTrackArray->GetLast() << std::endl;

	// loop over tracks
	// =================
	for (Int_t nTrack = 0; nTrack < 1 + fTrackArray->GetLast(); nTrack++)
	{
		WCSimRootTrack *myTrack = (WCSimRootTrack *)(fTrackArray->At(nTrack));

		// get track properties
		// ====================
		ipdg = myTrack->GetIpnu();
		ipdgparent = myTrack->GetParenttype();
		iflag = myTrack->GetFlag();

		trkP = myTrack->GetP();
		trkE = myTrack->GetE();
		trkM = myTrack->GetM();
		trkK = trkE - trkM;

		g4vx = myTrack->GetStart(0);
		g4vy = myTrack->GetStart(1);
		g4vz = myTrack->GetStart(2);
		g4ex = myTrack->GetStop(0);
		g4ey = myTrack->GetStop(1);
		g4ez = myTrack->GetStop(2);
		px = myTrack->GetDir(0);
		py = myTrack->GetDir(1);
		pz = myTrack->GetDir(2);

		vx = -99999.9;
		vy = -99999.9;
		vz = -99999.9;
		ex = -99999.9;
		ey = -99999.9;
		ez = -99999.9;

		vtxReg = 0;
		endReg = 0;

		ds = 0.0;
		beta = 0.0;

		qx = g4ex - g4vx;
		qy = g4ey - g4vy;
		qz = g4ez - g4vz;
		g4ds = sqrt(qx * qx + qy * qy + qz * qz);

		if (g4ds > 0.0)
		{
			qx /= g4ds;
			qy /= g4ds;
			qz /= g4ds;
		}

		// entry point
		if (endReg >= 0)
		{
			if (WCSimGeometry::Instance()->InsideDetector(g4vx, g4vy, g4vz))
			{
				vx = g4vx;
				vy = g4vy;
				vz = g4vz;
			}
			else
			{
				WCSimGeometry::Instance()->ProjectToNearEdge(g4vx, g4vy, g4vz, px, py, pz, vx, vy, vz, vtxReg);
				if (vtxReg < 0)
				{
					WCSimGeometry::Instance()->ProjectToNearEdge(g4vx, g4vy, g4vz, qx, qy, qz, vx, vy, vz, vtxReg);
				}
			}
			if (vtxReg < 0)
				endReg = -1;
		}

		// exit point
		if (vtxReg >= 0)
		{
			if (WCSimGeometry::Instance()->InsideDetector(g4ex, g4ey, g4ez))
			{
				ex = g4ex;
				ey = g4ey;
				ez = g4ez;
			}
			else
			{
				WCSimGeometry::Instance()->ProjectToNearEdge(g4ex, g4ey, g4ez, -px, -py, -pz, ex, ey, ez, endReg);
				if (endReg < 0)
				{
					WCSimGeometry::Instance()->ProjectToNearEdge(g4ex, g4ey, g4ez, -qx, -qy, -qz, ex, ey, ez, endReg);
				}
			}
			if (endReg < 0)
				vtxReg = -1;
		}

		// distance travelled through detector
		if (vtxReg >= 0 && endReg >= 0)
		{
			ds = sqrt((ex - vx) * (ex - vx) + (ey - vy) * (ey - vy) + (ez - vz) * (ez - vz));
		}

		// figure of merit
		fom = trkK * ds; // <- no idea if this is right...

		// possible track
		possible_track = 0;
		if (fabs(ipdg) == 13 || fabs(ipdg) == 11 || fabs(ipdg) == 211 || fabs(ipdg) == 321 // <- no idea if this is right...
			|| fabs(ipdg) == 2212 || ipdg == 22)
		{
			if (trkK > 0.0)
			{
				beta = trkP / trkE;
				if (fom > 0.0 && trkP / trkE > 0.75)
				{
					possible_track = 1;
				}
			}
		}

		// new track
		new_track = 0;
		if (possible_track)
		{
			if (nTrack >= 2 && trkK > fEnergyThreshold && ds > fRangeThreshold)
			{
				new_track = 1;
			}
		}

		// make new track
		// ==============
		if (new_track)
		{
			WCSimTrueTrack *trueTrack = new WCSimTrueTrack(ipdg, ipdgparent, g4vx, g4vy, g4vz, g4ex, g4ey, g4ez, vx, vy,
														   vz, ex, ey, ez, px, py, pz, trkE, trkP);
			fTrueEvent->AddTrack(trueTrack);

			fTrackList->push_back(trueTrack);

			// assign primary track
			if (fom > maxfom)
			{
				primeIpdg = ipdg;
				primeG4Vx = g4vx;
				primeG4Vy = g4vy;
				primeG4Vz = g4vz;
				primeG4Ex = g4ex;
				primeG4Ey = g4ey;
				primeG4Ez = g4ez;
				primeVx = vx;
				primeVy = vy;
				primeVz = vz;
				primeEx = ex;
				primeEy = ey;
				primeEz = ez;
				primePx = px;
				primePy = py;
				primePz = pz;
				primeTrkE = trkE;
				primeTrkP = trkP;
				maxfom = fom;
			}
		}

		// print out track info
		//    std::cout << "   [" << nTrack; if( new_track) std::cout << "*"; std::cout << "] iflag=" << iflag << " pdg=" << ipdg << ", p=" << trkP << " E=" << trkE << " K=" << trkK << std::endl;
		//    std::cout << "     g4vtx=(" << g4vx << "," << g4vy << "," << g4vz << "), g4end=(" << g4ex << "," << g4ey << "," << g4ez << ") " << std::endl;
		//    std::cout << "     g4dir=(" << px << "," << py << "," << pz << "), g4trk=(" << qx << "," << qy << "," << qz << ")" << std::endl;
		//    std::cout << "     vtx=(" << vx << "," << vy << "," << vz << "), end=(" << ex << "," << ey << "," << ez << ") " << std::endl;
		//    std::cout << "     ds=" << ds << " beta=" << beta << " fom=" << fom << " ";
		//    if( new_track ) std::cout << " [TRACK] "; // PRINT: MAKE TRACK
		//    std::cout << std::endl; // PRINT: NEW LINE

		// print warning messages
		//if (vtxReg < 0) {
		//	std::cout << "      <warning> couldn't find detector entry point! " << std::endl;
		//}

		//if (endReg < 0) {
		//	std::cout << "      <warning> couldn't find detector exit point! " << std::endl;
		//}
	}

	// set true event header
	fTrueEvent->SetHeader(primeIpdg, primeG4Vx, primeG4Vy, primeG4Vz, primeG4Ex, primeG4Ey, primeG4Ez, primeVx, primeVy,
						  primeVz, primeEx, primeEy, primeEz, primePx, primePy, primePz, primeTrkE, primeTrkP);

	std::cout << "  true event info: " << std::endl
			  << "   Ipdg = " << primeIpdg << std::endl
			  << "   G4 Vertex = ("
			  << primeG4Vx << "," << primeG4Vy << "," << primeG4Vz << ") " << std::endl
			  << "   G4 End = (" << primeG4Ex
			  << "," << primeG4Ey << "," << primeG4Ez << ") " << std::endl
			  << "   Vertex = (" << primeVx << "," << primeVy
			  << "," << primeVz << ") " << std::endl
			  << "   End = (" << primeEx << "," << primeEy << "," << primeEz
			  << ") " << std::endl
			  << "   Direction = (" << primePx << "," << primePy << "," << primePz << ") "
			  << std::endl
			  << "   Energy = " << primeTrkE << " MeV" << std::endl;

	return;
}

void WCSimInterface::BuildRecoEvent(WCSimRootTrigger *myTrigger)
{
	// Build Reco Event
	// ================
	std::cout << " *** WCSimInterface::BuildRecoEvent(...) *** " << std::endl;

	// write header
	// ============
	Int_t runNum = myTrigger->GetHeader()->GetRun();
	Int_t eventNum = myTrigger->GetHeader()->GetEvtNum();
	Int_t triggerNum = myTrigger->GetHeader()->GetSubEvtNumber();

	fRecoEvent->SetHeader(runNum, eventNum, triggerNum);

	// get array of digitized hits
	// ===========================
	TClonesArray *fDigiHitArray = (TClonesArray *)(myTrigger->GetCherenkovDigiHits());

	////////////////////
	double lowCut = WCSimParameters::Instance()->GetIgnoreLowCut();
	double highCut = WCSimParameters::Instance()->GetIgnoreHighCut();
	double length = WCSimGeometry::Instance()->GetCylLength();
	std::cout << "lowCut -> " << lowCut << ", highCut -> " << highCut << ", length -> " << length << std::endl;
	double capCut = (length / 2) - 20;
	std::cout << "capCut -> " << capCut << std::endl;
	////////////////////

	// loop over digits
	// ================
	for (Int_t nDigit = 0; nDigit < 1 + fDigiHitArray->GetLast(); nDigit++)
	{

		WCSimRootCherenkovDigiHit *myDigit = (WCSimRootCherenkovDigiHit *)(fDigiHitArray->At(nDigit));
		Int_t tube = myDigit->GetTubeId();
		Double_t rawQ = myDigit->GetQ();
		Double_t rawT = myDigit->GetT();

		Double_t calQ = rawQ;
		Double_t calT = rawT - WCSimParameters::TimeSlew(rawQ);

		Int_t region = WCSimGeometry::Instance()->GetRegion(tube);
		Double_t x = WCSimGeometry::Instance()->GetX(tube);
		Double_t y = WCSimGeometry::Instance()->GetY(tube);
		Double_t z = WCSimGeometry::Instance()->GetZ(tube);

		////////////////////
		if (lowCut != 0.0 && z < lowCut && z > -capCut)
		{
			continue;
		}
		else if (highCut != 0.0 && z > highCut && z < capCut)
		{
			continue;
		}
		else
		{
			WCSimRecoDigit *recoDigit = new WCSimRecoDigit(region, tube, x, y, z, rawT, rawQ, calT, calQ);
			if (region != WCSimGeometry::kVeto)
			{
				fRecoEvent->AddDigit(recoDigit);
				fDigitList->push_back(recoDigit);
			}
			else
			{
				fRecoEvent->AddVetoDigit(recoDigit);
				fVetoDigitList->push_back(recoDigit);
			}
		}
		////////////////////

		/*
		WCSimRecoDigit* recoDigit = new WCSimRecoDigit(region, tube, x, y, z, rawT, rawQ, calT, calQ);
		if (region != WCSimGeometry::kVeto) {
			fRecoEvent->AddDigit(recoDigit);
			fDigitList->push_back(recoDigit);
		} else {
			fRecoEvent->AddVetoDigit(recoDigit);
			fVetoDigitList->push_back(recoDigit);
		}
		*/
	}

	std::cout << "   Number of Digits = " << fRecoEvent->GetNDigits() << std::endl;
	std::cout << "   Number of Veto Digits = " << fRecoEvent->GetNVetoDigits() << std::endl;
	std::cout << "   Tracks found = " << fTrackList->size() << std::endl;

	return;
}

WCSimLikelihoodDigitArray WCSimInterface::GetWCSimLikelihoodDigitArray(int ievent)
{
	WCSimRootEvent *myEvent = (WCSimRootEvent *)(GetWCSimEvent(ievent));
	fLikelihoodDigitArray = WCSimLikelihoodDigitArray(myEvent);
	return fLikelihoodDigitArray;
}

std::vector<WCSimLikelihoodTrackBase *> WCSimInterface::GetPi0PhotonTracks(const WCSimTruthSummary &sum, double energy,
																		   TClonesArray *trajCont)
{
	std::vector<WCSimLikelihoodTrackBase *> photons;
	// Iterate through the TClonesArray looking for what we want...
	bool gotPi0 = false;
	int gotPhotons = 0;
	int pi0ID = -1; // This is the geant4 track ID for the parent pi0

	// Temporarily store the information
	TVector3 pi0Vtx;
	TVector3 pi0Dir;
	double photonEn1, photonEn2;
	TVector3 photonDir1, photonDir2;

	TVector3 mainVtx = sum.GetVertex();
	for (int e = 0; e < trajCont->GetEntries(); ++e)
	{
		// Skip the first entry for particle gun events as it is repeated
		if (e == 0 && sum.IsParticleGunEvent())
			continue;

		// Get the track
		WCSimRootTrack *trk = (WCSimRootTrack *)(*trajCont)[e];

		if (!gotPi0)
		{
			// Is it a pizero?
			if (trk->GetIpnu() != 111)
				continue;
			if (fabs(trk->GetE() - energy) > 1e-2)
				continue;
			// This is our pi0
			gotPi0 = true;
			pi0ID = trk->GetId();

			pi0Vtx = TVector3(trk->GetStart(0), trk->GetStart(1), trk->GetStart(2));
			pi0Dir = TVector3(trk->GetDir(0), trk->GetDir(1), trk->GetDir(2));
		}
		else
		{
			// Presumably the photons are after the pi0? Let's look for them, then
			if (trk->GetIpnu() != 22)
				continue; // Is it a photon?
			if (trk->GetParenttype() != 111)
				continue; // Was the parent a pi0?
			if (trk->GetParentId() != pi0ID)
				continue; // Pi0 stop and photon start at the same point?
			// This is a photon that we are looking for!
			if (gotPhotons == 0)
			{
				photonEn1 = trk->GetE();
				photonDir1 = TVector3(trk->GetDir(0), trk->GetDir(1), trk->GetDir(2));
			}
			else
			{
				photonEn2 = trk->GetE();
				photonDir2 = TVector3(trk->GetDir(0), trk->GetDir(1), trk->GetDir(2));
			}
			++gotPhotons;
		}
	}
	if (gotPi0 && gotPhotons == 2)
	{
		WCSimLikelihoodTrackBase *track;
		std::map<FitterParameterType::Type, double> extraPars;
		FitterParameterType::Type type = FitterParameterType::kConversionDistance;
		extraPars[type] = 0.0;

		track = WCSimLikelihoodTrackFactory::MakeTrack(TrackType::PhotonLike, pi0Vtx.X(), pi0Vtx.Y(), pi0Vtx.Z(),
													   sum.GetVertexT(), TMath::ACos(photonDir1.Z()), TMath::ATan2(photonDir1.Y(), photonDir1.X()), photonEn1,
													   extraPars);
		photons.push_back(track);
		track = WCSimLikelihoodTrackFactory::MakeTrack(TrackType::PhotonLike, pi0Vtx.X(), pi0Vtx.Y(), pi0Vtx.Z(),
													   sum.GetVertexT(), TMath::ACos(photonDir2.Z()), TMath::ATan2(photonDir2.Y(), photonDir2.X()), photonEn2,
													   extraPars);
		photons.push_back(track);
	}
	return photons;
}

bool WCSimInterface::EventIsVetoed()
{
	return AnyTrueLeptonEscaped();
}

bool WCSimInterface::AnyTrueLeptonEscaped()
{
	WCSimTruthSummary ts(fEvent->GetTruthSummary());
	TVector3 vtx = ts.GetVertex() * 0.1;
	TVector3 dir = ts.GetPrimaryDir(0);

	if (fEvent->GetNumberOfEvents() <= 0)
	{
		return false;
	}
	bool escapes = false;

	WCSimRootTrigger *fTrigger = fEvent->GetTrigger(0);
	TObjArray *allTracks = fTrigger->GetTracks();
	int nTracks = fTrigger->GetNtrack();

	std::set<int> startVols;
	std::set<int> endVols;

	for (int i = 0; i < nTracks; ++i)
	{
		// Only look at electron and muon tracks
		WCSimRootTrack *myTrack = (WCSimRootTrack *)allTracks->At(i);

		if (myTrack->GetIpnu() != 11 && myTrack->GetIpnu() != 13)
		{
			continue;
		}

		startVols.insert(myTrack->GetStartvol());
		endVols.insert(myTrack->GetStopvol());
	}
	escapes = escapes || ((startVols.size() != 1) || (endVols.size() != 1));
	return escapes;
}
