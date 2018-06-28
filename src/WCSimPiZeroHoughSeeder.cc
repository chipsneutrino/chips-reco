/*
 * WCSimPiZeroHoughSeeder.cc
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */
#include "WCSimFitterConfig.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimPiZeroHoughSeeder.hh"
#include "WCSimPiZeroSeed.hh"
#include "WCSimPiZeroSeeder.hh"
#include "WCSimRecoSeed.hh"

#include "TLorentzVector.h"

#include <algorithm>

#ifndef REFLEX_DICTIONARY
ClassImp (WCSimPiZeroHoughSeeder)
#endif

WCSimPiZeroHoughSeeder::WCSimPiZeroHoughSeeder(WCSimFitterConfig * config) :
		WCSimPiZeroSeeder(config) {
	fMadeSeeds = false;
}

WCSimPiZeroHoughSeeder::~WCSimPiZeroHoughSeeder() {
	for (unsigned int iHoughResult = 0; iHoughResult < fHoughResults.size(); ++iHoughResult) {
		delete fHoughResults.at(iHoughResult).first;
		delete fHoughResults.at(iHoughResult).second;
	}
	fHoughResults.clear();
}

void WCSimPiZeroHoughSeeder::MakeSeeds() {
	if (fMadeSeeds) {
		return;
	}

	// Run the Hough transform
	RunHough();

	// Do we get two separate vertices or one?
	FitUsingHoughResults();
	std::cout << "There are " << fLikelihoodDigitArray->GetNDigits() << " digits" << std::endl;
	std::cout << "Digit 0 charge = ";
	std::cout << fLikelihoodDigitArray->GetDigit(0)->GetQ() << std::endl;

	fMadeSeeds = true;
	return;
}

void WCSimPiZeroHoughSeeder::RunHough() {
	std::cout << " *** WCSimPiZeroHoughSeeder::RunHough() *** " << std::endl;
	fHoughResults.clear();

	// Run the old Hough transform reco
	//WCSimRecoSeed * myReco = dynamic_cast<WCSimRecoSeed*>(WCSimRecoFactory::Instance()->MakeReco("seed")); // This calls new() - delete when done
	WCSimRecoSeed * myReco = new WCSimRecoSeed(); // This calls new() - delete when done
	myReco->SetNumberOfTracks(2); // We only run this on pi0 -> gamma gamma hypotheses

	WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
	std::vector<WCSimRecoEvent*> slicedEvents = myReco->RunSeed(recoEvent);

	// Make a vector of all of the available rings we have.
	std::vector < std::pair<WCSimRecoRing*, double> > primaryRings;
	std::vector < std::pair<WCSimRecoRing*, double> > otherRings;

	// Get the two highest rings from different slices:
	for (unsigned int se = 0; se < slicedEvents.size(); ++se) {
		std::pair<WCSimRecoRing*, double> tempPair;
		tempPair = std::make_pair<WCSimRecoRing*, double>(slicedEvents[se]->GetRing(0), slicedEvents[se]->GetVtxTime());
		primaryRings.push_back(tempPair);
	}
	std::sort(primaryRings.begin(), primaryRings.end(), WCSimLikelihoodFitter::RingSort);
//	for(int i = 0; i < primaryRings.size(); ++i)
//	{
//		std::cout << "Primary ring " << i << "/" << primaryRings.size() << " has height " << primaryRings.at(i).first->GetHeight() << "  ";
//		std::cout << "Vtx (" << primaryRings.at(i).first->GetVtxX() << ", " << primaryRings.at(i).first->GetVtxY() << ", " << primaryRings.at(i).first->GetVtxZ() << ")";
//		std::cout << "   Dir(" << primaryRings.at(i).first->GetDirX() << ", " << primaryRings.at(i).first->GetDirY() << ", " << primaryRings.at(i).first->GetDirZ() << ")" << std::endl;
//	}

	// Get the two highest rings from any slice:
	for (unsigned int se = 0; se < slicedEvents.size(); ++se) {
		for (int r = 0; r < slicedEvents[se]->GetNRings(); ++r) {
			std::pair<WCSimRecoRing*, double> tempPair;
			tempPair = std::make_pair<WCSimRecoRing*, double>(slicedEvents[se]->GetRing(r),
					slicedEvents[se]->GetVtxTime());
			otherRings.push_back(tempPair);
		}
	}
	std::sort(otherRings.begin(), otherRings.end(), RingSort);
	for (int i = 0; (size_t) i < otherRings.size(); ++i) {
		std::cout << "Primary ring " << i << "/" << otherRings.size() << " has height "
				<< otherRings.at(i).first->GetHeight() << "  ";
		std::cout << "Vtx (" << otherRings.at(i).first->GetVtxX() << ", " << otherRings.at(i).first->GetVtxY() << ", "
				<< otherRings.at(i).first->GetVtxZ() << ")";
		std::cout << "   Dir(" << otherRings.at(i).first->GetDirX() << ", " << otherRings.at(i).first->GetDirY() << ", "
				<< otherRings.at(i).first->GetDirZ() << ")" << std::endl;
	}

	// Take the leading two ring pairs from each list
	// Loop over an array of pointers to the two sets of rings to avoid duplicating code
	std::vector<std::pair<WCSimRecoRing*, double> >* houghRings[2] = {&primaryRings, &otherRings};
	for (int i = 0; i <= 1; ++i) {
		std::vector < std::pair<WCSimRecoRing*, double> > ringList = *(houghRings[i]);

		double seedX, seedY, seedZ, seedT;
		double dirX, dirY, dirZ;
		double seedTheta, seedPhi;
		double seedE = 1000; //< Default value; will fit this later
		if (ringList.size() >= 2) {
			// Ignore if the first ring is much higher than the second
			if (ringList.at(0).first->GetHeight() > 2.0 * ringList.at(1).first->GetHeight()) {
				break;
			}

			seedX = ringList.at(0).first->GetVtxX();
			seedY = ringList.at(0).first->GetVtxY();
			seedZ = ringList.at(0).first->GetVtxZ();
			seedT = ringList.at(0).second;
			dirX = ringList.at(0).first->GetDirX();
			dirY = ringList.at(0).first->GetDirY();
			dirZ = ringList.at(0).first->GetDirZ();

			seedTheta = TMath::ACos(dirZ);
			if (dirY != 0.0) {
				seedPhi = TMath::ATan2(dirY, dirX);
			} // Ensure range is -pi to pi
			else {
				seedPhi = (dirX < 0.0) ? 0.5 * TMath::Pi() : -0.5 * TMath::Pi();
			}
			WCSimLikelihoodTrackBase * track1 = WCSimLikelihoodTrackFactory::MakeTrack(TrackType::PhotonLike, seedX,
					seedY, seedZ, seedT, seedTheta, seedPhi, seedE, 0.0);

			seedX = ringList.at(1).first->GetVtxX();
			seedY = ringList.at(1).first->GetVtxY();
			seedZ = ringList.at(1).first->GetVtxZ();
			seedT = ringList.at(1).second;
			dirX = ringList.at(1).first->GetDirX();
			dirY = ringList.at(1).first->GetDirY();
			dirZ = ringList.at(1).first->GetDirZ();
			seedTheta = TMath::ACos(dirZ);
			if (dirY != 0.0) {
				seedPhi = TMath::ATan2(dirY, dirX);
			} // Ensure range is -pi to pi
			else {
				seedPhi = (dirX < 0.0) ? 0.5 * TMath::Pi() : -0.5 * TMath::Pi();
			}
			WCSimLikelihoodTrackBase * track2 = WCSimLikelihoodTrackFactory::MakeTrack(TrackType::PhotonLike, seedX,
					seedY, seedZ, seedT, seedTheta, seedPhi, seedE, 0.0);
			SetStartingEnergies(ringList.at(0).first, ringList.at(1).first, track1, track2);

			fHoughResults.push_back(std::make_pair(track1, track2));
		}
	}
	if (fHoughResults.size() == 0) {
		std::cout << "Didn't find two rings with the Hough transform" << std::endl;
	}

	// Tidy up
	delete myReco;

	// Need to delete the elements of slicedEvents as they are not destroyed by WCSimRecoSlicer
	for (unsigned int v = 0; v < slicedEvents.size(); ++v) {
		delete slicedEvents[v];
	}
}

void WCSimPiZeroHoughSeeder::FitUsingHoughResults() {
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> >::iterator houghResultItr;
	std::cout << "Size of fHoughResults = " << fHoughResults.size() << std::endl;
	houghResultItr = fHoughResults.begin();

	WCSimFitterConfig * backupConfig = fFitterConfig;
	while (houghResultItr != fHoughResults.end()) {
		double conv1, conv2;
		TVector3 vertex = FindSharedVertex(*houghResultItr, conv1, conv2);

		WCSimFitterConfig config;

		config.SetNumTracks(2);

		std::cout << "Fit using Hough results" << std::endl;
		config.SetTrackType(0, "PhotonLike");
		config.SetTrackType(1, "PhotonLike");
		config.SetJoinParametersTogether(0, 1, "kVtxX");
		config.SetJoinParametersTogether(0, 1, "kVtxY");
		config.SetJoinParametersTogether(0, 1, "kVtxZ");
		config.SetJoinParametersTogether(0, 1, "kVtxT");
		if (fFitterConfig->GetIsPiZeroFit()) {
			config.SetIsPiZeroFit(fFitterConfig->GetIsPiZeroFit());
			config.SetForcePiZeroMass(fFitterConfig->GetForcePiZeroMass());
		}

		config.SetParameter(0, "kVtxX", fFitterConfig->GetParMin(0, "kVtxX"), fFitterConfig->GetParMax(0, "kVtxX"),
				vertex.X(), fFitterConfig->GetParStep(0, "kVtxX"), false);
		config.SetParameter(0, "kVtxY", fFitterConfig->GetParMin(0, "kVtxY"), fFitterConfig->GetParMax(0, "kVtxY"),
				vertex.Y(), fFitterConfig->GetParStep(0, "kVtxY"), false);
		config.SetParameter(0, "kVtxZ", fFitterConfig->GetParMin(0, "kVtxZ"), fFitterConfig->GetParMax(0, "kVtxZ"),
				vertex.Z(), fFitterConfig->GetParStep(0, "kVtxZ"), false);
		config.SetParameter(0, "kVtxT", fFitterConfig->GetParMin(0, "kVtxT"), fFitterConfig->GetParMax(0, "kVtxT"),
				houghResultItr->first->GetT(), fFitterConfig->GetParStep(0, "kVtxT"), false);
		config.SetParameter(0, "kDirTh", fFitterConfig->GetParMin(0, "kDirTh"), fFitterConfig->GetParMax(0, "kDirTh"),
				houghResultItr->first->GetTheta(), fFitterConfig->GetParStep(0, "kDirTh"), false);
		config.SetParameter(0, "kDirPhi", fFitterConfig->GetParMin(0, "kDirPhi"),
				fFitterConfig->GetParMax(0, "kDirPhi"), houghResultItr->first->GetPhi(),
				fFitterConfig->GetParStep(0, "kDirPhi"), false);
		config.SetParameter(0, "kEnergy", fFitterConfig->GetParMin(0, "kEnergy"),
				fFitterConfig->GetParMax(0, "kEnergy"), houghResultItr->first->GetE(),
				fFitterConfig->GetParStep(0, "kEnergy"), false);
		config.SetParameter(0, "kConversionDistance", fFitterConfig->GetParMin(0, "kConversionDistance"),
				fFitterConfig->GetParMax(0, "kConversionDistance"), houghResultItr->first->GetConversionDistance(),
				fFitterConfig->GetParStep(0, "kConversionDistance"), false);
		config.SetParameter(1, "kDirTh", fFitterConfig->GetParMin(1, "kDirTh"), fFitterConfig->GetParMax(1, "kDirTh"),
				houghResultItr->second->GetTheta(), fFitterConfig->GetParStep(1, "kDirTh"), false);
		config.SetParameter(1, "kDirPhi", fFitterConfig->GetParMin(1, "kDirPhi"),
				fFitterConfig->GetParMax(1, "kDirPhi"), houghResultItr->second->GetPhi(),
				fFitterConfig->GetParStep(1, "kDirPhi"), false);
		config.SetParameter(1, "kEnergy", fFitterConfig->GetParMin(1, "kEnergy"),
				fFitterConfig->GetParMax(1, "kEnergy"), houghResultItr->second->GetE(),
				fFitterConfig->GetParStep(1, "kEnergy"), false);
		config.SetParameter(1, "kConversionDistance", fFitterConfig->GetParMin(1, "kConversionDistance"),
				fFitterConfig->GetParMax(1, "kConversionDistance"), houghResultItr->second->GetConversionDistance(),
				fFitterConfig->GetParStep(1, "kConversionDistance"), false);

		fFitterConfig = &config;
		fFitterTrackParMap = WCSimFitterTrackParMap(&config);
		fFitterTrackParMap.Set();

		fBestFit.clear();
		for (unsigned int i = 0; i < fBestFit.size(); ++i) {
			delete fBestFit.at(i);
			fBestFit.at(i) = 0x0;
		}
		double twoLnL = FitTwoTracks();
		WCSimPiZeroSeed * seed = new WCSimPiZeroSeed(fBestFit.at(0), fBestFit.at(1), twoLnL);
		fPiZeroSeeds.push_back(seed);
		houghResultItr++;
		fFitterConfig = backupConfig;
	}

	// Reset the TrackParMap
	std::sort(fPiZeroSeeds.begin(), fPiZeroSeeds.end());
	fFitterTrackParMap = WCSimFitterTrackParMap(fFitterConfig);
	fFitterTrackParMap.Set();
	std::cout << "Finished Hough fitting" << std::endl;
}

double WCSimPiZeroHoughSeeder::FitTwoTracks() {
	FixVertex();
	FixDirection();
	FreeConversionLength();
	FreeEnergy();
	double twoLnL = FitAndGetLikelihood("Simplex");
	return twoLnL;
}

void WCSimPiZeroHoughSeeder::SetStartingEnergies(WCSimRecoRing* ring1, WCSimRecoRing* ring2,
		WCSimLikelihoodTrackBase* track1, WCSimLikelihoodTrackBase* track2) {
	// Set the energies of the two tracks we found with the Hough transform
	// so that the invariant mass of the tracks is that of a pi0, and the
	// ratio of the two energies is the same as the ratio of the two Hough
	// peak heights

	double heightRatio = ring1->GetHeight() / ring2->GetHeight();
	TVector3 dir1(ring1->GetDirX(), ring1->GetDirY(), ring1->GetDirZ());
	TVector3 dir2(ring2->GetDirX(), ring2->GetDirY(), ring2->GetDirZ());
	double oneMinusCosTheta = 1.0 - dir1.Dot(dir2);
	double mPiZero = 134.9766;

	double e1 = mPiZero * sqrt(heightRatio / (2 * oneMinusCosTheta));
	double e2 = e1 / heightRatio;

	track1->SetE(e1);
	track2->SetE(e2);
	return;
}

TVector3 WCSimPiZeroHoughSeeder::FindSharedVertex(
		std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> tracks, double &convDistance1,
		double &convDistance2) {
	convDistance1 = 0.0;
	convDistance2 = 0.0;

	WCSimLikelihoodTrackBase * track1 = tracks.first;
	WCSimLikelihoodTrackBase * track2 = tracks.second;

	TVector3 vtxDiff = track1->GetVtx() - track2->GetVtx();

	// Vertices are effectively the same
	if (vtxDiff.Mag() < 5) {
		return track1->GetVtx();
	}

	// Solving for the points of closest approach of two vectors uses the following quantities
	// Write the trajectories as: l1 = vtx1 + alpha * momentum1 unit vector
	// 							  l2 = vtx2 + beta  * momentum2 unit vector
	// Then define the vector d = l1 - l2 and solve for the alpha and beta that mimimise d^2
	// using partial derivatives
	TVector3 dirDiff = track1->GetDir() - track2->GetDir();
	TVector3 dirSum = track1->GetDir() + track2->GetDir();
	double cosTheta = track1->GetDir().Dot(track2->GetDir());

	// Directions are the same: 1 +/- cosTheta will blow up
	if (cosTheta <= -(1 - 1e6) || cosTheta >= (1 - 1e6)) {
		return track1->GetVtx();
	}

	// The values of alpha and beta at the point of closest approach
	double distTrack1 = (-vtxDiff.Dot(dirDiff)) / (2 - 2 * cosTheta) - (dirSum.Dot(vtxDiff)) / (2 + 2 * cosTheta);
	double distTrack2 = (-vtxDiff.Dot(dirDiff)) / (2 - 2 * cosTheta) + (dirSum.Dot(vtxDiff)) / (2 + 2 * cosTheta);

	// Need to be sure we're extrapolating _backwards_ from the vertex
	if (distTrack1 <= 0 && distTrack2 <= 0) {
		convDistance1 = fabs(distTrack1);
		convDistance2 = fabs(distTrack2);
		return 0.5 * (track1->GetPropagatedPos(distTrack1) + track2->GetPropagatedPos(distTrack2));
	}

	// Failsafe
	return track1->GetVtx();
}
