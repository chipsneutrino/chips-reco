/*
 * WCSimTimePredictor.cc
 *
 *  Created on: 7 Jul 2015
 *      Author: andy
 */

#include "WCSimTimePredictor.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimEmissionProfiles.hh"
#include "WCSimEmissionProfileManager.hh"
#include "WCSimLikelihoodTuner.hh"

#ifndef REFLEX_DICTIONARY
ClassImp (WCSimTimePredictor)
#endif

WCSimTimePredictor::WCSimTimePredictor() {
	fDigitArray = 0x0;
	fTuner = 0x0;
}

WCSimTimePredictor::WCSimTimePredictor(WCSimLikelihoodDigitArray* myArray,
		WCSimEmissionProfileManager * myEmissionProfileManager) {
	fDigitArray = myArray;
	fEmissionProfileManager = myEmissionProfileManager;
	fTuner = new WCSimLikelihoodTuner(myArray, myEmissionProfileManager);
}

WCSimTimePredictor::~WCSimTimePredictor() {
	// TODO Auto-generated destructor stub
	if (fTuner != 0x0) {
		delete fTuner;
		fTuner = 0x0;
	}
}

void WCSimTimePredictor::AddTrack(WCSimLikelihoodTrackBase* myTrack) {
	// std::cout << "Time predictor: adding track" << std::endl;
	fTracks.push_back(myTrack);
}

void WCSimTimePredictor::SetTracks(std::vector<WCSimLikelihoodTrackBase*> myTracks) {
	// std::cout << "Time predictor: setting tracks" << std::endl;
	fTracks = myTracks;
	// std::cout << "Time predictor has " << fTracks.size() << " tracks" << std::endl;
}

void WCSimTimePredictor::ClearTracks() {
	// std::cout << "Time predictor clearing tracks" << std::endl;
	fTracks.clear();
}

void WCSimTimePredictor::UpdateDigitArray(WCSimLikelihoodDigitArray* myDigitArray) {
	fDigitArray = myDigitArray;
}

double WCSimTimePredictor::GetPredictedTime(unsigned int iDigit) {
	WCSimLikelihoodDigit * digit = fDigitArray->GetDigit(iDigit);
	return GetPredictedTime(digit);
}

double WCSimTimePredictor::GetPredictedTime(WCSimLikelihoodDigit* myDigit) {
	double earliestArrival = -999.9;
	std::vector<WCSimLikelihoodTrackBase*>::iterator trackItr = fTracks.begin();
	while (trackItr != fTracks.end()) {
		bool isHit = false;
		double arrivalTime = GetHitTime(*trackItr, myDigit, isHit);
		if (isHit and (arrivalTime < earliestArrival or earliestArrival == -999.9)) {
			earliestArrival = arrivalTime;
		}
		++trackItr;
	}
	return earliestArrival;
}

double WCSimTimePredictor::GetHitTime(WCSimLikelihoodTrackBase* myTrack, WCSimLikelihoodDigit* myDigit, bool &isHit) {

	// Get survival distance for the particle
	double survivalDistance = GetSurvivalDistance(myTrack);

	// Get escape distance
	double escapeDistance = GetEscapeDistance(myTrack);

	return GetHitTime(myTrack, myDigit, isHit, survivalDistance, escapeDistance);
}

double WCSimTimePredictor::GetSurvivalDistance(WCSimLikelihoodTrackBase* myTrack) {
	return fEmissionProfileManager->GetStoppingDistance(myTrack);
}

// Get the distance the charge particle travels before photons released at the Cherenkov angle start
double WCSimTimePredictor::GetTravelDistance(WCSimLikelihoodTrackBase* myTrack, WCSimLikelihoodDigit* myDigit) {

	// We need the distance and angle from vertex to PMT
	TVector3 pmt = myDigit->GetPos(); // Get it in metres
	TVector3 vertex = myTrack->GetVtx(); // Get it into metres
	double dist = (pmt - vertex).Mag();

	TVector3 direction = myTrack->GetDir();
	double cosAngToDir = (pmt - vertex).Dot(direction) / dist; // We need the cosine later so work it out first to save a cos() call
	double angToDir = TMath::ACos(cosAngToDir);

	// Cherenkov angle needs to be < 90 so we can give up here if angToDir isn't
	if (angToDir > 0.5 * TMath::Pi()) {
		return -999;
	}
	double sinAngToDir = TMath::Sin(angToDir);  // We'll also need the sine later

	// std::cout << "PMT:" << std::endl;
	// pmt.Print();
	// std::cout << "Vtx: " << std::endl;
	// vertex.Print();
	// std::cout << "Direction: " << std::endl;
	// (myTrack->GetDir() * (1.0 / (myTrack->GetDir().Mag()))).Print();

	// Now get some values relating to the Cherenkov angle
	double refInd = myDigit->GetAverageRefIndex();
	double cosCherenkovAngle = 1.0 / refInd;
	double sinCherenkovAngle = sqrt(1.0 - cosCherenkovAngle * cosCherenkovAngle); // Is about 40 degrees so no floating point worries

	// We want the s distance at which a photon emitted at the Cherenkov angle hits the PMT
	// This eventually gives us a quadratic in s to solve, of the form as^2 + bs + c = 0
	double a = 1.0;
	double b = -2.0 * dist * cosAngToDir;
	double c = dist * dist * (1 - ((sinAngToDir * sinAngToDir) / (sinCherenkovAngle * sinCherenkovAngle)));

	double toReturn = -999;
	double discrim = b * b - 4 * a * c;
	double sqrtDiscrim = sqrt(discrim);
	if (discrim < 0) {
		return -999;
	} else {
		double s1 = -b + sqrtDiscrim;
		double s2 = -b - sqrtDiscrim;
		if (s1 <= 0 && s2 <= 0) {
			return -999;
		} else if (s1 > 0 && s2 <= 0) { /*std::cout << "s1 = " << s1 << std::endl;*/
			toReturn = s1;
		} else if (s2 > 0 && s1 <= 0) { /*std::cout << "s2 = " << s2 << std::endl;*/
			toReturn = s2;
		} else if (s1 > 0 && s2 > 0 && s1 > s2) { /*std::cout << "s2 = " << s2 << std::endl;*/
			toReturn = s2;
		} else if (s2 > 0 && s1 > 0 && s2 < s1) { /*std::cout << "s1 = " << s1 << std::endl;*/
			toReturn = s1;
		}
		// if(myDigit->GetZ() < -900 && toReturn)
		// {
		//   std::cout << "Digit " << myDigit->GetTubeId() << "  Z = " <<myDigit->GetZ() << "   distance = " << toReturn << "    s1 = " << s1 << "    s2 = " << s2 << std::endl;
		// }
	}
	return toReturn;
}

std::vector<double> WCSimTimePredictor::GetAllPredictedTimes() {
	// std::cout << "GetAllPredictedTimes" << std::endl;

	std::vector<double> allPredictedTimes;

	// Get survival distance for the particle
	std::vector<WCSimLikelihoodTrackBase*>::iterator trackItr;
	std::vector<double> survivalDistances;
	std::vector<double> escapeDistances;

	for (trackItr = fTracks.begin(); trackItr != fTracks.end(); ++trackItr) {
		// (*trackItr)->Print();
		survivalDistances.push_back(GetSurvivalDistance(*trackItr));
		escapeDistances.push_back(GetEscapeDistance(*trackItr));
	}
	// std::cout << "Size of fTracks = " << fTracks.size() << std::endl;
	int nonDefaultCounter = 0;
	for (int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit) {
		WCSimLikelihoodDigit * myDigit = fDigitArray->GetDigit(iDigit);
		double earliestArrival = -999.9;

		for (unsigned int jTrack = 0; jTrack < fTracks.size(); ++jTrack) {
			bool isHit = false;
			double arrivalTime = GetHitTime(fTracks.at(jTrack), myDigit, isHit, survivalDistances.at(jTrack),
					escapeDistances.at(jTrack));
			if (isHit && (arrivalTime < earliestArrival or earliestArrival == -999.9)) {
				earliestArrival = arrivalTime;
			}
		}
		allPredictedTimes.push_back(earliestArrival);
		if (earliestArrival > -999.9) {
			nonDefaultCounter++;
		}
	}
	// std::cout << "Non -999.9 times = " << nonDefaultCounter << std::endl;
	return allPredictedTimes;
}

double WCSimTimePredictor::GetHitTime(WCSimLikelihoodTrackBase* myTrack, WCSimLikelihoodDigit* myDigit, bool& isHit,
		const double& survivalDistance, const double& escapeDistance) {

	double arrivalTime = -999;
	// How far does particle have to travel before light given out at Cherenkov angle hits the PMT?
	double travelDistance = GetTravelDistance(myTrack, myDigit);

	if ((travelDistance > 0) && (travelDistance < survivalDistance) && (travelDistance < escapeDistance)) {
		double photonDistance = (myDigit->GetPos() - myTrack->GetPropagatedPos(travelDistance)).Mag();

		double n = myDigit->GetAverageRefIndex();
		arrivalTime = myTrack->GetT()
				+ (0.01 * travelDistance / TMath::C() + 0.01 * photonDistance / (TMath::C() / n)) * 1e9;
		//std::cout << "TravelDistance = " << travelDistance << " and photonDistance = " << photonDistance << " and c = " << TMath::C() << " so arrives as " << arrivalTime << std::endl;
	}

	isHit = (arrivalTime != -999);
	return arrivalTime;
}

double WCSimTimePredictor::GetEscapeDistance(WCSimLikelihoodTrackBase* myTrack) {
	if (fTuner == 0x0) {
		fTuner = new WCSimLikelihoodTuner(fDigitArray, fEmissionProfileManager);
	}
	double cutoff = fTuner->GetCutoff(myTrack);
	return cutoff;
}
