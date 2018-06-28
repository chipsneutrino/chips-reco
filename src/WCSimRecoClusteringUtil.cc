#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>

#include "WCSimRecoClusteringUtil.hh"
#include "WCSimRecoDigit.hh"

#include <TDirectory.h>
#include <TH3D.h>
#include <TFile.h>

// A couple of little sort functions for digits
bool SortDigitsByCharge(WCSimRecoDigit* a, WCSimRecoDigit *b);
bool SortDigitsByTime(WCSimRecoDigit* a, WCSimRecoDigit *b);

WCSimRecoClusteringUtil::WCSimRecoClusteringUtil() {

}

WCSimRecoClusteringUtil::WCSimRecoClusteringUtil(std::vector<WCSimRecoDigit*> digitVec, double qCut, double tCut,
		double dCut, unsigned int minSize) {

	this->SetDigitVector(digitVec);
	this->SetChargeCut(qCut);
	this->SetTimeCut(tCut);
	this->SetDistanceCut(dCut);
	this->SetMinSize(minSize);
}

WCSimRecoClusteringUtil::~WCSimRecoClusteringUtil() {

}

void WCSimRecoClusteringUtil::Reset() {

	fChargeCutDigits.clear();
	for (unsigned int i = 0; i < fTimeSlicedDigits.size(); ++i) {
		fTimeSlicedDigits[i].clear();
	}
	fTimeSlicedDigits.clear();
	for (unsigned int i = 0; i < fFullSlicedDigits.size(); ++i) {
		fFullSlicedDigits[i].clear();
	}
	fFullSlicedDigits.clear();

}

void WCSimRecoClusteringUtil::PerformClustering() {

	if (fDigitVec.size() == 0) {
		std::cerr << "WCSimRecoClusteringUtil: Digit vector empty, exiting." << std::endl;
		return;
	}

	this->Reset();
	this->ApplyChargeCut();
	this->SliceByTime();
	this->SliceInSpace();

}

void WCSimRecoClusteringUtil::ApplyChargeCut() {

	for (unsigned int i = 0; i < fDigitVec.size(); ++i) {
		WCSimRecoDigit* digit = fDigitVec.at(i);
		if (digit->GetRawQPEs() > fChargeCut) {
			fChargeCutDigits.push_back(digit);
		}
	}

}

void WCSimRecoClusteringUtil::SliceByTime() {

	if (fChargeCutDigits.size() == 0)
		return;

	// Firstly, sort the digit vector.
	std::sort(fChargeCutDigits.begin(), fChargeCutDigits.end(), SortDigitsByTime);

	// Look for gaps in the time distribution of the events
	std::vector<unsigned int> gapIndexVec; // Stores the last index in the current window
	gapIndexVec.push_back(0); // An initial 0 comes in helpful.
	for (unsigned int i = 1; i < fChargeCutDigits.size(); ++i) {
		double diff = fChargeCutDigits[i]->GetRawTime() - fChargeCutDigits[i - 1]->GetRawTime();
		if (diff > fTimeCut) {
			std::cout << "Found time gap = " << diff << ", " << fTimeCut << std::endl;
			gapIndexVec.push_back(i - 1);
		}
	}
	gapIndexVec.push_back(fChargeCutDigits.size() - 1); // Last element must be the end of the slice

	// Now make a series of vectors containing just those time sliced elements.
	for (unsigned int v = 1; v < gapIndexVec.size(); ++v) {
		std::vector<WCSimRecoDigit*>::iterator first = fChargeCutDigits.begin() + gapIndexVec[v - 1];
		std::vector<WCSimRecoDigit*>::iterator last = fChargeCutDigits.begin() + gapIndexVec[v];
		std::vector<WCSimRecoDigit*> tempVec(first, last);
		fTimeSlicedDigits.push_back(tempVec);
	}

}

void WCSimRecoClusteringUtil::SliceInSpace() {

	if (fTimeSlicedDigits.size() == 0)
		return;

	// LEIGH: Let's have a look at what we're making...
	TDirectory* tmpd = gDirectory;
	TFile *slicePlotFile = new TFile("slicePlotFile.root", "recreate");

	for (unsigned int t = 0; t < fTimeSlicedDigits.size(); ++t) {

		std::vector<WCSimRecoDigit*> timeSlice = fTimeSlicedDigits[t];
		std::sort(timeSlice.begin(), timeSlice.end(), SortDigitsByCharge);
		//  std::cout << " Slicing " << digitVector.size() << " digits from " << digitVector[0]->GetRawQPEs()
		//            << " to " << digitVector[digitVector.size()-1]->GetRawQPEs() << " PE." << std::endl;

		fIsDigitClustered.clear();
		for (unsigned int d = 0; d < timeSlice.size(); ++d) {
			fIsDigitClustered.push_back(0);
		}

		std::vector < std::vector<WCSimRecoDigit*> > slicedDigits;
		// Don't do any slicing if we don't have enough digits
		if (timeSlice.size() < 3) {
			std::cerr << "Not performing slicing, too few digits" << std::endl;
		} else {
			unsigned int clusterNum = 0;
			std::vector<WCSimRecoDigit*> tempVec;

			while (!AreAllDigitsClustered()) {
				slicedDigits.push_back(tempVec);
				unsigned int nextDigitSeed = GetFirstUnclusteredDigit();
				// Seed the cluster with the next unclustered digit.
				slicedDigits[clusterNum].push_back(timeSlice[nextDigitSeed]);
				fIsDigitClustered[nextDigitSeed] = true;

				// Loop over the filtered digits in the event
				bool addedHit = true;
				while (addedHit) {
					int nAdded = 0;
					for (unsigned int d = 0; d < timeSlice.size(); ++d) {
						if (fIsDigitClustered[d])
							continue;

						// Check if we should add the hit.
						fIsDigitClustered[d] = AddDigitIfClose(timeSlice[d], slicedDigits[clusterNum]);
						if (fIsDigitClustered[d])
							++nAdded;
					} // End the for loop.

					// If we added any hits then we should iterate over the list again.
					if (nAdded > 0) {
						addedHit = true;
					} else {
						addedHit = false;
					}
				} // End clustering while loop

				++clusterNum;
			} // End the while loop
		}

		// Now we should have all of our clusters from this time slice.
		// Want to copy those that are big enough into the main slice vector.
		for (unsigned int v = 0; v < slicedDigits.size(); ++v) {
			if (slicedDigits[v].size() >= fMinSliceSize) {

				fFullSlicedDigits.push_back(slicedDigits[v]);
				std::cout << "Found cluster in time slice " << t << " with " << slicedDigits[v].size() << " hits."
						<< std::endl;

				// For each slice, make a TH3D we can look at.
				std::stringstream plotName;
				plotName << "slice_" << (t * 10) + v;
				TH3D *hHist = new TH3D(plotName.str().c_str(), "", 150, -1500, 1500, 150, -1500, 1500, 150, -1500,
						1500);
				for (unsigned int z = 0; z < slicedDigits[v].size(); ++z) {
					hHist->Fill(slicedDigits[v][z]->GetX(), slicedDigits[v][z]->GetY(), slicedDigits[v][z]->GetZ());
				}
				hHist->Write();
			}
		}
	}
	// Close the file
	slicePlotFile->Close();
	delete slicePlotFile;
	gDirectory = tmpd;

}

bool WCSimRecoClusteringUtil::AddDigitIfClose(WCSimRecoDigit* digit, std::vector<WCSimRecoDigit*>& cluster) {

	for (unsigned int d = 0; d < cluster.size(); ++d) {
		if (GetDistanceBetweenDigits(digit, cluster[d]) < fDistanceCut) {
			cluster.push_back(digit);
			return true;
		}
	}
	return false;
}

double WCSimRecoClusteringUtil::GetDistanceBetweenDigits(WCSimRecoDigit* digit1, WCSimRecoDigit* digit2) {
	double dist = (digit1->GetX() - digit2->GetX()) * (digit1->GetX() - digit2->GetX());
	dist += (digit1->GetY() - digit2->GetY()) * (digit1->GetY() - digit2->GetY());
	dist += (digit1->GetZ() - digit2->GetZ()) * (digit1->GetZ() - digit2->GetZ());
	return sqrt(dist);
}

unsigned int WCSimRecoClusteringUtil::GetFirstUnclusteredDigit() {
	for (unsigned int i = 0; i < fIsDigitClustered.size(); ++i) {
		if (fIsDigitClustered[i] == false) {
			return i;
		}
	}
	return 999999;
}

bool WCSimRecoClusteringUtil::AreAllDigitsClustered() {

	if (GetFirstUnclusteredDigit() == 999999)
		return true;
	else
		return false;

}

// Sort the digits charge (largest first)
bool SortDigitsByCharge(WCSimRecoDigit *a, WCSimRecoDigit *b) {
	return a->GetRawQPEs() > b->GetRawQPEs();
}

// Sort the digits time (smallest first)
bool SortDigitsByTime(WCSimRecoDigit *a, WCSimRecoDigit *b) {
	return a->GetRawTime() < b->GetRawTime();
}

