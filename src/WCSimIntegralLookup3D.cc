/*
 * WCSimIntegralLookup3D.cc
 *
 *  Created on: 23 Mar 2015
 *      Author: ajperch
 */

#include <stdexcept>
#include "WCSimEmissionProfiles.hh"
#include "WCSimIntegralLookup3D.hh"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TString.h>

#ifndef REFLEX_DICTIONARY
ClassImp (WCSimIntegralLookup3D)
ClassImp (WCSimIntegralLookupHistArray)
#endif

//////////////////////////////////////////////////////////////////////////////////////
// Now the main class to do the looking-up                                          //
//////////////////////////////////////////////////////////////////////////////////////

WCSimIntegralLookup3D::WCSimIntegralLookup3D(TString fileName) {

	fHistFile = 0x0;
	fHistFile = new TFile(fileName.Data(), "READ");
	if (!fHistFile) {
		throw std::runtime_error("Could not open file");
	}
	if (fHistFile->IsZombie()) {
		throw std::runtime_error("Integral lookup TFile is a zombie");
	}

	// This line sometimes fails. When it tries to load the VERY LARGE WCSimIntegralLookupHistArray from
	// the file it ends up quiting the execution of the whole program with no error. Some wierd ROOT thing.
	// TODO: Figure this out!!!
	fHistFile->GetObject("WCSimIntegralLookupHistArray", fIntegrals);

	if (fIntegrals == NULL) {
		throw std::runtime_error("Could not get WCSimIntegralLookupHistArray from file");
	}
}

WCSimIntegralLookup3D::~WCSimIntegralLookup3D() {
	fHistFile->Close();
	delete fHistFile;
}

double WCSimIntegralLookup3D::GetRhoIntegral(const double& E, const double &s) {
	return fIntegrals->GetRhoIntegral(E, 0);
}

double WCSimIntegralLookup3D::GetRhoSIntegral(const double& E, const double &s) {
	return fIntegrals->GetRhoIntegral(E, 1);
}

double WCSimIntegralLookup3D::GetRhoSSIntegral(const double& E, const double &s) {
	return fIntegrals->GetRhoIntegral(E, 2);
}

double WCSimIntegralLookup3D::GetRhoGIntegral(const double& E, const double &s, const double& R0,
		const double& cosTh0) {
	return fIntegrals->GetRhoGIntegral(E, R0, cosTh0, 0);
}

double WCSimIntegralLookup3D::GetRhoGSIntegral(const double& E, const double &s, const double& R0,
		const double& cosTh0) {
	return fIntegrals->GetRhoGIntegral(E, R0, cosTh0, 1);
}

double WCSimIntegralLookup3D::GetRhoGSSIntegral(const double& E, const double &s, const double& R0,
		const double& cosTh0) {
	return fIntegrals->GetRhoGIntegral(E, R0, cosTh0, 2);
}

//////////////////////////////////////////////////////////////////////
// WCSimIntegralLookupHistArray is a class to hold an array of      //
// tabulated integrals as a function of energy, binned in distance  //
// and angle from the track to the PMT.  							//
//////////////////////////////////////////////////////////////////////

void WCSimIntegralLookupHistArray::SetRhoG(double R0, double cosTh0, TGraph gr) {
	SetGraph(GetArrayIndex(R0, cosTh0), fRhoGArr, gr);
}

void WCSimIntegralLookupHistArray::SetRhoGS(double R0, double cosTh0, TGraph gr) {
	SetGraph(GetArrayIndex(R0, cosTh0), fRhoGSArr, gr);
}

void WCSimIntegralLookupHistArray::SetRhoGSS(double R0, double cosTh0, TGraph gr) {
	SetGraph(GetArrayIndex(R0, cosTh0), fRhoGSSArr, gr);
}

void WCSimIntegralLookupHistArray::SetGraph(int bin, std::vector<TGraph*>& graphVec, TGraph gr) {
	if (graphVec.at(bin) != NULL) {
		delete graphVec.at(bin);
	}
	graphVec.at(bin) = new TGraph(gr);
}

void WCSimIntegralLookupHistArray::SetRhoGSpline(double R0, double cosTh0, TSpline3 * spline) {
	SetSpline(GetArrayIndex(R0, cosTh0), fRhoGSplineArr, spline);
}

void WCSimIntegralLookupHistArray::SetRhoGSSpline(double R0, double cosTh0, TSpline3 * spline) {
	SetSpline(GetArrayIndex(R0, cosTh0), fRhoGSSplineArr, spline);
}

void WCSimIntegralLookupHistArray::SetRhoGSSSpline(double R0, double cosTh0, TSpline3 * spline) {
	SetSpline(GetArrayIndex(R0, cosTh0), fRhoGSSSplineArr, spline);
}

void WCSimIntegralLookupHistArray::SetRhoSSpline(TSpline3* spline) {
	if (fRhoSSpline != NULL) {
		delete fRhoSSpline;
	}
	fRhoSSpline = (TSpline3*) spline->Clone();
}

void WCSimIntegralLookupHistArray::SetRhoSSSpline(TSpline3* spline) {
	if (fRhoSSSpline != NULL) {
		delete fRhoSSSpline;
	}
	fRhoSSSpline = (TSpline3*) spline->Clone();
}

void WCSimIntegralLookupHistArray::SetSpline(int bin, std::vector<TSpline3*>& splineVec, TSpline3* spline) {
	if (splineVec.at(bin) != NULL) {
		delete splineVec.at(bin);
	}
	if (spline != NULL) {
		std::cout << "Setting spline!" << std::endl;
		splineVec.at(bin) = (TSpline3*) (spline->Clone());
	} else {
		splineVec.at(bin) = NULL;
	}
	return;
}

double WCSimIntegralLookupHistArray::GetRhoIntegral(double E, const int sPower) {
	if (sPower == 0) {
		return 1;
	}
	if (E < fEMin) {
		E = fEMin;
	}
	if (E > fEMax) {
		E = fEMax;
	} else if (sPower == 1 && fRhoSSpline != NULL) {
		return fRhoSSpline->Eval(E);
	} else if (sPower == 2 && fRhoSSSpline != NULL) {
		return fRhoSSSpline->Eval(E);
	}

	assert(sPower == 0 || sPower == 1 || sPower == 2);
	return 0;
}

double WCSimIntegralLookupHistArray::GetRhoGIntegral(double E, const double R0, const double cosTh0, const int sPower) {
	unsigned long int bin = GetArrayIndex(R0, cosTh0);
	if (bin >= fRhoGSplineArr.size()) {
		std::cout << "R0 = " << R0 << " and cosTh0 = " << cosTh0 << " and sPower = " << sPower << " so bin is " << bin
				<< " and array has " << fRhoGSplineArr.size() << " bins" << std::endl;
	}

	if (E < fEMin) {
		E = fEMin;
	}
	if (E > fEMax) {
		E = fEMax;
	}
	if (sPower == 0 && fRhoGSplineArr.at(bin) != NULL) {
		return fRhoGSplineArr.at(bin)->Eval(E);
	}
	if (sPower == 1 && fRhoGSSplineArr.at(bin) != NULL) {
		return fRhoGSSplineArr.at(bin)->Eval(E);
	}
	if (sPower == 2 && fRhoGSSSplineArr.at(bin) != NULL) {
		return fRhoGSSSplineArr.at(bin)->Eval(E);
	}

	assert(sPower == 0 || sPower == 1 || sPower == 2);
	return 0;
}

WCSimIntegralLookupHistArray::WCSimIntegralLookupHistArray() :
		fEMin(-999), fEMax(-999), fR0Min(-999), fR0Max(-999), fR0Bins(-999), fCosTheta0Min(-999), fCosTheta0Max(-999), fCosTheta0Bins(
				-999), fRhoInt(NULL), fRhoSInt(NULL), fRhoSSInt(NULL), fRhoSSpline(NULL), fRhoSSSpline(NULL) {
}

WCSimIntegralLookupHistArray::WCSimIntegralLookupHistArray(double EMin, double EMax, double R0Max, double R0Min,
		int R0Bins, double cosTh0Max, double cosTh0Min, int cosTh0Bins) :
		fEMin(EMin), fEMax(EMax), fR0Min(R0Min), fR0Max(R0Max), fR0Bins(R0Bins), fCosTheta0Min(cosTh0Min), fCosTheta0Max(
				cosTh0Max), fCosTheta0Bins(cosTh0Bins), fRhoInt(NULL), fRhoSInt(NULL), fRhoSSInt(NULL), fRhoSSpline(
				NULL), fRhoSSSpline(NULL) {
	std::cout << "fEMin = " << fEMin << std::endl;
	std::cout << "fEMax = " << fEMax << std::endl;
	ResetArrays();
}

WCSimIntegralLookupHistArray::WCSimIntegralLookupHistArray(const WCSimIntegralLookupHistArray& other) :
		fEMin(other.fEMin), fEMax(other.fEMax), fR0Min(other.fR0Min), fR0Max(other.fR0Max), fR0Bins(other.fR0Bins), fCosTheta0Min(
				other.fCosTheta0Min), fCosTheta0Max(other.fCosTheta0Max), fCosTheta0Bins(other.fCosTheta0Bins), fRhoInt(
				NULL), fRhoSInt(NULL), fRhoSSInt(NULL), fRhoSSpline(NULL), fRhoSSSpline(NULL) {
	fRhoInt = (TGraph*) other.fRhoInt->Clone();
	fRhoSInt = (TGraph*) other.fRhoSInt->Clone();
	fRhoSSInt = (TGraph*) other.fRhoSSInt->Clone();
	fRhoSSpline = (TSpline3*) other.fRhoSSpline->Clone();
	fRhoSSSpline = (TSpline3*) other.fRhoSSSpline->Clone();

	std::vector<TGraph*>::const_iterator grIt;
	for (grIt = other.fRhoGArr.begin(); grIt != fRhoGArr.end(); ++grIt) {
		fRhoGArr.push_back((TGraph*) ((*grIt)->Clone()));
	}
	for (grIt = other.fRhoGSArr.begin(); grIt != fRhoGSArr.end(); ++grIt) {
		fRhoGSArr.push_back((TGraph*) ((*grIt)->Clone()));
	}
	for (grIt = other.fRhoGSSArr.begin(); grIt != fRhoGSSArr.end(); ++grIt) {
		fRhoGSSArr.push_back((TGraph*) ((*grIt)->Clone()));
	}

	std::vector<TSpline3*>::const_iterator splineIt;
	for (splineIt = other.fRhoGSplineArr.begin(); splineIt != other.fRhoGSplineArr.end(); ++splineIt) {
		fRhoGSplineArr.push_back((TSpline3*) ((*splineIt)->Clone()));
	}
	for (splineIt = other.fRhoGSSplineArr.begin(); splineIt != other.fRhoGSSplineArr.end(); ++splineIt) {
		fRhoGSSplineArr.push_back((TSpline3*) ((*splineIt)->Clone()));
	}
	for (splineIt = other.fRhoGSSSplineArr.begin(); splineIt != other.fRhoGSSSplineArr.end(); ++splineIt) {
		fRhoGSSSplineArr.push_back((TSpline3*) ((*splineIt)->Clone()));
	}
}

WCSimIntegralLookupHistArray& WCSimIntegralLookupHistArray::operator =(const WCSimIntegralLookupHistArray& rhs) {
	// Assign all the simple variables
	fEMin = rhs.fEMin;
	fEMax = rhs.fEMax;
	fR0Min = rhs.fR0Min;
	fR0Max = rhs.fR0Max;
	fR0Bins = rhs.fR0Bins;

	fCosTheta0Min = rhs.fCosTheta0Min;
	fCosTheta0Max = rhs.fCosTheta0Max;
	fCosTheta0Bins = rhs.fCosTheta0Bins;

	// We need to be a bit careful with the TGraphs so that this is safe against self-assignment
	// So make vectors with copies of all histograms from the vectors belonging to rhs
	// Then clear our own vectors and delete everything they hold
	// Then assign the copied vectors to be this class's new members
	std::vector<TGraph*>::const_iterator grIt;
	std::vector<TGraph*> tmpRhoGArr, tmpRhoGSArr, tmpRhoGSSArr; // Temporary home for the copies

	// Copy all the graphs into the temporary vectors
	for (grIt = rhs.fRhoGArr.begin(); grIt != rhs.fRhoGArr.end(); ++grIt) {
		tmpRhoGArr.push_back((TGraph*) ((*grIt)->Clone()));
	}
	for (grIt = rhs.fRhoGSArr.begin(); grIt != rhs.fRhoGSArr.end(); ++grIt) {
		tmpRhoGSArr.push_back((TGraph*) ((*grIt)->Clone()));
	}
	for (grIt = rhs.fRhoGSSArr.begin(); grIt != rhs.fRhoGSSArr.end(); ++grIt) {
		tmpRhoGSSArr.push_back((TGraph*) ((*grIt)->Clone()));
	}

	// Copy all the splines into the temporary vectors
	std::vector<TSpline3*>::const_iterator splineIt;
	std::vector<TSpline3*> tmpRhoGSplineArr, tmpRhoGSSplineArr, tmpRhoGSSSplineArr; // Temporary home for the copies
	for (splineIt = rhs.fRhoGSplineArr.begin(); splineIt != rhs.fRhoGSplineArr.end(); ++splineIt) {
		if ((*splineIt) != NULL) {
			tmpRhoGSplineArr.push_back((TSpline3*) ((*splineIt)->Clone()));
		} else {
			tmpRhoGSplineArr.push_back(NULL);
		}
	}
	for (splineIt = rhs.fRhoGSSplineArr.begin(); splineIt != rhs.fRhoGSSplineArr.end(); ++splineIt) {
		if ((*splineIt) != NULL) {
			tmpRhoGSSplineArr.push_back((TSpline3*) ((*splineIt)->Clone()));
		} else {
			tmpRhoGSSplineArr.push_back(NULL);
		}
	}
	for (splineIt = rhs.fRhoGSSSplineArr.begin(); splineIt != rhs.fRhoGSSSplineArr.end(); ++splineIt) {
		if ((*splineIt) != NULL) {
			tmpRhoGSSSplineArr.push_back((TSpline3*) ((*splineIt)->Clone()));
		} else {
			tmpRhoGSSSplineArr.push_back(NULL);
		}
	}

	// Clear out our own vectors
	for (grIt = fRhoGArr.begin(); grIt != fRhoGArr.end(); ++grIt) {
		if (*grIt != NULL) {
			delete *grIt;
		}
	}
	for (grIt = fRhoGSArr.begin(); grIt != fRhoGSArr.end(); ++grIt) {
		if (*grIt != NULL) {
			delete *grIt;
		}
	}
	for (grIt = fRhoGSSArr.begin(); grIt != fRhoGSSArr.end(); ++grIt) {
		if (*grIt != NULL) {
			delete *grIt;
		}
	}
	for (splineIt = fRhoGSplineArr.begin(); splineIt != fRhoGSplineArr.end(); ++splineIt) {
		if (*splineIt != NULL) {
			delete *splineIt;
		}
	}
	for (splineIt = fRhoGSSplineArr.begin(); splineIt != fRhoGSSplineArr.end(); ++splineIt) {
		if (*splineIt != NULL) {
			delete *splineIt;
		}
	}
	for (splineIt = fRhoGSSSplineArr.begin(); splineIt != fRhoGSSSplineArr.end(); ++splineIt) {
		if (*splineIt != NULL) {
			delete *splineIt;
		}
	}

	// Now assign the tempoary vectors to be our new members
	fRhoGArr = tmpRhoGArr;
	fRhoGSArr = tmpRhoGSArr;
	fRhoGSSArr = tmpRhoGSSArr;
	fRhoGSplineArr = tmpRhoGSplineArr;
	fRhoGSSplineArr = tmpRhoGSSplineArr;
	fRhoGSSSplineArr = tmpRhoGSSSplineArr;

	// Lastly we have to play the same trick with the raw TGraphs and splines
	TGraph * tmpRhoInt = (TGraph*) rhs.fRhoInt->Clone();
	TGraph * tmpRhoSInt = (TGraph*) rhs.fRhoSInt->Clone();
	TGraph * tmpRhoSSInt = (TGraph*) rhs.fRhoSSInt->Clone();
	TSpline3 * tmpRhoSSpline = NULL;
	if (rhs.fRhoSSpline != NULL) {
		tmpRhoSSpline = (TSpline3*) rhs.fRhoSSpline->Clone();
	}
	TSpline3 * tmpRhoSSSpline = NULL;
	if (rhs.fRhoSSSpline != NULL) {
		tmpRhoSSSpline = (TSpline3*) rhs.fRhoSSSpline->Clone();
	}
	if (fRhoInt != NULL) {
		delete fRhoInt;
	}
	if (fRhoSInt != NULL) {
		delete fRhoSInt;
	}
	if (fRhoSSInt != NULL) {
		delete fRhoSSInt;
	}
	if (fRhoSSpline != NULL) {
		delete fRhoSSpline;
	}
	if (fRhoSSSpline != NULL) {
		delete fRhoSSSpline;
	}
	fRhoInt = tmpRhoInt;
	fRhoSInt = tmpRhoSInt;
	fRhoSSInt = tmpRhoSSInt;
	fRhoSSpline = tmpRhoSSpline;
	fRhoSSSpline = tmpRhoSSSpline;

	return *this;
}

WCSimIntegralLookupHistArray::~WCSimIntegralLookupHistArray() {

	// Delete all our splines
	std::vector<TSpline3*>::iterator splineIt;
	for (splineIt = fRhoGSplineArr.begin(); splineIt != fRhoGSplineArr.end(); ++splineIt) {
		if (*splineIt != NULL) {
			delete *splineIt;
		}
	}
	fRhoGSplineArr.clear();
	for (splineIt = fRhoGSSplineArr.begin(); splineIt != fRhoGSSplineArr.end(); ++splineIt) {
		if (*splineIt != NULL) {
			delete *splineIt;
		}
	}
	fRhoGSSplineArr.clear();
	for (splineIt = fRhoGSSSplineArr.begin(); splineIt != fRhoGSSSplineArr.end(); ++splineIt) {
		if (*splineIt != NULL) {
			delete *splineIt;
		}
	}
	fRhoGSSSplineArr.clear();

	// Delete all our TGraphs
	std::vector<TGraph*>::iterator grIt;
	for (grIt = fRhoGArr.begin(); grIt != fRhoGArr.end(); ++grIt) {
		if (*grIt != NULL) {
			delete *grIt;
		}
	}
	fRhoGArr.clear();
	for (grIt = fRhoGSArr.begin(); grIt != fRhoGSArr.end(); ++grIt) {
		if (*grIt != NULL) {
			delete *grIt;
		}
	}
	fRhoGSArr.clear();
	for (grIt = fRhoGSSArr.begin(); grIt != fRhoGSSArr.end(); ++grIt) {
		if (*grIt != NULL) {
			delete *grIt;
		}
	}
	fRhoGSSArr.clear();

	if (fRhoInt != NULL) {
		delete fRhoInt;
		fRhoInt = NULL;
	}
	if (fRhoSInt != NULL) {
		delete fRhoSInt;
		fRhoSInt = NULL;
	}
	if (fRhoSSInt != NULL) {
		delete fRhoSSInt;
		fRhoSSInt = NULL;
	}
	if (fRhoSSpline != NULL) {
		delete fRhoSSpline;
		fRhoSSpline = NULL;
	}
	if (fRhoSSSpline != NULL) {
		delete fRhoSSSpline;
		fRhoSSSpline = NULL;
	}
}

void WCSimIntegralLookupHistArray::SetBins(double R0Max, double R0Min, int R0Bins, double cosTh0Max, double cosTh0Min,
		int cosTh0Bins) {
	fR0Max = R0Max;
	fR0Min = R0Min;
	fR0Bins = R0Bins;

	fCosTheta0Min = cosTh0Min;
	fCosTheta0Max = cosTh0Max;
	fCosTheta0Bins = cosTh0Bins;

	ResetArrays();
}

TGraph* WCSimIntegralLookupHistArray::GetRhoGInt(double R0, double cosTh0) const {
	return fRhoGArr.at(this->GetArrayIndex(R0, cosTh0));
}

TGraph* WCSimIntegralLookupHistArray::GetRhoGSInt(double R0, double cosTh0) const {
	return fRhoGSArr.at(GetArrayIndex(R0, cosTh0));
}

TGraph* WCSimIntegralLookupHistArray::GetRhoGSSInt(double R0, double cosTh0) const {
	return fRhoGSSArr.at(GetArrayIndex(R0, cosTh0));
}

TSpline3* WCSimIntegralLookupHistArray::GetRhoGSpline(double R0, double cosTh0) const {
	return fRhoGSplineArr.at(this->GetArrayIndex(R0, cosTh0));
}

TSpline3* WCSimIntegralLookupHistArray::GetRhoGSSpline(double R0, double cosTh0) const {
	return fRhoGSSplineArr.at(GetArrayIndex(R0, cosTh0));
}

TSpline3* WCSimIntegralLookupHistArray::GetRhoGSSSpline(double R0, double cosTh0) const {
	return fRhoGSSSplineArr.at(GetArrayIndex(R0, cosTh0));
}

void WCSimIntegralLookupHistArray::ResetArrays() {
	std::cout << " *** WCSimIntegralLookupHistArray::ResetArrays() *** " << std::endl;

	if (fRhoSInt != NULL) {
		delete fRhoSInt;
	}
	if (fRhoSSInt != NULL) {
		delete fRhoSSInt;
	}

	std::cout << " -- Making Rho graphs and splines -- " << std::endl;
	fRhoInt = new TGraph();
	fRhoInt->SetName("fRhoInt");
	fRhoSInt = new TGraph();
	fRhoSInt->SetName("fRhoSInt");
	fRhoSSInt = new TGraph();
	fRhoSSInt->SetName("fRhoSSInt");
	fRhoSSpline = new TSpline3();
	fRhoSSpline->SetName("fRhoSSpline");
	fRhoSSSpline = new TSpline3();
	fRhoSSSpline->SetName("fRhoSSSpline");

	int nBins = fR0Bins * fCosTheta0Bins;

	std::vector<TGraph*> * allGraphs[3] = { &fRhoGArr, &fRhoGSArr, &fRhoGSSArr };
	std::vector<TSpline3*> * allSplines[3] = { &fRhoGSplineArr, &fRhoGSSplineArr, &fRhoGSSSplineArr };

	for (int i = 0; i < 3; ++i) {
		std::cout << " -- Nulling TGraphs and Splines -- " << i << "/3" << std::endl;
		for (size_t j = 0; j < allGraphs[i]->size(); ++j) {
			TGraph * gr = allGraphs[i]->at(j);
			if (gr != NULL) {
				delete gr;
				gr = NULL;
			}
		}
		for (size_t j = 0; j < allSplines[i]->size(); ++j) {
			TSpline3 * spline = allSplines[i]->at(j);
			if (spline != NULL) {
				delete spline;
				spline = NULL;
			}
		}
	}

	std::cout << " -- Resizing arrays -- " << std::endl;
	fRhoGArr.resize(nBins, NULL);
	fRhoGSArr.resize(nBins, NULL);
	fRhoGSSArr.resize(nBins, NULL);
	fRhoGSplineArr.resize(nBins, NULL);
	fRhoGSSplineArr.resize(nBins, NULL);
	fRhoGSSSplineArr.resize(nBins, NULL);

	double R0Width = (fR0Max - fR0Min) / fR0Bins;
	double th0Width = (fCosTheta0Max - fCosTheta0Min) / fCosTheta0Bins;

	int lastPercent = -1;
	for (int iR = 0; iR < fR0Bins; ++iR) {
		double R0 = fR0Min + R0Width * iR;
		for (int iTh = 0; iTh < fCosTheta0Bins; ++iTh) {

			int bin = iR * fCosTheta0Bins + iTh;
			int percent = (iTh - 1) * (iR - 1) * 100 / nBins;
			if (percent > lastPercent) {
				std::cout << "\rMaking graphs, " << percent << "% complete" << std::flush;
				lastPercent = percent;
			}
			double cosTheta0 = fCosTheta0Min + iTh * th0Width;
			TString nameRhoG = Form("gRhoG_%f_%f", R0, cosTheta0);
			TString nameRhoGS = Form("gRhoGS_%f_%f", R0, cosTheta0);
			TString nameRhoGSS = Form("gRhoGSS_%f_%f", R0, cosTheta0);
			TString r0Title = Form("%d", (int) R0);
			fRhoGArr.at(bin) = new TGraph();
			fRhoGSArr.at(bin) = new TGraph();
			fRhoGSSArr.at(bin) = new TGraph();

			fRhoGArr.at(bin)->SetName(nameRhoG.Data());
			fRhoGSArr.at(bin)->SetName(nameRhoGS.Data());
			fRhoGSSArr.at(bin)->SetName(nameRhoGSS.Data());

			fRhoGArr.at(bin)->SetTitle(r0Title.Data());
			fRhoGSArr.at(bin)->SetTitle(r0Title.Data());
			fRhoGSSArr.at(bin)->SetTitle(r0Title.Data());

			// If you make a TSpline3 with the default (blank) constructor
			// and then don't put a graph into it, it segfaults when you try
			// to Draw()/Eval() etc. it so we'll only make the splines
			// once we've filled the graphs
		}
	}

	std::cout << "Done" << std::endl;
}

unsigned long int WCSimIntegralLookupHistArray::GetArrayIndex(double R0, double cosTheta0) const {

	double widthR = (fR0Max - fR0Min) / fR0Bins;
	double widthTh = (fCosTheta0Max - fCosTheta0Min) / fCosTheta0Bins;

	int R0Bin = (R0 - fR0Min) / widthR;
	int th0Bin = (cosTheta0 - fCosTheta0Min) / widthTh;
	return (R0Bin * fCosTheta0Bins) + th0Bin;
}

void WCSimIntegralLookupHistArray::ClearGraphs() {
	std::vector<TGraph*>::iterator grItr;
	std::vector<TGraph*> * arrays[3] = { &fRhoGArr, &fRhoGSArr, &fRhoGSSArr };
	for (int i = 0; i < 3; ++i) {
		for (grItr = arrays[i]->begin(); grItr != arrays[i]->end(); ++grItr) {
			delete *grItr;
			*grItr = NULL;
		}
		arrays[i]->clear();
	}
	return;

}

void WCSimIntegralLookupHistArray::Verify() {
	return;

	TFile * test = new TFile("test.root", "RECREATE");
	test->mkdir("gRhoG");
	test->mkdir("gRhoGS");
	test->mkdir("gRhoGSS");
	test->mkdir("splines");

	std::vector<TGraph*>::iterator grItr;
	std::vector<TGraph*> * arrays[3] = { &fRhoGArr, &fRhoGSArr, &fRhoGSSArr };

	for (int i = 0; i < 3; ++i) {
		for (grItr = arrays[i]->begin(); grItr != arrays[i]->end(); ++grItr) {
			if ((*grItr)->GetN() == 0) {
				continue;
			}
			TString name((*grItr)->GetName());
			if (name.Contains("gRhoG_")) {
				test->cd("gRhoG");
			} else if (name.Contains("gRhoGS_")) {
				test->cd("gRhoGS");
			} else if (name.Contains("gRhoGSS_")) {
				test->cd("gRhoGSS");
			}
			TDirectory * dir = (TDirectory*) gDirectory->Get((*grItr)->GetTitle());
			if (!dir)
				dir = gDirectory->mkdir((*grItr)->GetTitle());
			gDirectory->cd((*grItr)->GetTitle());

			(*grItr)->SetMarkerSize(1.0);
			(*grItr)->SetMarkerStyle(20);
			(*grItr)->SetMarkerColor(i + 1);
			(*grItr)->Write();
			test->cd();
		}
	}

	std::vector<TSpline3*>::iterator splineItr;
	std::vector<TSpline3*> * splineArrays[3] = { &fRhoGSplineArr, &fRhoGSSplineArr, &fRhoGSSSplineArr };
	for (int i = 0; i < 3; ++i) {
		test->cd("splines");
		for (splineItr = splineArrays[i]->begin(); splineItr != splineArrays[i]->end(); ++splineItr) {
			if ((*splineItr) == NULL) {
				continue;
			}
			(*splineItr)->SetLineColor(i + 1);
			std::cout << (*splineItr) << std::endl;
			(*splineItr)->Write();
		}
		test->cd();
	}
	if (fRhoInt != NULL) {
		fRhoInt->Write();
	}
	if (fRhoSInt != NULL) {
		fRhoSInt->Write();
	}
	if (fRhoSSInt != NULL) {
		fRhoSSInt->Write();
	}
	if (fRhoSSpline != NULL) {
		fRhoSSpline->Write();
	}
	if (fRhoSSSpline != NULL) {
		fRhoSSSpline->Write();
	}

	test->Close();

}

void WCSimIntegralLookup3D::SaveIntegrals(const double E, const double s, const double R0, const double cosTh0) {
	TCanvas * can = new TCanvas("can", "", 1500, 500);
	TGraph * rhoG = fIntegrals->GetRhoGInt(R0, cosTh0);
	TGraph * rhoGS = fIntegrals->GetRhoGSInt(R0, cosTh0);
	TGraph * rhoGSS = fIntegrals->GetRhoGSSInt(R0, cosTh0);
	TSpline3 * rhoGSpline = fIntegrals->GetRhoGSpline(R0, cosTh0);
	TSpline3 * rhoGSSpline = fIntegrals->GetRhoGSSpline(R0, cosTh0);
	TSpline3 * rhoGSSSpline = fIntegrals->GetRhoGSSSpline(R0, cosTh0);
	rhoG->SetMarkerSize(0.7);
	rhoGS->SetMarkerSize(0.7);
	rhoGSS->SetMarkerSize(0.7);
	can->Divide(3, 1);
	can->cd(1);
	rhoG->Draw("AP");
	rhoGSpline->Draw("SAME");
	can->cd(2);
	rhoGS->Draw("AP");
	rhoGSSpline->Draw("SAME");
	can->cd(3);
	rhoGSS->Draw("AP");
	rhoGSSSpline->Draw("SAME");
	can->cd(3);
	can->SaveAs(TString::Format("integrals_%.02f_%.04f.png", R0, cosTh0).Data());
	can->SaveAs(TString::Format("integrals_%.02f_%.04f.root", R0, cosTh0).Data());
	delete can;

}
