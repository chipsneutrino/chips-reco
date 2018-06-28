#include "WCSimHoughTransform.hh"

#include "TH2.h"
#include "TMath.h"
#include "TSpectrum2.h"

#include <cmath>

ClassImp (WCSimHoughTransform)

// Creator.  fHoughArray has one point for each (x,y) combination
// so make this the same size as # phi bins * # costheta bins
WCSimHoughTransform::WCSimHoughTransform(Int_t phiBins, Int_t costhetaBins) {
	fXbins = phiBins;
	fYbins = costhetaBins;
	fHoughArray = new Double_t[fXbins * fYbins];

	this->Reset();
}

// Destructor
WCSimHoughTransform::~WCSimHoughTransform() {
	delete[] fHoughArray;
}

// Set every entry in fHoughArray to 0
void WCSimHoughTransform::Reset() {
	for (Int_t xbin = 0; xbin < fXbins; xbin++) {
		for (Int_t ybin = 0; ybin < fYbins; ybin++) {
			fHoughArray[xbin * fYbins + ybin] = 0.0;
		}
	}

	fEntries = 0;
	fIntegral = 0.0;

	return;
}

// Add point to fHoughArray for a given point on a ring around a PMT hit 
// - convert into phi & costheta first
void WCSimHoughTransform::Fill(Double_t nx, Double_t ny, Double_t nz, Double_t weight) {
	Double_t phi = 0.0;
	Double_t costheta = 0.0;

	phi = TMath::ATan2(ny, nx);

	phi *= 180.0 / TMath::Pi();
	costheta = nz;

	this->Fill(phi, costheta, weight);
}

// Add entry to fHoughArray for given phi & costheta
void WCSimHoughTransform::Fill(Double_t phi, Double_t costheta, Double_t weight) {
	Int_t xbin = this->GetBinX(phi);
	Int_t ybin = this->GetBinY(costheta);

	if (xbin >= 0 && xbin < fXbins && ybin >= 0 && ybin < fYbins) {
		fHoughArray[xbin * fYbins + ybin] += weight;
		fEntries += 1;
		fIntegral += weight;
	}
}

// Loops over all bins to find the Hough peak.  Note (x,y) are not points in the
// detector in Cartesian space, but points (phi, costheta) in Hough space
void WCSimHoughTransform::FindPeak(Double_t& phi, Double_t& costheta, Double_t& height) {
	phi = 0.0;
	costheta = 0.0;
	// Handle empties
	if (fEntries <= 0 || fIntegral <= 0.0) {
		return;
	}

	Int_t xbinmax = -1;
	Int_t ybinmax = -1;
	Double_t max = 0.0;

	// Locate and record the peak of fHoughArray
	for (Int_t xbin = 0; xbin < fXbins; xbin++) {
		for (Int_t ybin = 0; ybin < fYbins; ybin++) {
			if (fHoughArray[xbin * fYbins + ybin] > max) {
				xbinmax = xbin;
				ybinmax = ybin;
				max = fHoughArray[xbin * fYbins + ybin];
			}
		}
	}

	// 5x5 sliding window to be moved around the peak to refine
	// its reconstruction based on neighbouring bins
	Int_t window = 2;  // window size either side of peak
	Double_t Swx = 0.0;
	Double_t Swy = 0.0;
	Double_t Sw = 0.0;
	Double_t Sn = 0.0;

	// Looks at each bin in a 5x5 window round the peak. Sums the contents of each of these
	// bins, weighting them by their (signed) distance from the peak in x or y, and divides by the
	// total contents of all 25 bins.  Converts this into units of phi or costheta and shifts the
	// peak by this much.  Hence no shift for a symmetric distribution, and a big shift for an
	// asymmetric one with high values far to one side of the peak.  Also averages the peak height.
	for (Int_t nx = -window; nx <= +window; nx++) {
		for (Int_t ny = -window; ny <= +window; ny++) {

			// fHoughArray squashes an m x n 2d array into a 1d array with m*n entries
			// This part handles the window sliding off the end of one row or column and loops
			// it back to the start
			Int_t xbin = xbinmax + nx;
			if (xbin < 0)
				xbin += fXbins;
			if (xbin >= fXbins)
				xbin -= fXbins;
			Int_t ybin = ybinmax + ny;

			if (xbin >= 0 && xbin < fXbins && ybin >= 0 && ybin < fYbins) {
				// Do the summations (phi goes from -180 to +180 and costheta from -1 to 1)
				Swx += 360.0 * (double(nx) / double(fXbins)) * fHoughArray[xbin * fYbins + ybin];
				Swy += 2.0 * (double(ny) / double(fYbins)) * fHoughArray[xbin * fYbins + ybin];
				Sw += fHoughArray[xbin * fYbins + ybin];
				Sn += 1.0;
			}

		}
	}
	// Shift the peak as above and scale its height to the average
	Double_t peakx = this->GetX(xbinmax) + Swx / Sw;
	Double_t peaky = this->GetY(ybinmax) + Swy / Sw;
	Double_t peaksize = Sw / Sn;

	if (peakx < -180.0)
		peakx += 360.0;
	if (peakx > +180.0)
		peakx -= 360.0;

	phi = peakx;
	costheta = peaky;
	height = peaksize;
}

// As above but called with, and converts back to, (x,y,z)
void WCSimHoughTransform::FindPeak(Double_t& hx, Double_t& hy, Double_t& hz, Double_t& height) {
	Double_t phi = 0.0;
	Double_t phiradians = 0.0;
	Double_t costheta = 0.0;
	Double_t sintheta = 0.0;
	Double_t tempheight = 0.0;

	this->FindPeak(phi, costheta, tempheight);

	phiradians = (TMath::Pi() / 180.0) * phi;
	sintheta = sqrt(1.0 - costheta * costheta);

	hx = sintheta * cos(phiradians);
	hy = sintheta * sin(phiradians);
	hz = costheta;
	height = tempheight;
}

// As above, overloaded to just accept and return the height
void WCSimHoughTransform::FindPeak(Double_t& height) {
	Double_t phi = 0.0;
	Double_t costheta = 0.0;
	Double_t tempheight = 0.0;

	this->FindPeak(phi, costheta, tempheight);

	height = tempheight;
}

// Histogram the phi-costheta Hough space that contains the peak
TH2D* WCSimHoughTransform::GetTH2D(const char* histname) {
	TH2D* hist = new TH2D(histname, "", fXbins, -180.0, +180.0, fYbins, -1.0, +1.0);

	for (Int_t nx = 0; nx < fXbins; nx++) {
		for (Int_t ny = 0; ny < fYbins; ny++) {
			hist->SetBinContent(nx + 1, ny + 1, fHoughArray[nx * fYbins + ny]);
		}
	}

	hist->SetEntries(fEntries);

	hist->GetXaxis()->SetTitle("#phi / degrees");
	hist->GetYaxis()->SetTitle("cos(#theta)");

	return hist;
}

// Shift histogram bins by a constant amount in the x-axis, representing a rotation in phi
TH2D* WCSimHoughTransform::GetRotatedTH2D(const char* histname, Double_t deltaPhi) {
	TH2D* hist = new TH2D(histname, "", fXbins, -180.0, +180.0, fYbins, -1.0, +1.0);
	Int_t deltaBinX = hist->GetXaxis()->FindBin(deltaPhi) - hist->GetXaxis()->FindBin(0.0);
	Int_t deltaBinY = 0;

	for (Int_t nx = 0; nx < fXbins; nx++) {
		Int_t nxbin = nx + deltaBinX;
		if (nxbin >= fXbins)
			nxbin -= fXbins;
		for (Int_t ny = 0; ny < fYbins; ny++) {
			Int_t nybin = ny + deltaBinY;
			if (nybin >= fYbins)
				nybin -= fYbins;
			hist->SetBinContent(nxbin + 1, nybin + 1, fHoughArray[nx * fYbins + ny]);
		}
	}

	hist->SetEntries(fEntries);

	hist->GetXaxis()->SetTitle("#phi / degrees");
	hist->GetYaxis()->SetTitle("cos(#theta)");

	return hist;
}
