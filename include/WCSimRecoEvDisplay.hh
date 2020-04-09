/*
 * WCSimRecoEvDisplay.hh
 *
 *  Created on: 13 Feb 2015
 *      Author: ajperch
 */

#pragma once

#include "WCSimEvDisplay.hh"
#include <vector>
class WCSimRecoSummary;
class TLine;
class TPad;
class TPolyMarker;
class TGWindow;
class TPaveText;
class TLegend;
class TH1D;

class WCSimRecoEvDisplay : public WCSimEvDisplay
{
public:
	WCSimRecoEvDisplay();
	WCSimRecoEvDisplay(const TGWindow *p, UInt_t w, UInt_t h);
	virtual ~WCSimRecoEvDisplay();

	// All "slot" commands, ie those that are called by buttons
	// need to be public.
	// Open the file to display events from
	void OpenFile();
	void OpenFile(std::string name);
	void OpenWCSimRecoFile(std::string name);

	// Button to show truth or reco
	void ShowFit();
	void ShowFitOverlay();
	void ShowTruthOverlay();

	// Functions to show the different likelihood contributions
	void SetViewTotalLnL();
	void SetViewQLnL();
	void SetViewTLnL();

	// Functions to show the different likelihood contributions from the correct track
	void SetViewCorrectTotalLnL();
	void SetViewCorrectQLnL();
	void SetViewCorrectTLnL();

	// Show the reco-true distributions
	void SetViewChargeRMT();
	void SetViewTimeRMT();

	// Show the predicted distributions at the best-fit
	void SetViewChargePrediction();
	void SetViewTimePrediction();

	// Show the predicted distributions for the correct tracks
	void SetViewCorrectChargePrediction();
	void SetViewCorrectTimePrediction();

	// Special display view for saving plots
	void ShowDisplayView();

protected:
	void MakeGraphColours();

	void FillPlots();
	void FillPlotsFromRecoFile();
	void FillPlotsFromLikelihood();
	void FillPlotsFromCorrectLikelihood();
	void FillPlotsFromRMT();
	void FillPlotsFromPrediction();
	void FillPlotsFromCorrectPrediction();
	// Draw the reco plots to their pads, but don't show yet.
	void UpdateRecoPads();
	// Update the truth TPaveText panel
	void UpdateFitPave();

	// Draw the fit information to its pad, but don't show yet.
	void UpdateFitPad();
	// Draw the fit overlay rings to their pads
	void UpdateFitOverlayPad();
	// The display pad will show both reco and true
	void UpdateDisplayPad();

	// Match the reco and true events
	void MatchRecoToTrue();
	void UpdateDisplayInfo();
	std::string GetTrackInfo(bool isReco, unsigned int trackNo);

	// Encapsulate some of the GUI generation
	void CreateDebugButtonBar();
	void CreateSubButtonBar();
	void CreateMainButtonBar();

	// Few functions to draw truth rings
	void DrawFitRing(unsigned int particleNo, int colour);
	void ClearFitMarkerVectors();
	int GetFitRingColour(int ring) const;
	void DrawFitOverlays();
	void HideFitOverlays();
	void UpdateCanvases();

	void ResizePads();
	TPad *fFitOverlayPad;
	TPad *fFitPad;
	TPad *fDisplayPad;
	void ResizeFitTexts(bool commonVertex);

	// Need to overload the SetPlotZAxes function
	void SetPlotZAxes();

	// Vectors of TPolyMarkers to store the truth rings for each detector region
	std::vector<TPolyMarker *> fFitMarkersTop;
	std::vector<TPolyMarker *> fFitMarkersBarrel;
	std::vector<TPolyMarker *> fFitMarkersBottom;
	std::vector<TLine *> fFitLines; // One line per truth ring to fill the legend
	TLegend *fFitLegend;

	WCSimRecoSummary *fRecoSummary;

	TChain *fRecoSummaryChain;

	// TChain for the HitComparison tree
	TChain *fHitComparisonChain;

	// The truth display is all contained within TPaveText objects
	TPaveText *fFitTextMain;
	TPaveText *fFitTextPrimaries;

	// Display view text
	TPaveText *fDisplayTextTrue;
	TPaveText *fDisplayTextReco;
	// Pairs of reconstructed and true tracks for the display view
	std::vector<std::pair<unsigned int, unsigned int>> fDisplayTrackPairs;
	TLegend *fDisplayLegend;

	// Needed for the likelihood plotting
	double fLnLMin;
	double fLnLMax;
	std::vector<double> fLnLBins; // Lower edges of the bins
	void CalculateLnLBins();
	unsigned int GetLnLBin(double lnl) const;

	// Needed for the reco-true plotting
	double fQRMTMin;
	double fQRMTMax;
	double fTRMTMin;
	double fTRMTMax;
	std::vector<double> fQRMTBins; // Lower edges of the bins
	std::vector<double> fTRMTBins; // Lower edges of the bins
	void CalculateRMTBins();
	unsigned int GetQRMTBin(double rmt) const;
	unsigned int GetTRMTBin(double rmt) const;
	TH1D *fChargeRMTHist;
	TH1D *fTimeRMTHist;
	TH1D *fChargePredHist;
	TH1D *fTimePredHist;
	TH1D *fCorrectChargePredHist;
	TH1D *fCorrectTimePredHist;

	std::vector<Int_t> fColoursNormal;
	std::vector<Int_t> fColoursDiff;

	ClassDef(WCSimRecoEvDisplay, 0)
};
