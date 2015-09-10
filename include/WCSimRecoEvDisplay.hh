/*
 * WCSimRecoEvDisplay.hh
 *
 *  Created on: 13 Feb 2015
 *      Author: ajperch
 */

#ifndef WCSIMRECOEVDISPLAY_HH_
#define WCSIMRECOEVDISPLAY_HH_

#include "WCSimEvDisplay.hh"
#include <vector>
class WCSimRecoSummary;
class TLine;
class TPad;
class TPolyMarker;
class TGWindow;
class TPaveText;
class TLegend;


class WCSimRecoEvDisplay: public WCSimEvDisplay {
public:
	WCSimRecoEvDisplay();
	WCSimRecoEvDisplay(const TGWindow *p,UInt_t w,UInt_t h);
	virtual ~WCSimRecoEvDisplay();


	// All "slot" commands, ie those that are called by buttons
	// need to be public.
	// Open the file to display events from
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

  // Show the reco-true distributions
  void SetViewChargeRMT();
  void SetViewTimeRMT();

protected:
	
	void FillPlots();
	void FillPlotsFromRecoFile();
  void FillPlotsFromLikelihood();
  void FillPlotsFromRMT();
  // Draw the reco plots to their pads, but don't show yet.
  void UpdateRecoPads();
  // Update the truth TPaveText panel
	void UpdateFitPave();
	
  // Draw the fit information to its pad, but don't show yet.
	void UpdateFitPad();
	// Draw the fit overlay rings to their pads
	void UpdateFitOverlayPad();

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
	TPad * fFitOverlayPad;
	TPad * fFitPad;

  // Need to overload the SetPlotZAxes function
  void SetPlotZAxes();

	// Vectors of TPolyMarkers to store the truth rings for each detector region
	std::vector<TPolyMarker*> fFitMarkersTop;
	std::vector<TPolyMarker*> fFitMarkersBarrel;
	std::vector<TPolyMarker*> fFitMarkersBottom;
	std::vector<TLine*> fFitLines; // One line per truth ring to fill the legend
	TLegend * fFitLegend;

	WCSimRecoSummary * fRecoSummary;

	TChain * fRecoSummaryChain;

  // TChain for the HitComparison tree
  TChain *fHitComparisonChain;

	// The truth display is all contained within TPaveText objects
	TPaveText *fFitTextMain;
	TPaveText *fFitTextPrimaries;

  // Needed for the likelihood plotting
  double fLnLMin;
  double fLnLMax;
  std::vector<double> fLnLBins; // Lower edges of the bins
  void CalculateLnLBins();
  unsigned int GetLnLBin(double lnl) const;

  // Needed for the reco-true plotting
  double fRMTMin;
  double fRMTMax;
  std::vector<double> fRMTBins; // Lower edges of the bins
  void CalculateRMTBins();
  unsigned int GetRMTBin(double rmt) const;


	ClassDef(WCSimRecoEvDisplay,0)
};

#endif /* WCSIMRECOEVDISPLAY_HH_ */
