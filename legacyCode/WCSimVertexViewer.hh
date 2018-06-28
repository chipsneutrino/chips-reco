#ifndef WCSIMVERTEXVIEWER_HH
#define WCSIMVERTEXVIEWER_HH

#include "WCSimDisplay.hh"

#include <vector>

#include "TString.h"
#include "TLatex.h"

class TCanvas;
class TPad;
class TH1D;
class TH2D;
class TMarker;
class TLegend;
class TF1;

class WCSimRecoVertex;
class WCSimRecoEvent;
class WCSimTrueEvent;

class WCSimVertexViewer: public WCSimDisplay {
	public:
		WCSimVertexViewer();
		~WCSimVertexViewer();

		void DrawDisplay(WCSimRecoEvent*) {
		}
		void DrawCleanDisplay(WCSimRecoEvent*) {
		}
		void DrawTrueEvent(WCSimTrueEvent*) {
		}

		void DrawRecoEvent(WCSimRecoEvent* recoevent);
		void DrawNewRecoEvent(WCSimRecoEvent* recoevent, WCSimRecoVertex* vertex = 0);

		void ResetDisplay();
		void PrintDisplay();

	private:
		void Initialize();
		void Reset();

		void BuildDisplay();
		void DrawNewDisplay();

		void PrintDisplayGIF();
		void PrintDisplayEPS();

		Int_t QPEtoCOL(Double_t Q);
		Double_t QPEtoSIZE(Double_t Q);

		Int_t fRbins;
		Double_t fRmin;
		Double_t fRmax;
		Double_t fNewRmin;
		Double_t fNewRmax;

		Int_t fLbins;
		Double_t fLmin;
		Double_t fLmax;
		Double_t fNewLmin;
		Double_t fNewLmax;

		Int_t fTbins;
		Double_t fTmin;
		Double_t fTmax;
		Double_t fNewTmin;
		Double_t fNewTmax;

		TString fTitleString;
		TLatex fTitleLatex;

		// run/event numbers
		Int_t fRunNumber;
		Int_t fEventNumber;
		Int_t fTriggerNumber;

		// event displays
		TCanvas* fVertexCanvas2D;
		TH2D* fVertexHistTimeFit;
		TMarker* fVertexMarkerDigits;
		TMarker* fVertexMarkerFilterDigits;
		TF1* fVertexTimeFit;
		TLegend* fVertexLegend2D;

		TCanvas* fVertexCanvas1D;
		TH1D* fVertexHistDigitsAll;
		TH1D* fVertexHistDigitsCone;
		TH1D* fVertexHistFilterDigitsCone;
		TLegend* fVertexLegend1D;

		TCanvas* fTimeCanvas;
		TH1D* fTimeHistDigitsAll;
		TH1D* fTimeHistDigitsCone;
		TH1D* fTimeHistFilterDigitsCone;
		TLegend* fTimeLegend;

		TCanvas* fTransCanvas;
		TH1D* fTransHistDigits;
		TH1D* fTransHistFilterDigits;
		TLegend* fTransLegend;

		TCanvas* fLongCanvas;
		TH1D* fLongHistDigits;
		TH1D* fLongHistFilterDigits;
		TLegend* fLongLegend;

		// list of digits
		std::vector<TMarker*> wcVertexDigits;
		std::vector<TMarker*> wcVertexFilterDigits;

		ClassDef(WCSimVertexViewer,0)

};

#endif

