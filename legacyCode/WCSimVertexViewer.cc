#include "WCSimVertexViewer.hh"
#include "WCSimVertexGeometry.hh"

#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TMath.h"
#include "TF1.h"

#include "WCSimRecoEvent.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoVertex.hh"

#include "WCSimParameters.hh"

#include <iostream>
#include <cmath>
#include <cassert>

ClassImp (WCSimVertexViewer)

WCSimVertexViewer::WCSimVertexViewer() :
		WCSimDisplay(), fVertexCanvas2D(0), fVertexHistTimeFit(0), fVertexMarkerDigits(0), fVertexMarkerFilterDigits(0), fVertexTimeFit(
				0), fVertexLegend2D(0), fVertexCanvas1D(0), fVertexHistDigitsAll(0), fVertexHistDigitsCone(0), fVertexHistFilterDigitsCone(
				0), fVertexLegend1D(0), fTimeCanvas(0), fTimeHistDigitsAll(0), fTimeHistDigitsCone(0), fTimeHistFilterDigitsCone(
				0), fTimeLegend(0), fTransCanvas(0), fTransHistDigits(0), fTransHistFilterDigits(0), fTransLegend(0), fLongCanvas(
				0), fLongHistDigits(0), fLongHistFilterDigits(0), fLongLegend(0) {
	std::cout << " *** WCSimVertexViewer::WCSimVertexViewer() *** " << std::endl;

	this->Initialize();
}

void WCSimVertexViewer::Initialize() {
	std::cout << " *** WCSimVertexViewer::Initialize() *** " << std::endl;

	this->BuildDisplay();
}

WCSimVertexViewer::~WCSimVertexViewer() {

}

void WCSimVertexViewer::Reset() {
	std::cout << " *** WCSimVertexViewer::Reset() *** " << std::endl;

	if (fVertexCanvas2D)
		delete fVertexCanvas2D;
	if (fVertexHistTimeFit)
		delete fVertexHistTimeFit;
	if (fVertexLegend2D)
		delete fVertexLegend2D;

	if (fVertexCanvas1D)
		delete fVertexCanvas1D;
	if (fVertexHistDigitsAll)
		delete fVertexHistDigitsAll;
	if (fVertexHistDigitsCone)
		delete fVertexHistDigitsCone;
	if (fVertexHistFilterDigitsCone)
		delete fVertexHistFilterDigitsCone;
	if (fVertexTimeFit)
		delete fVertexTimeFit;
	if (fVertexLegend1D)
		delete fVertexLegend1D;

	if (fTimeCanvas)
		delete fTimeCanvas;
	if (fTimeHistDigitsAll)
		delete fTimeHistDigitsAll;
	if (fTimeHistDigitsCone)
		delete fTimeHistDigitsCone;
	if (fTimeHistFilterDigitsCone)
		delete fTimeHistFilterDigitsCone;
	if (fTimeLegend)
		delete fTimeLegend;

	if (fTransCanvas)
		delete fTransCanvas;
	if (fTransHistDigits)
		delete fTransHistDigits;
	if (fTransHistFilterDigits)
		delete fTransHistFilterDigits;
	if (fTransLegend)
		delete fTransLegend;

	if (fLongCanvas)
		delete fLongCanvas;
	if (fLongHistDigits)
		delete fLongHistDigits;
	if (fLongHistFilterDigits)
		delete fLongHistFilterDigits;
	if (fLongLegend)
		delete fLongLegend;

	return;
}

void WCSimVertexViewer::BuildDisplay() {
	std::cout << " *** WCSimVertexViewer::BuildDisplay() *** " << std::endl;

	// common binning
	// ==============
	fRbins = 300;
	fRmin = -3000.0;
	fRmax = +33000.0;
	fNewRmin = fRmin;
	fNewRmax = fRmax;

	fLbins = 240;
	fLmin = -1000.0;
	fLmax = +11000.0;
	fNewLmin = fLmin;
	fNewLmax = fLmax;

	fTbins = 720;
	fTmin = -300.0;
	fTmax = +3300.0;
	fNewTmin = fTmin;
	fNewTmax = fTmax;

	// speed of light
	// ==============
	Double_t fC = WCSimParameters::SpeedOfLight();

	// set up 1D display
	// =================
	fVertexCanvas1D = new TCanvas("wcsim_vertex_timefit_1d", "WCSim Vertex TimeFit Display (1D)", 750, 800);
	fVertexCanvas1D->SetLeftMargin(0.2);

	fVertexHistDigitsAll = new TH1D("hVertexTimeDigitsAll", "", 60, -20.0, +40.0);
	fVertexHistDigitsAll->GetXaxis()->SetTitle(" Residual / ns");
	fVertexHistDigitsAll->GetYaxis()->SetTitle(" Digits ");
	fVertexHistDigitsAll->GetXaxis()->SetLabelSize(0.05);
	fVertexHistDigitsAll->GetYaxis()->SetLabelSize(0.05);
	fVertexHistDigitsAll->GetXaxis()->SetTitleSize(0.06);
	fVertexHistDigitsAll->GetYaxis()->SetTitleSize(0.06);
	fVertexHistDigitsAll->GetYaxis()->SetTitleOffset(1.5);
	fVertexHistDigitsAll->GetXaxis()->CenterTitle();
	fVertexHistDigitsAll->GetYaxis()->CenterTitle();
	fVertexHistDigitsAll->SetLineWidth(2);
	fVertexHistDigitsAll->SetLineStyle(7);

	fVertexHistDigitsCone = new TH1D("hVertexTimeDigitsCone", "", 60, -20.0, +40.0);
	fVertexHistDigitsCone->SetLineWidth(2);

	fVertexHistFilterDigitsCone = new TH1D("hVertexTimeFilterDigitsCone", "", 60, -20.0, +40.0);
	fVertexHistFilterDigitsCone->SetLineWidth(1);
	fVertexHistFilterDigitsCone->SetLineColor(kBlue + 1);
	fVertexHistFilterDigitsCone->SetFillColor(kBlue + 1);

	fVertexLegend1D = new TLegend(0.51, 0.74, 0.92, 0.88);
	fVertexLegend1D->AddEntry(fVertexHistDigitsAll, "All Digits", "l");
	fVertexLegend1D->AddEntry(fVertexHistDigitsCone, "All Digits (Cone)", "l");
	fVertexLegend1D->AddEntry(fVertexHistFilterDigitsCone, "Filtered Digits (Cone)", "lf");
	fVertexLegend1D->SetFillColor(kWhite);
	fVertexLegend1D->SetFillStyle(0);
	fVertexLegend1D->SetLineColor(kWhite);

	// set up 2D display
	// =================
	fVertexCanvas2D = new TCanvas("wcsim_vertex_timefit_2d", "WCSim Vertex TimeFit Display (2D)", 900, 800);

	fVertexHistTimeFit = new TH2D("hVertexHistTimeFit", "", fRbins, fRmin, fRmax, fTbins, fTmin, fTmax);

	fVertexHistTimeFit->GetXaxis()->SetTitle(" Optical Path From Vertex / cm");
	fVertexHistTimeFit->GetYaxis()->SetTitle(" Time From Vertex / ns");
	fVertexHistTimeFit->GetXaxis()->SetLabelSize(0.05);
	fVertexHistTimeFit->GetYaxis()->SetLabelSize(0.05);
	fVertexHistTimeFit->GetXaxis()->SetTitleSize(0.06);
	fVertexHistTimeFit->GetYaxis()->SetTitleSize(0.06);
	fVertexHistTimeFit->GetYaxis()->SetTitleOffset(1.4);
	fVertexHistTimeFit->GetXaxis()->CenterTitle();
	fVertexHistTimeFit->GetYaxis()->CenterTitle();

	fVertexTimeFit = new TF1("VertexTimeFit", "[0]*x", -1000.0, +10000.0);
	fVertexTimeFit->SetParameter(0, 1.0 / fC);
	fVertexTimeFit->SetLineWidth(2);
	fVertexTimeFit->SetLineColor(kRed);

	fVertexMarkerDigits = new TMarker(-99999.9, -99999.9, 4);
	fVertexMarkerFilterDigits = new TMarker(-99999.9, -99999.9, 8);

	fVertexMarkerDigits->SetMarkerColor(kGray);
	fVertexMarkerDigits->SetMarkerSize(2.0);

	fVertexMarkerFilterDigits->SetMarkerColor(kBlack);
	fVertexMarkerFilterDigits->SetMarkerSize(2.0);

	fVertexLegend2D = new TLegend(0.21, 0.73, 0.56, 0.88);
	fVertexLegend2D->AddEntry(fVertexMarkerDigits, "All Digits", "p");
	fVertexLegend2D->AddEntry(fVertexMarkerFilterDigits, "Filtered Digits", "p");
	fVertexLegend2D->AddEntry(fVertexTimeFit, "Timing Fit", "l");
	fVertexLegend2D->SetLineColor(kWhite);
	fVertexLegend2D->SetFillColor(kWhite);
	fVertexLegend2D->SetFillStyle(0);

	// set up time display
	// ====================
	fTimeCanvas = new TCanvas("wcsim_vertex_timeprofile", "WCSim Vertex TimeProfile Display", 1050, 800);
	fTimeCanvas->SetLeftMargin(0.2);

	fTimeHistDigitsAll = new TH1D("hTimeDigitsAll", "", fTbins, fTmin, fTmax);
	fTimeHistDigitsAll->GetXaxis()->SetTitle(" Time / ns");
	fTimeHistDigitsAll->GetYaxis()->SetTitle(" Digits ");
	fTimeHistDigitsAll->GetXaxis()->SetLabelSize(0.05);
	fTimeHistDigitsAll->GetYaxis()->SetLabelSize(0.05);
	fTimeHistDigitsAll->GetXaxis()->SetTitleSize(0.06);
	fTimeHistDigitsAll->GetYaxis()->SetTitleSize(0.06);
	fTimeHistDigitsAll->GetYaxis()->SetTitleOffset(1.1);
	fTimeHistDigitsAll->GetXaxis()->CenterTitle();
	fTimeHistDigitsAll->GetYaxis()->CenterTitle();
	fTimeHistDigitsAll->SetLineWidth(2);
	fTimeHistDigitsAll->SetLineStyle(7);

	fTimeHistDigitsCone = new TH1D("hTimeDigitsCone", "", fTbins, fTmin, fTmax);
	fTimeHistDigitsCone->SetLineWidth(2);

	fTimeHistFilterDigitsCone = new TH1D("hTimeFilterDigitsCone", "", fTbins, fTmin, fTmax);
	fTimeHistFilterDigitsCone->SetLineWidth(1);
	fTimeHistFilterDigitsCone->SetLineColor(kBlue + 1);
	fTimeHistFilterDigitsCone->SetFillColor(kBlue + 1);

	fTimeLegend = new TLegend(0.51, 0.74, 0.92, 0.88);
	fTimeLegend->AddEntry(fTimeHistDigitsAll, "All Digits", "l");
	fTimeLegend->AddEntry(fTimeHistDigitsCone, "All Digits (Cone)", "l");
	fTimeLegend->AddEntry(fTimeHistFilterDigitsCone, "Filtered Digits (Cone)", "lf");
	fTimeLegend->SetFillColor(kWhite);
	fTimeLegend->SetFillStyle(0);
	fTimeLegend->SetLineColor(kWhite);

	// set up transverse display
	// =========================
	fTransCanvas = new TCanvas("wcsim_transverse_display", "WCSim Transverse Display", 750, 800);
	fTransCanvas->SetLeftMargin(0.2);

	fTransHistDigits = new TH1D("hTransHistDigits ", "", 65, 0, 130.0);
	fTransHistDigits->GetXaxis()->SetTitle(" Angle / Degrees ");
	fTransHistDigits->GetYaxis()->SetTitle(" Digits ");
	fTransHistDigits->GetXaxis()->SetLabelSize(0.05);
	fTransHistDigits->GetYaxis()->SetLabelSize(0.05);
	fTransHistDigits->GetXaxis()->SetTitleSize(0.06);
	fTransHistDigits->GetYaxis()->SetTitleSize(0.06);
	fTransHistDigits->GetYaxis()->SetTitleOffset(1.5);
	fTransHistDigits->GetXaxis()->CenterTitle();
	fTransHistDigits->GetYaxis()->CenterTitle();
	fTransHistDigits->SetLineWidth(2);

	fTransHistFilterDigits = new TH1D("hTransHistFilterDigits", "", 65, 0, 130.0);
	fTransHistFilterDigits->SetLineWidth(1);
	fTransHistFilterDigits->SetLineColor(kBlue + 1);
	fTransHistFilterDigits->SetFillColor(kBlue + 1);

	fTransLegend = new TLegend(0.51, 0.78, 0.88, 0.88);
	fTransLegend->AddEntry(fTransHistDigits, "All Digits", "l");
	fTransLegend->AddEntry(fTransHistFilterDigits, "Filtered Digits", "lf");
	fTransLegend->SetFillColor(kWhite);
	fTransLegend->SetFillStyle(0);
	fTransLegend->SetLineColor(kWhite);

	// set up longitudinal display
	// ===========================
	fLongCanvas = new TCanvas("wcsim_longitudinal_display", "WCSim Longitudinal Display", 1050, 800);

	fLongHistDigits = new TH1D("hLongHistDigits ", "", fLbins, fLmin, fLmax);
	fLongHistDigits->GetXaxis()->SetTitle(" Distance Along Track / cm ");
	fLongHistDigits->GetYaxis()->SetTitle(" Digits ");
	fLongHistDigits->GetXaxis()->SetLabelSize(0.05);
	fLongHistDigits->GetYaxis()->SetLabelSize(0.05);
	fLongHistDigits->GetXaxis()->SetTitleSize(0.06);
	fLongHistDigits->GetYaxis()->SetTitleSize(0.06);
	fLongHistDigits->GetYaxis()->SetTitleOffset(1.1);
	fLongHistDigits->GetXaxis()->CenterTitle();
	fLongHistDigits->GetYaxis()->CenterTitle();
	fLongHistDigits->SetLineWidth(2);

	fLongHistFilterDigits = new TH1D("hLongHistFilterDigits", "", fLbins, fLmin, fLmax);
	fLongHistFilterDigits->SetLineWidth(1);
	fLongHistFilterDigits->SetLineColor(kBlue + 1);
	fLongHistFilterDigits->SetFillColor(kBlue + 1);

	fLongLegend = new TLegend(0.51, 0.78, 0.88, 0.88);
	fLongLegend->AddEntry(fLongHistDigits, "All Digits", "l");
	fLongLegend->AddEntry(fLongHistFilterDigits, "Filtered Digits", "lf");
	fLongLegend->SetFillColor(kWhite);
	fLongLegend->SetFillStyle(0);
	fLongLegend->SetLineColor(kWhite);

	// Draw everything
	// ===============
	this->DrawNewDisplay();

	return;
}

void WCSimVertexViewer::DrawRecoEvent(WCSimRecoEvent* myEvent) {
	// check event
	// ===========
	if (myEvent == 0) {
		return;
	}

	// draw display
	// ============
	this->DrawNewRecoEvent(myEvent);

}

void WCSimVertexViewer::ResetDisplay() {
	// all digits
	// ==========
	for (UInt_t i = 0; i < wcVertexDigits.size(); i++) {
		delete (TMarker*) (wcVertexDigits.at(i));
	}

	wcVertexDigits.clear();

	// filtered digits
	// ===============
	for (UInt_t i = 0; i < wcVertexFilterDigits.size(); i++) {
		delete (TMarker*) (wcVertexFilterDigits.at(i));
	}

	wcVertexFilterDigits.clear();

	// reset histograms
	// ================
	if (fVertexHistDigitsAll)
		fVertexHistDigitsAll->Reset();
	if (fVertexHistDigitsCone)
		fVertexHistDigitsCone->Reset();
	if (fVertexHistFilterDigitsCone)
		fVertexHistFilterDigitsCone->Reset();
	if (fTimeHistDigitsAll)
		fTimeHistDigitsAll->Reset();
	if (fTimeHistDigitsCone)
		fTimeHistDigitsCone->Reset();
	if (fTimeHistFilterDigitsCone)
		fTimeHistFilterDigitsCone->Reset();
	if (fTransHistDigits)
		fTransHistDigits->Reset();
	if (fTransHistFilterDigits)
		fTransHistFilterDigits->Reset();
	if (fLongHistDigits)
		fLongHistDigits->Reset();
	if (fLongHistFilterDigits)
		fLongHistFilterDigits->Reset();

	// reset event info
	// ================
	fRunNumber = -1;
	fEventNumber = -1;
	fTriggerNumber = -1;

	// draw new display
	// ================
	this->DrawNewDisplay();

	return;
}

void WCSimVertexViewer::DrawNewDisplay() {
	// Generate Title String
	// =====================
	fTitleString = "";

	if (fRunNumber >= 0) {
		fTitleString.Append("Run: ");
		fTitleString += fRunNumber;
		fTitleString.Append("  ");
		fTitleString.Append("Event: ");
		fTitleString += fEventNumber;
		fTitleString.Append("  ");
		fTitleString.Append("SubEvent: ");
		fTitleString += fTriggerNumber;

		fTitleLatex.SetTextFont(42);
		fTitleLatex.SetTextSize(0.04);
		fTitleLatex.SetTextAlign(11);
	}

	// Draw 1D Display
	// ===============
	if (fVertexCanvas1D) {
		fVertexCanvas1D->cd();
		fVertexCanvas1D->Clear();
		fVertexHistDigitsAll->Draw();
		fVertexHistDigitsCone->Draw("same");
		fVertexHistFilterDigitsCone->Draw("same");
		fVertexLegend1D->Draw();
		fVertexCanvas1D->Update();
	}

	// Draw 2D Display
	// ===============
	if (fVertexCanvas2D) {
		fVertexCanvas2D->cd();
		fVertexHistTimeFit->Draw();
		fVertexLegend2D->Draw();
		fVertexCanvas2D->Update();
	}

	// Draw TimeProfile Display
	// ========================
	if (fTimeCanvas) {
		fTimeCanvas->cd();
		fTimeCanvas->Clear();
		fTimeHistDigitsAll->Draw();
		fTimeHistDigitsCone->Draw("same");
		fTimeHistFilterDigitsCone->Draw("same");
		fTimeLegend->Draw();
		fTimeCanvas->Update();
	}

	// Draw Transverse Display
	// =======================
	if (fTransCanvas) {
		fTransCanvas->cd();
		fTransHistDigits->Draw();
		fTransHistFilterDigits->Draw("same");
		fTransLegend->Draw();
		fTransCanvas->Update();
	}

	// Draw Longitidinal Display
	// =========================
	if (fLongCanvas) {
		fLongCanvas->cd();
		fLongHistDigits->Draw();
		fLongHistFilterDigits->Draw("same");
		fLongLegend->Draw();
		fLongCanvas->Update();
	}

	return;
}

void WCSimVertexViewer::PrintDisplay() {
	// print GIF
	if (fPrintGIF)
		this->PrintDisplayGIF();

	// print EPS
	if (fPrintEPS)
		this->PrintDisplayEPS();

	return;
}

void WCSimVertexViewer::PrintDisplayGIF() {
	if (fVertexCanvas1D) {
		TString outfile("wcresiduals");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".gif");
		fVertexCanvas1D->SaveAs(outfile.Data());
	}

	if (fVertexCanvas2D) {
		TString outfile("wctimefit");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".gif");
		fVertexCanvas2D->SaveAs(outfile.Data());
	}

	if (fTimeCanvas) {
		TString outfile("wctimeprofile");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".gif");
		fTimeCanvas->SaveAs(outfile.Data());
	}

	if (fTransCanvas) {
		TString outfile("wctransverse");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".gif");
		fTransCanvas->SaveAs(outfile.Data());
	}

	if (fLongCanvas) {
		TString outfile("wclongitudinal");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".gif");
		fLongCanvas->SaveAs(outfile.Data());
	}

	return;
}

void WCSimVertexViewer::PrintDisplayEPS() {
	if (fVertexCanvas1D) {
		TString outfile("wcresiduals");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".eps");
		fVertexCanvas1D->SaveAs(outfile.Data());
	}

	if (fVertexCanvas2D) {
		TString outfile("wctimefit");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".eps");
		fVertexCanvas2D->SaveAs(outfile.Data());
	}

	if (fTimeCanvas) {
		TString outfile("wctimeprofile");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".eps");
		fTimeCanvas->SaveAs(outfile.Data());
	}

	if (fTransCanvas) {
		TString outfile("wctransverse");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".eps");
		fTransCanvas->SaveAs(outfile.Data());
	}

	if (fLongCanvas) {
		TString outfile("wclongitudinal");
		outfile.Append(".");
		outfile += fRunNumber;
		outfile.Append(".");
		outfile += fEventNumber;
		outfile.Append(".eps");
		fLongCanvas->SaveAs(outfile.Data());
	}

	return;
}

void WCSimVertexViewer::DrawNewRecoEvent(WCSimRecoEvent* myEvent, WCSimRecoVertex* myVertex) {
	std::cout << " *** WCSimVertexViewer::DrawEvent(...) *** " << std::endl;

	// reset display
	// =============
	this->ResetDisplay();

	// Check for Event and Digits
	// ==========================
	if (myEvent == 0 || myEvent->GetNDigits() == 0)
		return;

	// Check for Vertex
	// ================
	WCSimRecoVertex* inputVertex = 0;

	if (myVertex) {
		inputVertex = myVertex;
	} else if (myEvent->FoundVertex()) {
		inputVertex = myEvent->GetVertex();
	}

	if (inputVertex == 0) {
		std::cout << " <warning>: running vertex viewer with no vertex! " << std::endl;
		return;
	}

	// Get Event Info
	// ==============
	fRunNumber = myEvent->GetRun();
	fEventNumber = myEvent->GetEvent();
	fTriggerNumber = myEvent->GetTrigger();

	fTitleString = "";

	fTitleString.Append("Run: ");
	fTitleString += fRunNumber;
	fTitleString.Append("  ");
	fTitleString.Append("Event: ");
	fTitleString += fEventNumber;
	fTitleString.Append("  ");
	fTitleString.Append("SubEvent: ");
	fTitleString += fTriggerNumber;

	fTitleLatex.SetTextFont(42);
	fTitleLatex.SetTextSize(0.04);
	fTitleLatex.SetTextAlign(11);

	// Reco Vertex
	// ===========
	Double_t vtxX = inputVertex->GetX();
	Double_t vtxY = inputVertex->GetY();
	Double_t vtxZ = inputVertex->GetZ();
	Double_t vtxTime = inputVertex->GetTime();

	Double_t dirX = inputVertex->GetDirX();
	Double_t dirY = inputVertex->GetDirY();
	Double_t dirZ = inputVertex->GetDirZ();

	std::cout << "  Displaying vertex: (vx,vy,vz,vt)=(" << vtxX << "," << vtxY << "," << vtxZ << "," << vtxTime << ") "
			<< std::endl;
	std::cout << "                     (px,py,pz,th)=(" << dirX << "," << dirY << "," << dirZ << ") " << std::endl;

	// Useful Variables
	// ================
	Double_t min = 0.0;
	Double_t max = 0.0;
	Double_t hmax = 0.0;

	Double_t fNewRmin = fRmax;
	Double_t fNewRmax = fRmin;
	Double_t fNewLmin = fLmax;
	Double_t fNewLmax = fLmin;
	Double_t fNewTmin = fTmax;
	Double_t fNewTmax = fTmin;

	TMarker* marker = 0;

	// Calculate Residuals
	// ===================
	WCSimVertexGeometry::Instance()->CalcResiduals(myEvent, inputVertex);

	// Loop over Digits
	// ================
	for (Int_t idigit = 0; idigit < WCSimVertexGeometry::Instance()->GetNDigits(); idigit++) {

		Bool_t Pass = WCSimVertexGeometry::Instance()->IsFiltered(idigit);
		Double_t Q = WCSimVertexGeometry::Instance()->GetDigitQ(idigit);

		if (Q < GetPulseHeightCut())
			continue;

		Double_t thetadeg = WCSimVertexGeometry::Instance()->GetConeAngle(idigit);
		Double_t phideg = WCSimVertexGeometry::Instance()->GetAngle(idigit);
		Double_t dt = WCSimVertexGeometry::Instance()->GetDeltaTime(idigit);
		Double_t Ltrack = WCSimVertexGeometry::Instance()->GetDistTrack(idigit);
		Double_t delta = WCSimVertexGeometry::Instance()->GetExtendedResidual(idigit);
		Double_t path = WCSimVertexGeometry::Instance()->GetExtendedPath(idigit);

		// all digits
		// ==========
		marker = new TMarker(path, dt, 4);
		marker->SetMarkerSize(1.0);
		marker->SetMarkerColor(kGray);

		fTransHistDigits->Fill(phideg);
		fVertexHistDigitsAll->Fill(delta);
		fTimeHistDigitsAll->Fill(dt);
		if (phideg < thetadeg) {
			fLongHistDigits->Fill(Ltrack);
			fVertexHistDigitsCone->Fill(delta);
			fTimeHistDigitsCone->Fill(dt);
		}

		wcVertexDigits.push_back(marker);

		// filter digits
		// =============
		if (Pass) {
			marker = new TMarker(path, dt, 8);
			marker->SetMarkerSize(1.0);
			marker->SetMarkerColor(kBlack);

			//
			// Int_t markercolor = this->QPEtoCOL(Q);
			// marker->SetMarkerColor(markercolor);
			// marker->SetMarkerSize(1.0);
			//
			// Double_t markersize = this->QPEtoSIZE(Q);
			// marker->SetMarkerSize(markersize);
			// marker->SetMarkerColor(kBlack);
			//

			fTransHistFilterDigits->Fill(phideg);
			if (phideg < thetadeg) {
				fLongHistFilterDigits->Fill(Ltrack);
				fVertexHistFilterDigitsCone->Fill(delta);
				fTimeHistFilterDigitsCone->Fill(dt);
			}

			if (dt < fNewTmin)
				fNewTmin = dt;
			if (dt > fNewTmax)
				fNewTmax = dt;
			if (path < fNewRmin)
				fNewRmin = path;
			if (path > fNewRmax)
				fNewRmax = path;
			if (Ltrack < fNewLmin)
				fNewLmin = Ltrack;
			if (Ltrack > fNewLmax)
				fNewLmax = Ltrack;

			wcVertexFilterDigits.push_back(marker);
		}
	}

	// set T limits
	// ============
	fNewTmin -= 10.0;
	fNewTmax += 100.0;

	if (fNewTmax < fNewTmin) {
		fNewTmin = fTmin;
		fNewTmax = fTmax;
	}

	fVertexHistTimeFit->GetYaxis()->SetLimits(fNewTmin, fNewTmax);
	fTimeHistDigitsAll->GetXaxis()->SetRangeUser(fNewTmin, fNewTmax);

	// set R limits
	// ============
	fNewRmin -= 100.0;
	fNewRmax += 100.0;

	if (fNewRmax < fNewRmin) {
		fNewRmin = fRmin;
		fNewRmax = fRmax;
	}

	fVertexHistTimeFit->GetXaxis()->SetLimits(fNewRmin, fNewRmax);

	// set L limits
	// ============
	fNewLmin -= 50.0;
	fNewLmax += 50.0;

	if (fNewLmax < fNewLmin) {
		fNewLmin = fLmin;
		fNewLmax = fLmax;
	}

	fLongHistDigits->GetXaxis()->SetRangeUser(fNewLmin, fNewLmax);

	// set digits min/max (timing residuals)
	// =====================================
	hmax = 0.0;
	for (Int_t i = 0; i < fVertexHistDigitsAll->GetXaxis()->GetNbins(); i++) {
		if (fVertexHistDigitsAll->GetBinContent(i + 1) > hmax) {
			hmax = fVertexHistDigitsAll->GetBinContent(i + 1);
		}
	}

	max = (Double_t)(100.0 * (1.5 + (Int_t)(hmax / 100.0)));
	min = 0.0;

	fVertexHistDigitsAll->SetMinimum(min);
	fVertexHistDigitsAll->SetMaximum(max);

	// set digits min/max (timing profile)
	// ===================================
	hmax = 0.0;
	for (Int_t i = 0; i < fTimeHistDigitsAll->GetXaxis()->GetNbins(); i++) {
		if (fTimeHistDigitsAll->GetBinContent(i + 1) > hmax) {
			hmax = fTimeHistDigitsAll->GetBinContent(i + 1);
		}
	}

	max = (Double_t)(100.0 * (1.5 + (Int_t)(hmax / 100.0)));
	min = 0.0;

	fTimeHistDigitsAll->SetMinimum(min);
	fTimeHistDigitsAll->SetMaximum(max);

	// set digits min/max (transverse)
	// ===============================
	hmax = 0.0;
	for (Int_t i = 0; i < fTransHistDigits->GetXaxis()->GetNbins(); i++) {
		if (fTransHistDigits->GetBinContent(i + 1) > hmax) {
			hmax = fTransHistDigits->GetBinContent(i + 1);
		}
	}

	max = (Double_t)(100.0 * (1.5 + (Int_t)(hmax / 100.0)));
	min = 0.0;

	fTransHistDigits->SetMinimum(min);
	fTransHistDigits->SetMaximum(max);

	// set digits min/max (longitudinal)
	// =================================
	hmax = 0.0;
	for (Int_t i = 0; i < fLongHistDigits->GetXaxis()->GetNbins(); i++) {
		if (fLongHistDigits->GetBinContent(i + 1) > hmax) {
			hmax = fLongHistDigits->GetBinContent(i + 1);
		}
	}

	max = (Double_t)(100.0 * (1.5 + (Int_t)(hmax / 100.0)));
	min = 0.0;

	fLongHistDigits->SetMinimum(min);
	fLongHistDigits->SetMaximum(max);

	// draw 1D canvas
	// ==============
	if (fVertexCanvas1D) {

		fVertexCanvas1D->cd();
		fVertexCanvas1D->Clear();
		fVertexHistDigitsAll->Draw();
		fVertexHistDigitsCone->Draw("same");
		fVertexHistFilterDigitsCone->Draw("same");
		fVertexLegend1D->Draw();

		min = fVertexHistDigitsAll->GetMinimum();
		max = fVertexHistDigitsAll->GetMaximum();
		fTitleLatex.DrawLatex(-20.0, max + 0.015 * (max - min), fTitleString.Data());

		fVertexCanvas1D->Update();
	}

	// draw 2D canvas
	// ==============
	if (fVertexCanvas2D) {

		fVertexCanvas2D->cd();
		fVertexCanvas2D->Clear();

		fVertexHistTimeFit->Draw();

		for (UInt_t idigit = 0; idigit < wcVertexDigits.size(); idigit++) {
			TMarker* marker = (TMarker*) (wcVertexDigits.at(idigit));
			if (marker->GetX() > fNewRmin && marker->GetX() < fNewRmax && marker->GetY() > fNewTmin
					&& marker->GetY() < fNewTmax) {
				marker->Draw();
			}
		}

		for (UInt_t idigit = 0; idigit < wcVertexFilterDigits.size(); idigit++) {
			TMarker* marker = (TMarker*) (wcVertexFilterDigits.at(idigit));
			if (marker->GetX() > fNewRmin && marker->GetX() < fNewRmax && marker->GetY() > fNewTmin
					&& marker->GetY() < fNewTmax) {
				marker->Draw();
			}
		}

		fVertexTimeFit->Draw("same");
		fVertexLegend2D->Draw();

		fTitleLatex.DrawLatex(fNewRmin, fNewTmax + 0.015 * (fNewTmax - fNewTmin), fTitleString.Data());

		fVertexCanvas2D->Update();
	}

	// draw canvas
	// ===========
	if (fTimeCanvas) {

		fTimeCanvas->cd();
		fTimeCanvas->Clear();

		fTimeHistDigitsAll->Draw();
		fTimeHistDigitsCone->Draw("same");
		fTimeHistFilterDigitsCone->Draw("same");
		fTimeLegend->Draw();

		min = fTimeHistDigitsAll->GetMinimum();
		max = fTimeHistDigitsAll->GetMaximum();

		fTitleLatex.DrawLatex(fNewTmin, max + 0.015 * (max - min), fTitleString.Data());

		fTimeCanvas->Update();
	}

	// draw canvas
	// ===========
	if (fTransCanvas) {

		fTransCanvas->cd();
		fTransCanvas->Clear();

		fTransHistDigits->Draw();
		fTransHistFilterDigits->Draw("same");
		fTransLegend->Draw();

		min = fTransHistDigits->GetMinimum();
		max = fTransHistDigits->GetMaximum();

		fTitleLatex.DrawLatex(0.0, max + 0.015 * (max - min), fTitleString.Data());

		fTransCanvas->Update();
	}

	// draw canvas
	// ===========
	if (fLongCanvas) {

		fLongCanvas->cd();
		fLongCanvas->Clear();

		fLongHistDigits->Draw();
		fLongHistFilterDigits->Draw("same");
		fLongLegend->Draw();

		min = fLongHistDigits->GetMinimum();
		max = fLongHistDigits->GetMaximum();

		fTitleLatex.DrawLatex(fNewLmin, max + 0.015 * (max - min), fTitleString.Data());

		fLongCanvas->Update();
	}

	return;
}

Int_t WCSimVertexViewer::QPEtoCOL(Double_t Q) {
	Int_t colourCode = kWhite;

	if (Q < 0.8)
		colourCode = kYellow - 7;
	if (Q >= 0.8 && Q < 1.5)
		colourCode = kCyan - 7;
	if (Q >= 1.5 && Q < 2.5)
		colourCode = kCyan + 1;
	if (Q >= 2.5 && Q < 5.0)
		colourCode = kBlue - 4;
	if (Q >= 5.0 && Q < 10.0)
		colourCode = kBlue + 1;
	if (Q >= 10.0 && Q < 15.0)
		colourCode = kMagenta + 1;
	if (Q >= 15.0 && Q < 20.0)
		colourCode = kMagenta + 2;
	if (Q >= 20.0 && Q < 30.0)
		colourCode = kRed - 4;
	if (Q >= 30.0)
		colourCode = kRed;

	return colourCode;
}

Double_t WCSimVertexViewer::QPEtoSIZE(Double_t Q) {
	Double_t markerSize = 0.2;

	if (Q >= 0.0 && Q < 1.0)
		markerSize = 0.2;
	if (Q >= 1.0 && Q < 2.5)
		markerSize = 0.5;
	if (Q >= 2.5 && Q < 5.0)
		markerSize = 0.8;
	if (Q >= 5.0 && Q < 7.5)
		markerSize = 1.0;
	if (Q >= 7.5 && Q < 10.0)
		markerSize = 1.2;
	if (Q >= 10.0 && Q < 15.0)
		markerSize = 1.4;
	if (Q >= 15.0 && Q < 20.0)
		markerSize = 1.7;
	if (Q >= 20.0 && Q < 25.0)
		markerSize = 2.0;
	if (Q >= 25.0 && Q < 50.0)
		markerSize = 2.3;
	if (Q >= 50.0)
		markerSize = 2.7;

	return markerSize;
}
