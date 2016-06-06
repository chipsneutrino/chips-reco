/*
 * WCSimRecoEvDisplay.cc
 *
 *  Created on: 13 Feb 2015
 *      Author: ajperch
 */

#include "WCSimEvDisplay.hh"
#include "WCSimRecoEvDisplay.hh"
#include "WCSimRecoSummary.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include <TGButton.h>
#include <TGFrame.h>
#include <TGFileDialog.h>
#include <TGMenu.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TPaveText.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TFile.h>
#include <TChain.h>
#include <TH2D.h>
#include <TPolyMarker.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TColor.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <utility>

ClassImp(WCSimRecoEvDisplay)

WCSimRecoEvDisplay::WCSimRecoEvDisplay(){
	// TODO Auto-generated constructor stub

}

WCSimRecoEvDisplay::WCSimRecoEvDisplay(const TGWindow* p, UInt_t w, UInt_t h)
{
	// Initialise some histogram pointers
	fBarrelHist = 0x0;
	fTopHist = 0x0;
	fBottomHist = 0x0;
	fChargeHist = 0x0;
	fTimeHist = 0x0;
	fChargeRMTHist = 0x0;
	fTimeRMTHist = 0x0;

	// Initialise the TChain pointers
	fChain = 0x0;
	fGeomTree = 0x0;

	// Initialise the truth object
	fTruthSummary = 0x0;
	fTruthLegend = 0x0;
	fWhichPads = 0; // Default is reco
	fTruthTextMain = new TPaveText(0.05, 0.45, 0.95, 0.90, "NDC");
	fTruthTextPrimaries = new TPaveText(0.05, 0.1, 0.95, 0.40, "NDC");
  // The above TPaveText objects will be resized should we need to show this.
  fTruthTextOverlay = new TPaveText(0.05,0.05,0.95,0.15,"NDC");
	fDatabasePDG = 0x0;

	// Initialise the TGNumberEntry
	fPEInput = 0x0;
	fChargeCut = 0;

  fLogZCharge = false;
  
  // Create the TGraph vectors with default TGraphs
  this->MakeGraphColours();

  for(unsigned int g = 0; g < fColours.size(); ++g){
    fTopGraphs.push_back(new TGraph());
    fBarrelGraphs.push_back(new TGraph());
    fBottomGraphs.push_back(new TGraph());
    // Initialise the graphs
    this->InitialiseGraph(fTopGraphs[g],g);
    this->InitialiseGraph(fBarrelGraphs[g],g);
    this->InitialiseGraph(fBottomGraphs[g],g);
  }

	// Set up some plot style
	this->SetStyle();

	std::cout << "Starting WCSim Event Display" << std::endl;

	// Create the menu bar that lives at the top of the window
	TGMenuBar *mainMenu = new TGMenuBar(this, 400, 20, kHorizontalFrame);
	//	mainMenu->AddPopup("&File",fMenuFile,new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
	this->AddFrame(mainMenu,
			new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0,
					0, 2));

	// Create canvas widget
	fHitMapCanvas = new TRootEmbeddedCanvas("barrelCanvas", this, 400, 200);
	TCanvas *can = fHitMapCanvas->GetCanvas();
	can->cd();
	fBarrelPad = new TPad("fBarrelPad", "", 0.0, 0.6, 1.0, 1.0);
	fTopPad = new TPad("fTopPad", "", 0.0, 0.2, 0.487, 0.6);
	fBottomPad = new TPad("fBottomPad", "", 0.487, 0.2, 1.0, 0.6);
	fChargePad = new TPad("fChargePad", "", 0.0, 0.0, 0.5, 0.2);
	fTimePad = new TPad("fTimePad", "", 0.5, 0.0, 1.0, 0.2);
  // Do a bit of beauty work
  fBarrelPad->SetLeftMargin(0.05);
  fBarrelPad->SetRightMargin(0.075);
  fBarrelPad->SetBottomMargin(0.12);
  fTopPad->SetRightMargin(0.0);
  fTopPad->SetBottomMargin(0.12);
  fBottomPad->SetLeftMargin(0.0);
  fBottomPad->SetRightMargin(0.147);
  fBottomPad->SetBottomMargin(0.12);
	fBarrelPad->Draw();
	fTopPad->Draw();
	fBottomPad->Draw();
	fChargePad->Draw();
	fTimePad->Draw();
  fBarrelPad->SetTicks(1,1);
  fTopPad->SetTicks(1,1);
  fBottomPad->SetTicks(1,1);

	// Create the truth information pad, but don't draw it
	fTruthPad = new TPad("fTruthPad", "", 0.0, 0.0, 1.0, 1.0);
	fTruthOverlayPad = new TPad("fTruthOverlayPad", "", 0.0, 0.0, 1.0, 0.2);

	this->AddFrame(fHitMapCanvas,
			new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));

	// These build from top to bottom, so add them in the correct order.
	this->CreateSubButtonBar();
  this->CreateDebugButtonBar();
	this->CreateMainButtonBar();

	// Initialise the histograms and display them.
	//	this->ResizePlotsFromNtuple();

	// Set a name to the main frame
	this->SetWindowName("CHIPS Event Display");
	// Map all subwindows of main frame
	this->MapSubwindows();
	// Set the window size
	this->Resize(w, h);
	// Map main frame
	this->MapWindow();

	// Create temporary dummy plots then update them later
	this->MakeDefaultPlots();
	this->UpdateCanvases();

	// Set the current event to the default value of 0.
	fCurrentEvent = 0;
	fFileType = -1;
	this->HideFrame(hWCSimButtons);
	fViewType = 0; // Look at charge by default

	// By default show the 1D plots
	fShow1DHists = 1;

	// Create the truth information pad, but don't draw it
	fFitPad = new TPad("fFitPad", "", 0.0, 0.0, 1.0, 1.0);
	fFitOverlayPad = new TPad("fFitOverlayPad", "", 0.0, 0.0, 1.0, 0.2);
	fFitTextMain = new TPaveText(0.05, 0.45, 0.95, 0.90, "NDC");
	fFitTextPrimaries = new TPaveText(0.05, 0.1, 0.95, 0.40, "NDC");

  // This display pad takes up the space where the 1D plots usually live.
  fDisplayPad = new TPad("fDisplayPad","",0.0,0.0,1.0,0.2);

	fFitLegend = 0x0;
	fDisplayLegend = 0x0;

	fRecoSummary = 0x0;

	fRecoSummaryChain = 0x0;
  fHitComparisonChain = 0x0;

}

WCSimRecoEvDisplay::~WCSimRecoEvDisplay() {
	// TODO Auto-generated destructor stub
}

void WCSimRecoEvDisplay::UpdateFitPad() {


	  fFitPad->cd();
	  fFitTextMain->Draw();
	  fFitTextPrimaries->Draw();

	  fHitMapCanvas->GetCanvas()->cd();
}

void WCSimRecoEvDisplay::UpdateFitOverlayPad() {
	  fFitOverlayPad->cd();
	  // Draw the TLegend

	  if(fFitLegend){
      fFitLegend->Draw();
	  }
	  else{
	    std::cout << "No fit rings found" << std::endl;
	  }

	  fHitMapCanvas->GetCanvas()->cd();
}

void WCSimRecoEvDisplay::CreateMainButtonBar() {
	// Create a horizontal frame widget with buttons
	TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 200, 40);
	TGTextButton *open = new TGTextButton(hframe, "&Open");
	open->Connect("Clicked()", "WCSimRecoEvDisplay", this, "OpenFile()");
	hframe->AddFrame(open, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
	TGTextButton *prev = new TGTextButton(hframe, "&Prev");
	prev->Connect("Clicked()", "WCSimRecoEvDisplay", this, "PrevEvent()");
	hframe->AddFrame(prev, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
	TGTextButton *next = new TGTextButton(hframe, "&Next");
	next->Connect("Clicked()", "WCSimRecoEvDisplay", this, "NextEvent()");
	hframe->AddFrame(next, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  // Numeric entry field for choosing events
  fEventInput = new TGNumberEntry(hframe,0,9,999,TGNumberFormat::kNESInteger,
                               TGNumberFormat::kNEANonNegative);
  fEventInput->Connect("ValueSet(Long_t)","WCSimRecoEvDisplay",this,"SetEvent()");
  (fEventInput->GetNumberEntry())->Connect("ReturnPressed()","WCSimEvDisplay",this,"SetEvent()");
  // Make a label to go along side it
  TGLabel *eventLabel = new TGLabel(hframe,"Event:");
	hframe->AddFrame(eventLabel, new TGLayoutHints(kLHintsCenterX&&kLHintsCenterY,5,5,3,4));
	hframe->AddFrame(fEventInput, new TGLayoutHints(kLHintsCenterX,5,5,3,4)); 

  TGTextButton *disp = new TGTextButton(hframe, "&Print View");
  disp->Connect("Clicked()","WCSimRecoEvDisplay",this,"ShowDisplayView()");
  hframe->AddFrame(disp,new TGLayoutHints(kLHintsCenterX,5,5,3,4));
	TGTextButton *save = new TGTextButton(hframe, "&Save");
	save->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SaveEvent()");
	hframe->AddFrame(save, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
	TGTextButton *exit = new TGTextButton(hframe, "&Exit",
			"gApplication->Terminate(0)");
	hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
	this->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

}

void WCSimRecoEvDisplay::CreateSubButtonBar() {
	std::cout << "Doing it in the reco EvDisplay" << std::endl;

	// Create a horizontal frame to store buttons specific to WCSim files
	hWCSimButtons = new TGHorizontalFrame(this, 200, 40);

  // Switch between veto and inner detector PMTs
  TGTextButton *togVeto = new TGTextButton(hWCSimButtons,"&Veto");
  togVeto->Connect("Clicked()","WCSimEvDisplay",this,"ShowVeto()");
  hWCSimButtons->AddFrame(togVeto, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

	// Show Charge on z axis
	TGTextButton *togCharge = new TGTextButton(hWCSimButtons, "&Charge");
	togCharge->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewCharge()");
	hWCSimButtons->AddFrame(togCharge, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show time on z axis
	TGTextButton *togTime = new TGTextButton(hWCSimButtons, "&Time");
	togTime->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewTime()");
	hWCSimButtons->AddFrame(togTime, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Toggle the 1D plots
	TGTextButton *togHists = new TGTextButton(hWCSimButtons, "Toggle 1D Plots");
	togHists->Connect("Clicked()", "WCSimRecoEvDisplay", this, "Toggle1DHists()");
	hWCSimButtons->AddFrame(togHists, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show the standard reco view
	TGTextButton *showReco = new TGTextButton(hWCSimButtons, "&Reco");
	showReco->Connect("Clicked()", "WCSimRecoEvDisplay", this, "ShowReco()");
	hWCSimButtons->AddFrame(showReco, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show the truth overlay view
	TGTextButton *showTruthOverlay = new TGTextButton(hWCSimButtons, "&Truth Overlay");
	showTruthOverlay->Connect("Clicked()", "WCSimRecoEvDisplay", this, "ShowTruthOverlay()");
	hWCSimButtons->AddFrame(showTruthOverlay, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show the truth summary information
	TGTextButton *showTruth = new TGTextButton(hWCSimButtons, "&Truth");
	showTruth->Connect("Clicked()", "WCSimRecoEvDisplay", this, "ShowTruth()");
	hWCSimButtons->AddFrame(showTruth, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show the fit overlay view
	TGTextButton *showFitOverlay = new TGTextButton(hWCSimButtons, "&Fit Overlay");
	showFitOverlay->Connect("Clicked()", "WCSimRecoEvDisplay", this, "ShowFitOverlay()");
	hWCSimButtons->AddFrame(showFitOverlay, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show the truth summary information
	TGTextButton *showFit = new TGTextButton(hWCSimButtons, "&Fit");
	showFit->Connect("Clicked()", "WCSimRecoEvDisplay", this, "ShowFit()");
	hWCSimButtons->AddFrame(showFit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Numeric entry field for the PE threshold
	fPEInput = new TGNumberEntry(hWCSimButtons, 0, 9, 999, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative);

	// TGNumberEntry has two ways to set numbers, so two connects
	fPEInput->Connect("ValueSet(Long_t)", "WCSimRecoEvDisplay", this,"SetChargeCut()");
	(fPEInput->GetNumberEntry())->Connect("ReturnPressed()", "WCSimRecoEvDisplay", this, "SetChargeCut()");

	// Make a label to go along side it
	TGLabel *peLabel = new TGLabel(hWCSimButtons, "Charge Cut:");
	hWCSimButtons->AddFrame(peLabel, new TGLayoutHints(kLHintsCenterX && kLHintsCenterY, 5, 5, 3, 4));
	hWCSimButtons->AddFrame(fPEInput, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  // Toggle button for logging the charge z axis
  TGTextButton *logZToggle = new TGTextButton(hWCSimButtons,"&Log Charge");
  logZToggle->Connect("Clicked()","WCSimRecoEvDisplay",this,"ToggleLogZ()");
  hWCSimButtons->AddFrame(logZToggle, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

	// Add the TGHorizontalFrame to the main layout
	this->AddFrame(hWCSimButtons, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

}

void WCSimRecoEvDisplay::CreateDebugButtonBar() {

	// Buttons for reconstruction debugging
	TGHorizontalFrame *hDebugButtons = new TGHorizontalFrame(this, 200, 40);

	// Show total lnl on z axis
	TGTextButton *togChargeRMT = new TGTextButton(hDebugButtons, "&Charge Reco-True");
	togChargeRMT->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewChargeRMT()");
	hDebugButtons->AddFrame(togChargeRMT, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show total lnl on z axis
	TGTextButton *togTimeRMT = new TGTextButton(hDebugButtons, "&Time Reco-True");
	togTimeRMT->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewTimeRMT()");
	hDebugButtons->AddFrame(togTimeRMT, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show total lnl on z axis
	TGTextButton *togTotLnL = new TGTextButton(hDebugButtons, "&Total LnL");
	togTotLnL->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewTotalLnL()");
	hDebugButtons->AddFrame(togTotLnL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show charge lnl on z axis
	TGTextButton *togQLnL = new TGTextButton(hDebugButtons, "&Charge LnL");
	togQLnL->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewQLnL()");
	hDebugButtons->AddFrame(togQLnL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show time lnl on z axis
	TGTextButton *togTLnL = new TGTextButton(hDebugButtons, "&Time LnL");
	togTLnL->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewTLnL()");
	hDebugButtons->AddFrame(togTLnL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Add the TGHorizontalFrame to the main layout
	this->AddFrame(hDebugButtons, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

}

void WCSimRecoEvDisplay::DrawFitRing(unsigned int particleNo, int colour) {

	  // Does the truth information exist
	  if(fRecoSummary==0x0) return;

	  // Get the track vertex and direction
	  TVector3 trkVtx = fRecoSummary->GetVertex(particleNo);
	  trkVtx = trkVtx * 0.001; // Convert to m
	  TVector3 trkDir;
	  trkDir = fRecoSummary->GetPrimaryDir(particleNo);

	  // Get the projection of the track vertex onto the wall
	  TVector3 trkVtxProj; // Point where the track would hit the wall, filled by ProjectToWall
	  unsigned int detRegion; // Region of the detector, filled by ProjectToWall
	  this->ProjectToWall(trkVtx,trkDir,trkVtxProj, detRegion);

	  // Now that we have the projected point on the wall, iterate around the Cherenkov cone and
	  // find where all the points of the cone intercept the wall
	  unsigned int nMarkers = 1440; // Number of points that will appear on the final plots
	  double dPhi = 360 / (double)nMarkers; // Angle step around the direction vector

	  // Get the particle mass
	  if(fDatabasePDG == 0x0){
	    fDatabasePDG = new TDatabasePDG();
	  }
	  int pdgCode;
	  double en;
	  pdgCode = fRecoSummary->GetPrimaryPDG(particleNo);
	  en = fRecoSummary->GetPrimaryEnergy(particleNo);

	  double mass = 1000*fDatabasePDG->GetParticle(pdgCode)->Mass();
	  double beta = sqrt(en*en - mass*mass) / en; // beta = p / E
	  double refrac = 1.33; // Refractive index of water
	  double thetaC = (180.0 / TMath::Pi()) * TMath::ACos(1./(refrac*beta)); // Cherenkov angle

	  // Also need 6 vectors to store the 2D-coordinates for the 3 regions
	  std::vector<double> topPos1; // For top, this is the y coord
	  std::vector<double> topPos2; // For top, this is the x coord
	  std::vector<double> barrelPos1; // This is the phi coord
	  std::vector<double> barrelPos2; // This is the z coord
	  std::vector<double> bottomPos1; // This is the y coord
	  std::vector<double> bottomPos2; // This is the x coord

	  for(unsigned int n = 0; n < nMarkers; ++n){

	    double phi = n * dPhi;

	    // Find a point on the Cherenkov cone and its direction
	    TVector3 circPos;
	    TVector3 circDir;
	    this->FindCircle(trkVtxProj,trkVtx,thetaC,phi,circPos,circDir);

	    // Now we have the point on the circle, project this onto the wall
	    TVector3 finalPos;
	    this->ProjectToWall(circPos,circDir,finalPos,detRegion);

	    if(detRegion == 0){
	      topPos1.push_back(finalPos.Y());
	      topPos2.push_back(finalPos.X());
	    }
	    else if(detRegion == 1){
	      barrelPos1.push_back(TMath::ATan2(finalPos.Y(),finalPos.X()));
	      barrelPos2.push_back(finalPos.Z());
	    }
	    else{
	      bottomPos1.push_back(finalPos.Y());
	      bottomPos2.push_back(finalPos.X());
	    }
	  }

	  // Now we have the vectors of coordinates that we need for each region. Now to make the TPolyMarkers
	  // and add an entry to the legend
	  std::stringstream legendText;
    std::string particleName = this->GetParticleName(pdgCode);
	  legendText << particleName << " :: vertex = (" << trkVtx.X() << ", " << trkVtx.Y() << ", " << trkVtx.Z() 
               << ") m, direction = (" << trkDir.X() << ", " << trkDir.Y() << ", " << trkDir.Z() 
               << ") and energy = " << en << " MeV";
	  TLine* line = new TLine();
	  line->SetLineColor(colour);
    line->SetLineWidth(2);
	  fFitLines.push_back(line);
	  fFitLegend->AddEntry(fFitLines[fFitLines.size()-1],legendText.str().c_str(),"l");
	  if(topPos1.size() != 0){
	    this->MakePolyMarker(topPos1,topPos2,fFitMarkersTop,colour);
	  }
	  if(barrelPos1.size() != 0){
	    this->MakePolyMarker(barrelPos1,barrelPos2,fFitMarkersBarrel,colour);
	  }
	  if(bottomPos1.size() != 0){
	    this->MakePolyMarker(bottomPos1,bottomPos2,fFitMarkersBottom,colour);
	  }
}

void WCSimRecoEvDisplay::ClearFitMarkerVectors() {
	this->DeleteAndClearElements(fFitMarkersTop);
	this->DeleteAndClearElements(fFitMarkersBarrel);
	this->DeleteAndClearElements(fFitMarkersBottom);
	// Also delete the vector of lines
	for(unsigned int l = 0; l < fFitLines.size(); ++l){
		delete (TLine*)fFitLines[l];
	    fFitLines[l] = 0x0;
	}
	fFitLines.clear();
}

int WCSimRecoEvDisplay::GetFitRingColour(int ring) const {
	if(ring == 1) { return kViolet; }
	if(ring == 2) { return kBlue; }
	if(ring == 3) { return kOrange+1; }
	if(ring == 4) { return kGreen+1; }
	return kBlack;
}

// Match the reconstruction to the truth for the display view.
// For now, match by pdg and then best vertex delta. Without a pdg match,
// just don't bother.
void WCSimRecoEvDisplay::MatchRecoToTrue(){

  fDisplayTrackPairs.clear();
  // Loop over the reconstructed objects.
  for(unsigned int r = 0; r < fRecoSummary->GetNPrimaries(); ++r){
    unsigned int bestMatch = 999;
    double vertexDelta = 1.e8;
    TVector3 recoVtx = fRecoSummary->GetVertex(r);
    int recoPDG = fRecoSummary->GetPrimaryPDG(r);
    // Is this a neutrino event?
    if(fTruthSummary->IsNeutrinoEvent()){
      // Loop over the true tracks 
      for(unsigned int t = 0; t < fTruthSummary->GetNPrimaries(); ++t){
        TVector3 trueVtx = fTruthSummary->GetVertex();
        if(recoPDG == fTruthSummary->GetPrimaryPDG(t)){
          if((recoVtx-trueVtx).Mag() < vertexDelta){
            bestMatch = t;
            vertexDelta = (recoVtx-trueVtx).Mag();
          }
        }
      }
      // Also check the true overlaid objects
      if(fTruthSummary->IsOverlayEvent()){
        for(unsigned int o = 0; o < fTruthSummary->GetNOverlays(); ++o){
          TVector3 trueVtx = fTruthSummary->GetOverlayVertex(); // Only one vertex for all overlaid tracks
          if(recoPDG == fTruthSummary->GetOverlayPDG(o)){
            if((recoVtx-trueVtx).Mag() < vertexDelta){
              bestMatch = o + fTruthSummary->GetNPrimaries();
              vertexDelta = (recoVtx-trueVtx).Mag();
            }
          }
        }
      }
    }
    else{
      // Use the particle gun information
      TVector3 trueVtx = fTruthSummary->GetVertex();
      if(recoPDG == fTruthSummary->GetBeamPDG()){
        if((recoVtx-trueVtx).Mag() < vertexDelta){
          bestMatch = 0; // Only one particle gun.
          vertexDelta = (recoVtx-trueVtx).Mag();
        }
      }
    }

   
    // Fill the match vector (999 means no match found).
    fDisplayTrackPairs.push_back(std::make_pair<unsigned int,unsigned int>(r,bestMatch));
  }

}

void WCSimRecoEvDisplay::UpdateDisplayInfo(){

  if(fDisplayLegend != 0x0){
	    if(fDisplayPad->GetListOfPrimitives()->FindObject(fDisplayLegend)){
	      fDisplayPad->GetListOfPrimitives()->Remove(fDisplayLegend);
	    }
    delete fDisplayLegend;
  }
  fDisplayPad->cd();
  fDisplayLegend = new TLegend(0.05,0.1,0.4,0.9);
//  fDisplayLegend = new TLegend(0.05,0.1,0.95,0.9);
  fDisplayLegend->SetBorderSize(0);
  fDisplayLegend->SetTextSize(0.1);
//  fDisplayLegend->SetTextSize(0.2);
//  fDisplayLegend->SetNColumns(2); 
  // Loop over the display object vector. No lines for the truth.
  for(unsigned int i = 0; i < fDisplayTrackPairs.size(); ++i){
    unsigned int recoTrack = fDisplayTrackPairs[i].first;
    fDisplayLegend->AddEntry(fFitLines[recoTrack],(this->GetTrackInfo(true,recoTrack)).c_str(),"l");

    unsigned int trueTrack = fDisplayTrackPairs[i].second;
//    if(fTruthSummary->IsNeutrinoEvent()){
          fDisplayLegend->AddEntry((TObject*)0x0,(this->GetTrackInfo(false,trueTrack)).c_str(),"");
//    }
    
  }

  fDisplayLegend->Draw();
  // Go back to the main pad.
	fHitMapCanvas->GetCanvas()->cd();

}

std::string WCSimRecoEvDisplay::GetTrackInfo(bool isReco, unsigned int trackNo){

  bool returnEmpty = true;

  std::stringstream info;
  
  std::cout << "Looking for track info: " << isReco << ", " << trackNo << std::endl;

  int pdg;
  TVector3 vtx;
  TVector3 dir;
  double energy;

  if(isReco){
    pdg = fRecoSummary->GetPrimaryPDG(trackNo);
    vtx = fRecoSummary->GetVertex(trackNo);
    dir = fRecoSummary->GetPrimaryDir(trackNo);
    energy = fRecoSummary->GetPrimaryEnergy(trackNo);
    returnEmpty = false;
  }
  else{
    if(fTruthSummary->IsNeutrinoEvent()){
      if(trackNo != 999){
        returnEmpty = false;
        // Is this a primary particle?
        if(trackNo < fTruthSummary->GetNPrimaries()){
          pdg = fTruthSummary->GetPrimaryPDG(trackNo);
          vtx = fTruthSummary->GetVertex();
          dir = fTruthSummary->GetPrimaryDir(trackNo);
          energy = fTruthSummary->GetPrimaryEnergy(trackNo);
        }
        else{
          int overlayNo = trackNo - fTruthSummary->GetNPrimaries();
          pdg = fTruthSummary->GetOverlayPDG(overlayNo);
          vtx = fTruthSummary->GetOverlayVertex();
          dir = fTruthSummary->GetOverlayDir(overlayNo);
          energy = fTruthSummary->GetOverlayEnergy(overlayNo);
        }
      }
    }
    else{
      // If this is a particle gun then we just use the beam particle... there aren't others.
      if(trackNo == 0){
        returnEmpty = false;
        pdg = fTruthSummary->GetBeamPDG();
        vtx = fTruthSummary->GetVertex();
        dir = fTruthSummary->GetBeamDir();
        energy = fTruthSummary->GetBeamEnergy();  
      }
    }
  }

  if(returnEmpty){
    info << "No truth match";
  }
  else{
    if(isReco){
      info << "Reconstructed: ";
    }
    else{
      info << "Truth: ";
    }
    // Remember to convert the vertex into metres
    info << this->GetParticleName(pdg) << std::fixed << std::setprecision(3) << " at (" << vtx.X()*0.001 << ", " << vtx.Y()*0.001 << ", " << vtx.Z()*0.001
         << ")m with direction = (" << dir.X() << ", " << dir.Y() << ", " << dir.Z()
         << ") and energy = " << energy*0.001 << " GeV";
  }



  return info.str();

}

void WCSimRecoEvDisplay::OpenFile(std::string name) {
	TGFileInfo fileInfo;
	const char *filetypes[] = {"ROOT files", "*.root", 0, 0};
	fileInfo.fFileTypes = filetypes;
	fileInfo.fIniDir = StrDup(".");

	// Open the browser if no file is supplied.
	if(name == ""){
		new TGFileDialog(gClient->GetDefaultRoot(), this, kFDOpen, &fileInfo);
		if(fileInfo.fFilename){
			name = fileInfo.fFilename;
		}
		std::cout << "'" << name << "' selected." << std::endl;
	}

	if(name != ""){
		// Quick test to see if the file contains what we need
		TFile tempFile(name.c_str(),"READ");
    if(tempFile.Get("wcsimT")){
    	std::cout << "Opening reco file " << name << std::endl;
			this->OpenWCSimRecoFile(name);

		}
		else{
			std::cout << "Selected file is not a WCSim file." << std::endl;
		}
		tempFile.Close();
	}
	std::cout << "Open!" << name << std::endl;
}

void WCSimRecoEvDisplay::ShowFit() {
	if (fRecoSummary != 0x0) {
		if (fWhichPads != 4) {

			fWhichPads = 4;

			std::cout << "== Event Fit Information ==" << std::endl;

			std::string intType = this->ConvertTrueEventType();
			TVector3 vtx = fRecoSummary->GetVertex(0);
			// The interaction
			std::cout << "= Interaction = " << std::endl;
			std::cout << " - Vertex   : (" << vtx.X() << ","
    					  << vtx.Y() << "," << vtx.Z() << ")" << std::endl;

			// Primaries
			std::cout << "= Particles = " << std::endl;
			for (unsigned int n = 0; n < fRecoSummary->GetNPrimaries();	++n) {
				TVector3 dir = fRecoSummary->GetPrimaryDir(n);
				std::cout << " - Particle: "
						  << fRecoSummary->GetPrimaryPDG(n);
				std::cout << " with energy "
						  << fRecoSummary->GetPrimaryEnergy(n);
				std::cout << " MeV and direction (" << dir.X() << ","
						  << dir.Y() << "," << dir.Z() << ")" << std::endl;
			}
			this->ResizePads();
		} else {
			std::cerr << "Already displaying the fit information"
					<< std::endl;
		}
	} else {
		std::cerr << "No fit information found, so can't display it"
				<< std::endl;
	}
}

void WCSimRecoEvDisplay::UpdateFitPave() {
	  fFitTextMain->Clear();
	  fFitTextMain->SetFillColor(kWhite);
	  fFitTextMain->SetBorderSize(0);

	  fFitTextPrimaries->Clear();
	  fFitTextPrimaries->SetFillColor(kWhite);
	  fFitTextPrimaries->SetBorderSize(0);

	  // Clear truth ring vectors
	  this->ClearFitMarkerVectors();

	  // Stream to parse things into strings
	  std::stringstream tmpS;

    if(fRecoSummary->HasCommonVertex()){
  	  TVector3 vtx = fRecoSummary->GetVertex(0);
  
  	  std::cout << "Vertex = " << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << std::endl;
  	  std::cout << "Radius = " << fWCRadius << ", " << "Height = " << fWCLength << std::endl;
  	  tmpS << vtx.X() << "," << vtx.Y() << "," << vtx.Z();
  	  fFitTextMain->AddText(("Vertex at ("+tmpS.str()+") mm").c_str());
      tmpS.str("");
      tmpS << "Vertex time = " << fRecoSummary->GetVertexT(0) << " ns";
      fFitTextMain->AddText(tmpS.str().c_str());
    }
	  int nFitRings = 0;
	  // Create the TLegend for the truth overlays
	  if(fFitLegend != 0x0){
	    // Remove it from the pad
	    if(fFitOverlayPad->GetListOfPrimitives()->FindObject(fFitLegend)){
	      fFitOverlayPad->GetListOfPrimitives()->Remove(fFitLegend);
	    }
	    delete fFitLegend;
	    fFitLegend = 0x0;
	  }
	  fFitLegend = new TLegend(0.0,0.1,0.5,0.9, "Reconstructed Particles");
	  fFitLegend->SetFillColor(kWhite);
	  fFitLegend->SetBorderSize(0);
    fFitLegend->SetTextSize(0.1);

	  // The primary particles list
	  fFitTextPrimaries->AddText("List of fitted particles");
	  TVector3 dir;
	  for(unsigned int n = 0; n < fRecoSummary->GetNPrimaries(); ++n){
		  int pdg = fRecoSummary->GetPrimaryPDG(n);
	    double energy = fRecoSummary->GetPrimaryEnergy(n);
	    std::string mod = "";
	    // Add a true ring to the display
	    ++nFitRings;
	    int ringColour = this->GetFitRingColour(nFitRings);
	    this->DrawFitRing(n,ringColour);

	    dir = fRecoSummary->GetPrimaryDir(n);
		  tmpS.str("");
      if(fRecoSummary->HasCommonVertex() == false){
  	    TVector3 vtx = fRecoSummary->GetVertex(n);
    	  tmpS << "Vertex = (" << vtx.X() << "," << vtx.Y() << "," << vtx.Z() << ") mm";
        tmpS << " at " << fRecoSummary->GetVertexT(n) << " ns";
		    fFitTextPrimaries->AddText(tmpS.str().c_str());
        tmpS.str("");
      }
		  tmpS << mod << " ";
		  tmpS << "Particle: " << pdg;
		  tmpS << " with energy " << energy;
		  tmpS << " MeV and direction (" << dir.X() << "," << dir.Y() << "," << dir.Z() << ")";
		  tmpS << " " << mod;
		  fFitTextPrimaries->AddText(tmpS.str().c_str());
	  }
}

void WCSimRecoEvDisplay::OpenWCSimRecoFile(std::string name) {
	fFileType = 1; // We have a saved fit file

	// Check what we should be showing.
	this->ResizePads();
    // Show the WCSim buttons if they are not visible.
	if(!this->IsVisible(hWCSimButtons)){
		this->ShowFrame(hWCSimButtons);
	}

	// Sort the main chain first
	if(fChain != 0x0){
	  delete fChain;
	}

	fChain = new TChain("wcsimT");
	fChain->Reset();
	fChain->Add(name.c_str());

	if(fRecoSummaryChain != 0x0)
	{
		delete fRecoSummaryChain;
	}

	fRecoSummaryChain = new TChain("fRecoSummaryTree");
	fRecoSummaryChain->Reset();
	fRecoSummaryChain->Add(name.c_str());

  if(fHitComparisonChain != 0x0){
    delete fHitComparisonChain;
  }
  fHitComparisonChain = new TChain("fHitComparisonTree");
  fHitComparisonChain->Reset();
  fHitComparisonChain->Add(name.c_str());

	// Now the geometry

	if(fGeomTree != 0x0){
	  delete fGeomTree;
	}
	fGeomTree = new TChain("wcsimGeoT");
	fGeomTree->Reset();
	fGeomTree->Add(name.c_str());
	this->ResizePlotsFromGeometry();

	// Each entry is an event, so use this to set the limits.
	fMinEvent = 0;
	fCurrentEvent = 0;
	fMaxEvent = fChain->GetEntries() - 1;

	std::cout << "Filling plots" << std::endl;
	this->FillPlots();

}

void WCSimRecoEvDisplay::FillPlots() {
  if(fViewType < 2){
	  FillPlotsFromRecoFile();
  }
  else if (fViewType < 5){
    FillPlotsFromLikelihood();
  }
  else{
    FillPlotsFromRMT();
  }
}

void WCSimRecoEvDisplay::FillPlotsFromRecoFile() {
	// First things first, clear the histograms.
	this->ClearPlots();

	// Need to load the events.	
	WCSimRootEvent *wcSimEvt = new WCSimRootEvent();
 	fChain->SetBranchAddress("wcsimrootevent",&wcSimEvt);
  // Force deletion to prevent memory leak 
  fChain->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	// Get the geometry
	WCSimRootGeom *geo = new WCSimRootGeom();
	fGeomTree->SetBranchAddress("wcsimrootgeom",&geo);
	fGeomTree->GetEntry(0);

	// Get the current event;
	fChain->GetEntry(fCurrentEvent);
	if(wcSimEvt==0x0) std::cout << "Null pointer :( " << std::endl;
	WCSimRootTrigger* wcSimTrigger = wcSimEvt->GetTrigger(0);
    
  // Get the fit information  
  fRecoSummaryChain->SetBranchAddress("WCSimRecoSummary",&fRecoSummary);
  fRecoSummaryChain->GetBranch("WCSimRecoSummary")->SetAutoDelete(kTRUE);
	fRecoSummaryChain->GetEntry(fCurrentEvent);
  this->UpdateFitPave();

  // Get the truth information
  if(fTruthSummary != 0x0){
    delete fTruthSummary;
    fTruthSummary = 0x0;
  }
  fTruthSummary = new WCSimTruthSummary(wcSimEvt->GetTruthSummary());
  this->ClearPi0Vector();

  // Quick check for pi zeroes and their decay photons
  if(fTruthSummary->IsPrimaryPiZero()){
    std::vector<double> pi0EnVec = fTruthSummary->GetPiZeroEnergies();
    for(unsigned int p = 0; p < pi0EnVec.size(); ++p){
      this->SearchForPi0Photons(pi0EnVec[p],wcSimTrigger->GetTracks());
    }
  }
  // If we have a pi-zero gun, make sure we treat it properly.
  else if (fTruthSummary->GetBeamPDG() == 111){
    this->SearchForPi0Photons(fTruthSummary->GetBeamEnergy(),wcSimTrigger->GetTracks());
  }

  // Update the truth view
  this->UpdateTruthPave();

  int nDigiHits = wcSimTrigger->GetNcherenkovdigihits();
	std::cout << "Number of PMTs hit: " << nDigiHits << std::endl;

  // Need to loop through the hits once to find the charge and time ranges
  fQMin = 1e10;
  fQMax = -1e10;
  fTMin = 1e10;
  fTMax = -1e10;
  for(int i = 0; i < nDigiHits; ++i){
		TObject *element = (wcSimTrigger->GetCherenkovDigiHits())->At(i);
		WCSimRootCherenkovDigiHit *hit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
		WCSimRootPMT pmt = geo->GetPMTFromTubeID(hit->GetTubeId());
    // Skip veto PMTs
    if(fViewVeto){
      if(pmt.GetCylLoc() != 3) continue;
    }
    else{
      if(pmt.GetCylLoc() == 3) continue;
    }

    double q = hit->GetQ();
    double t = hit->GetT();
    if(q < fQMin) fQMin = q;
    if(q > fQMax) fQMax = q;
    if(t < fTMin) fTMin = t;
    if(t > fTMax) fTMax = t;
  }
  // For now, always want charge to start at 0.
  fQMin = 0;

  this->CalculateChargeAndTimeBins();
  this->ResetGraphs();

  // Now loop through again and fill things
	for (int i=0;i<nDigiHits;i++)
	{
		// Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
		TObject *element = (wcSimTrigger->GetCherenkovDigiHits())->At(i);

		WCSimRootCherenkovDigiHit *wcSimDigiHit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
		WCSimRootPMT pmt = geo->GetPMTFromTubeID(wcSimDigiHit->GetTubeId());

    // Skip veto PMTs
    if(fViewVeto){
      if(pmt.GetCylLoc() != 3) continue;
    }
    else{
      if(pmt.GetCylLoc() == 3) continue;
    }

		double pmtX = pmt.GetPosition(0) * 0.01; // convert to m
		double pmtY = pmt.GetPosition(1) * 0.01; // convert to m
		double pmtZ = pmt.GetPosition(2) * 0.01; // convert to m
		double pmtPhi = TMath::ATan2(pmtY,pmtX);
		double pmtQ = wcSimDigiHit->GetQ();
		double pmtT = wcSimDigiHit->GetT();
		// Set the z-axis to be charge or time.
		double colourAxis = pmtQ;
    fBarrelHist->GetZaxis()->SetTitle("Charge (p.e.)");
		if(fViewType == 1){
      colourAxis = pmtT;
      fBarrelHist->GetZaxis()->SetTitle("Time (ns)");
    }
    fTopHist->GetZaxis()->SetTitle("Charge (p.e.)");
		if(fViewType == 1){
      fTopHist->GetZaxis()->SetTitle("Time (ns)");
    }
    fBottomHist->GetZaxis()->SetTitle("Charge (p.e.)");
		if(fViewType == 1){
      fBottomHist->GetZaxis()->SetTitle("Time (ns)");
    }

    // Make sure we pass the charge cut
    if(pmtQ > fChargeCut){
      unsigned int bin;
      if(fViewType == 0){ bin = this->GetChargeBin(pmtQ);}
      else{ bin = this->GetTimeBin(pmtT);}

      // Set the underflow bin of the histograms to make sure the colour axis shows
      fTopHist->SetBinContent(0,1);
      fBarrelHist->SetBinContent(0,1);
      fBottomHist->SetBinContent(0,1);
    
	  	// Top cap
  		if(pmt.GetCylLoc() == 0){
//  			fTopHist->Fill(pmtY,pmtX,colourAxis);
          fTopGraphs[bin]->SetPoint(fTopGraphs[bin]->GetN(),pmtY,pmtX);
  		}
  		// Bottom cap
  		else if(pmt.GetCylLoc() == 2){
//  			fBottomHist->Fill(pmtY,pmtX,colourAxis);
          fBottomGraphs[bin]->SetPoint(fBottomGraphs[bin]->GetN(),pmtY,pmtX);
      
  		}
  		// Barrel
  		else if(pmt.GetCylLoc() == 1){
//  			fBarrelHist->Fill(pmtPhi,pmtZ,colourAxis);
          fBarrelGraphs[bin]->SetPoint(fBarrelGraphs[bin]->GetN(),pmtPhi,pmtZ);
  		}
      // Veto
      else{
//        std::cout << "Veto PMT hit: " <<  pmtX << ", " << pmtY << ", " << pmtZ << " :: " << pmtPhi << ", " << pmtQ << ", " << pmtT << std::endl;
        if(pmt.GetOrientation(2) > 0.99){
          fTopGraphs[bin]->SetPoint(fTopGraphs[bin]->GetN(),pmtY,pmtX);
        }
        else if(pmt.GetOrientation(2) < -0.99){
          fBottomGraphs[bin]->SetPoint(fBottomGraphs[bin]->GetN(),pmtY,pmtX);
        }
        else{
          fBarrelGraphs[bin]->SetPoint(fBarrelGraphs[bin]->GetN(),pmtPhi,pmtZ);
        }
      }
		  // Now fill the 1D histograms
		  fChargeHist->Fill(pmtQ);	 
		  fTimeHist->Fill(pmtT);	 
    }


	} // End of loop over Cherenkov digihits

	delete geo;
	geo = 0x0;

  // Update the pads
  this->UpdateRecoPads();
  this->UpdateTruthPad();
  this->UpdateTruthOverlayPad();
	this->UpdateFitPad();
  this->UpdateFitOverlayPad();
  // Attempt to match the reco and true objects
  this->MatchRecoToTrue();    
  // Update the display
  this->UpdateDisplayInfo();
  // Now draw whichever pad we need
	this->UpdateCanvases();
}

void WCSimRecoEvDisplay::ResizePads() {
	// Get list of objects attached to the main canvas
	TList *list = fHitMapCanvas->GetCanvas()->GetListOfPrimitives();
	// UpdateCanvases draw the pads we want, so remove them all here
	if (list->FindObject(fTruthPad)) {
		list->Remove(fTruthPad);
	}
	if (list->FindObject(fTruthOverlayPad)) {
		list->Remove(fTruthOverlayPad);
	}
	if (list->FindObject(fFitPad)) {
		list->Remove(fFitPad);
	}
	if (list->FindObject(fFitOverlayPad)) {
		list->Remove(fFitOverlayPad);
	}
	if (list->FindObject(fBarrelPad)) {
		list->Remove(fBarrelPad);
	}
	if (list->FindObject(fTopPad)) {
		list->Remove(fTopPad);
	}
	if (list->FindObject(fBottomPad)) {
		list->Remove(fBottomPad);
	}
	if (list->FindObject(fChargePad)) {
		list->Remove(fChargePad);
	}
	if (list->FindObject(fTimePad)) {
		list->Remove(fTimePad);
	}
  if (list->FindObject(fDisplayPad)){
    list->Remove(fDisplayPad);
  }

	// If we want to show truth
	if (fWhichPads == 1) {
		fTruthPad->SetPad(0.0, 0.0, 1.0, 1.0);
	}
	// Or else show the reco
	else if (fWhichPads == 0) {
		// Resize the reco pads
		if (fShow1DHists) {
			fBarrelPad->SetPad(0.0, 0.6, 1.0, 1.0);
			fTopPad->SetPad(0.0, 0.2, 0.487, 0.6);
			fBottomPad->SetPad(0.487, 0.2, 1.0, 0.6);
			fChargePad->SetPad(0.0, 0.0, 0.5, 0.2);
			fTimePad->SetPad(0.5, 0.0, 1.0, 0.2);
		} else {
			// Resize the reco pads
			fBarrelPad->SetPad(0.0, 0.5, 1.0, 1.0);
			fTopPad->SetPad(0.0, 0.0, 0.487, 0.5);
			fBottomPad->SetPad(0.487, 0.0, 1.0, 0.5);
		}
	}
	// Else show the truth overlays
	else if (fWhichPads == 2) {
		// Make sure to leave space for the truth overlay pad at the bottom
		fBarrelPad->SetPad(0.0, 0.6, 1.0, 1.0);
		fTopPad->SetPad(0.0, 0.2, 0.487, 0.6);
		fBottomPad->SetPad(0.487, 0.2, 1.0, 0.6);
		fTruthOverlayPad->SetPad(0.0, 0.0, 1.0, 0.2);
	}
	// Else show the fit overlays
	else if(fWhichPads == 3) {
		// Make sure to leave space for the fit overlay pad at the bottom
		fBarrelPad->SetPad(0.0, 0.6, 1.0, 1.0);
		fTopPad->SetPad(0.0, 0.2, 0.487, 0.6);
		fBottomPad->SetPad(0.487, 0.2, 1.0, 0.6);
		fFitOverlayPad->SetPad(0.0, 0.0, 1.0, 0.2);
	}
	// Else show the fit text
	else if(fWhichPads == 4){
		fFitPad->SetPad(0.0, 0.0, 1.0, 1.0);
	}
  // Otherwise show the display view
  else{
		fBarrelPad->SetPad(0.0, 0.6, 1.0, 1.0);
		fTopPad->SetPad(0.0, 0.2, 0.487, 0.6);
		fBottomPad->SetPad(0.487, 0.2, 1.0, 0.6);
		fDisplayPad->SetPad(0.0, 0.0, 1.0, 0.2);
  }
	std::cout << "Done!  Update canvases..." << std::endl;
	this->UpdateCanvases();
}

void WCSimRecoEvDisplay::UpdateCanvases() {

	TCanvas *canvas = fHitMapCanvas->GetCanvas();
	canvas->cd();

	// Remove the truth overlays
	this->HideTruthOverlays();
  this->HideFitOverlays();

  fBarrelPad->SetLogz(fLogZCharge);
  fTopPad->SetLogz(fLogZCharge);
  fBottomPad->SetLogz(fLogZCharge);

	if (fWhichPads == 0) {
		// Now draw the pads
		fBarrelPad->Draw();
		fTopPad->Draw();
		fBottomPad->Draw();
		if (fShow1DHists) {
			fChargePad->Draw();
			fTimePad->Draw();
		}
	} else if (fWhichPads == 1) {
		fTruthPad->Draw();
	} else if (fWhichPads == 2) {
    this->DrawTruthOverlays();
		canvas->cd(); // Need to cd back here since the above changes directory
		fBarrelPad->Draw();
		fTopPad->Draw();
		fBottomPad->Draw();
		fTruthOverlayPad->Draw();
	} else if (fWhichPads == 3) {
    this->DrawFitOverlays();
		canvas->cd(); // Need to cd back here since the above changes directory
		fBarrelPad->Draw();
		fTopPad->Draw();
		fBottomPad->Draw();
		fFitOverlayPad->Draw();
	}
  else if (fWhichPads == 4){
    fFitPad->Draw();
  }
  // Else show the display view
  else{
    this->DrawFitOverlays();
		canvas->cd(); // Need to cd back here since the above changes directory
		fBarrelPad->Draw();
		fTopPad->Draw();
		fBottomPad->Draw();
		fDisplayPad->Draw();
  }

	canvas->Modified();
	canvas->Update();
}

void WCSimRecoEvDisplay::DrawFitOverlays() {
  std::cout << "Drawing fit overlays" << std::endl;
	// Take the plots one by one and draw them.
	fBarrelPad->cd();
	fBarrelHist->Draw("colz");
  fBarrelTitle->Draw();
  this->DrawHitGraphs(fBarrelGraphs);
	// Draw the fit rings
  std::cout << "There are " << fFitMarkersBarrel.size() << " barrel markers" << std::endl;
	for (unsigned int r = 0; r < fFitMarkersBarrel.size(); ++r) {
		fFitMarkersBarrel[r]->Draw("C");
	}
	fBarrelPad->Modified();
	fBarrelPad->Update();

	fTopPad->cd();
//	fTopHist->Draw("colz");
	fTopHist->Draw();
  fTopTitle->Draw();
  this->DrawHitGraphs(fTopGraphs);
	// Draw the fit rings
  std::cout << "There are " << fFitMarkersTop.size() << " top markers" << std::endl;
	for (unsigned int r = 0; r < fFitMarkersTop.size(); ++r) {
		fFitMarkersTop.at(r)->Draw("C");
	}
	fTopPad->Modified();
	fTopPad->Update();

	fBottomPad->cd();
	fBottomHist->Draw("colz");
  fBottomTitle->Draw();
  this->DrawHitGraphs(fBottomGraphs);
	// Draw the truth rings
	for (unsigned int r = 0; r < fFitMarkersBottom.size(); ++r) {
		fFitMarkersBottom[r]->Draw("C");
	}
	fBottomPad->Modified();
	fBottomPad->Update();

  this->AdjustBarrelZAxis();
}

void WCSimRecoEvDisplay::HideFitOverlays() {
	TList *list = fBarrelPad->GetListOfPrimitives();
	for (unsigned int r = 0; r < fFitMarkersBarrel.size(); ++r) {
		if (list->FindObject(fFitMarkersBarrel[r])) {
			list->Remove(fFitMarkersBarrel[r]);
		}
	}

	list = fTopPad->GetListOfPrimitives();
	for (unsigned int r = 0; r < fFitMarkersTop.size(); ++r) {
		if (list->FindObject(fFitMarkersTop[r])) {
			list->Remove(fFitMarkersTop[r]);
		}
	}

	list = fBottomPad->GetListOfPrimitives();
	for (unsigned int r = 0; r < fFitMarkersBottom.size(); ++r) {
		if (list->FindObject(fFitMarkersBottom[r])) {
			list->Remove(fFitMarkersBottom[r]);
		}
	}
}

void WCSimRecoEvDisplay::ShowFitOverlay() {
	if (fWhichPads != 3) {
		fWhichPads = 3;
    // If we aren't showing the reco plots, then show them.
    if(fViewType > 1){
      SetViewCharge();
    }
		this->ResizePads();
	} else {
		std::cerr << "Already displaying fit overlays" << std::endl;
	}
}

void WCSimRecoEvDisplay::ShowTruthOverlay(){
  if(fWhichPads != 2){
    fWhichPads = 2;
    // If we aren't showing the reco plots, then show them.
    if(fViewType > 1){
      SetViewCharge();
    }
    this->ResizePads();
  }
  else{
    std::cerr << "Already displaying truth overlays" << std::endl;
  }
}

void WCSimRecoEvDisplay::ShowDisplayView(){

  if(fWhichPads != 5){
    fWhichPads = 5;

    // Now we have all the information, redraw things.
    this->ResizePads();
  } 

}

void WCSimRecoEvDisplay::SetPlotZAxes(){

  double min = fQMin;
  double max = fQMax;
  if(max > 100.) max = 100.;

  if(fViewType == 1){
    min = fTMin;
    max = fTMax;
  }
  if(fViewType >= 2 && fViewType < 5){
    min = fLnLMin;
    max = fLnLMax;
  }
  if(fViewType == 5){
    min = fQRMTMin;
    max = fQRMTMax;
  }
  if(fViewType == 6){
    min = fTRMTMin;
    max = fTRMTMax;
  }

  // Make sure we've updated the palette
  gStyle->SetPalette(fColours.size(), &(fColours[0]));

	// Make sure the histogram max / min are correct
	fBarrelHist->SetMaximum(max);
	fBarrelHist->SetMinimum(min);
	fTopHist->SetMaximum(max);
	fTopHist->SetMinimum(min);
	fBottomHist->SetMaximum(max);
	fBottomHist->SetMinimum(min);

}

// Switch the z-axis scale to show to total likelihood.
void WCSimRecoEvDisplay::SetViewTotalLnL(){
	if(fViewType != 2){
		fViewType = 2;
    if(fWhichPads != 0){
      this->ShowReco();
    }
   	if(fShow1DHists == 1){
      this->Toggle1DHists();
    }
		std::cout << "Setting colour axis to total likelihood" << std::endl;
		this->FillPlots();
	}
	else{
		std::cerr << "Already viewing total likelihood." << std::endl;
	}
}

// Switch the z-axis scale to show the charge likelihood.
void WCSimRecoEvDisplay::SetViewQLnL(){
	if(fViewType != 3){
		fViewType = 3;
    if(fWhichPads != 0){
      this->ShowReco();
    }
   	if(fShow1DHists == 1){
      this->Toggle1DHists();
    }
		std::cout << "Setting colour axis to charge likelihood" << std::endl;
		this->FillPlots();
	}
	else{
		std::cerr << "Already viewing charge likelihood." << std::endl;
	}
}

// Switch the z-axis scale to show the time likelihood.
void WCSimRecoEvDisplay::SetViewTLnL(){
	if(fViewType != 4){
		fViewType = 4;
    if(fWhichPads != 0){
      this->ShowReco();
    }
   	if(fShow1DHists == 1){
      this->Toggle1DHists();
    }
		std::cout << "Setting colour axis to time likelihood" << std::endl;
		this->FillPlots();
	}
	else{
		std::cout << "Already viewing time likelihood." << std::endl;
	}
}

// Switch the z-axis scale to show the reco-true for charge.
void WCSimRecoEvDisplay::SetViewChargeRMT(){
	if(fViewType != 5){
		fViewType = 5;
    if(fWhichPads != 0){
      this->ShowReco();
    }
		std::cout << "Setting colour axis to reco-true for charge" << std::endl;
		this->FillPlots();
	}
	else{
		std::cout << "Already viewing time reco-true for charge" << std::endl;
	}
}

// Switch the z-axis scale to show the reco-true for time.
void WCSimRecoEvDisplay::SetViewTimeRMT(){
	if(fViewType != 6){
		fViewType = 6;
    if(fWhichPads != 0){
      this->ShowReco();
    }
		std::cout << "Setting colour axis to reco-true for time" << std::endl;
		this->FillPlots();
	}
	else{
		std::cout << "Already viewing time reco-true for time" << std::endl;
	}
}

// Function to fill the standard plots using the likelihood values
void WCSimRecoEvDisplay::FillPlotsFromLikelihood(){
	// First things first, clear the histograms.
	this->ClearPlots();

  // Set up access to the chain
  int pmtID = -1;
  int eventID = -1;
  fHitComparisonChain->SetBranchAddress("eventID",&eventID);
  fHitComparisonChain->SetBranchAddress("pmtID",&pmtID);
  
  // Only load the LnL values that we need.
  double lnlVal;
  std::string lnlVarName = "minus2LnL";
  if(fViewType == 3) lnlVarName = "charge2LnL";
  if(fViewType == 4) lnlVarName = "time2LnL";
  fHitComparisonChain->SetBranchAddress(lnlVarName.c_str(),&lnlVal);

	// Get the geometry
	WCSimRootGeom *geo = new WCSimRootGeom();
	fGeomTree->SetBranchAddress("wcsimrootgeom",&geo);
	fGeomTree->GetEntry(0);

  // Loop over the chain to find the min and max values of the likelihood
  fLnLMin = 1e10;
  fLnLMax = -1e10;
  // Loop over the chain to fill the plots
  for(int i = 0; i < fHitComparisonChain->GetEntries(); ++i){
    fHitComparisonChain->GetEntry(i);
    // Check if we are on the event that we care about
    if(eventID < fCurrentEvent) continue;
    if(eventID > fCurrentEvent) break;
   
    if(lnlVal < fLnLMin) fLnLMin = lnlVal; 
    if(lnlVal > fLnLMax) fLnLMax = lnlVal; 
  }
  // Let's keep the min at zero for now.
  fLnLMin = 0;
  this->CalculateLnLBins();
  this->ResetGraphs();

  // Loop over the chain to fill the plots
  for(int i = 0; i < fHitComparisonChain->GetEntries(); ++i){
    fHitComparisonChain->GetEntry(i);

    // Check if we are on the event that we care about
    if(eventID < fCurrentEvent) continue;
    if(eventID > fCurrentEvent) break;

    // Check which bin this value of LnL should go in
    unsigned int bin = this->GetLnLBin(lnlVal);

    pmtID++; // Looks like we have the off-by-one problem again

    // Get the PMT from the geometry. The coordinates are stored in the tree
    // but it is convenient to use the pmt to find out the region.
		WCSimRootPMT pmt = geo->GetPMTFromTubeID(pmtID);
		double pmtX = pmt.GetPosition(0);
		double pmtY = pmt.GetPosition(1);
		double pmtZ = pmt.GetPosition(2);
		double pmtPhi = TMath::ATan2(pmtY,pmtX);
    // Top cap
		if(pmt.GetCylLoc() == 0){
        fTopGraphs[bin]->SetPoint(fTopGraphs[bin]->GetN(),pmtY,pmtX);
		}
		// Bottom cap
		else if(pmt.GetCylLoc() == 2){
        fBottomGraphs[bin]->SetPoint(fBottomGraphs[bin]->GetN(),pmtY,pmtX);    
		}
		// Barrel
		else{
        fBarrelGraphs[bin]->SetPoint(fBarrelGraphs[bin]->GetN(),pmtPhi,pmtZ);
		}
  }  
  delete geo;
  geo = 0x0;

  // Set the underflow bin of the histograms to make sure the colour axis shows
  fTopHist->SetBinContent(0,1);
  fBarrelHist->SetBinContent(0,1);
  fBottomHist->SetBinContent(0,1);

  // Reset all the addresses linked to this chain
  fHitComparisonChain->ResetBranchAddresses();

  // Update the pads
  this->UpdateRecoPads();
  // Now draw whichever pad we need
	this->UpdateCanvases();
}

// Function to fill the standard plots using the likelihood values
void WCSimRecoEvDisplay::FillPlotsFromRMT(){
	// First things first, clear the histograms.
	this->ClearPlots();

  // Set up access to the chain
  int pmtID = -1;
  int eventID = -1;
  fHitComparisonChain->SetBranchAddress("eventID",&eventID);
  fHitComparisonChain->SetBranchAddress("pmtID",&pmtID);
  
  double recoQ = -999;
  double trueQ = -999;
  double recoT = -999;
  double trueT = -999;
  fHitComparisonChain->SetBranchAddress("predQ",&recoQ);
  fHitComparisonChain->SetBranchAddress("trueQ",&trueQ);
  fHitComparisonChain->SetBranchAddress("predT",&recoT);
  fHitComparisonChain->SetBranchAddress("trueT",&trueT);

	// Get the geometry
	WCSimRootGeom *geo = new WCSimRootGeom();
	fGeomTree->SetBranchAddress("wcsimrootgeom",&geo);
	fGeomTree->GetEntry(0);

  // Loop over the chain to find the min and max values of the likelihood
  fQRMTMin = 1e10;
  fQRMTMax = -1e10;
  fTRMTMin = 1e10;
  fTRMTMax = -1e10;
  // Loop over the chain to fill the plots
  for(int i = 0; i < fHitComparisonChain->GetEntries(); ++i){
    fHitComparisonChain->GetEntry(i);
    // Check if we are on the event that we care about
    if(eventID < fCurrentEvent) continue;
    if(eventID > fCurrentEvent) break;
    double rmtValQ = recoQ - trueQ;
    double rmtValT = recoT - trueT;
    // Time can have default values
    if(trueT <= 0 || recoT <= 0) rmtValT = 0;

    if(rmtValQ < fQRMTMin) fQRMTMin = rmtValQ; 
    if(rmtValQ > fQRMTMax) fQRMTMax = rmtValQ; 
    if(rmtValT < fTRMTMin) fTRMTMin = rmtValT; 
    if(rmtValT > fTRMTMax) fTRMTMax = rmtValT; 
  }
  // Nice to have these plots symmetric about 0:
  if(fabs(fQRMTMin) > fQRMTMax){
    fQRMTMax = -fQRMTMin;
  }
  else{
    fQRMTMin = -fQRMTMax;
  }
  if(fabs(fTRMTMin) > fTRMTMax){
    fTRMTMax = -fTRMTMin;
  }
  else{
    fTRMTMin = -fTRMTMax;
  }
  this->CalculateRMTBins();
  this->ResetGraphs();

  // Make the two 1D histograms
  if(fChargeRMTHist != 0x0){
    delete fChargeRMTHist;
  }
  fChargeRMTHist = new TH1D("hChargeRMT","Charge Reco - True; Charge Reco - True / PE",100,fQRMTMin,fQRMTMax);
  if(fTimeRMTHist != 0x0){
    delete fTimeRMTHist;
  }
  fTimeRMTHist = new TH1D("hTimeRMT","Time Reco - True; Time Reco - True / ns",100,fTRMTMin,fTRMTMax);

  // Loop over the chain to fill the plots
  for(int i = 0; i < fHitComparisonChain->GetEntries(); ++i){
    fHitComparisonChain->GetEntry(i);

    // Check if we are on the event that we care about
    if(eventID < fCurrentEvent) continue;
    if(eventID > fCurrentEvent) break;

    // Make the 1D plots
    fChargeRMTHist->Fill(recoQ-trueQ);
    if(!(trueT <= 0 || recoT <= 0)) fTimeRMTHist->Fill(recoT-trueT);

    // If we are looking at time, ignore those that shouldn't be hit
    if(fViewType == 6 && (trueT <= 0 || recoT <= 0)) continue;

    double rmtVal =  -999;
    if(fViewType == 5) rmtVal = recoQ - trueQ;
    if(fViewType == 6) rmtVal = recoT - trueT;
    // Check which bin this value of LnL should go in
    unsigned int bin = -1;
    if(fViewType == 5) bin = this->GetQRMTBin(rmtVal);
    if(fViewType == 6) bin = this->GetTRMTBin(rmtVal);

    pmtID++; // Looks like we have the off-by-one problem again

    // Get the PMT from the geometry. The coordinates are stored in the tree
    // but it is convenient to use the pmt to find out the region.
		WCSimRootPMT pmt = geo->GetPMTFromTubeID(pmtID);
		double pmtX = pmt.GetPosition(0);
		double pmtY = pmt.GetPosition(1);
		double pmtZ = pmt.GetPosition(2);
		double pmtPhi = TMath::ATan2(pmtY,pmtX);

    // Top cap
		if(pmt.GetCylLoc() == 0){
        fTopGraphs[bin]->SetPoint(fTopGraphs[bin]->GetN(),pmtY,pmtX);
		}
		// Bottom cap
		else if(pmt.GetCylLoc() == 2){
        fBottomGraphs[bin]->SetPoint(fBottomGraphs[bin]->GetN(),pmtY,pmtX);    
		}
		// Barrel
		else{
        fBarrelGraphs[bin]->SetPoint(fBarrelGraphs[bin]->GetN(),pmtPhi,pmtZ);
		}
  }  
  delete geo;
  geo = 0x0;

  // Set the underflow bin of the histograms to make sure the colour axis shows
  fTopHist->SetBinContent(0,1);
  fBarrelHist->SetBinContent(0,1);
  fBottomHist->SetBinContent(0,1);

  // Reset all the addresses linked to this chain
  fHitComparisonChain->ResetBranchAddresses();

  // Update the pads
  this->UpdateRecoPads();
  // Now draw whichever pad we need
	this->UpdateCanvases();
}


// Calculate the bins for the likelihood plots
void WCSimRecoEvDisplay::CalculateLnLBins(){
  // Firstly, clear the existing vector
  fLnLBins.clear();
  
  double delta = (fLnLMax - fLnLMin) / static_cast<double>(fColours.size());

  for(size_t i = 0; i < fColours.size(); ++i){
    fLnLBins.push_back(fLnLMin+i*delta);
  }
}

unsigned int WCSimRecoEvDisplay::GetLnLBin(double lnl) const{
  unsigned int bin = fColours.size()-1;
  for(size_t i = 0; i < fColours.size(); ++i){
    if(lnl < fLnLBins[i]){
      bin = i - 1;
      break;
    }
  }
  return bin;
}

// Calculate the bins for the likelihood plots
void WCSimRecoEvDisplay::CalculateRMTBins(){
  // Firstly, clear the existing vector
  fQRMTBins.clear();
  fTRMTBins.clear();
  
  double deltaQ = (fQRMTMax - fQRMTMin) / static_cast<double>(fColours.size());
  double deltaT = (fTRMTMax - fTRMTMin) / static_cast<double>(fColours.size());

  for(size_t i = 0; i < fColours.size(); ++i){
    fQRMTBins.push_back(fQRMTMin+i*deltaQ);
    fTRMTBins.push_back(fTRMTMin+i*deltaT);
  }
}

unsigned int WCSimRecoEvDisplay::GetQRMTBin(double rmt) const{
  unsigned int bin = fColours.size()-1;
  for(unsigned int i = 1; i < fQRMTBins.size(); ++i){
    if(rmt < fQRMTBins[i]){
      bin = i - 1;
      break;
    }
  }
  return bin;
}

unsigned int WCSimRecoEvDisplay::GetTRMTBin(double rmt) const{
  unsigned int bin = fTRMTBins.size()-1;
  for(unsigned int i = 1; i < fTRMTBins.size(); ++i){
    if(rmt < fTRMTBins[i]){
      bin = i - 1;
      break;
    }
  }
  return bin;
}

// Change the size of the reco text displays
void WCSimRecoEvDisplay::ResizeFitTexts(bool commonVertex){
 if(commonVertex){
    fFitTextMain->SetY1NDC(0.45);
    fFitTextPrimaries->SetY1NDC(0.10);
    fFitTextPrimaries->SetY2NDC(0.40);
  }
  else{
    fFitTextMain->SetY1NDC(0.89);
    fFitTextPrimaries->SetY1NDC(0.10);
    fFitTextPrimaries->SetY2NDC(0.90);
  }
  fFitPad->cd();
  fFitTextMain->Draw();
  fFitTextPrimaries->Draw();
  fFitPad->Modified();
  fFitPad->Update();
  fHitMapCanvas->GetCanvas()->cd();
  fHitMapCanvas->GetCanvas()->Modified();
  fHitMapCanvas->GetCanvas()->Update();
}

// Draw the reco plots to the reco pads
void WCSimRecoEvDisplay::UpdateRecoPads(){

  this->ResizeFitTexts(fRecoSummary->HasCommonVertex());

  this->MakeGraphColours();

	this->SetPlotZAxes();
  // Set the styles how we want them
  this->MakePlotsPretty(fBarrelHist);
  this->MakePlotsPretty(fTopHist);
  this->MakePlotsPretty(fBottomHist);
  this->MakePlotsPretty(fChargeHist);
  this->MakePlotsPretty(fTimeHist);
  if(fChargeRMTHist != 0x0) this->MakePlotsPretty(fChargeRMTHist);
  if(fTimeRMTHist != 0x0) this->MakePlotsPretty(fTimeRMTHist);

  // Take the plots one by one and draw them.
  fBarrelPad->cd();
  fBarrelHist->Draw("colz");
  fBarrelTitle->Draw();
  this->DrawHitGraphs(fBarrelGraphs);
  fBarrelPad->Modified();
  fBarrelPad->Update();

  fTopPad->cd();
//  fTopHist->Draw("colz");
  fTopHist->Draw();
  fTopTitle->Draw();
  this->DrawHitGraphs(fTopGraphs);
  fTopPad->Modified();
  fTopPad->Update();

  fBottomPad->cd();
  fBottomHist->Draw("colz");
  fBottomTitle->Draw();
  this->DrawHitGraphs(fBottomGraphs);
  fBottomPad->Modified();
  fBottomPad->Update();

  this->AdjustBarrelZAxis();

  fChargePad->cd();
  if(fViewType < 5){ 
    fChargeHist->Draw();
  }
  else{
    fChargeRMTHist->Draw();
  }
  fChargePad->Modified();
  fChargePad->Update();

  fTimePad->cd();
  if(fViewType < 5){
    fTimeHist->Draw();
  }
  else{
    fTimeRMTHist->Draw();
  }
  fTimePad->Modified();
  fTimePad->Update();

  fHitMapCanvas->GetCanvas()->cd();
}

void WCSimRecoEvDisplay::MakeGraphColours(){

  Int_t startColour = 0;
  const UInt_t nContours = 100;
  if(fViewType < 5 && fColoursNormal.size() < nContours){
    fColoursNormal.clear();
    // Make a palette
    const Int_t nRGBs = 9;
    // This is the kCool scheme taken from ROOT 6.
    Double_t stops[nRGBs] = {0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
    Double_t red[nRGBs] = { 33./255., 31./255., 42./255., 68./255., 86./255., 111./255., 141./255., 172./255., 227./255.};
    Double_t green[nRGBs] = { 255./255., 175./255., 145./255., 106./255., 88./255., 55./255., 15./255., 0./255., 0./255.};
    Double_t blue[nRGBs] = { 255./255., 205./255., 202./255., 203./255., 208./255., 205./255., 203./255., 206./255., 231./255.};
    startColour = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nContours);
    gStyle->SetNumberContours(nContours);
    for(size_t i = 0; i < nContours; ++i)
    {
      fColoursNormal.push_back(startColour+i);
    }
  }
  else if( fColoursDiff.size() < nContours){
    // Make a palette
    fColoursDiff.clear();
    const Int_t nRGBs = 9;
    Double_t stops[nRGBs] =  { 0.000, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000 };
    Double_t red[nRGBs] =    { 0.000, 0.250, 0.500, 0.750, 1.000, 1.000, 1.000, 1.000, 1.000 };
    Double_t green[nRGBs] =  { 0.000, 0.250, 0.500, 0.750, 0.750, 0.750, 0.500, 0.250, 0.000 };
    Double_t blue[nRGBs] =   { 1.000, 1.000, 1.000, 1.000, 1.000, 0.750, 0.500, 0.250, 0.000 };
    startColour = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nContours);
    gStyle->SetNumberContours(nContours);
    for(size_t i = 0; i < nContours; ++i)
    {
      fColoursDiff.push_back(startColour+i);
    }
  }
  
  // Update the palette
  if( fViewType < 5 ){
    fColours = fColoursNormal;
  }
  else{
    fColours = fColoursDiff;
  }

  for(unsigned int g = 0; g < fTopGraphs.size(); ++g){
    // Initialise the graphs
    this->InitialiseGraph(fTopGraphs[g],g);
    this->InitialiseGraph(fBarrelGraphs[g],g);
    this->InitialiseGraph(fBottomGraphs[g],g);
  }
}

