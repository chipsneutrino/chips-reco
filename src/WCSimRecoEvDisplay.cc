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

#include <iostream>
#include <string>
#include <sstream>

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

	// Initialise the TChain pointers
	fChain = 0x0;
	fGeomTree = 0x0;

	// Initialise the truth object
	fTruthSummary = 0x0;
	fTruthLegend = 0x0;
	fWhichPads = 0; // Default is reco
	fTruthTextMain = new TPaveText(0.05, 0.45, 0.95, 0.90, "NDC");
	fTruthTextPrimaries = new TPaveText(0.05, 0.1, 0.95, 0.40, "NDC");
	fDatabasePDG = 0x0;

	// Initialise the TGNumberEntry
	fPEInput = 0x0;
	fChargeCut = 0;
  
  // Create the TGraph vectors with default TGraphs
  this->MakeGraphColours();
  gStyle->SetPalette(10,fColours);
  for(unsigned int g = 0; g < 10; ++g){
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
	fTopPad = new TPad("fTopPad", "", 0.0, 0.2, 0.5, 0.6);
	fBottomPad = new TPad("fBottomPad", "", 0.5, 0.2, 1.0, 0.6);
	fChargePad = new TPad("fChargePad", "", 0.0, 0.0, 0.5, 0.2);
	fTimePad = new TPad("fTimePad", "", 0.5, 0.0, 1.0, 0.2);
	fBarrelPad->Draw();
	fTopPad->Draw();
	fBottomPad->Draw();
	fChargePad->Draw();
	fTimePad->Draw();
	// Create the truth information pad, but don't draw it
	fTruthPad = new TPad("fTruthPad", "", 0.0, 0.0, 1.0, 1.0);
	fTruthOverlayPad = new TPad("fTruthOverlayPad", "", 0.0, 0.0, 1.0, 0.2);

	this->AddFrame(fHitMapCanvas,
			new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));

	// These build from top to bottom, so add them in the correct order.
	CreateSubButtonBar();
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

	fFitLegend = 0x0;

	fRecoSummary = 0x0;

	fRecoSummaryChain = 0x0;
  fHitComparisonChain = 0x0;

	std::cout << "Constructor finished" << std::endl;
}

WCSimRecoEvDisplay::~WCSimRecoEvDisplay() {
	// TODO Auto-generated destructor stub
}

void WCSimRecoEvDisplay::UpdateFitPad() {

    std::cout << "Changing to " << fFitPad << std::endl;
	  fFitPad->cd();
    std::cout << "Drawing text " << fFitTextMain << std::endl;
	  fFitTextMain->Draw();
    std::cout << "Drawing primaries text " << fFitTextPrimaries << std::endl;
	  fFitTextPrimaries->Draw();

	  fHitMapCanvas->GetCanvas()->cd();
}

void WCSimRecoEvDisplay::UpdateFitOverlayPad() {
    std::cout << "Changing to fit overlay pad " << fFitOverlayPad << std::endl;
	  fFitOverlayPad->cd();
	  // Draw the TLegend

	  std::cout << "Drawing fit legend" << fFitLegend << std::endl;
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

	// Show Charge on z axis
	TGTextButton *togCharge = new TGTextButton(hWCSimButtons, "&Charge");
	togCharge->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewCharge()");
	hWCSimButtons->AddFrame(togCharge, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show time on z axis
	TGTextButton *togTime = new TGTextButton(hWCSimButtons, "&Time");
	togTime->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewTime()");
	hWCSimButtons->AddFrame(togTime, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show total lnl on z axis
	TGTextButton *togTotLnL = new TGTextButton(hWCSimButtons, "&Total LnL");
	togTotLnL->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewTotalLnL()");
	hWCSimButtons->AddFrame(togTotLnL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show charge lnl on z axis
	TGTextButton *togQLnL = new TGTextButton(hWCSimButtons, "&Charge LnL");
	togQLnL->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewQLnL()");
	hWCSimButtons->AddFrame(togQLnL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	// Show time lnl on z axis
	TGTextButton *togTLnL = new TGTextButton(hWCSimButtons, "&Time LnL");
	togTLnL->Connect("Clicked()", "WCSimRecoEvDisplay", this, "SetViewTLnL()");
	hWCSimButtons->AddFrame(togTLnL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

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

	// Add the TGHorizontalFrame to the main layout
	this->AddFrame(hWCSimButtons, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

}

void WCSimRecoEvDisplay::DrawFitRing(unsigned int particleNo, int colour) {

	  // Does the truth information exist
	  if(fRecoSummary==0x0) return;

	  // Get the track vertex and direction
	  TVector3 trkVtx = fRecoSummary->GetVertex();
	  trkVtx = trkVtx * 0.1; // Convert to cm
	  TVector3 trkDir;
	  trkDir = fRecoSummary->GetPrimaryDir(particleNo);

	  // Get the projection of the track vertex onto the wall
	  TVector3 trkVtxProj; // Point where the track would hit the wall, filled by ProjectToWall
	  unsigned int detRegion; // Region of the detector, filled by ProjectToWall
	  this->ProjectToWall(trkVtx,trkDir,trkVtxProj, detRegion);

	  // Now that we have the projected point on the wall, iterate around the Cherenkov cone and
	  // find where all the points of the cone intercept the wall
	  unsigned int nMarkers = 360; // Number of points that will appear on the final plots
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
	  legendText << "PDG code = " << pdgCode << " and total energy = " << en << " MeV";
	  TLine* line = new TLine();
	  line->SetLineColor(colour);
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

void WCSimRecoEvDisplay::OpenFile(std::string name) {
	std::cout << "Opening" << name << std::endl;
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
			TVector3 vtx = fRecoSummary->GetVertex();
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
    std::cout << "Updating fFitTextMain = " << fFitTextMain << std::endl;
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
	  TVector3 vtx = fRecoSummary->GetVertex();

	  std::cout << "Vertex = " << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << std::endl;
	  std::cout << "Radius = " << fWCRadius << ", " << "Height = " << fWCLength << std::endl;
	  tmpS << vtx.X() << "," << vtx.Y() << "," << vtx.Z();
	  fFitTextMain->AddText(("Vertex at ("+tmpS.str()+") mm").c_str());
    tmpS.str("");
    tmpS << "Vertex time = " << fRecoSummary->GetVertexT() << " ns";
    fFitTextMain->AddText(tmpS.str().c_str());
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
	  fFitLegend = new TLegend(0.2,0.2,0.8,0.8, "Fit results");
	  fFitLegend->SetFillColor(kWhite);
	  fFitLegend->SetBorderSize(0);

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
		  tmpS << mod << " ";
		  tmpS << "Particle: " << pdg;
		  tmpS << " with energy " << energy;
		  tmpS << " MeV and direction (" << dir.X() << "," << dir.Y() << "," << dir.Z() << ")";
		  tmpS << " " << mod;
		  fFitTextPrimaries->AddText(tmpS.str().c_str());
	  }
}

void WCSimRecoEvDisplay::OpenWCSimRecoFile(std::string name) {
	std::cout << "Opening reco file from " << name << std::endl;
	fFileType = 1; // We have a saved fit file

	// Check what we should be showing.
	this->ResizePads();
    // Show the WCSim buttons if they are not visible.
	if(!this->IsVisible(hWCSimButtons)){
		this->ShowFrame(hWCSimButtons);
	}

	std::cout << "First the main chain" << std::endl;
	// Sort the main chain first
	if(fChain != 0x0){
	  delete fChain;
	}

	fChain = new TChain("wcsimT");
	fChain->Reset();
	fChain->Add(name.c_str());

	std::cout << "Then the reco summary chain" << std::endl;
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

	std::cout << "Then the geometry" << std::endl;
	if(fGeomTree != 0x0){
	  delete fGeomTree;
	}
	fGeomTree = new TChain("wcsimGeoT");
	fGeomTree->Reset();
	fGeomTree->Add(name.c_str());
	this->ResizePlotsFromGeometry();


	std::cout << "Huzzahs are in order" << std::endl;
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
  else{
    FillPlotsFromLikelihood();
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
  std::cout << "Getting reco summary" << std::endl;
  std::cout << "Chain has " << fRecoSummaryChain->GetEntries() << std::endl;
  fRecoSummaryChain->SetBranchAddress("WCSimRecoSummary",&fRecoSummary);
  fRecoSummaryChain->GetBranch("WCSimRecoSummary")->SetAutoDelete(kTRUE);
	fRecoSummaryChain->GetEntry(fCurrentEvent);
  std::cout << "Got it: " << fRecoSummary << std::endl;
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
		double pmtX = pmt.GetPosition(0);
		double pmtY = pmt.GetPosition(1);
		double pmtZ = pmt.GetPosition(2);
		double pmtPhi = TMath::ATan2(pmtY,pmtX);
		double pmtQ = wcSimDigiHit->GetQ();
		double pmtT = wcSimDigiHit->GetT();
		// Set the z-axis to be charge or time.
		double colourAxis = pmtQ;
		if(fViewType == 1) colourAxis = pmtT;

    // Make sure we pass the charge cut
    if(pmtQ > fChargeCut){
    unsigned int bin;
    if(fViewType == 0) bin = this->GetChargeBin(pmtQ);
    else bin = this->GetTimeBin(pmtT);

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
  		else{
//  			fBarrelHist->Fill(pmtPhi,pmtZ,colourAxis);
          fBarrelGraphs[bin]->SetPoint(fBarrelGraphs[bin]->GetN(),pmtPhi,pmtZ);
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
  // Now draw whichever pad we need
	this->UpdateCanvases();
}

void WCSimRecoEvDisplay::ResizePads() {
	std::cout << "Resizing plots - pads" << std::endl;
	// Get list of objects attached to the main canvas
	TList *list = fHitMapCanvas->GetCanvas()->GetListOfPrimitives();
	std::cout << "Clearing primitives" << std::endl;
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

	std::cout << "Setting pads" << std::endl;
	// If we want to show truth
	if (fWhichPads == 1) {
		fTruthPad->SetPad(0.0, 0.0, 1.0, 1.0);
	}
	// Or else show the reco
	else if (fWhichPads == 0) {
		// Resize the reco pads
		if (fShow1DHists) {
			fBarrelPad->SetPad(0.0, 0.6, 1.0, 1.0);
			fTopPad->SetPad(0.0, 0.2, 0.5, 0.6);
			fBottomPad->SetPad(0.5, 0.2, 1.0, 0.6);
			fChargePad->SetPad(0.0, 0.0, 0.5, 0.2);
			fTimePad->SetPad(0.5, 0.0, 1.0, 0.2);
		} else {
			// Resize the reco pads
			fBarrelPad->SetPad(0.0, 0.5, 1.0, 1.0);
			fTopPad->SetPad(0.0, 0.0, 0.5, 0.5);
			fBottomPad->SetPad(0.5, 0.0, 1.0, 0.5);
		}
	}
	// Else show the truth overlays
	else if (fWhichPads == 2) {
		// Make sure to leave space for the truth overlay pad at the bottom
		fBarrelPad->SetPad(0.0, 0.6, 1.0, 1.0);
		fTopPad->SetPad(0.0, 0.2, 0.5, 0.6);
		fBottomPad->SetPad(0.5, 0.2, 1.0, 0.6);
		fTruthOverlayPad->SetPad(0.0, 0.0, 1.0, 0.2);
	}
	// Else show the fit overlays
	else if(fWhichPads == 3) {
		// Make sure to leave space for the fit overlay pad at the bottom
		fBarrelPad->SetPad(0.0, 0.6, 1.0, 1.0);
		fTopPad->SetPad(0.0, 0.2, 0.5, 0.6);
		fBottomPad->SetPad(0.5, 0.2, 1.0, 0.6);
		fFitOverlayPad->SetPad(0.0, 0.0, 1.0, 0.2);
	}
	// Otherwise show the fit text
	else
	{
		fFitPad->SetPad(0.0, 0.0, 1.0, 1.0);
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
  else{
    fFitPad->Draw();
  }
	canvas->Modified();
	canvas->Update();
}

void WCSimRecoEvDisplay::DrawFitOverlays() {
  std::cout << "Drawing fit overlays" << std::endl;
	// Take the plots one by one and draw them.
	fBarrelPad->cd();
	fBarrelHist->Draw("colz");
  this->DrawHitGraphs(fBarrelGraphs);
	// Draw the fit rings
  std::cout << "There are " << fFitMarkersBarrel.size() << " barrel markers" << std::endl;
	for (unsigned int r = 0; r < fFitMarkersBarrel.size(); ++r) {
    std::cout << "Drawing fit marker " << r << "/" << fFitMarkersBarrel.size() << std::endl;
    std::cout << "Address " << fFitMarkersBarrel.at(r) << std::endl;
		fFitMarkersBarrel[r]->Draw("C");
	}
	fBarrelPad->Modified();
	fBarrelPad->Update();
  std::cout << "Done barrel" << std::endl;

	fTopPad->cd();
	fTopHist->Draw("colz");
  this->DrawHitGraphs(fTopGraphs);
	// Draw the fit rings
  std::cout << "There are " << fFitMarkersTop.size() << " top markers" << std::endl;
	for (unsigned int r = 0; r < fFitMarkersTop.size(); ++r) {
    std::cout << "Drawing fit marker " << r << "/" << fFitMarkersTop.size() << std::endl;
    std::cout << "Address " << fFitMarkersTop.at(r) << std::endl;
		fFitMarkersTop.at(r)->Draw("C");
	}
	fTopPad->Modified();
	fTopPad->Update();
  std::cout << "Done le top" << std::endl;

	fBottomPad->cd();
	fBottomHist->Draw("colz");
  this->DrawHitGraphs(fBottomGraphs);
	// Draw the truth rings
	for (unsigned int r = 0; r < fFitMarkersBottom.size(); ++r) {
		fFitMarkersBottom[r]->Draw("C");
	}
	fBottomPad->Modified();
	fBottomPad->Update();
  std::cout << "Done le bottom" << std::endl;
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

void WCSimRecoEvDisplay::SetPlotZAxes(){

  double min = fQMin;
  double max = fQMax;

  if(fViewType == 1){
    min = fTMin;
    max = fTMax;
  }
  if(fViewType >= 2){
    min = fLnLMin;
    max = fLnLMax;
  }

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
		std::cout << "Setting colour axis to total likelihood" << std::endl;
		this->FillPlots();
	}
	else{
		std::cout << "Already viewing total likelihood." << std::endl;
	}
}

// Switch the z-axis scale to show the charge likelihood.
void WCSimRecoEvDisplay::SetViewQLnL(){
	if(fViewType != 3){
		fViewType = 3;
    if(fWhichPads != 0){
      this->ShowReco();
    }
		std::cout << "Setting colour axis to charge likelihood" << std::endl;
		this->FillPlots();
	}
	else{
		std::cout << "Already viewing charge likelihood." << std::endl;
	}
}

// Switch the z-axis scale to show the time likelihood.
void WCSimRecoEvDisplay::SetViewTLnL(){
	if(fViewType != 4){
		fViewType = 4;
    if(fWhichPads != 0){
      this->ShowReco();
    }
		std::cout << "Setting colour axis to time likelihood" << std::endl;
		this->FillPlots();
	}
	else{
		std::cout << "Already viewing time likelihood." << std::endl;
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

// Calculate the bins for the likelihood plots
void WCSimRecoEvDisplay::CalculateLnLBins(){
  // Firstly, clear the existing vector
  fLnLBins.clear();
  
  double delta = (fLnLMax - fLnLMin) / 10.;

  for(int i = 0; i < 10; ++i){
    fLnLBins.push_back(fLnLMin+i*delta);
  }
}

unsigned int WCSimRecoEvDisplay::GetLnLBin(double lnl) const{
  unsigned int bin = 9;
  for(unsigned int i = 1; i < fLnLBins.size(); ++i){
    if(lnl < fLnLBins[i]){
      bin = i - 1;
      break;
    }
  }
  return bin;
}

