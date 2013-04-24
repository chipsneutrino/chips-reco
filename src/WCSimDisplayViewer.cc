#include "WCSimDisplayViewer.hh"

#include "WCSimDisplayFactory.hh"
#include "WCSimDisplay.hh"

#include "WCSimRecoFactory.hh"
#include "WCSimReco.hh"
#include "WCSimRecoEvent.hh"

#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"

#include <iostream>
#include <string>
#include <cassert>

ClassImp(WCSimDisplayViewer)

WCSimDisplayViewer::WCSimDisplayViewer() :  
  fRootStyle(0),
  fControlCanvas(0),
  buttonStart(0),
  buttonNext(0),
  buttonPrevious(0),
  buttonGoTo(0),
  buttonPrint(0),
  buttonQuit(0),
  fCounter(-1),
  fReco(0),
  fRecoEvent(0),
  fTrueEvent(0),
  fDisplayRecoClusters(0),
  fDisplayRecoEvent(0),
  fDisplayTrueEvent(0),
  fPulseHeightCut(0.0),
  fPrintGIF(0),
  fPrintEPS(0),
  fDisplayReady(0)
{
  this->Initialize(); 
}

WCSimDisplayViewer::~WCSimDisplayViewer()
{  
  if(fControlCanvas) delete fControlCanvas;
  if(buttonStart) delete buttonStart;
  if(buttonNext) delete buttonNext;
  if(buttonPrevious) delete buttonPrevious;
  if(buttonGoTo) delete buttonGoTo;
  if(buttonPrint) delete buttonPrint;
  if(buttonQuit) delete buttonQuit;

  if(fReco) delete fReco;
}

void WCSimDisplayViewer::RootStyle()
{
  if( fRootStyle ) return;

  std::cout << " *** WCSimDisplayViewer::ForceRootStyle() *** " << std::endl;

  // create new TStyle
  fRootStyle = new  TStyle("My Style", "");

  // set the background color to white
  fRootStyle->SetFillColor(10);
  fRootStyle->SetFrameFillColor(10);
  fRootStyle->SetCanvasColor(10);
  fRootStyle->SetPadColor(10);
  fRootStyle->SetTitleFillColor(0);
  fRootStyle->SetStatColor(10);

  // don't put a colored frame around the plots
  fRootStyle->SetFrameBorderMode(0);
  fRootStyle->SetCanvasBorderMode(0);
  fRootStyle->SetPadBorderMode(0);
  fRootStyle->SetLegendBorderSize(0);

  // use the primary color palette
  fRootStyle->SetPalette(1,0);

  // set the default line color for a histogram to be black
  fRootStyle->SetHistLineColor(kBlack);

  // set the default line color for a fit function to be red
  fRootStyle->SetFuncColor(kRed);

  // make the axis labels black
  fRootStyle->SetLabelColor(kBlack,"xyz");

  // set the default title color to be black
  fRootStyle->SetTitleColor(kBlack);

  // set the margins
  fRootStyle->SetPadBottomMargin(0.18);
  fRootStyle->SetPadTopMargin(0.08);
  fRootStyle->SetPadRightMargin(0.08);
  fRootStyle->SetPadLeftMargin(0.17);

  // set axis label and title text sizes
  fRootStyle->SetLabelFont(42,"xyz");
  fRootStyle->SetLabelSize(0.06,"xyz");
  fRootStyle->SetLabelOffset(0.015,"xyz");
  fRootStyle->SetTitleFont(42,"xyz");
  fRootStyle->SetTitleSize(0.07,"xyz");
  fRootStyle->SetTitleOffset(1.1,"yz");
  fRootStyle->SetTitleOffset(1.0,"x");
  fRootStyle->SetStatFont(42);
  fRootStyle->SetStatFontSize(0.07);
  fRootStyle->SetTitleBorderSize(0);
  fRootStyle->SetStatBorderSize(0);
  fRootStyle->SetTextFont(42);

  // set line widths
  fRootStyle->SetFrameLineWidth(2);
  fRootStyle->SetFuncWidth(2);

  // set the number of divisions to show
  fRootStyle->SetNdivisions(506, "xy");

  // turn off xy grids
  fRootStyle->SetPadGridX(0);
  fRootStyle->SetPadGridY(0);

  // set the tick mark style
  fRootStyle->SetPadTickX(1);
  fRootStyle->SetPadTickY(1);

  // turn off stats
  fRootStyle->SetOptStat(0);

  // marker settings
  fRootStyle->SetMarkerStyle(20);
  fRootStyle->SetMarkerSize(0.6);
  fRootStyle->SetLineWidth(1);

  // done
  fRootStyle->cd();
  gROOT->ForceStyle();
  gStyle->ls();

  return;
}

void WCSimDisplayViewer::Initialize()
{
  std::cout << " *** WCSimDisplayViewer::Initialize() *** " << std::endl;

  // ROOT style
  this->RootStyle();

  // initializing control canvas
  if( !fControlCanvas ){
    fControlCanvas = new TCanvas("ControlCanvas",
                                 "ControlCanvas",0,0,150,600);
    fControlCanvas->cd();

    Int_t butNumber = 6;
    Double_t butLowerY,butUpperY,butSpaceY;

    butSpaceY = 1.0/((double)(butNumber));

    butUpperY = 1.0;
    butLowerY = butUpperY-butSpaceY;
    buttonStart = new TButton("Start","viewer->controlStart()",0.0,butLowerY,1.0,butUpperY);
    buttonStart->SetTextSize(0.33);
    buttonStart->SetTextColor(1);
    buttonStart->SetFillColor(9);
    buttonStart->Draw();

    butUpperY = butLowerY;
    butLowerY = butUpperY-butSpaceY;
    buttonNext = new TButton("Next","viewer->controlNext()",0.0,butLowerY,1.0,butUpperY);
    buttonNext->SetTextSize(0.33);
    buttonNext->SetTextColor(1);
    buttonNext->SetFillColor(9);
    buttonNext->Draw();

    butUpperY = butLowerY;
    butLowerY = butUpperY-butSpaceY;
    buttonPrevious = new TButton("Previous","viewer->controlPrevious()",0.0,butLowerY,1.0,butUpperY);
    buttonPrevious->SetTextSize(0.33);
    buttonPrevious->SetTextColor(1);
    buttonPrevious->SetFillColor(9);
    buttonPrevious->Draw();

    butUpperY = butLowerY;
    butLowerY = butUpperY-butSpaceY;
    buttonGoTo = new TButton("GoTo","viewer->controlGoTo()",0.0,butLowerY,1.0,butUpperY);
    buttonGoTo->SetTextSize(0.33);
    buttonGoTo->SetTextColor(1);
    buttonGoTo->SetFillColor(9);
    buttonGoTo->Draw();

    butUpperY = butLowerY;
    butLowerY = butUpperY-butSpaceY;
    buttonPrint = new TButton("Print","viewer->controlPrint()",0.0,butLowerY,1.0,butUpperY);
    buttonPrint->SetTextSize(0.33);
    buttonPrint->SetTextColor(1);
    buttonPrint->SetFillColor(9);
    buttonPrint->Draw();

    butUpperY = butLowerY;
    butLowerY = butUpperY-butSpaceY;
    buttonQuit = new TButton("Quit",".q",0.0,butLowerY,1.0,butUpperY);
    buttonQuit->SetTextSize(0.33);
    buttonQuit->SetTextColor(1);
    buttonQuit->SetFillColor(9);
    buttonQuit->Draw();
  }
  
  fRunNum = -1;
  fEventNum = -1;
  fTriggerNum = -1;
  
  fPulseHeightCut = 0.0;

  fPrintGIF = 1;
  fPrintEPS = 0;

  fDisplayReady = 0;

  return;
}

void WCSimDisplayViewer::UseDisplay(const char* newdisplay)
{
  Bool_t foundit = 0;

  for( UInt_t n=0; n<fDisplayNameList.size(); n++ ){
    const char* olddisplay = fDisplayNameList.at(n);
    if( strcmp(olddisplay,newdisplay)==0 ) foundit = 1;
  }

  if( foundit==0 ){
    fDisplayNameList.push_back(newdisplay);
  }

  return;
}

void WCSimDisplayViewer::PrintDisplayList()
{
  std::cout << " *** WCSimDisplayViewer::PrintDisplayList() *** " << std::endl;  
  for( UInt_t n=0; n<fDisplayNameList.size(); n++ ){
    const char* displayname = fDisplayNameList.at(n);
    std::cout << "  [" << n << "] " << displayname << std::endl;
  }
}

void WCSimDisplayViewer::SetPulseHeightCut(Double_t cut)
{
  fPulseHeightCut = cut;

  for( UInt_t n=0; n<fDisplayList.size(); n++ ){
    WCSimDisplay* myDisplay = (WCSimDisplay*)(fDisplayList.at(n));
    myDisplay->SetPulseHeightCut(fPulseHeightCut);
  }

  return;
}
  
void WCSimDisplayViewer::PrintGIF(Bool_t yesno)
{ 
  fPrintGIF = yesno; 

  for( UInt_t n=0; n<fDisplayList.size(); n++ ){
    WCSimDisplay* myDisplay = (WCSimDisplay*)(fDisplayList.at(n));
    myDisplay->PrintGIF(fPrintGIF);
  }

  return;
}
  
void WCSimDisplayViewer::PrintEPS(Bool_t yesno)
{ 
  fPrintEPS = yesno; 

  for( UInt_t n=0; n<fDisplayList.size(); n++ ){
    WCSimDisplay* myDisplay = (WCSimDisplay*)(fDisplayList.at(n));
    myDisplay->PrintEPS(fPrintEPS);
  }

  return;
}

void WCSimDisplayViewer::BuildDisplay()
{
  std::cout << " *** WCSimDisplayViewer::BuildDisplay() *** " << std::endl;

  // Building Displays
  if( !fDisplayReady ){

    if( fDisplayNameList.size()==0 ){
      std::cout << " --- using default display --- " << std::endl;
      this->UseDisplay("default");
    }

    for( UInt_t n=0; n<fDisplayNameList.size(); n++ ){
      const char* displayname = fDisplayNameList.at(n);

      std::cout << "   Building display: " << displayname << std::endl;
      WCSimDisplay* myDisplay = WCSimDisplayFactory::Instance()->MakeDisplay(displayname);

      if( myDisplay ){
        myDisplay->SetPulseHeightCut(fPulseHeightCut);
        myDisplay->PrintGIF(fPrintGIF);
        myDisplay->PrintEPS(fPrintEPS);
        fDisplayList.push_back(myDisplay);
        fDisplayReady = 1;
      }
      else{
        std::cout << " --- no displays available: EXIT --- " << std::endl;
        assert(myDisplay);
      }
    }
  }

  // Building Reconstruction
  if( !fReco ){
    fReco = WCSimRecoFactory::Instance()->MakeReco();  
  }

  return;
}

void WCSimDisplayViewer::DisplayEvent(Int_t entry)
{ 
  if( !fDisplayReady ){
    std::cout << " --- starting display --- " << std::endl;
    this->BuildDisplay();
  }

  if( fDisplayReady ){
    this->UpdateDisplay(entry);
  }
}

void WCSimDisplayViewer::UpdateDisplay(Int_t entry)
{
  std::cout << " --- new event --- " << std::endl;

  WCSimInterface::LoadEvent(entry);

  if( WCSimInterface::TouchEvent() ){
    
    // create new event    
    fRecoEvent = WCSimInterface::RecoEvent();
    fTrueEvent = WCSimInterface::TrueEvent();

    // get run/event/trigger
    fRunNum = fRecoEvent->GetRun();
    fEventNum = fRecoEvent->GetEvent();
    fTriggerNum = fRecoEvent->GetTrigger();

    // run reconstruction
    if( fDisplayRecoClusters
     || fDisplayRecoEvent ){
      fReco->RunFilter(fRecoEvent);
    }

    if( fDisplayRecoEvent ){
      fReco->RunRecoVertex(fRecoEvent);
      fReco->RunRecoRings(fRecoEvent);
    } 

    // run displays  
    std::cout << " *** WCSimDisplayViewer::UpdateDisplay(...) *** " << std::endl;
    std::cout << "  Displaying Entry [" << entry << "] " << std::endl;
    std::cout << "  Displaying Event: " << fRunNum << "/" << fEventNum << "/" << fTriggerNum << std::endl;
    for( UInt_t n=0; n<fDisplayList.size(); n++ ){
      WCSimDisplay* myDisplay = (WCSimDisplay*)(fDisplayList.at(n));

      if( fDisplayRecoClusters ){
        myDisplay->DrawCleanDisplay(fRecoEvent);
      }
      else{
        myDisplay->DrawDisplay(fRecoEvent);
      }

      if( fDisplayTrueEvent ){
        myDisplay->DrawTrueEvent(fTrueEvent);
      }     

      if( fDisplayRecoEvent ){
        myDisplay->DrawRecoEvent(fRecoEvent);
      }
    }

    std::cout << "  Displaying Entry [" << entry << "]  DONE! " << std::endl; 
  }

  return;

}

void WCSimDisplayViewer::controlStart()
{
  if( fDisplayReady ){
    std::cout << " --- display running --- " << std::endl;
    return;
  }

  if( WCSimInterface::TouchData()==0 ){
    return;
  }

  if( WCSimGeometry::TouchGeometry()==0 ){
    return;
  }

  if( !fDisplayReady ){
    std::cout << " --- starting display --- " << std::endl;
    this->BuildDisplay();
    if( fDisplayReady ) this->controlNext();
  }

  return;
}

void WCSimDisplayViewer::controlNext()
{
  Bool_t display = 0;
  Bool_t displaynext = 0;

  if( fDisplayReady ){
    if( fCounter<WCSimInterface::GetNumEvents()-1 ){
      fCounter++;
      display = 1;
      displaynext = 1;
    }
    if( display ){
      this->UpdateDisplay(fCounter);
    }
    else{
      if( displaynext ) this->controlNext();
      else std::cout << " --- end of file --- " << std::endl;
    }
  }

  else{
    this->controlStart();
  }

  return;
}

void WCSimDisplayViewer::controlPrevious()
{ 
  Bool_t display = 0;
  Bool_t displayprev = 0;

  if( fDisplayReady ){
    if( fCounter>0){
      fCounter--;
      display = 1;
      displayprev = 1;
    }
    if( display ){
      this->UpdateDisplay(fCounter);
    }
    else{
      if( displayprev ) this->controlPrevious();
      else std::cout << " --- beginning of file --- " << std::endl;
    }
  }
  return;
}

void WCSimDisplayViewer::controlGoTo()
{
  if( fDisplayReady ){
    std::cout << " --- GOTO --- " << std::endl;

    // goto: run
    // --- run numbers not implemented yet ---
    /*
    Int_t run = 0;
    std::cout << "  Enter Run : "; 
    std::cin >> run; 
    std::cout << "  [confirm:" << run << "] " << std::endl;
    */

    // go to: event
    Int_t event = -1;
    std::cout << "  Enter Event : "; 
    std::cin >> event; 
    std::cout << "  [confirm:" << event << "] " << std::endl;

    // search for event
    std::cout << " Searching for event: " << event << std::endl;

    Bool_t foundit = 0;

    Int_t oldCounter = fCounter;
    Int_t newCounter = fCounter;

    Int_t tempCounter = 0;
    Int_t tempEvent = -1;

    // go forwards
    if( !foundit
     && event>fEventNum ){
      tempCounter=oldCounter;
      tempEvent=fEventNum;
      while( tempCounter<WCSimInterface::GetNumEvents()-1
       && (tempEvent<0||tempEvent<event)
       && !(tempEvent==event) ){
        tempCounter++;
        WCSimInterface::LoadEvent(tempCounter);
        if( WCSimInterface::TouchEvent() ){
          tempEvent = WCSimInterface::GetEventNumber();
	}
      }
      if( tempEvent==event ){
        newCounter=tempCounter; foundit=1;
      }
    }

    // go backwards
    if( !foundit
     && event<fEventNum ){
      tempCounter=oldCounter;
      tempEvent=fEventNum;
      while( tempCounter>0
       && (tempEvent<0||tempEvent>event)
       && !(tempEvent==event) ){
        tempCounter--;
        WCSimInterface::LoadEvent(tempCounter);
        if( WCSimInterface::TouchEvent() ){
          tempEvent = WCSimInterface::GetEventNumber();
	}
      }
      if( tempEvent==event ){
        newCounter=tempCounter; foundit=1;
      }
    }

    // display new event
    if( foundit ){
      fCounter = newCounter;
      this->UpdateDisplay(fCounter);
    }
    else{
      std::cout << "  <warning> failed to find event " << std::endl;
    }
  }

  return;
}

void WCSimDisplayViewer::controlPrint()
{
  if( fDisplayReady ){
    std::cout << "  saving images of event display " << std::endl;
    for( UInt_t n=0; n<fDisplayList.size(); n++ ){
      WCSimDisplay* myDisplay = (WCSimDisplay*)(fDisplayList.at(n));
      myDisplay->PrintDisplay();
    }
  }

  return;
}
