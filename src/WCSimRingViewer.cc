#include "WCSimRingViewer.hh"

#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMarker.h"
#include "TLegend.h"

#include "WCSimRecoEvent.hh"

ClassImp(WCSimRingViewer)

WCSimRingViewer::WCSimRingViewer()
{
  fCanvas = new TCanvas("wcsim_ringfit_canvas","WCSim RingViewer Display");
}

WCSimRingViewer::~WCSimRingViewer()
{
  if( fCanvas ) delete fCanvas;
}

void WCSimRingViewer::Initialize()
{

}

void WCSimRingViewer::Display(WCSimRecoEvent* myEvent)
{


  return;
}

