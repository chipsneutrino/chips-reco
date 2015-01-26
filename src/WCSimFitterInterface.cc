/*
 * WCSimFitterInterface.cc
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#include "WCSimFitterInterface.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimFitterTree.hh"
#include <TString.h>
#include <cassert>

ClassImp(WCSimFitterInterface)
static WCSimFitterInterface * fgFitterInterface = NULL;



WCSimFitterInterface::WCSimFitterInterface() :
		fFileName(""), fNumTracks(0), fFitter(0),
		fFitterPlots(0), fMakeFits(true), fMakeSurfaces(true){
      fFitterPlots = new WCSimFitterPlots();
      fFitterTree = new WCSimFitterTree();
      fFitter = NULL;
      Init();
	// TODO Auto-generated constructor stub
}

WCSimFitterInterface::~WCSimFitterInterface() {
	// TODO Auto-generated destructor stub
  if( fFitter != NULL) { delete fFitter; }
  if( fFitterPlots != NULL) { delete fFitterPlots; }
  if( fFitterTree != NULL ) { delete fFitterTree; }
}

void WCSimFitterInterface::Init()
{
  if( fFitter == NULL ) { fFitter = new WCSimLikelihoodFitter() ; }
  fFitter->SetFitterPlots( fFitterPlots );
  fFitter->SetFitterTree( fFitterTree );
}

WCSimFitterInterface* WCSimFitterInterface::Instance()
{
  if( !fgFitterInterface ){
    fgFitterInterface = new WCSimFitterInterface();
  }
  return fgFitterInterface;
}


//  void WCSimFitterInterface::SetFile(const char * fileName) {
//  	fFileName = TString(fileName);
//  }

void WCSimFitterInterface::FixParameter(unsigned int numTrack, const char* name, bool doIt) {
	WCSimFitterConfig::Instance()->FixTrackParameter(numTrack, name, true);
}

void WCSimFitterInterface::FreeParameter(unsigned int numTrack, const char* name, bool doIt) {
	WCSimFitterConfig::Instance()->FreeTrackParameter(numTrack, name, true);
}

void WCSimFitterInterface::SetNumTracks(unsigned int numTracks) {
	WCSimFitterConfig::Instance()->SetNumTracks(numTracks);
}

void WCSimFitterInterface::SetTrackType( unsigned int numTrack, const char * typeName)
{
  std::cout << "WCSimFitterInterface::SetTrackType(" << numTrack << ", " << typeName << ")" << std::endl;
	WCSimFitterConfig::Instance()->SetTrackType(numTrack, typeName);
}

void WCSimFitterInterface::JoinParametersTogether(unsigned int numTrack1, unsigned int numTrack2, const char * name)
{
	assert(numTrack1 < fNumTracks && numTrack2 < fNumTracks && numTrack1 != numTrack2);
	if( numTrack1 < numTrack2 ){ WCSimFitterConfig::Instance()->SetJoinParametersTogether(numTrack1, numTrack2, name);}
	else{ WCSimFitterConfig::Instance()->SetJoinParametersTogether(numTrack2, numTrack1, name); }


}
unsigned int WCSimFitterInterface::GetNumTracks() const {
	return WCSimFitterConfig::Instance()->GetNumTracks();
}

void WCSimFitterInterface::SetParameter(unsigned int numTrack, const char * name, double min, double max, double start, bool fixIt)
{
  WCSimFitterConfig::Instance()->SetParameter(numTrack, name, min, max, start, fixIt);
}

void WCSimFitterInterface::SetParMin(unsigned int numTrack, const char* name,
		double min) {
	WCSimFitterConfig::Instance()->SetParMin(numTrack, name, min);
}

void WCSimFitterInterface::SetParMax(unsigned int numTrack, const char* name,
		double max) {
	WCSimFitterConfig::Instance()->SetParMax(numTrack, name, max);
}

void WCSimFitterInterface::SetParStart(unsigned int numTrack, const char* name,
		double start) {
	WCSimFitterConfig::Instance()->SetParStart(numTrack, name, start);
}

void WCSimFitterInterface::SetParRange(unsigned int numTrack, const char* name,
		double min, double max) {
	WCSimFitterConfig::Instance()->SetParRange(numTrack, name, min, max);
}

Double_t WCSimFitterInterface::GetParMin(unsigned int numTrack, const char* name) {
	return WCSimFitterConfig::Instance()->GetParMin(numTrack, name);
}

Double_t WCSimFitterInterface::GetParMax(unsigned int numTrack, const char* name) {
	return WCSimFitterConfig::Instance()->GetParMax(numTrack, name);
}

Double_t WCSimFitterInterface::GetParStart(unsigned int numTrack,
		const char* name) {
	return WCSimFitterConfig::Instance()->GetParStart(numTrack, name);
}

std::pair<Double_t, Double_t> WCSimFitterInterface::GetParRange(unsigned int numTrack,
		const char* name) {
	return WCSimFitterConfig::Instance()->GetParRange(numTrack, name);
}

void WCSimFitterInterface::SetMakeFits(bool doIt) {
	fMakeFits = doIt;
}

bool WCSimFitterInterface::GetMakeFits() {
	return fMakeFits;
}

void WCSimFitterInterface::SetNumEventsToFit(int nEvts) {
	WCSimFitterConfig::Instance()->SetNumEventsToFit(nEvts);
	return;
}

int WCSimFitterInterface::GetNumEventsToFit() {
	return WCSimFitterConfig::Instance()->GetNumEventsToFit();
}

void WCSimFitterInterface::PlotForEachEvent(const char* name, Bool_t doIt) {
	fFitterPlots->SetPlotForEachEvent(name, doIt);
}

bool WCSimFitterInterface::GetPlotForEachEvent(const char* name) {
	return fFitterPlots->GetPlotForEachEvent(name);
}

void WCSimFitterInterface::PlotRecoMinusTrue(const char* name, Bool_t doIt) {
	fFitterPlots->SetPlotRecoMinusTrue(name);
}

Bool_t WCSimFitterInterface::GetPlotRecoMinusTrue(const char* name) {
	return fFitterPlots->GetPlotRecoMinusTrue(name);
}

void WCSimFitterInterface::SetMakeSurfaces( Bool_t doIt ) {
	fMakeSurfaces = doIt;
}

bool WCSimFitterInterface::GetMakeSurfaces() {
	return fMakeSurfaces;
}

void WCSimFitterInterface::SetNumSurfaceBins( unsigned int nBins ) {
  fFitterPlots->SetNumSurfaceBins( nBins );
}

unsigned int WCSimFitterInterface::GetNumSurfaceBins() const {
  return fFitterPlots->GetNumSurfaceBins();
}

void WCSimFitterInterface::Make1DSurface(unsigned int nTrack, const char* name, Bool_t doIt) {
	fFitterPlots->Make1DSurface(name, doIt, nTrack);
}

bool WCSimFitterInterface::GetMake1DSurface(unsigned int nTrack, const char* name) {
	return fFitterPlots->GetMake1DSurface(name, nTrack);
}

void WCSimFitterInterface::Make2DSurface(unsigned int nTrack, const char* name, unsigned int nTrack2, const char* name2, Bool_t doIt) {
	fFitterPlots->Make2DSurface(name, name2, doIt, nTrack, nTrack2 );
}

bool WCSimFitterInterface::GetMake2DSurface(unsigned int nTrack, const char* name,
		unsigned int nTrack2, const char* name2) {
	return fFitterPlots->GetMake2DSurface(name, name2, nTrack, nTrack2);
}

void WCSimFitterInterface::Print() {
  PrintFitConfiguration();
  PrintPlotsConfiguration();
  PrintSurfaceConfiguration();
}

void WCSimFitterInterface::PrintFitConfiguration() {
	WCSimFitterConfig::Instance()->Print();
}

void WCSimFitterInterface::PrintPlotsConfiguration() {
	fFitterPlots->Print();
}

void WCSimFitterInterface::PrintSurfaceConfiguration() {
	fFitterPlots->PrintSurfaces();
}

void WCSimFitterInterface::Run() {
  std::cout << " *** WCSimFitterInterface::Run() *** " << std::endl;
  Init();
  std::cout << "  Making histograms" << std::endl;
  fFitterPlots->MakeHistograms(WCSimFitterConfig::Instance());
  std::cout << "  Making tree" << std::endl;
  fFitterTree->MakeTree();

  std::cout << "  Running fits " << std::endl;
  if(fMakeFits) { fFitter->RunFits(); }
  std::cout << "  Running surfaces " << std::endl;
  if(fMakeSurfaces) { fFitter->RunSurfaces(); }
  std::cout << "  Saving plots " << std::endl;
  fFitterPlots->SavePlots();
  std::cout << "  Saving tree " << std::endl;
  fFitterTree->SaveTree(fFitterPlots->GetSaveFileName());
  std::cout << " *********************************** " << std::endl;
}
