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
#include "WCSimPiZeroFitter.hh"
#include <TString.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <cassert>

ClassImp(WCSimFitterInterface)

WCSimFitterInterface::WCSimFitterInterface() :
		fFitterConfig(NULL), fFileName(""), fNumTracks(0), fFitter(0x0), fPiZeroFitter(0x0),
		fFitterPlots(0x0), fFitterTree(0x0), fMakeFits(true), fMakeSurfaces(true){
	TTimeStamp ts;
	unsigned int year, month, day, hour, minute, second;
	ts.GetDate(true, 0, &year, &month, &day);
	ts.GetTime(true, 0, &hour, &minute, &second);
	TString saveName = Form("fitterPlots_%02d%02d%02d_%02d%02d%02d", year, month, day, hour, minute, second);

	fFitterConfig = new WCSimFitterConfig();
	fFitterPlots = new WCSimFitterPlots(saveName);
	fFitterTree = new WCSimFitterTree(saveName);
	// TODO Auto-generated constructor stub
}

WCSimFitterInterface::~WCSimFitterInterface() {
	// TODO Auto-generated destructor stub
  if( fFitter != NULL) { delete fFitter; }
  if( fFitterPlots != NULL) { delete fFitterPlots; }
  if( fFitterTree != NULL ) { delete fFitterTree; }
  if( fFitterConfig != NULL) { delete fFitterConfig; }
}

void WCSimFitterInterface::InitFitter()
{
  if( fFitter == 0x0 && !GetIsPiZeroFit()) { 
    fFitter = new WCSimLikelihoodFitter(fFitterConfig) ; 
    fFitter->SetFitterPlots(fFitterPlots);
    fFitter->SetFitterTree(fFitterTree);
  }
  else if( fPiZeroFitter == 0x0 && GetIsPiZeroFit() ) { 
    fPiZeroFitter = new WCSimPiZeroFitter(fFitterConfig) ; 
    fPiZeroFitter->SetFitterPlots(fFitterPlots);
    fPiZeroFitter->SetFitterTree(fFitterTree);
  }
}

void WCSimFitterInterface::FixParameter(unsigned int numTrack, const char* name, bool doIt) { 
  fFitterConfig->FixTrackParameter(numTrack, name, true);
}

void WCSimFitterInterface::FreeParameter(unsigned int numTrack, const char* name, bool doIt) {
	fFitterConfig->FreeTrackParameter(numTrack, name, true);
}

void WCSimFitterInterface::SetNumTracks(unsigned int numTracks) {
	fFitterConfig->SetNumTracks(numTracks);
}

void WCSimFitterInterface::SetTrackType( unsigned int numTrack, const char * typeName)
{
  //std::cout << "WCSimFitterInterface::SetTrackType(" << numTrack << ", " << typeName << ")" << std::endl;
	fFitterConfig->SetTrackType(numTrack, typeName);
}

void WCSimFitterInterface::JoinParametersTogether(unsigned int numTrack1, unsigned int numTrack2, const char * name)
{
	assert(numTrack1 != numTrack2);
	if( numTrack1 < numTrack2 ){ fFitterConfig->SetJoinParametersTogether(numTrack1, numTrack2, name);}
	else{ fFitterConfig->SetJoinParametersTogether(numTrack2, numTrack1, name); }


}
unsigned int WCSimFitterInterface::GetNumTracks() const {
	return fFitterConfig->GetNumTracks();
}

TrackType::Type WCSimFitterInterface::GetTrackType(const unsigned int &numTrack)
{
  return fFitterConfig->GetTrackType(numTrack);
}

void WCSimFitterInterface::SetParameter(unsigned int numTrack, const char * name, double min, double max, double start, bool fixIt)
{
  fFitterConfig->SetParameter(numTrack, name, min, max, start, fixIt);
}

void WCSimFitterInterface::SetParMin(unsigned int numTrack, const char* name,
		double min) {
	fFitterConfig->SetParMin(numTrack, name, min);
}

void WCSimFitterInterface::SetParMax(unsigned int numTrack, const char* name,
		double max) {
	fFitterConfig->SetParMax(numTrack, name, max);
}

void WCSimFitterInterface::SetParStart(unsigned int numTrack, const char* name,
		double start) {
	fFitterConfig->SetParStart(numTrack, name, start);
}

void WCSimFitterInterface::SetParRange(unsigned int numTrack, const char* name,
		double min, double max) {
	fFitterConfig->SetParRange(numTrack, name, min, max);
}

Double_t WCSimFitterInterface::GetParMin(unsigned int numTrack, const char* name) {
	return fFitterConfig->GetParMin(numTrack, name);
}

Double_t WCSimFitterInterface::GetParMax(unsigned int numTrack, const char* name) {
	return fFitterConfig->GetParMax(numTrack, name);
}

Double_t WCSimFitterInterface::GetParStart(unsigned int numTrack,
		const char* name) {
	return fFitterConfig->GetParStart(numTrack, name);
}

std::pair<Double_t, Double_t> WCSimFitterInterface::GetParRange(unsigned int numTrack,
		const char* name) {
	return fFitterConfig->GetParRange(numTrack, name);
}

void WCSimFitterInterface::SetMakeFits(bool doIt) {
	fFitterConfig->SetMakeFits(doIt);
}

bool WCSimFitterInterface::GetMakeFits() {
	return fFitterConfig->GetMakeFits();
}

void WCSimFitterInterface::SetNumEventsToFit(int nEvts) {
	fFitterConfig->SetNumEventsToFit(nEvts);
	return;
}

int WCSimFitterInterface::GetNumEventsToFit() {
	return fFitterConfig->GetNumEventsToFit();
}

void WCSimFitterInterface::SetFirstEventToFit(int iEvt)
{
	fFitterConfig->SetFirstEventToFit(iEvt);
}

int WCSimFitterInterface::GetFirstEventToFit() const
{
	return fFitterConfig->GetFirstEventToFit();
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
	fFitterConfig->Print();
}

void WCSimFitterInterface::PrintPlotsConfiguration() {
	fFitterPlots->Print();
}

void WCSimFitterInterface::PrintSurfaceConfiguration() {
	fFitterPlots->PrintSurfaces();
}

void WCSimFitterInterface::SetInputFileName(const char * inputfile)
{
	fFileName = inputfile;
}

void WCSimFitterInterface::SaveResults()
{
  //std::cout << "  Saving tree " << std::endl;
  fFitterTree->SaveTree();
  //std::cout << "  Saving plots " << std::endl;
  fFitterPlots->SavePlots();
}

void WCSimFitterInterface::SaveProfiles()
{
  //std::cout << "  Saving profiles " << std::endl;
  fFitterPlots->SaveProfiles();
  //std::cout << "  Saved profiles" << std::endl;
}

void WCSimFitterInterface::Run() {
  std::cout << " *** WCSimFitterInterface::Run() *** " << std::endl;
  std::cout << " *** InitOutputFiles *** " << std::endl;
  InitOutputFiles();
  std::cout << " *** InitFitter *** " << std::endl;
  InitFitter();

  std::cout << " *** Making histograms *** " << std::endl;
  fFitterPlots->MakeHistograms(fFitterConfig);
  std::cout << " *** Making tree *** " << std::endl;
  fFitterTree->MakeTree();

  std::cout << "  Running fits " << std::endl;
  if(fMakeFits) 
  { 
    if(fFitterConfig->GetIsPiZeroFit())
    {
      std::cout << "Pi zero" << std::endl;
      fPiZeroFitter->RunFits();
    }
    else
    {
      fFitter->RunFits(); 
    }
  }
  std::cout << "  Running surfaces " << std::endl;
  if(fMakeSurfaces) { fFitter->RunSurfaces(); }
  SaveProfiles();
  std::cout << "  Saving tree " << std::endl;
  fFitterTree->SaveTree();
  std::cout << "  Saving plots " << std::endl;
  fFitterPlots->SavePlots();
  std::cout << " *********************************** " << std::endl;
}

void WCSimFitterInterface::SetIsPiZeroFit(const bool &isPiZero)
{
  fFitterConfig->SetIsPiZeroFit(isPiZero);
  return;
}

bool WCSimFitterInterface::GetIsPiZeroFit() const
{
  return fFitterConfig->GetIsPiZeroFit();
}

void WCSimFitterInterface::SetForcePiZeroMass(const bool& doIt)
{
	fFitterConfig->SetForcePiZeroMass(doIt);
}

bool WCSimFitterInterface::GetForcePiZeroMass() const
{
	return fFitterConfig->GetForcePiZeroMass();
}

void WCSimFitterInterface::InitOutputFiles()
{
	TTimeStamp ts;
    unsigned int year, month, day, hour, minute, second;
	ts.GetDate(true, 0, &year, &month, &day);
	ts.GetTime(true, 0, &hour, &minute, &second);
    TString time = Form("%02d%02d%02d", hour, minute, second);

    TString basename = gSystem->BaseName(fFileName.Data());
    basename.ReplaceAll(".root","");

	TString saveNamePlots = Form("fit_%s_%04d_to_%04d_%s_plots.root", 
                                 basename.Data(), 
                                 fFitterConfig->GetFirstEventToFit(),
                                 fFitterConfig->GetNumEventsToFit() + fFitterConfig->GetFirstEventToFit(),
                                 time.Data());
	TString saveNameTree = Form("fit_%s_%04d_to_%04d_%s_tree.root", 
                                 basename.Data(), 
                                 fFitterConfig->GetFirstEventToFit(),
                                 fFitterConfig->GetNumEventsToFit() + fFitterConfig->GetFirstEventToFit(),
                                 time.Data());
	fFitterPlots->SetSaveFileName(saveNamePlots);
	fFitterTree->SetSaveFileName(saveNameTree);
    fFitterPlots->MakeSaveFile();
    fFitterTree->MakeSaveFile();
    return; 
}
