/*
 * WCSimDetectorParameters.cc
 *
 *  Created on: 29 Jul 2015
 *      Author: andy
 */

#include "WCSimDetectorParameters.hh"
#include "WCSimPMTConfig.hh"
#include "WCSimPMTManager.hh"
#include "WCSimTrackParameterEnums.hh"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimDetectorParameters)
#endif
/// Static pointer to self (so we can make object a singleton)
static WCSimDetectorParameters * fgDetectorParameters = NULL;

WCSimDetectorParameters * WCSimDetectorParameters::Instance()
{
	if( fgDetectorParameters == NULL){
		fgDetectorParameters = new WCSimDetectorParameters();
	}
	assert(fgDetectorParameters != NULL);
	return fgDetectorParameters;

}

WCSimDetectorParameters::WCSimDetectorParameters() : fPMTManager(){
	// TODO Auto-generated constructor stub
  fSpectrumFile = 0x0;
  OpenFile();
}

WCSimDetectorParameters::~WCSimDetectorParameters() {
	// TODO Auto-generated destructor stub
  fSpectrumFile->Close();
  delete fSpectrumFile;
}

double WCSimDetectorParameters::GetWavelengthAveragedQE(const std::string& pmtName) {
	if(!IsInMap(pmtName, &fAverageQEMap))
	{
		double averageQE = WorkOutAverageQE(pmtName);
		fAverageQEMap[pmtName] = averageQE;
	}
	return fAverageQEMap[pmtName];
}

double WCSimDetectorParameters::GetQEAveragedRefIndex(const std::string& pmtName) {
	if(!IsInMap(pmtName, &fQEAveragedRefIndexMap))
	{
		double averageRefIndex = WorkOutAverageRefIndex(pmtName);
		fQEAveragedRefIndexMap[pmtName] = averageRefIndex;
	}
	return fQEAveragedRefIndexMap[pmtName];
}

bool WCSimDetectorParameters::IsInMap(const std::string& pmtName,
		std::map<std::string, double>* map) {
	return(map->find(pmtName) != map->end());
}

double WCSimDetectorParameters::WorkOutAverageRefIndex(
		const std::string& pmtName) {

	TH1F * spectrumHist = 0x0;
	TGraph * qeGraph = new TGraph();

  TString spectrumName = Form("hSpectrum");
  fSpectrumFile->GetObject(spectrumName.Data(),spectrumHist);
	if(spectrumHist == 0x0)
  {
    std::cerr << "Couldn't get " << spectrumName << std::endl;
    assert(spectrumHist != 0x0);
  }

	WCSimPMTConfig config = fPMTManager.GetPMTByName(pmtName);
	std::vector<std::pair<double,double> > effVector = config.GetEfficiencyVector();

	std::vector<std::pair<double, double> >::const_iterator effItr = effVector.begin();
	while(effItr != effVector.end())
	{
		qeGraph->SetPoint(qeGraph->GetN(), effItr->first, effItr->second);
		effItr++;
	}

	TGraph * refIndGraph = 0x0;
	fSpectrumFile->GetObject("gRefIndex",refIndGraph);


	double refInd = AverageHistWithGraph(MultiplyHistByGraph(spectrumHist, qeGraph), refIndGraph);
  delete qeGraph;
	return refInd;
}

double WCSimDetectorParameters::WorkOutAverageQE(const std::string& pmtName) {
	
  TH1F * spectrumHist = 0x0;
	TGraph * qeGraph = new TGraph();
  TString spectrumName = Form("hSpectrum");
  fSpectrumFile->GetObject(spectrumName.Data(), spectrumHist);
	assert(spectrumHist != 0x0);

	WCSimPMTConfig config = fPMTManager.GetPMTByName(pmtName);
	std::vector<std::pair<double,double> > effVector = config.GetEfficiencyVector();

	std::vector<std::pair<double, double> >::const_iterator effItr = effVector.begin();
	while(effItr != effVector.end())
	{
		qeGraph->SetPoint(qeGraph->GetN(), effItr->first, effItr->second);
		effItr++;
	}

  double wavelengthAveragedQE = AverageHistWithGraph(spectrumHist, qeGraph, 280.0, 660.0);
  // WCSim has a hardcoded PMT QE cutoff for wavelengths less than 280nm and greater than 660nm
  // in WCSimStackingAction.cc and WCSimWCSD.cc
  // n.b. WCSim calls GetPMTQE with a lower cutoff of 240nm, but the function itself has 280nm hardcoded
  // as well

  delete qeGraph;


  return wavelengthAveragedQE;
}

void WCSimDetectorParameters::OpenFile() {
  if(fSpectrumFile == 0x0)
  {
    TString filename = TString(getenv("WCSIMANAHOME")) + "/config/spectra.root";
    fSpectrumFile = new TFile(filename.Data(), "READ");
	}
  return;
}

double WCSimDetectorParameters::AverageHistWithGraph(TH1F* hist,
		TGraph* graph) {
  double minX = TMath::MinElement(graph->GetN(),graph->GetX()); 
  double maxX = TMath::MaxElement(graph->GetN(),graph->GetX()); 
  return AverageHistWithGraph(hist, graph, minX, maxX);
}


double WCSimDetectorParameters::AverageHistWithGraph(TH1F* hist,
		TGraph* graph, const double &min, const double &max) {
	// Take a TH1 and a TGraph
	// For every bin in the histogram, interpolate x between the two nearest points
	// on the graph to find the value at the bin centre, then multiply the bin contents
	// by that interpolated value.  Sum these, and divide by the number of hist entries.
  // 
  // This function doesn't allow negative values (because it's for refrative indices and
  // QE) so will set negative values to zero
  //
  // Min and max are used to consider only a range - if x is less than min or greater
  // than max, the histogram value for this bin is multiplied by zero even if the graph
  // is nonzero at this point (the WCSim QE is forced to zero belo 280nm and above 660nm
  // even though the QE graph goes further than this)
  // If min and max are not set then they will be set to the minimum and maximum x value
  // of the graph (ie. won't have any effect)
  // 
	// If the histogram goes off the front or back of the graph, set it to zero for this bin


  double integral = hist->Integral();
	double weightedAverage = 0.0;
	for(int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin)
	{
    double x = hist->GetXaxis()->GetBinCenter(iBin);
    double graphVal = 0.0;
    if(x >= min && x <= max) 
    {
      graphVal = graph->Eval(hist->GetXaxis()->GetBinCenter(iBin));
    }
    if(graphVal < 0){ graphVal = 0.0; }
		weightedAverage += hist->GetBinContent(iBin) * graphVal;
	}
	weightedAverage /= integral;
	return weightedAverage;
}

TH1F* WCSimDetectorParameters::MultiplyHistByGraph(TH1F* hist, TGraph* graph) {
	for(int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin)
	{
		double newContent = hist->GetBinContent(iBin);
    double graphVal = graph->Eval(hist->GetXaxis()->GetBinCenter(iBin));
    if(graphVal < 0){ graphVal = 0.0; }
    newContent *= graphVal;
		hist->SetBinContent(iBin, newContent);
	}
	return hist;
}

double WCSimDetectorParameters::WavelengthAveragedQE(
		const std::string& pmtName) {
	return WCSimDetectorParameters::Instance()->GetWavelengthAveragedQE(pmtName);
}

double WCSimDetectorParameters::QEAveragedRefIndex(const std::string& pmtName) {
	return WCSimDetectorParameters::Instance()->GetQEAveragedRefIndex(pmtName);
}

double WCSimDetectorParameters::PMTExposeHeight(const std::string &pmtName){
  return WCSimDetectorParameters::Instance()->GetPMTExposeHeight(pmtName);
}

double WCSimDetectorParameters::GetPMTExposeHeight(const std::string &pmtName){
	if(!IsInMap(pmtName, &fExposeHeightMap))
	{
		double exposeHeight = WorkOutExposeHeight(pmtName);
		fExposeHeightMap[pmtName] = exposeHeight;
	}
	return fExposeHeightMap[pmtName];
}

double WCSimDetectorParameters::WorkOutExposeHeight(const std::string &pmtName){
	WCSimPMTConfig config = fPMTManager.GetPMTByName(pmtName);
  return config.GetExposeHeight();
}
