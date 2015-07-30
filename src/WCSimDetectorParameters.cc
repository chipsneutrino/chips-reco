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
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"

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

  double wavelengthAveragedQE = AverageHistWithGraph(spectrumHist, qeGraph);
  delete qeGraph;


  // If I run WCSim with QE set to Stacking_Only and to Sensitive_Detector_Only
  // the difference in number of tracks differs from my QE weighting
  // calculations by this factor - need to try and pin this down, but for now
  // just put in this constant
  wavelengthAveragedQE *= 2.528;

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
	// Take a TH1 and a TGraph
	// For every bin in the histogram, interpolate x between the two nearest points
	// on the graph to find the value at the bin centre, then multiply the bin contents
	// by that interpolated value.  Sum these, and divide by the number of hist entries.

	// If the histogram goes off the front or back of the graph, extrapolate graph
	double integral = hist->Integral();
	double weightedAverage = 0.0;
	for(int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin)
	{
		weightedAverage += hist->GetBinContent(iBin) * graph->Eval(hist->GetXaxis()->GetBinCenter(iBin));
	}
	weightedAverage /= integral;
	return weightedAverage;
}

TH1F* WCSimDetectorParameters::MultiplyHistByGraph(TH1F* hist, TGraph* graph) {
	for(int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin)
	{
		double newContent = hist->GetBinContent(iBin);
		newContent *= graph->Eval(hist->GetXaxis()->GetBinCenter(iBin));
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

