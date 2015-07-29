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
#include "TH1D.h"

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

}

WCSimDetectorParameters::~WCSimDetectorParameters() {
	// TODO Auto-generated destructor stub
  std::map<TrackType::Type, TFile *>::iterator fileMapItr = fFileMap.begin();
  while(fileMapItr != fFileMap.end())
  {
    (fileMapItr->second)->Close();
    delete fileMapItr->second;
    fileMapItr->second = 0x0;
  }
}

double WCSimDetectorParameters::GetWavelengthAveragedQE(
		const TrackType::Type& type, const std::string& pmtName) {
	if(!IsInMap(type, pmtName, &fAverageQEMap))
	{
		double averageQE = WorkOutAverageQE(type, pmtName);
		fAverageQEMap[type][pmtName] = averageQE;
	}
	return fAverageQEMap[type][pmtName];
}

double WCSimDetectorParameters::GetQEAveragedRefIndex(
		const TrackType::Type& type, const std::string& pmtName) {
	if(!IsInMap(type, pmtName, &fQEAveragedRefIndexMap))
	{
		double averageRefIndex = WorkOutAverageRefIndex(type, pmtName);
		fQEAveragedRefIndexMap[type][pmtName] = averageRefIndex;
	}
	return fQEAveragedRefIndexMap[type][pmtName];
}

bool WCSimDetectorParameters::IsInMap(const TrackType::Type& type,
		const std::string& pmtName,
		std::map<TrackType::Type, std::map<std::string, double> >* map) {
	if(map->find(type) != map->end())
	{
		if((*map)[type].find(pmtName) == (*map)[type].end())
		{
			return true;
		}
	}
	return false;
}

bool WCSimDetectorParameters::IsInMap(const TrackType::Type& type,
		std::map<TrackType::Type, TFile*> * map) {
	return (map->find(type) != map->end());
}

double WCSimDetectorParameters::WorkOutAverageRefIndex(
		const TrackType::Type& type, const std::string& pmtName) {

	if(!IsInMap(type, &fFileMap))
	{
		OpenFile(type);
	}

	TH1D * spectrumHist = 0x0;
	TGraph * qeGraph = 0x0;

	fFileMap[type]->GetObject("hSpectrum",spectrumHist);
	assert(spectrumHist != 0x0);

	WCSimPMTConfig config = fPMTManager.GetPMTByName(pmtName);
	std::vector<std::pair<double,double> > effVector = config.GetEfficiencyVector();

	std::vector<std::pair<double, double> >::const_iterator effItr = effVector.begin();
	while(effItr != effVector.end())
	{
		qeGraph->SetPoint(qeGraph->GetN(), effItr->first, effItr->second);
		effItr++;
	}

	TGraph * refIndGraph = 0x0;
	fFileMap[type]->GetObject("gRefIndex",refIndGraph);


	double refInd = AverageHistWithGraph(MultiplyHistByGraph(spectrumHist, qeGraph), refIndGraph);
	return refInd;
}

double WCSimDetectorParameters::WorkOutAverageQE(const TrackType::Type& type,
		const std::string& pmtName) {

	if(!IsInMap(type, &fFileMap))
	{
		OpenFile(type);
	}
	TH1D * spectrumHist = 0x0;
	TGraph * qeGraph = 0x0;
	fFileMap[type]->GetObject("hSpectrum",spectrumHist);
	assert(spectrumHist != 0x0);

	WCSimPMTConfig config = fPMTManager.GetPMTByName(pmtName);
	std::vector<std::pair<double,double> > effVector = config.GetEfficiencyVector();

	std::vector<std::pair<double, double> >::const_iterator effItr = effVector.begin();
	while(effItr != effVector.end())
	{
		qeGraph->SetPoint(qeGraph->GetN(), effItr->first, effItr->second);
		effItr++;
	}

	return AverageHistWithGraph(spectrumHist, qeGraph);
}

void WCSimDetectorParameters::OpenFile(const TrackType::Type& type) {
	TString filename("");
	switch(type)
	{
		case TrackType::MuonLike:
			filename = "spectrumMuon.root";
      break;
		case TrackType::ElectronLike:
			filename = "spectrumElectron.root";
      break;
		default:
			assert(type == TrackType::MuonLike || type == TrackType::ElectronLike);
	}
  TFile * f = new TFile(filename.Data(), "READ");
	fFileMap.insert(std::make_pair(type, f));
	return;
}

double WCSimDetectorParameters::AverageHistWithGraph(TH1D* hist,
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

TH1D* WCSimDetectorParameters::MultiplyHistByGraph(TH1D* hist, TGraph* graph) {
	for(int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin)
	{
		double newContent = hist->GetBinContent(iBin);
		newContent *= graph->Eval(hist->GetXaxis()->GetBinCenter(iBin));
		hist->SetBinContent(iBin, newContent);
	}
	return hist;
}

double WCSimDetectorParameters::WavelengthAveragedQE(
		const TrackType::Type& type, const std::string& pmtName) {
	return WCSimDetectorParameters::Instance()->GetWavelengthAveragedQE(type, pmtName);
}

double WCSimDetectorParameters::QEAveragedRefIndex(const TrackType::Type& type,
		const std::string& pmtName) {
	return WCSimDetectorParameters::Instance()->GetQEAveragedRefIndex(type, pmtName);
}

