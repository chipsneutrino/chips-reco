/*
 * WCSimFitterPlots.cc
 *
 *  Created on: 3 Nov 2014
 *      Author: andy
 */

#include "WCSimFitterConfig.hh"
#include "WCSimFitterPlots.hh"
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTimeStamp.h>

WCSimFitterPlots::WCSimFitterPlots() : fSaveFile(NULL), fNumSurfaceBins(20){

	// TODO Auto-generated constructor stub
	TTimeStamp ts;
	unsigned int year, month, day, hour, minute, second;
	ts.GetDate(true, 0, &year, &month, &day);
	ts.GetTime(true, 0, &hour, &minute, &second);
	fSaveFileName.Form("fitterPlots_%d%d%d_%d%d.root", year, month, day, hour, minute);

}

WCSimFitterPlots::~WCSimFitterPlots() {
	// TODO Auto-generated destructor stub
}

void WCSimFitterPlots::SetPlotForEachEvent(const char* name, bool doIt) {
  if(doIt){ std::cout << " *** SetPlotForEachEvent *** ... plotting " << name << std::endl; }
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	std::vector<TH1D*> plotVec;
	fForEachEvent.insert(std::make_pair(type,plotVec));
}

bool WCSimFitterPlots::GetPlotForEachEvent(const char* name) const {
	return (fForEachEvent.find(FitterParameterType::FromName(name)) != fForEachEvent.end());
}

void WCSimFitterPlots::SetPlotRecoMinusTrue(const char* name, bool doIt) {
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	std::vector<TH1D*> plotVec;
	fRecoMinusTrue.insert(std::make_pair(type, plotVec));
}

bool WCSimFitterPlots::GetPlotRecoMinusTrue(const char* name) const {
	return (fRecoMinusTrue.find(FitterParameterType::FromName(name)) != fRecoMinusTrue.end());
}

void WCSimFitterPlots::Make1DSurface(const char* name, bool doIt, unsigned int trackNum) {
	std::pair<unsigned int, FitterParameterType::Type> toProfile = std::make_pair(trackNum, FitterParameterType::FromName(name));
	TH1D* plotVec = 0;
	fSurfaces1D.insert(std::make_pair(toProfile, plotVec));
}

bool WCSimFitterPlots::GetMake1DSurface(const char* name, unsigned int trackNum) const {
	return (fSurfaces1D.find(std::make_pair(trackNum, FitterParameterType::FromName(name))) != fSurfaces1D.end());
}

void WCSimFitterPlots::Make2DSurface(const char* name, const char* name2,
		bool doIt, unsigned int trackNum, unsigned int trackNum2) {

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	FitterParameterType::Type type2 = FitterParameterType::FromName(name2);

	std::pair<unsigned int, FitterParameterType::Type> toProfile1 = std::make_pair(trackNum, type);
	std::pair<unsigned int, FitterParameterType::Type> toProfile2 = std::make_pair(trackNum2, type2);

	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > types;
	if( trackNum < trackNum2){ types = std::make_pair(toProfile1, toProfile2); }
	else if( trackNum2 < trackNum ) { types = std::make_pair(toProfile2, toProfile1);}
	else if( type < type2) { types = std::make_pair(toProfile1, toProfile2); }
	else{ types = std::make_pair(toProfile2, toProfile1); }

	TH2D* plotHist = 0;
	fSurfaces2D.insert(std::make_pair(types, plotHist));
}

bool WCSimFitterPlots::GetMake2DSurface(const char* name,
		const char* name2, unsigned int trackNum, unsigned int trackNum2) const {
	    std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > types;
      std::pair<unsigned int, FitterParameterType::Type> type1, type2;
      type1 = std::make_pair( trackNum, FitterParameterType::FromName(name));
      type2 = std::make_pair( trackNum2, FitterParameterType::FromName(name2) );
      types = std::make_pair(type1, type2);
	return (fSurfaces2D.find(types) != fSurfaces2D.end());
}

void WCSimFitterPlots::SetNumSurfaceBins(unsigned int numBins) {
	fNumSurfaceBins = numBins;
}

unsigned int WCSimFitterPlots::GetNumSurfaceBins() {
	return fNumSurfaceBins;
}

void WCSimFitterPlots::MakeHistograms(WCSimFitterConfig * fitterConfig) {
  std::cout << " *** WCSimFitterPlots::MakeHistograms() *** " << std::endl;
  std::cout << "  MakePlotsForEachEvent " << std::endl;
	MakePlotsForEachEvent(fitterConfig);
  std::cout << "  MakeRecoMinusTrue " << std::endl;
	MakeRecoMinusTrue(fitterConfig);
  std::cout << "  Make1DSurfaces " << std::endl;
	Make1DSurfaces(fitterConfig);
  std::cout << "  Make2DSurfaces " << std::endl;
	Make2DSurface(fitterConfig);
  std::cout << "  Done!" << std::endl;
	return;
}

void WCSimFitterPlots::MakePlotsForEachEvent(WCSimFitterConfig* fitterConfig) {
	std::map<FitterParameterType::Type, std::vector<TH1D *> >::iterator plotItr = fForEachEvent.begin();
	for( ; plotItr != fForEachEvent.end(); ++plotItr)
	{
    (*plotItr).second.clear();
		for(unsigned int iTrack = 0; iTrack < fitterConfig->GetNumTracks() ; ++iTrack)
		{
			TH1D * hist = NULL;
			FitterParameterType::Type type = ((*plotItr).first);
			const char * parName = FitterParameterType::AsString(type).c_str();
			if( hist != NULL ){ delete hist; }

			TString name, title;
			name.Form("h%s_track%d", parName, iTrack);
			title.Form("Value of %s for track %d, for each of %d fitted events",
					   parName, iTrack, fitterConfig->GetNumTracks());

			int numBins = 20;
			if( fitterConfig->GetNumEventsToFit() > 400)
			{
				numBins = fitterConfig->GetNumEventsToFit() / 20;
				numBins += (10 - (numBins % 10));
			}
			hist = new TH1D(title.Data(), name.Data(), numBins,
							fitterConfig->GetParMin(iTrack, parName),
							fitterConfig->GetParMax(iTrack, parName));

			hist->GetXaxis()->SetTitle(FitterParameterType::AsString(type).c_str());
			hist->GetXaxis()->CenterTitle();
			hist->GetXaxis()->SetLabelSize(0.05);
			hist->GetXaxis()->SetTitleSize(0.05);
			hist->GetYaxis()->SetTitle("Events");
			hist->GetYaxis()->SetLabelSize(0.05);
			hist->GetYaxis()->SetTitleSize(0.05);
			hist->GetYaxis()->CenterTitle();
			hist->SetLineWidth(2);
			hist->SetLineColor(kRed);
      (*plotItr).second.push_back(hist);
		}

	}
}

void WCSimFitterPlots::MakeRecoMinusTrue(WCSimFitterConfig* fitterConfig) {
	std::map<FitterParameterType::Type, std::vector<TH1D *> >::iterator plotItr = fRecoMinusTrue.begin();
	for( ; plotItr != fRecoMinusTrue.end(); ++plotItr)
	{
    (*plotItr).second.clear();
		for(unsigned int iTrack = 0; iTrack < fitterConfig->GetNumTracks() ; ++iTrack)
		{
			TH1D * hist = NULL;;
			FitterParameterType::Type type = ((*plotItr).first);
			std::string parName = FitterParameterType::AsString(type).c_str();
			if( hist != NULL ){ delete hist; }

			TString name, title;
			name.Form("h_%s_track%d_reco_minus_true", parName.c_str(), iTrack);
			title.Form("Reco - true value of %s for track %d, for each of %d fitted events",
					   parName.c_str(), iTrack, fitterConfig->GetNumTracks());

			int numBins = 20;
			if( fitterConfig->GetNumEventsToFit() > 400)
			{
				numBins = fitterConfig->GetNumEventsToFit() / 20;
				numBins += (10 - (numBins % 10));
			}
			hist = new TH1D(title.Data(), name.Data(), numBins,
							2.0 * fitterConfig->GetParMin(iTrack, parName.c_str()),
							2.0 * fitterConfig->GetParMax(iTrack, parName.c_str()));

			TString xTitle;
			xTitle.Form("Reco - true: %s", FitterParameterType::AsString(type).c_str());
			hist->GetXaxis()->SetTitle(xTitle.Data());
			hist->GetXaxis()->CenterTitle();
			hist->GetXaxis()->SetLabelSize(0.05);
			hist->GetXaxis()->SetTitleSize(0.05);
			hist->GetYaxis()->SetTitle("Events");
			hist->GetYaxis()->SetLabelSize(0.05);
			hist->GetYaxis()->SetTitleSize(0.05);
			hist->GetYaxis()->CenterTitle();
			hist->SetLineWidth(2);
			hist->SetLineColor(kBlue);
      (*plotItr).second.push_back(hist);
		}
	}
}

void WCSimFitterPlots::Make1DSurfaces(WCSimFitterConfig* fitterConfig) {
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D *>::iterator plotItr = fSurfaces1D.begin();
	for( ; plotItr != fSurfaces1D.end(); ++plotItr)
	{
		TH1D * hist = NULL;
		FitterParameterType::Type type = ((*plotItr).first).second;
		const char * parName = FitterParameterType::AsString(type).c_str();
    unsigned int trackNum = ((*plotItr).first).first;


		TString name, title;
		name.Form("h_%s_likelihood_profile", parName);
		title.Form("Profile of likelihood when %s is varied for track %d, with other parameters fixed",
				   parName, trackNum);

		hist = new TH1D(title.Data(), name.Data(), fNumSurfaceBins,
						fitterConfig->GetParMin(trackNum, parName),
						fitterConfig->GetParMax(trackNum, parName));
		TString xTitle;
		xTitle.Form("%s", FitterParameterType::AsString(type).c_str());
		hist->GetXaxis()->SetTitle(xTitle.Data());
		hist->GetXaxis()->CenterTitle();
		hist->GetXaxis()->SetLabelSize(0.05);
		hist->GetXaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitle("-2 Ln(L)");
		hist->GetYaxis()->SetLabelSize(0.05);
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->CenterTitle();
		hist->SetLineWidth(2);
		hist->SetLineColor(kMagenta+2);

		if( (*plotItr).second != NULL ){ delete (*plotItr).second; }
    (*plotItr).second = hist;

	}
}

void WCSimFitterPlots::Make2DSurface(WCSimFitterConfig* fitterConfig) {

	std::map< std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> >, TH2D *>::iterator plotItr = fSurfaces2D.begin();
	for( ; plotItr != fSurfaces2D.end(); ++plotItr)
	{
		TH2D * hist = NULL;
    
	  std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > types = ((*plotItr).first);

		const char * parName1 = FitterParameterType::AsString(types.first.second).c_str();
		const char * parName2 = FitterParameterType::AsString(types.second.second).c_str();
		unsigned int trackNum1 = (types.first).first;
		unsigned int trackNum2 = (types.second).first;

		TString name, title;
		name.Form("h_%s_vs_%s_likelihood_profile", parName2, parName1);
		title.Form("Profile of likelihood when %s & %s are varied for tracks %d and %d, with other parameters fixed",
				   parName2, parName1, trackNum2, trackNum1);

		hist = new TH2D(title.Data(), name.Data(),
						fNumSurfaceBins,
						fitterConfig->GetParMin(trackNum1, parName1),
						fitterConfig->GetParMax(trackNum1, parName1),
						fNumSurfaceBins,
						fitterConfig->GetParMin(trackNum2, parName2),
						fitterConfig->GetParMax(trackNum2, parName2));

		hist->GetXaxis()->SetTitle(parName1);
		hist->GetXaxis()->CenterTitle();
		hist->GetXaxis()->SetLabelSize(0.05);
		hist->GetXaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitle(parName2);
		hist->GetYaxis()->SetLabelSize(0.05);
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->CenterTitle();
		hist->GetZaxis()->SetTitle("-2 Ln(L)");
		hist->SetOption("COLZ");
		
    if( (*plotItr).second != NULL ){ delete (*plotItr).second; }
    (*plotItr).second = hist;
	}
}

void WCSimFitterPlots::Print() {
	std::cout << " -- WCSimFitterPlots: Listing all plots -- " << std::endl;
	PrintVariables();
	PrintRecoMinusTrue();
	Print1DSurfaces();
	Print2DSurfaces();
}

void WCSimFitterPlots::PrintVariables() {
	std::cout << " -- WCSimFitterPlots: Variables to be plotted -- " << std::endl;
	std::map<FitterParameterType::Type, std::vector<TH1D*> >::const_iterator plotItr = fForEachEvent.begin();
	for( ; plotItr != fForEachEvent.end(); ++plotItr)
	{
		FitterParameterType::Type type = ((*plotItr).first);
		for(unsigned int iHist = 0; iHist < ((*plotItr).second).size(); ++ iHist)
		{
			TH1D * hist = ((*plotItr).second).at(iHist);
			std::cout << "\t Plotting " << FitterParameterType::AsString(type) << " in "
					<< hist->GetNbinsX() << " bins from " << hist->GetXaxis()->GetXmin() << " to "
					<< hist->GetXaxis()->GetXmax() << " in " << hist->GetName() << std::endl;
		}
	}
}

void WCSimFitterPlots::PrintRecoMinusTrue()
{
	std::cout << " -- WCSimFitterPlots: Reco - true variables to be plotted -- " << std::endl;
	std::map<FitterParameterType::Type, std::vector<TH1D*> >::const_iterator plotItr = fRecoMinusTrue.begin();
	for( ; plotItr != fRecoMinusTrue.end(); ++plotItr)
	{
		FitterParameterType::Type type = ((*plotItr).first);
		for(unsigned int iHist = 0; iHist < ((*plotItr).second).size(); ++ iHist)
		{
			TH1D * hist = ((*plotItr).second).at(iHist);
			std::cout << "\t Plotting " << FitterParameterType::AsString(type) << " in "
					<< hist->GetNbinsX() << " bins from " << hist->GetXaxis()->GetXmin() << " to "
					<< hist->GetXaxis()->GetXmax() << " in " << hist->GetName() << std::endl;
		}
	}
}

void WCSimFitterPlots::PrintSurfaces() {
	Print1DSurfaces();
	Print2DSurfaces();
}

void WCSimFitterPlots::Print1DSurfaces()
{
	std::cout << " -- WCSimFitterPlots: 1D likelihood profiles to be plotted -- " << std::endl;
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D* >::iterator plotItr = fSurfaces1D.begin();
	for( ; plotItr != fSurfaces1D.end(); ++plotItr)
	{
		std::pair<unsigned int, FitterParameterType::Type> type = ((*plotItr).first);
		TH1D * hist = ((*plotItr).second);
		std::cout << "\t Plotting " << FitterParameterType::AsString(type.second) << " for track " << type.first;

    if(hist != NULL){
				  std::cout << "in " << hist->GetNbinsX() << " bins from " << hist->GetXaxis()->GetXmin() << " to "
				  << hist->GetXaxis()->GetXmax() << std::endl;
    }
    else{ std::cout << std::endl; }
	}

}

void WCSimFitterPlots::CreateNtuple(WCSimFitterConfig* fitterConfig) {
}

void WCSimFitterPlots::FillNtuple(WCSimFitterConfig* fitterConfig,
		std::vector<WCSimLikelihoodTrack> bestFits) {
}

void WCSimFitterPlots::Fill1DProfile(std::pair<unsigned int, FitterParameterType::Type> theProfile,
		double stepVal, double minus2LnL) {
  std::cout << std::endl << "Fill 1D profile" << std::endl;
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D*>::iterator mapItr = fSurfaces1D.find(theProfile);
	if( mapItr != fSurfaces1D.end() )
	{
	
  	TH1D * profile = (*mapItr).second;
    std::cout << "Profile = " << profile << std::endl; 
    std::cout << "Profile title = " << profile->GetTitle() << std::endl;
		profile->SetBinContent(profile->GetXaxis()->FindBin(stepVal), minus2LnL);
	}
}

void WCSimFitterPlots::Fill2DProfile(
		std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > theProfile, double stepValX,
		double stepValY, double minus2LnL) {
	std::map<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> >, TH2D* >::iterator mapItr;
	mapItr = fSurfaces2D.find(theProfile);
	if( mapItr != fSurfaces2D.end() )
	{
		TH2D * profile = (*mapItr).second;
		profile->SetBinContent(profile->GetXaxis()->FindBin(stepValX), profile->GetYaxis()->FindBin(stepValY), minus2LnL);
	}

}

std::vector<std::pair<unsigned int, FitterParameterType::Type> > WCSimFitterPlots::GetAll1DSurfaceKeys() const {
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D*>::const_iterator itr1D = fSurfaces1D.begin();
	std::vector<std::pair<unsigned int, FitterParameterType::Type> > myVec;
	for( ; itr1D != fSurfaces1D.end() ; ++itr1D)
	{
		myVec.push_back((*itr1D).first);
	}
	return myVec;
}

std::vector<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > > WCSimFitterPlots::GetAll2DSurfaceKeys() const {
	std::map<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> >, TH2D*>::const_iterator itr2D = fSurfaces2D.begin();
	std::vector<std::pair <std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > > myVec;
	for( ; itr2D != fSurfaces2D.end() ; ++itr2D)
	{
		myVec.push_back((*itr2D).first);
	}
	return myVec;
}

void WCSimFitterPlots::SetSaveFileName(const char* filename) {
	fSaveFileName = TString(filename);
}

void WCSimFitterPlots::SetInputFileNamesToWrite( TObjArray nameArr )
{
  fInputFileNames = nameArr;
}

TObjArray WCSimFitterPlots::GetInputFileNamesToWrite() const
{
  return fInputFileNames;
}


void WCSimFitterPlots::SavePlots() {
	TDirectory* tmpd = 0;

	if( fSaveFile == NULL )
	{
		tmpd = gDirectory;
	    std::cout << " *** WCSimFitter::SavePlots() *** " << std::endl;
	    std::cout << "  opening file: " << fSaveFileName.Data() << std::endl;
	    fSaveFile = new TFile(fSaveFileName.Data(), "RECREATE");
	}

	fSaveFile->cd();
	fSaveFile->mkdir("Profiles");
	fSaveFile->cd("Profiles");

	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D*>::iterator surf1DItr = fSurfaces1D.begin();
	for( ; surf1DItr != fSurfaces1D.end() ; ++surf1DItr)
	{
    TH1D * hist = dynamic_cast<TH1D*>( ((*surf1DItr).second)->Clone() );
		hist->Write();
	}
	std::map<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> >, TH2D*>::iterator surf2DItr = fSurfaces2D.begin();
	for( ; surf2DItr != fSurfaces2D.end() ; ++surf2DItr)
	{
    TH2D * hist = dynamic_cast<TH2D*>( ((*surf2DItr).second)->Clone() );
    hist->Write();
	}

	fSaveFile->cd();
	fSaveFile->mkdir("FitResults");
	fSaveFile->cd("FitResults");
	std::map<FitterParameterType::Type, std::vector<TH1D*> >::iterator eachEventItr = fForEachEvent.begin();
	for( ; eachEventItr != fForEachEvent.end() ; ++eachEventItr)
	{
		std::vector<TH1D*>::iterator plotItr = (*eachEventItr).second.begin();
		for( ; plotItr != (*eachEventItr).second.end() ; ++plotItr)
		{
			TH1D * hist = dynamic_cast<TH1D*>( (*plotItr)->Clone() );
      hist->Write();
		}
	}


	fSaveFile->cd();
	fSaveFile->mkdir("RecoMinusTrue");
	fSaveFile->cd("RecoMinusTrue");
	std::map<FitterParameterType::Type, std::vector<TH1D*> >::iterator rmtItr = fRecoMinusTrue.begin();
	for( ; rmtItr != fRecoMinusTrue.end() ; ++rmtItr)
	{
		std::vector<TH1D*>::iterator plotItr = (*rmtItr).second.begin();
		for( ; plotItr != ((*rmtItr).second).end() ; ++plotItr)
		{
      TH1D * hist = dynamic_cast<TH1D*>((*plotItr)->Clone());
      hist->Write();
		}
	}

  fSaveFile->cd();
  fInputFileNames.Write();
  
  fSaveFile->Close();
	tmpd->cd();
	return;
}

void WCSimFitterPlots::Print2DSurfaces()
{
	std::cout << " -- WCSimFitterPlots: 2D likelihood profiles to be plotted -- " << std::endl;
	std::map<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> >, TH2D*>::iterator plotItr;
	for( plotItr = fSurfaces2D.begin(); plotItr != fSurfaces2D.end(); ++plotItr)
	{
		std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > types = ((*plotItr).first);
		FitterParameterType::Type type1 = (types.first).second;
		FitterParameterType::Type type2 = (types.second).second;
    TH2D * hist = ((*plotItr).second);
    
		std::cout << "\t Plotting " << FitterParameterType::AsString(type2);

    if(hist != NULL)
    {
      std::cout << " in " << hist->GetNbinsY() << " bins from " << hist->GetYaxis()->GetXmin() << " to ";
	    std::cout << hist->GetYaxis()->GetXmax() << std::endl;
    }
    else{ std::cout << std::endl; }
    std::cout << "\t\t Against" << FitterParameterType::AsString(type1);
    if(hist != NULL)
    {
       std::cout <<  " in "
  				       << hist->GetNbinsX() << " bins from " << hist->GetXaxis()->GetXmin() << " to "
  				       << hist->GetXaxis()->GetXmax() << std::endl;
    }
    else{ std::cout << std::endl; }
	}

  std::cout << " -------------------------------------------------- " << std::endl;
}

void WCSimFitterPlots::FillPlots(std::vector<WCSimLikelihoodTrack> bestFits) {
	std::map<FitterParameterType::Type, std::vector<TH1D*> >::iterator plotVecItr = fForEachEvent.begin();
	for(unsigned int iTrack = 0; iTrack < bestFits.size(); ++iTrack)
	{
    std::cout << "Filling fit result plots for track" << iTrack << std::endl;
    std::cout << "Size of best fits vector = " << bestFits.size() << std::endl;
		for( ; plotVecItr != fForEachEvent.end(); ++plotVecItr)
		{
			FitterParameterType::Type type = (*plotVecItr).first;
			TH1D * hist = (*plotVecItr).second.at(iTrack);

			double bestFitParam = bestFits.at(iTrack).GetTrackParameter(type);
			hist->Fill(bestFitParam);
		}
	}
	return;
}

void WCSimFitterPlots::FillRecoMinusTrue(
		std::vector<WCSimLikelihoodTrack> bestFits, std::vector<WCSimLikelihoodTrack*> * trueTracks) {
	std::vector<WCSimLikelihoodTrack>::iterator bfItr = bestFits.begin();
	std::vector<WCSimLikelihoodTrack*>::iterator truItr = (*trueTracks).begin();

	unsigned int iTrack = 0;
  for( ; bfItr != bestFits.end() && truItr != (*trueTracks).end() ; ++bfItr, ++truItr)
	{
		WCSimLikelihoodTrack bf = (*bfItr);
		WCSimLikelihoodTrack tt = *(*truItr);

		std::map<FitterParameterType::Type, std::vector<TH1D*> >::iterator plotVecItr = fRecoMinusTrue.begin();
		for( ; plotVecItr != fRecoMinusTrue.end() ; ++plotVecItr )
		{
			FitterParameterType::Type type = (*plotVecItr).first;
			TH1D * hist = (*plotVecItr).second.at(iTrack);

			double rmt = bf.GetTrackParameter(type) - tt.GetTrackParameter(type);
			hist->Fill(rmt);
		}
      iTrack++;
	}
	return;
}



