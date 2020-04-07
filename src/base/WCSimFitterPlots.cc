/*
 * WCSimFitterPlots.cc
 *
 *  Created on: 3 Nov 2014
 *      Author: andy
 */
#include "TCanvas.h"

#include "WCSimFitterConfig.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimParameters.hh"
#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TString.h>
#include <TSystem.h>
#include <TTimeStamp.h>

WCSimFitterPlots::WCSimFitterPlots() : fSaveFile(NULL), fNumSurfaceBins(20)
{
	// TODO Auto-generated constructor stub
}

WCSimFitterPlots::~WCSimFitterPlots()
{
	// TODO Auto-generated destructor stub
}

void WCSimFitterPlots::SetPlotForEachEvent(const char *name, bool doIt)
{
	if (doIt)
	{
		std::cout << " *** SetPlotForEachEvent *** ... plotting " << name << std::endl;
		FitterParameterType::Type type = FitterParameterType::FromName(name);
		//  std::cout << " *** Set it *** " << std::endl;
		std::vector<TH1D *> plotVec;
		fForEachEvent.insert(std::make_pair(type, plotVec));
	}
}

bool WCSimFitterPlots::GetPlotForEachEvent(const char *name) const
{
	return (fForEachEvent.find(FitterParameterType::FromName(name)) != fForEachEvent.end());
}

void WCSimFitterPlots::SetPlotRecoMinusTrue(const char *name, bool doIt)
{
	std::cout << " *** SetPlotRecoMinusTrue ***  ... plotting " << name << std::endl;
	FitterParameterType::Type type = FitterParameterType::FromName(name);
	std::vector<TH1D *> plotVec;
	fRecoMinusTrue.insert(std::make_pair(type, plotVec));
}

bool WCSimFitterPlots::GetPlotRecoMinusTrue(const char *name) const
{
	return (fRecoMinusTrue.find(FitterParameterType::FromName(name)) != fRecoMinusTrue.end());
}

void WCSimFitterPlots::Make1DSurface(const char *name, bool doIt, unsigned int trackNum)
{
	std::cout << " *** Make 1D surface ***  ... plotting " << name << std::endl;
	FitterParameterType type = FitterParameterType::FromName(name);
	//  std::cout << "Made type" << std::endl;
	std::pair<unsigned int, FitterParameterType::Type> toProfile = std::make_pair(trackNum, type);
	//  std::cout << "Made pair" << std::endl;
	TH1D *plotVec = 0;
	fSurfaces1D.insert(std::make_pair(toProfile, plotVec));
	//  std::cout << "Inserted to map" << std::endl;
}

bool WCSimFitterPlots::GetMake1DSurface(const char *name, unsigned int trackNum) const
{
	return (fSurfaces1D.find(std::make_pair(trackNum, FitterParameterType::FromName(name))) != fSurfaces1D.end());
}

void WCSimFitterPlots::Make2DSurface(const char *name, const char *name2, bool doIt, unsigned int trackNum,
									 unsigned int trackNum2)
{

	FitterParameterType::Type type = FitterParameterType::FromName(name);
	FitterParameterType::Type type2 = FitterParameterType::FromName(name2);

	std::pair<unsigned int, FitterParameterType::Type> toProfile1 = std::make_pair(trackNum, type);
	std::pair<unsigned int, FitterParameterType::Type> toProfile2 = std::make_pair(trackNum2, type2);

	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> types;
	if (trackNum < trackNum2)
	{
		types = std::make_pair(toProfile1, toProfile2);
	}
	else if (trackNum2 < trackNum)
	{
		types = std::make_pair(toProfile2, toProfile1);
	}
	else if (type < type2)
	{
		types = std::make_pair(toProfile1, toProfile2);
	}
	else
	{
		types = std::make_pair(toProfile2, toProfile1);
	}

	TH2D *plotHist = 0;
	fSurfaces2D.insert(std::make_pair(types, plotHist));
}

bool WCSimFitterPlots::GetMake2DSurface(const char *name, const char *name2, unsigned int trackNum,
										unsigned int trackNum2) const
{
	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> types;
	std::pair<unsigned int, FitterParameterType::Type> type1, type2;
	type1 = std::make_pair(trackNum, FitterParameterType::FromName(name));
	type2 = std::make_pair(trackNum2, FitterParameterType::FromName(name2));
	types = std::make_pair(type1, type2);
	return (fSurfaces2D.find(types) != fSurfaces2D.end());
}

void WCSimFitterPlots::SetNumSurfaceBins(unsigned int numBins)
{
	fNumSurfaceBins = numBins;
}

unsigned int WCSimFitterPlots::GetNumSurfaceBins()
{
	return fNumSurfaceBins;
}

void WCSimFitterPlots::MakeHistograms(WCSimFitterConfig *fitterConfig)
{
	std::cout << " *** WCSimFitterPlots::MakeHistograms() *** " << std::endl;
	MakePlotsForEachEvent(fitterConfig);
	MakeRecoMinusTrue(fitterConfig);
	Make1DSurfaces(fitterConfig);
	Make2DSurface(fitterConfig);
	return;
}

void WCSimFitterPlots::MakePlotsForEachEvent(WCSimFitterConfig *fitterConfig)
{
	std::cout << " *** WCSimFitterPlots::MakePlotsForEachEvent() *** " << std::endl;
	std::map<FitterParameterType::Type, std::vector<TH1D *>>::iterator plotItr = fForEachEvent.begin();
	for (; plotItr != fForEachEvent.end(); ++plotItr)
	{
		(*plotItr).second.clear();
		for (unsigned int iTrack = 0; iTrack < fitterConfig->GetNumTracks(); ++iTrack)
		{
			TH1D *hist = NULL;
			FitterParameterType::Type type = ((*plotItr).first);
			std::string parName = FitterParameterType::AsString(type).c_str();
			if (hist != NULL)
			{
				delete hist;
			}

			TString name, title;
			name.Form("h_%s_track%d", parName.c_str(), iTrack);
			title.Form("Value of %s for track %d, for each of %d fitted tracks", parName.c_str(), iTrack,
					   fitterConfig->GetNumTracks());

			int numBins = 200;
			if (fitterConfig->GetNumEventsToFit() > 400)
			{
				std::cout << "You shouldn't see this" << std::endl;
				numBins = fitterConfig->GetNumEventsToFit() / 20;
				numBins += (10 - (numBins % 10));
			}
			hist = new TH1D(name.Data(), title.Data(), numBins, fitterConfig->GetParMin(iTrack, parName.c_str()),
							fitterConfig->GetParMax(iTrack, parName.c_str()));
			std::cout << hist->GetName() << " hist has " << numBins << " from " << hist->GetXaxis()->GetXmin() << " to "
					  << hist->GetXaxis()->GetXmax() << std::endl;

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

void WCSimFitterPlots::MakeRecoMinusTrue(WCSimFitterConfig *fitterConfig)
{
	std::cout << " *** WCSimFitterPlots::MakeRecoMinusTrue() *** " << std::endl;
	std::map<FitterParameterType::Type, std::vector<TH1D *>>::iterator plotItr = fRecoMinusTrue.begin();
	for (; plotItr != fRecoMinusTrue.end(); ++plotItr)
	{
		(*plotItr).second.clear();
		for (unsigned int iTrack = 0; iTrack < fitterConfig->GetNumTracks(); ++iTrack)
		{
			TH1D *hist = NULL;
			;
			FitterParameterType::Type type = ((*plotItr).first);
			std::string parName = FitterParameterType::AsString(type).c_str();
			if (hist != NULL)
			{
				delete hist;
			}

			TString name, title;
			name.Form("h_%s_track%d_reco_minus_true", parName.c_str(), iTrack);
			title.Form("Reco - true value of %s for track %d, for each of %d fitted events", parName.c_str(), iTrack,
					   fitterConfig->GetNumTracks());

			int numBins = 100;
			double paramMax = WCSimParameters::Instance()->GetPlotMax(parName);

			/*
			if (fitterConfig->GetNumEventsToFit() > 400) {
				numBins = fitterConfig->GetNumEventsToFit() / 20;
				numBins += (10 - (numBins % 10));
			}

			double lowBin = fitterConfig->GetParMax(iTrack, parName.c_str())
					- 2
							* (fitterConfig->GetParMax(iTrack, parName.c_str())
									- fitterConfig->GetParMin(iTrack, parName.c_str()));

			hist = new TH1D(name.Data(), title.Data(), numBins, lowBin,
					2.0 * fitterConfig->GetParMax(iTrack, parName.c_str()));
			*/
			hist = new TH1D(name.Data(), title.Data(), numBins, -paramMax, paramMax);
			std::cout << hist->GetName() << " hist has " << numBins << " from " << hist->GetXaxis()->GetXmin() << " to "
					  << hist->GetXaxis()->GetXmax() << std::endl;

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

void WCSimFitterPlots::Make1DSurfaces(WCSimFitterConfig *fitterConfig)
{
	std::cout << " *** WCSimFitterPlots::Make1DSurfaces() *** " << std::endl;
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D *>::iterator plotItr = fSurfaces1D.begin();
	for (; plotItr != fSurfaces1D.end(); ++plotItr)
	{
		TH1D *hist = NULL;
		FitterParameterType::Type type = ((*plotItr).first).second;
		std::cout << FitterParameterType::AsString(type) << std::endl;
		std::string parName = FitterParameterType::AsString(type).c_str();
		unsigned int trackNum = ((*plotItr).first).first;

		TString name, title;
		name.Form("h_%s_likelihood_profile", parName.c_str());
		title.Form("Profile of likelihood when %s is varied for track %d, with other parameters fixed", parName.c_str(),
				   trackNum);

		hist = new TH1D(name.Data(), title.Data(), fNumSurfaceBins, fitterConfig->GetParMin(trackNum, parName.c_str()),
						fitterConfig->GetParMax(trackNum, parName.c_str()));
		std::cout << hist->GetName() << " hist has " << fNumSurfaceBins << " from " << hist->GetXaxis()->GetXmin()
				  << " to " << hist->GetXaxis()->GetXmax() << std::endl;

		TString xTitle(parName.c_str());
		hist->GetXaxis()->SetTitle(xTitle.Data());
		hist->GetXaxis()->CenterTitle();
		hist->GetXaxis()->SetLabelSize(0.05);
		hist->GetXaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitle("-2 Ln(L)");
		hist->GetYaxis()->SetLabelSize(0.05);
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->CenterTitle();
		hist->SetLineWidth(2);
		hist->SetLineColor(kMagenta + 2);

		if ((*plotItr).second != NULL)
		{
			delete (*plotItr).second;
		}
		(*plotItr).second = hist;
	}
}

void WCSimFitterPlots::Make2DSurface(WCSimFitterConfig *fitterConfig)
{
	std::cout << " *** WCSimFitterPlots::Make2DSurface() *** " << std::endl;
	std::map<
		std::pair<std::pair<unsigned int, FitterParameterType::Type>,
				  std::pair<unsigned int, FitterParameterType::Type>>,
		TH2D *>::iterator plotItr =
		fSurfaces2D.begin();
	for (; plotItr != fSurfaces2D.end(); ++plotItr)
	{
		TH2D *hist = NULL;

		std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> types =
			((*plotItr).first);

		std::string parName1 = FitterParameterType::AsString(types.first.second).c_str();
		std::string parName2 = FitterParameterType::AsString(types.second.second).c_str();
		unsigned int trackNum1 = (types.first).first;
		unsigned int trackNum2 = (types.second).first;

		TString name, title;
		name.Form("h_%s_vs_%s_likelihood_profile", parName2.c_str(), parName1.c_str());
		title.Form("Profile of likelihood when %s & %s are varied for tracks %d and %d, with other parameters fixed",
				   parName2.c_str(), parName1.c_str(), trackNum2, trackNum1);

		hist = new TH2D(name.Data(), title.Data(), fNumSurfaceBins,
						fitterConfig->GetParMin(trackNum1, parName1.c_str()),
						fitterConfig->GetParMax(trackNum1, parName1.c_str()), fNumSurfaceBins,
						fitterConfig->GetParMin(trackNum2, parName2.c_str()),
						fitterConfig->GetParMax(trackNum2, parName2.c_str()));
		std::cout << hist->GetName() << " hist has " << fNumSurfaceBins << " x-from " << hist->GetXaxis()->GetXmin()
				  << " to " << hist->GetXaxis()->GetXmax() << " y-from " << hist->GetYaxis()->GetXmin() << " to "
				  << hist->GetYaxis()->GetXmax() << std::endl;

		hist->GetXaxis()->SetTitle(parName1.c_str());
		hist->GetXaxis()->CenterTitle();
		hist->GetXaxis()->SetLabelSize(0.05);
		hist->GetXaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitle(parName2.c_str());
		hist->GetYaxis()->SetLabelSize(0.05);
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->CenterTitle();
		hist->GetZaxis()->SetTitle("-2 Ln(L)");
		hist->SetOption("COLZ");

		if ((*plotItr).second != NULL)
		{
			delete (*plotItr).second;
		}
		(*plotItr).second = hist;
	}
}

void WCSimFitterPlots::Print()
{
	std::cout << " -------------------------------------------------- " << std::endl;
	PrintVariables();
	PrintRecoMinusTrue();
	Print1DSurfaces();
	Print2DSurfaces();
	std::cout << " -------------------------------------------------- " << std::endl;
}

void WCSimFitterPlots::PrintVariables()
{
	std::cout << " -- WCSimFitterPlots: Variables to be plotted -- " << std::endl;
	std::map<FitterParameterType::Type, std::vector<TH1D *>>::const_iterator plotItr = fForEachEvent.begin();
	for (; plotItr != fForEachEvent.end(); ++plotItr)
	{
		FitterParameterType::Type type = ((*plotItr).first);
		std::cout << "\t Plotting " << FitterParameterType::AsString(type) << std::endl;
	}
}

void WCSimFitterPlots::PrintRecoMinusTrue()
{
	std::cout << " -- WCSimFitterPlots: Reco - true variables to be plotted -- " << std::endl;
	std::map<FitterParameterType::Type, std::vector<TH1D *>>::const_iterator plotItr = fRecoMinusTrue.begin();
	for (; plotItr != fRecoMinusTrue.end(); ++plotItr)
	{
		FitterParameterType::Type type = ((*plotItr).first);
		std::cout << "\t Plotting " << FitterParameterType::AsString(type) << std::endl;
	}
}

void WCSimFitterPlots::PrintSurfaces()
{
	Print1DSurfaces();
	Print2DSurfaces();
}

void WCSimFitterPlots::Print1DSurfaces()
{
	std::cout << " -- WCSimFitterPlots: 1D likelihood profiles to be plotted -- " << std::endl;
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D *>::iterator plotItr = fSurfaces1D.begin();
	for (; plotItr != fSurfaces1D.end(); ++plotItr)
	{
		std::pair<unsigned int, FitterParameterType::Type> type = ((*plotItr).first);
		TH1D *hist = ((*plotItr).second);
		std::cout << "\t Plotting " << FitterParameterType::AsString(type.second) << " for track " << type.first;

		if (hist != NULL)
		{
			std::cout << "in " << hist->GetNbinsX() << " bins from " << hist->GetXaxis()->GetXmin() << " to "
					  << hist->GetXaxis()->GetXmax() << std::endl;
		}
		else
		{
			std::cout << std::endl;
		}
	}
}

void WCSimFitterPlots::Print2DSurfaces()
{
	std::cout << " -- WCSimFitterPlots: 2D likelihood profiles to be plotted -- " << std::endl;
	std::map<
		std::pair<std::pair<unsigned int, FitterParameterType::Type>,
				  std::pair<unsigned int, FitterParameterType::Type>>,
		TH2D *>::iterator plotItr;
	for (plotItr = fSurfaces2D.begin(); plotItr != fSurfaces2D.end(); ++plotItr)
	{
		std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> types =
			((*plotItr).first);
		FitterParameterType::Type type1 = (types.first).second;
		FitterParameterType::Type type2 = (types.second).second;
		TH2D *hist = ((*plotItr).second);

		std::cout << "\t Plotting " << FitterParameterType::AsString(type2);

		if (hist != NULL)
		{
			std::cout << " in " << hist->GetNbinsY() << " bins from " << hist->GetYaxis()->GetXmin() << " to ";
			std::cout << hist->GetYaxis()->GetXmax() << std::endl;
		}
		else
		{
			std::cout << std::endl;
		}
		std::cout << "\t\t Against" << FitterParameterType::AsString(type1);
		if (hist != NULL)
		{
			std::cout << " in " << hist->GetNbinsX() << " bins from " << hist->GetXaxis()->GetXmin() << " to "
					  << hist->GetXaxis()->GetXmax() << std::endl;
		}
		else
		{
			std::cout << std::endl;
		}
	}
}

void WCSimFitterPlots::FillPlots(std::vector<WCSimLikelihoodTrackBase *> bestFits)
{
	std::cout << " *** WCSimFitterPlots::FillPlots() *** " << std::endl;
	std::map<FitterParameterType::Type, std::vector<TH1D *>>::iterator plotVecItr = fForEachEvent.begin();
	for (unsigned int iTrack = 0; iTrack < bestFits.size(); ++iTrack)
	{
		for (; plotVecItr != fForEachEvent.end(); ++plotVecItr)
		{
			FitterParameterType::Type type = (*plotVecItr).first;
			TH1D *hist = (*plotVecItr).second.at(iTrack);

			double bestFitParam = bestFits.at(iTrack)->GetTrackParameter(type);
			hist->Fill(bestFitParam);
		}
	}
}

void WCSimFitterPlots::FillRecoMinusTrue(std::vector<WCSimLikelihoodTrackBase *> bestFits,
										 std::vector<WCSimLikelihoodTrackBase *> *trueTracks)
{
	std::cout << " *** WCSimFitterPlots::FillRecoMinusTrue() *** " << std::endl;
	// First sort all the tracks by energy
	std::sort(bestFits.begin(), bestFits.end(), WCSimLikelihoodTrackBase::EnergyGreaterThanOrEqualPtrs);
	// Make a copy so we don't mess with the pointed vector
	std::vector<WCSimLikelihoodTrackBase *> trueTracksSorted = *trueTracks;
	std::sort(trueTracksSorted.begin(), trueTracksSorted.end(), WCSimLikelihoodTrackBase::EnergyGreaterThanOrEqualPtrs);

	unsigned int maxTracksToPlot = trueTracksSorted.size();
	if (trueTracksSorted.size() > bestFits.size())
	{
		std::cerr << "Warning: more true tracks than best-fit tracks.  Picking the highest KE one" << std::endl;
		maxTracksToPlot = bestFits.size();
	}
	else if (trueTracksSorted.size() < bestFits.size())
	{
		std::cerr << "Warning: more fitted tracks than true tracks.  Picking the highest KE one" << std::endl;
	}

	for (unsigned int iTrack = 0; iTrack < maxTracksToPlot; ++iTrack)
	{
		WCSimLikelihoodTrackBase *bft = (bestFits.at(iTrack));
		WCSimLikelihoodTrackBase *trt = (trueTracksSorted.at(iTrack));

		std::map<FitterParameterType::Type, std::vector<TH1D *>>::iterator plotVecItr = fRecoMinusTrue.begin();
		for (; plotVecItr != fRecoMinusTrue.end(); ++plotVecItr)
		{
			FitterParameterType::Type type = (*plotVecItr).first;
			TH1D *hist = (*plotVecItr).second.at(iTrack);

			double rmt = bft->GetTrackParameter(type) - trt->GetTrackParameter(type);
			hist->Fill(rmt);
		}
	}
}

std::vector<std::pair<unsigned int, FitterParameterType::Type>> WCSimFitterPlots::GetAll1DSurfaceKeys() const
{
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D *>::const_iterator itr1D = fSurfaces1D.begin();
	std::vector<std::pair<unsigned int, FitterParameterType::Type>> myVec;
	for (; itr1D != fSurfaces1D.end(); ++itr1D)
	{
		myVec.push_back((*itr1D).first);
	}
	return myVec;
}

double WCSimFitterPlots::Get1DSurfaceBinCenter(std::pair<unsigned int, FitterParameterType::Type> theProfile,
											   int iBin)
{
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D *>::iterator mapItr = fSurfaces1D.find(theProfile);
	if (mapItr != fSurfaces1D.end())
	{
		TH1D *profile = (*mapItr).second;
		return profile->GetXaxis()->GetBinCenter(iBin);
	}
	else
	{
		assert(0);
	}
}

void WCSimFitterPlots::Fill1DProfile(std::pair<unsigned int, FitterParameterType::Type> theProfile, int binNum,
									 double minus2LnL)
{
	std::cout << " *** WCSimFitterPlots::Fill1DProfile() *** " << std::endl;
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D *>::iterator mapItr = fSurfaces1D.find(theProfile);
	if (mapItr != fSurfaces1D.end())
	{
		TH1D *profile = (*mapItr).second;
		//std::cout << "Profile title = " << profile->GetTitle() << std::endl;
		//std::cout << "binNum -> " << binNum << ", minus2lnl -> " << minus2LnL << std::endl;
		profile->SetBinContent(binNum, minus2LnL);
		TCanvas *can1 = new TCanvas("can1", "can1", 800, 600);
		int bin = profile->GetMinimumBin();
		TMarker *marker = new TMarker(profile->GetXaxis()->GetBinCenter(bin), profile->GetBinContent(bin), 29);
		marker->SetMarkerSize(1.8);
		marker->SetMarkerColor(kGreen - 2);

		profile->Draw();
		marker->Draw();
		can1->SaveAs("can1.png");

		profile = 0x0;
		delete marker;
		delete can1;
	}
}

std::vector<
	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>>>
WCSimFitterPlots::GetAll2DSurfaceKeys() const
{
	std::map<
		std::pair<std::pair<unsigned int, FitterParameterType::Type>,
				  std::pair<unsigned int, FitterParameterType::Type>>,
		TH2D *>::const_iterator itr2D =
		fSurfaces2D.begin();
	std::vector<std::pair<std::pair<unsigned int, FitterParameterType::Type>,
						  std::pair<unsigned int, FitterParameterType::Type>>>
		myVec;
	for (; itr2D != fSurfaces2D.end(); ++itr2D)
	{
		myVec.push_back((*itr2D).first);
	}
	return myVec;
}

void WCSimFitterPlots::Get2DSurfaceBinCenters(
	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> theProfile,
	int binNumX, int binNumY, double &x, double &y)
{
	std::map<
		std::pair<std::pair<unsigned int, FitterParameterType::Type>,
				  std::pair<unsigned int, FitterParameterType::Type>>,
		TH2D *>::iterator mapItr;
	mapItr = fSurfaces2D.find(theProfile);
	if (mapItr != fSurfaces2D.end())
	{
		TH2D *profile = (*mapItr).second;
		x = profile->GetXaxis()->GetBinCenter(binNumX);
		y = profile->GetYaxis()->GetBinCenter(binNumY);
	}
	else
	{
		std::cerr << "Error: 2D profile not found!" << std::endl;
		assert(mapItr != fSurfaces2D.end());
	}
}

double WCSimFitterPlots::Get2DSurfaceBinCenterX(
	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> theProfile,
	int binNumX)
{
	double x = 0;
	double y = 0;
	int binNumY = 1;
	Get2DSurfaceBinCenters(theProfile, binNumX, binNumY, x, y);
	return x;
}

double WCSimFitterPlots::Get2DSurfaceBinCenterY(
	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> theProfile,
	int binNumY)
{
	double x = 0.0;
	double y = 0.0;
	int binNumX = 1;
	Get2DSurfaceBinCenters(theProfile, binNumX, binNumY, x, y);
	return y;
}

void WCSimFitterPlots::Fill2DProfile(
	std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type>> theProfile,
	int binNumX, int binNumY, double minus2LnL)
{
	std::cout << " *** WCSimFitterPlots::Fill2DProfile() *** " << std::endl;
	std::map<
		std::pair<std::pair<unsigned int, FitterParameterType::Type>,
				  std::pair<unsigned int, FitterParameterType::Type>>,
		TH2D *>::iterator mapItr;
	mapItr = fSurfaces2D.find(theProfile);
	if (mapItr != fSurfaces2D.end())
	{
		TH2D *profile = (*mapItr).second;
		profile->SetBinContent(binNumX, binNumY, minus2LnL);
	}
}

void WCSimFitterPlots::SaveProfiles()
{
	std::cout << " *** WCSimFitterPlots::SaveProfiles() *** " << std::endl;
	TDirectory *tmpd = 0;
	tmpd = gDirectory;
	if (fSaveFile == NULL)
	{
		fSaveFile = new TFile(fSaveFileName.Data(), "UPDATE");
	}
	std::cout << "Saving plots to " << fSaveFileName << std::endl;
	fSaveFile->cd();

	if (fSaveFile->GetDirectory("Profiles") == 0)
	{
		fSaveFile->mkdir("Profiles");
	}
	fSaveFile->cd("Profiles");

	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D *>::iterator surf1DItr = fSurfaces1D.begin();
	for (; surf1DItr != fSurfaces1D.end(); ++surf1DItr)
	{
		TH1D *hist = (TH1D *)((*surf1DItr).second);
		std::cout << "Saving plot " << hist->GetName() << " has " << hist->GetNbinsX() << " bins" << std::endl;
		hist->Write();
	}
	std::map<
		std::pair<std::pair<unsigned int, FitterParameterType::Type>,
				  std::pair<unsigned int, FitterParameterType::Type>>,
		TH2D *>::iterator surf2DItr =
		fSurfaces2D.begin();
	for (; surf2DItr != fSurfaces2D.end(); ++surf2DItr)
	{
		TH2D *hist = (TH2D *)((*surf2DItr).second);
		std::cout << "Saving plot " << hist->GetName() << " has " << hist->GetNbinsX() * hist->GetNbinsY() << " bins"
				  << std::endl;
		hist->Write();
	}

	fSaveFile->cd();

	tmpd->cd();
}

void WCSimFitterPlots::SaveLikelihoodPlots(TH1D *predictedChargeHist, TH1D *measuredHist, TH1D *pred_measured_Hist, TH1D *best2LnLsHist,
										   TH1D *charge2LnLsHist, TH1D *time2LnLsHist, TH1D *predictedTimesHist)
{
	std::cout << " *** WCSimFitterPlots::SaveLikelihoodPlots() *** " << std::endl;
	TDirectory *tmpd = 0;
	tmpd = gDirectory;
	if (fSaveFile == NULL)
	{
		fSaveFile = new TFile(fSaveFileName.Data(), "UPDATE");
	}
	std::cout << "Saving plots to " << fSaveFileName << std::endl;
	fSaveFile->cd();

	if (fSaveFile->GetDirectory("Likelihood") == 0)
	{
		fSaveFile->mkdir("Likelihood");
	}
	fSaveFile->cd("Likelihood");

	predictedChargeHist->Write();
	measuredHist->Write();
	pred_measured_Hist->Write();
	best2LnLsHist->Write();
	charge2LnLsHist->Write();
	time2LnLsHist->Write();
	predictedTimesHist->Write();

	fSaveFile->cd();
	tmpd->cd();
}

void WCSimFitterPlots::SavePlots()
{
	std::cout << " *** WCSimFitterPlots::SavePlots() *** " << std::endl;
	TDirectory *tmpd = 0;
	tmpd = gDirectory;
	if (fSaveFile == NULL)
	{
		std::cout << "File " << fSaveFileName.Data() << " exists... updating it" << std::endl;
		fSaveFile = new TFile(fSaveFileName.Data(), "UPDATE");
	}
	std::cout << "Saving plots to " << fSaveFileName << std::endl;
	fSaveFile->cd();

	TNamed inputFile("inputFile", fInputFileName.Data());
	inputFile.Write();

	fSaveFile->cd();
	if (fSaveFile->GetDirectory("FitResults") == 0)
	{
		fSaveFile->mkdir("FitResults");
	}
	fSaveFile->cd("FitResults");

	std::map<FitterParameterType::Type, std::vector<TH1D *>>::iterator forEachVecItr = fForEachEvent.begin();
	for (unsigned int iTrack = 0; iTrack < (*forEachVecItr).second.size(); ++iTrack)
	{
		for (; forEachVecItr != fForEachEvent.end(); ++forEachVecItr)
		{
			TH1D *hist = (*forEachVecItr).second.at(iTrack);
			std::cout << "Saving Plot " << hist->GetName() << std::endl;
			hist->Write("", TObject::kOverwrite);
		}
	}

	fSaveFile->cd();
	if (fSaveFile->GetDirectory("RecoMinusTrue") == 0)
	{
		fSaveFile->mkdir("RecoMinusTrue");
	}
	fSaveFile->cd("RecoMinusTrue");

	std::map<FitterParameterType::Type, std::vector<TH1D *>>::iterator recoMinusVecItr = fRecoMinusTrue.begin();
	for (unsigned int iTrack = 0; iTrack < (*recoMinusVecItr).second.size(); ++iTrack)
	{
		for (; recoMinusVecItr != fRecoMinusTrue.end(); ++recoMinusVecItr)
		{
			TH1D *hist = (*recoMinusVecItr).second.at(iTrack);
			std::cout << "Saving Plot " << hist->GetName() << std::endl;
			hist->Write("", TObject::kOverwrite);
		}
	}

	fSaveFile->cd();
	fSaveFile->Close();
	fSaveFile = NULL;
	tmpd->cd();
	return;
}

void WCSimFitterPlots::SetSaveFileName(const char *filename)
{
	fSaveFileName = TString(filename);
}

TString WCSimFitterPlots::GetSaveFileName() const
{
	return fSaveFileName;
}

void WCSimFitterPlots::SetInputFileName(const char *inputfile)
{
	fInputFileName = TString(inputfile);
}

void WCSimFitterPlots::MakeSaveFile()
{
	std::cout << " *** WCSimFitter::MakeSaveFile() *** " << std::endl;
	if (fSaveFile != NULL)
	{
		std::cerr << "Save file already exists, so you can't make it again" << std::endl;
		return;
	}
	TDirectory *tmpd = gDirectory;

	TString fileNameCompare = fSaveFileName;
	fSaveFile = new TFile(fSaveFileName.Data(), "CREATE");

	std::cout << "  Plots to be saved in file: " << fSaveFileName.Data() << std::endl;

	tmpd->cd();
}
