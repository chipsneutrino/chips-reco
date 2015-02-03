/*
 * WCSimFitterPlots.hh
 *
 *  Created on: 3 Nov 2014
 *      Author: andy
 */

#ifndef WCSIMFITTERPLOTS_HH_
#define WCSIMFITTERPLOTS_HH_
#include "WCSimFitterParameters.hh"
#include "WCSimLikelihoodTrack.hh"
#include <map>
#include <vector>
#include <TObjArray.h>
class TFile;
class TH1D;
class TH2D;
class WCSimFitterConfig;
class WCSimRootEvent;

class WCSimFitterPlots {


public:
	WCSimFitterPlots();
	virtual ~WCSimFitterPlots();

	void SetSaveFileName(const char * filename);
	TString GetSaveFileName() const;

	void SetNumSurfaceBins(unsigned int numBins);
	unsigned int GetNumSurfaceBins();

	void SetPlotForEachEvent(const char * name, bool doIt = true);
	bool GetPlotForEachEvent(const char * name) const;

	void SetPlotRecoMinusTrue(const char * name, bool doIt = true);
	bool GetPlotRecoMinusTrue(const char * name) const;

	void Make1DSurface(const char * name, bool doIt = true, unsigned int trackNum = 0);
	bool GetMake1DSurface(const char * name, unsigned int trackNum = 0) const;

	void Make2DSurface(const char * name, const char * name2, bool doIt = true, unsigned int trackNum = 0, unsigned int trackNum2 = 0);
	bool GetMake2DSurface(const char * name, const char * name2, unsigned int trackNum = 0, unsigned int trackNum2 = 0) const;
	std::vector<std::pair<unsigned int, FitterParameterType::Type> > GetAll1DSurfaceKeys() const;
	std::vector<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > > GetAll2DSurfaceKeys() const;

	void MakeHistograms(WCSimFitterConfig * fitterConfig);
	void MakePlotsForEachEvent( WCSimFitterConfig * fitterConfig);
	void MakeRecoMinusTrue(WCSimFitterConfig * fitterConfig);
	void Make1DSurfaces(WCSimFitterConfig * fitterConfig);
	void Make2DSurface(WCSimFitterConfig * fitterConfig);


	void CreateNtuple(WCSimFitterConfig * fitterConfig);
	void FillNtuple(WCSimFitterConfig * fitterConfig, std::vector<WCSimLikelihoodTrack> bestFits);
	void FillPlots(std::vector<WCSimLikelihoodTrack> bestFits);
	void FillRecoMinusTrue(std::vector<WCSimLikelihoodTrack> bestFits, std::vector<WCSimLikelihoodTrack*> * trueTracks);
  double Get1DSurfaceBinCenter(std::pair<unsigned int, FitterParameterType::Type> theProfile, int binNum);
	void Fill1DProfile( std::pair<unsigned int, FitterParameterType::Type> theProfile, int binNum, double minus2LnL);
  void Get2DSurfaceBinCenters(std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > theProfile, int binNumX, int binNumY, double &x, double &y);
  double Get2DSurfaceBinCenterX(std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > theProfile, int binNumX);
  double Get2DSurfaceBinCenterY(std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > theProfile, int binNumY);
	void Fill2DProfile( std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > theProfile, int binNumX, int binNumY, double minus2LnL);
	void SavePlots();

	void Print();
	void PrintVariables();
	void PrintRecoMinusTrue();
	void PrintSurfaces();
	void Print1DSurfaces();
	void Print2DSurfaces();

private:
	TString fSaveFileName;
	TFile * fSaveFile;
	std::map<std::pair<unsigned int, FitterParameterType::Type>, TH1D*> fSurfaces1D;
	std::map<std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> >, TH2D*> fSurfaces2D;
	std::map<FitterParameterType::Type, std::vector<TH1D*> > fForEachEvent;
	std::map<FitterParameterType::Type, std::vector<TH1D*> > fRecoMinusTrue;
	int fNumSurfaceBins;
};

#endif /* WCSIMFITTERPLOTS_HH_ */