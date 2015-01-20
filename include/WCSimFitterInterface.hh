/*
 * WCSimFitterInterface.hh
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#ifndef WCSIMFITTERINTERFACE_HH_
#define WCSIMFITTERINTERFACE_HH_

#include <TString.h>
class WCSimLikelihoodFitter;
class WCSimFitterConfig;
class WCSimFitterPlots;


class WCSimFitterInterface {
public:
	static WCSimFitterInterface * Instance();
	virtual ~WCSimFitterInterface();
  void Init();

	// void SetFile(const char * fileName);
	void SetNumTracks( unsigned int numTracks );
	unsigned int GetNumTracks() const;

	void SetTrackType( unsigned int numTrack, const char * typeName);


	void FixParameter(  unsigned int numTrack, const char * name, bool doIt = true );
	void FreeParameter( unsigned int numTrack, const char * name, bool doIt = true );
	void JoinParametersTogether(unsigned int numTrack1, unsigned int numTrack2, const char * name);

	void SetParMin(   unsigned int numTrack, const char * name, double min);
	void SetParMax(   unsigned int numTrack, const char * name, double max);
	void SetParStart( unsigned int numTrack, const char * name, double start);
	void SetParRange( unsigned int numTrack, const char * name, double min, double max);
	void SetParameter( unsigned int numTrack, const char * name, double min, double max, double start, bool fixed);

	Double_t GetParMin(   unsigned int numTrack, const char * name );
	Double_t GetParMax(   unsigned int numTrack, const char * name );
	Double_t GetParStart( unsigned int numTrack, const char * name );
	std::pair<Double_t, Double_t> GetParRange( unsigned int numTrack, const char * name );

	void SetMakeFits( bool doIt = true);
	bool GetMakeFits();

	void SetNumEventsToFit(int nEvts);
	int  GetNumEventsToFit();

	void PlotForEachEvent(const char * name, bool doIt = true);
	bool GetPlotForEachEvent(const char * name);

	void PlotRecoMinusTrue(const char * name, bool doIt = true);
	bool GetPlotRecoMinusTrue(const char * name);

	void SetMakeSurfaces(bool doIt = true);
	bool GetMakeSurfaces();

  void SetNumSurfaceBins( unsigned int nBins );
  unsigned int GetNumSurfaceBins() const;

	void Make1DSurface(unsigned int nTrack, const char * name, bool doIt = true );
	bool GetMake1DSurface(unsigned int nTrack, const char * name );

	void Make2DSurface(unsigned int nTrack, const char * name, unsigned int nTrack2, const char * name2, bool doIt = true );
	bool GetMake2DSurface(unsigned int nTrack, const char * name, unsigned int nTrack2, const char * name2 );

	void PrintFitConfiguration();
	void PrintPlotsConfiguration();
	void PrintSurfaceConfiguration();
  void Print();

	void Run();

private:
	WCSimFitterInterface();

	TString fFileName;
	unsigned int fNumTracks;
	WCSimLikelihoodFitter * fFitter;
	WCSimFitterPlots * fFitterPlots;

	Bool_t fMakeFits;
	Bool_t fMakeSurfaces;

  ClassDef(WCSimFitterInterface,0)
};

#endif /* WCSIMFITTERINTERFACE_HH_ */
