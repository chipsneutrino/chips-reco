/*
 * WCSimFitterInterface.hh
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#ifndef WCSIMFITTERINTERFACE_HH_
#define WCSIMFITTERINTERFACE_HH_

#include <TString.h>
#include "WCSimTrackParameterEnums.hh"
class WCSimLikelihoodFitter;
class WCSimFitterConfig;
class WCSimFitterPlots;
class WCSimFitterTree;
class WCSimPiZeroFitter;


class WCSimFitterInterface {
public:
	virtual ~WCSimFitterInterface();
	WCSimFitterInterface();

    void InitOutputFiles(); 
	void InitFitter();

	// void SetFile(const char * fileName);
	void SetNumTracks( unsigned int numTracks );
	unsigned int GetNumTracks() const;

	void SetTrackType( unsigned int numTrack, const char * typeName);
  TrackType::Type GetTrackType(const unsigned int &numTrack);


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

	void SetFirstEventToFit(int iEvt);
	int GetFirstEventToFit() const;

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

	void SetIsPiZeroFit(const bool &isPiZero);
	bool GetIsPiZeroFit() const;

	void SetForcePiZeroMass(const bool &doIt = true);
	bool GetForcePiZeroMass() const;

	void PrintFitConfiguration();
	void PrintPlotsConfiguration();
	void PrintSurfaceConfiguration();
  void Print();

	void Run();

  void SaveResults();
  void SaveProfiles();
	void SetInputFileName(const char * inputfile);

private:

	WCSimFitterConfig * fFitterConfig;
	TString fFileName;
	unsigned int fNumTracks;
	WCSimLikelihoodFitter * fFitter;
	WCSimPiZeroFitter * fPiZeroFitter;
	WCSimFitterPlots * fFitterPlots;
	WCSimFitterTree * fFitterTree;

	Bool_t fMakeFits;
	Bool_t fMakeSurfaces;

  ClassDef(WCSimFitterInterface,0)
};

#endif /* WCSIMFITTERINTERFACE_HH_ */
