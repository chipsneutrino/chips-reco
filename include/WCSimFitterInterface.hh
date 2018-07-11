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
class WCSimOutputTree;
class WCSimPiZeroFitter;
class WCSimCosmicFitter;

class WCSimFitterInterface {
	public:
		virtual ~WCSimFitterInterface();
		WCSimFitterInterface();

		void SetInputFileName(const char * inputfile, bool modifyFile);

		void SetMakeFits(bool makeFits);

		void AddFitterConfig(WCSimFitterConfig * config);

		void AddFitterPlots(WCSimFitterPlots * plots);

		void InitFitter(WCSimFitterConfig * config);

		void LoadWCSimData();

		void InitFitterPlots(TString outputName, WCSimFitterConfig * config);

		void Run();

	private:
		// Fit configuration parameters
		int fNumFits;
		std::vector<WCSimFitterConfig *> fFitterConfigs;

		// Input file parameters
		bool fModifyInputFile;
		TString fFileName;
		TString fOutputName;

		// Fitters, selections, particle identification, output etc...
		WCSimLikelihoodFitter * fFitter;
		WCSimPiZeroFitter * fPiZeroFitter;
		WCSimCosmicFitter * fCosmicFitter;
		WCSimOutputTree * fOutputTree;

		bool fMakeFits;

		bool fMakePlots;
		TString fPlotsName;
		WCSimFitterPlots * fFitterPlots;

		ClassDef(WCSimFitterInterface,0)
};

#endif /* WCSIMFITTERINTERFACE_HH_ */
