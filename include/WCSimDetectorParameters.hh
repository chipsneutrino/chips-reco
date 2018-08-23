/*
 * WCSimDetectorParameters.hh
 *
 *  Created on: 29 Jul 2015
 *      Author: andy
 */

#ifndef INCLUDE_WCSIMDETECTORPARAMETERS_HH_
#define INCLUDE_WCSIMDETECTORPARAMETERS_HH_

#include <string>
#include <map>
#include "WCSimTrackParameterEnums.hh"
#include "WCSimPMTManager.hh"
#include "TFile.h"
#include "TObject.h"
class TGraph;
class TH1D;

class WCSimDetectorParameters: public TObject {
	public:
		static WCSimDetectorParameters * Instance();
		static TGraph * WavelengthAveragedQE(const std::string &pmtName);
		static double QEAveragedRefIndex(const std::string &pmtName);
		static double PMTExposeHeight(const std::string &pmtName);

		TGraph * GetWavelengthAveragedQE(const std::string &pmtName);
		double GetQEAveragedRefIndex(const std::string &pmtName);
		double GetPMTExposeHeight(const std::string &pmtName);

	private:
		WCSimDetectorParameters();
		virtual ~WCSimDetectorParameters();

		void OpenFile();
		bool IsInMap(const std::string &pmtName, std::map<std::string, double> * map);
		bool IsInMap(const std::string &pmtName, std::map<std::string, TGraph*> * map);
		double WorkOutAverageRefIndex(const std::string &pmtName);
		TGraph * WorkOutAverageQE(const std::string &pmtName);
		double WorkOutExposeHeight(const std::string &pmtName);
		double AverageHistWithGraph(TH1D * hist, TGraph * graph);
		double AverageHistWithGraph(TH1D * hist, TGraph * graph, const double &min, const double &max);
		TH1D * MultiplyHistByGraph(TH1D * hist, TGraph * graph);

		TH1D * GetWavelengthSpectrum(const std::string &pmtName);

		TFile* fSpectrumFile;
		TFile* fSpectrumVsDistanceFile;

		WCSimPMTManager fPMTManager;

		std::map<std::string, TGraph*> fAverageQEMap;
		std::map<std::string, double> fQEAveragedRefIndexMap;
		std::map<std::string, double> fExposeHeightMap;

		ClassDef(WCSimDetectorParameters,1)
};

#endif /* INCLUDE_WCSIMDETECTORPARAMETERS_HH_ */
