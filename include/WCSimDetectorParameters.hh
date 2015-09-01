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
class TH1F;

class WCSimDetectorParameters : public TObject {
public:
	static WCSimDetectorParameters * Instance();
	static double WavelengthAveragedQE(const std::string &pmtName);
	static double QEAveragedRefIndex(const std::string &pmtName);

	double GetWavelengthAveragedQE(const std::string &pmtName);
	double GetQEAveragedRefIndex(const std::string &pmtName);



private:
	WCSimDetectorParameters();
	virtual ~WCSimDetectorParameters();

	void OpenFile();
	bool IsInMap(const std::string &pmtName, std::map<std::string, double> * map);
	double WorkOutAverageRefIndex(const std::string &pmtName);
	double WorkOutAverageQE(const std::string &pmtName);
	double AverageHistWithGraph(TH1F * hist, TGraph * graph);
	double AverageHistWithGraph(TH1F * hist, TGraph * graph, const double &min, const double &max);
	TH1F * MultiplyHistByGraph(TH1F * hist, TGraph * graph);

	TH1F * GetWavelengthSpectrum(const std::string &pmtName);



  TFile* fSpectrumFile;

	WCSimPMTManager fPMTManager;

	std::map<std::string, double> fAverageQEMap;
	std::map<std::string, double> fQEAveragedRefIndexMap;

  ClassDef(WCSimDetectorParameters,1)
};

#endif /* INCLUDE_WCSIMDETECTORPARAMETERS_HH_ */
