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
	static double WavelengthAveragedQE(const TrackType::Type &type, const std::string &pmtName);
	static double QEAveragedRefIndex(const TrackType::Type &type, const std::string &pmtName);

	double GetWavelengthAveragedQE(const TrackType::Type &type, const std::string &pmtName);
	double GetQEAveragedRefIndex(const TrackType::Type &type, const std::string &pmtName);



private:
	WCSimDetectorParameters();
	virtual ~WCSimDetectorParameters();

	void OpenFile(const TrackType::Type &type);
	bool IsInMap(const TrackType::Type &type, const std::string &pmtName, std::map<TrackType::Type, std::map<std::string, double> > * map);
	bool IsInMap(const TrackType::Type &type, std::map<TrackType::Type, TFile*> *map);
	double WorkOutAverageRefIndex(const TrackType::Type &type, const std::string &pmtName);
	double WorkOutAverageQE(const TrackType::Type &type, const std::string &pmtName);
	double AverageHistWithGraph(TH1F * hist, TGraph * graph);
	TH1F * MultiplyHistByGraph(TH1F * hist, TGraph * graph);

	TH1F * GetWavelengthSpectrum(const TrackType::Type &type, const std::string &pmtName);



	std::map<TrackType::Type, TFile*> fFileMap;

	WCSimPMTManager fPMTManager;

	std::map<TrackType::Type, std::map<std::string, double> > fAverageQEMap;
	std::map<TrackType::Type, std::map<std::string, double> > fQEAveragedRefIndexMap;

  ClassDef(WCSimDetectorParameters,1)
};

#endif /* INCLUDE_WCSIMDETECTORPARAMETERS_HH_ */
