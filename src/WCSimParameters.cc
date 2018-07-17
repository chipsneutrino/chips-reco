#include "WCSimParameters.hh"

#include "WCSimGeometry.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

ClassImp (WCSimParameters)

static WCSimParameters* fgParameters = 0;

WCSimParameters* WCSimParameters::Instance() {
	if (!fgParameters) {
		fgParameters = new WCSimParameters();
	}

	if (!fgParameters) {
		assert (fgParameters);
	}

	if (fgParameters) {

	}

	return fgParameters;
}

WCSimParameters::WCSimParameters() {

	// Slicer parameters
	fSlicerClusterDistance = 250.0;
	fSlicerMinSize = 25;
	fSlicerChargeCut = 1.0;
	fSlicerTimeCut = 30.0;
	fIterateSlicing = false;

	// Veto slicing parameters
	fVetoClusterDistance = 500.0;
	fVetoMinSize = 5;
	fVetoMinChargeCut = 2.0;
	fVetoMaxChargeCut = 20.0;
	fVetoTimeCut = 100.0;

	// Integral parameters
	fCalculateIntegrals = false;
	fTruncateIntegrals = false;
	fConstrainExtent = true;

	// Likelihood tuner parameters
	fUseTransmission = true;
	fUseAngularEfficiency = true;
	fUseGlassCathodeReflection = true;
	fUseScatteringTable = true;

	// Fitter parameters
	fUseTime = true;
	fUseCharge = true;
	fEqualiseChargeAndTime = false;
	fSaveWCSimRootEvent = false;
	fDigiType = "kPoisson";
	fSaveSeedInfo = false;
	fSaveStageInfo = false;
	fSaveHitComparison = false;

	// Speed of light parameters
	fUseCustomParticleSpeed = false;
	fUseCustomSpeedOfLight = false;
	fUseFittedSpeedOfLight = true;
	fCustomParticleSpeed = 1.0;
	fCustomSpeedOfLight = 1.0;
	fFittedSpeedOfLight = 0.7219;

	// Not currently used parameters
	fUseSimpleTimeResolution = false;
	fUseSimpleTimeSlew = false;
	fUseSimpleRefractiveIndex = false;

	// Load the file with the default configuration
	// and update the member variables above as necessary
	fConfName = getenv("WCSIMANAHOME");
	fConfName.append("/config/default.cfg");
	this->SetConfigFile(fConfName.c_str());
	return;
}

WCSimParameters::~WCSimParameters() {
	// empty
}

Double_t WCSimParameters::TimeResolution(Double_t Q) {
	if (WCSimParameters::Instance()->GetUseSimpleTimeResolution()) {
		return WCSimParameters::Instance()->GetSimpleTimeResolution(Q);
	} else {
		return WCSimParameters::Instance()->GetTimeResolution(Q);
	}
}

Double_t WCSimParameters::TimeSlew(Double_t Q) {
	if (WCSimParameters::Instance()->GetUseSimpleTimeSlew()) {
		return WCSimParameters::Instance()->GetSimpleTimeSlew();
	} else {
		return WCSimParameters::Instance()->GetTimeSlew(Q);
	}
}

Double_t WCSimParameters::RefractiveIndex(Double_t L) {
	if (WCSimParameters::Instance()->GetUseSimpleRefractiveIndex()) {
		return WCSimParameters::Instance()->GetSimpleRefractiveIndex();
	} else {
		return WCSimParameters::Instance()->GetRefractiveIndex(L);
	}
}

void WCSimParameters::PrintParameters() {
	WCSimParameters::Instance()->RunPrintParameters();
}

void WCSimParameters::RunPrintParameters() {
	std::cout << std::endl << " *** WCSimParameters::PrintParameters() *** " << std::endl;

	std::cout << "SlicerClusterDistance = " << fSlicerClusterDistance << std::endl;
	std::cout << "SlicerMinSize = " << fSlicerMinSize << std::endl;
	std::cout << "SlicerChargeCut = " << fSlicerChargeCut << std::endl;
	std::cout << "SlicerTimeCut = " << fSlicerTimeCut << std::endl;
	std::cout << "IterateSlicing = " << fIterateSlicing << std::endl;
	std::cout << "VetoClusterDistance = " << fVetoClusterDistance << std::endl;
	std::cout << "VetoMinSize = " << fVetoMinSize << std::endl;
	std::cout << "VetoMinChargeCut = " << fVetoMinChargeCut << std::endl;
	std::cout << "VetoMaxChargeCut = " << fVetoMaxChargeCut << std::endl;
	std::cout << "VetoTimeCut = " << fVetoTimeCut << std::endl;
	std::cout << "CalculateIntegrals = " << fCalculateIntegrals << std::endl;
	std::cout << "TruncateIntegrals = " << fTruncateIntegrals << std::endl;
	std::cout << "ConstrainExtent = " << fConstrainExtent << std::endl;
	std::cout << "UseTransmission = " << fUseTransmission << std::endl;
	std::cout << "UseAngularEfficiency = " << fUseAngularEfficiency << std::endl;
	std::cout << "UseGlassCathodeReflection = " << fUseGlassCathodeReflection << std::endl;
	std::cout << "UseScatteringTable = " << fUseScatteringTable << std::endl;
	std::cout << "UseTime = " << fUseTime << std::endl;
	std::cout << "UseCharge = " << fUseCharge << std::endl;
	std::cout << "EqualiseChargeAndTime = " << fEqualiseChargeAndTime << std::endl;
	std::cout << "SaveWCSimRootEvent = " << fSaveWCSimRootEvent << std::endl;
	std::cout << "DigiType = " << fDigiType << std::endl;
	std::cout << "SaveSeedInfo = " << fSaveSeedInfo << std::endl;
	std::cout << "SaveStageInfo = " << fSaveStageInfo << std::endl;
	std::cout << "SaveHitComparison = " << fSaveHitComparison << std::endl;
	std::cout << "UseCustomParticleSpeed = " << fUseCustomParticleSpeed << std::endl;
	std::cout << "UseCustomSpeedOfLight = " << fUseCustomSpeedOfLight << std::endl;
	std::cout << "UseFittedSpeedOfLight = " << fUseFittedSpeedOfLight << std::endl;
	std::cout << "CustomParticleSpeed = " << fCustomParticleSpeed << std::endl;
	std::cout << "CustomSpeedOfLight = " << fCustomSpeedOfLight << std::endl;
	std::cout << "FittedSpeedOfLight = " << fFittedSpeedOfLight << std::endl;
	std::cout << "UseSimpleTimeResolution = " << fUseSimpleTimeResolution << std::endl;
	std::cout << "UseSimpleTimeSlew = " << fUseSimpleTimeSlew << std::endl;
	std::cout << "UseSimpleRefractiveIndex = " << fUseSimpleRefractiveIndex << std::endl;

	std::cout << " ******************************************* " << std::endl << std::endl;

	return;
}

void WCSimParameters::SetUseTimeOnly(Bool_t doIt) {
	if (doIt) {
		std::cout << " *** WCSimParameters::SetUseTimeOnly *** " << std::endl;
		fUseTime = true;
		fUseCharge = false;
	}
	return;
}

void WCSimParameters::SetUseChargeOnly(Bool_t doIt) {
	if (doIt) {
		std::cout << " *** WCSimParameters::SetUseChargeOnly *** " << std::endl;
		fUseCharge = true;
		fUseTime = false;
	}
	return;
}

void WCSimParameters::SetUseChargeAndTime(Bool_t doIt) {
	if (doIt) {
		std::cout << " *** WCSimParameters::SetUseChargeAndTime *** " << std::endl;
		fUseCharge = true;
		fUseTime = true;
	}
	return;
}

Double_t WCSimParameters::SpeedOfLight() {
	return 29.9792458;  // velocity of light [cm/ns]
}

Double_t WCSimParameters::CherenkovAngle() {
	return 42.0;  // degrees
}

Double_t WCSimParameters::ThetaC() {
	return 42.0;  // degrees
}

Double_t WCSimParameters::CosThetaC() {
	return 0.743144825477394244;  // return TMath::Cos(42.0*TMath::Pi()/180.0);
}

Double_t WCSimParameters::GetTimeResolution(Double_t Q) {
	/*
	 // Old Parameterisation (lifted from WCSim)
	 // ========================================
	 Double_t qpes = Q;
	 if( qpes<0.5 ) qpes = 0.5;
	 if( qpes>32.0 ) qpes = 32.0;
	 Double_t res = 0.33 + sqrt(2.0/qpes);
	 */

	/*
	 // Sep'2010: parameterisation, including scattered light:
	 // ======================================================
	 Double_t c0 = +0.271, c1 = +3.037, c2 = +2.543;

	 // Aug'2011: re-parameterisation, excluding scattered light:
	 // ========================================================
	 Double_t c0 = +0.013, c1 = +3.592, c2 = -1.635; // [c2 is -ve, take care]

	 // Nov'2011: re-parameterisation, for 200 kton geometry:
	 // =====================================================
	 Double_t c0 = -0.005, c1 = +3.634, c2 = -1.458; // [c2 is -ve, take care]
	 */

	Double_t qpes = Q;
	Double_t qpesLow = 0.0;
	if (qpes < 0.25)
		qpes = 0.25;
	if (qpes > 40.0)
		qpes = 40.0;

	Double_t c0 = -0.005;
	Double_t c1 = +3.634;
	Double_t c2 = -1.458;

	if (c2 < 0.0) {
		qpesLow = (2.0 * c2 / c1) * (2.0 * c2 / c1);
		if (qpes < qpesLow)
			qpes = qpesLow;
	}

	Double_t res = c0 + c1 / sqrt(qpes) + c2 / qpes;

	return res;
}

Double_t WCSimParameters::WCSimTimeResolution(Double_t q, Double_t timeConst) {

	// Copied from WCSim
	Double_t Q = (q > 0.5) ? q : 0.5;
	Double_t timingResolution = 0.33 + sqrt(timeConst / Q);
	// looking at SK's jitter function for 20" tubes
	if (timingResolution < 0.58)
		timingResolution = 0.58;

	return timingResolution;
}

Double_t WCSimParameters::GetTimeSlew(Double_t Q) {
	/*
	 // Sep'2010: parameterisation, including scattered light:
	 // ======================================================
	 Double_t c0 = +3.406, c1 = -2.423, c2 = +0.335;

	 // Aug'2011: re-parameterisation, excluding scattered light:
	 // =========================================================
	 Double_t c0 = +2.234, c1 = -1.362, c2 = +0.125;

	 // Nov'2011: re-parameterisation, for 200 kton geometry:
	 // =====================================================
	 Double_t c0 = +2.436, c1 = -1.291, c2 = +0.089;
	 */

	Double_t qpes = Q;
	if (qpes < 0.25)
		qpes = 0.25;
	if (qpes > 40.0)
		qpes = 40.0;

	Double_t c0 = +2.436;
	Double_t c1 = -1.291;
	Double_t c2 = +0.089;

	Double_t dt = c0 + c1 * log(qpes) + c2 * log(qpes) * log(qpes);

	return dt;
}

Double_t WCSimParameters::GetRefractiveIndex(Double_t r) {
	Double_t c = 29.98;
	Double_t n0 = 1.33;        // Old Attempt:
	Double_t L0 = 0.0;         // 40.0
	Double_t dndx = 0.000123;  // 0.00015

	Double_t L = r / c;

	Double_t n = n0 * (1.0 + dndx * (L - L0));

	return n;
}

Double_t WCSimParameters::GetSimpleTimeResolution(Double_t Q) {
	Double_t qpes = Q;
	if (qpes < 0.25)
		qpes = 0.25;
	if (qpes > 64.0)
		qpes = 64.0;

	Double_t res = 2.0 / sqrt(qpes);

	return res;
}

void WCSimParameters::SetConfigFile(const char* config) {
	std::cout << "Setting config file" << std::endl;
	fConfName = config;
	this->LoadConfig();
	return;
}

void WCSimParameters::LoadConfig() {
	// Open the file and test opening:
	std::cout << "Loading configuration file" << std::endl;
	std::ifstream inFile;
	inFile.open(fConfName.c_str());

	if (!inFile.is_open()) {
		std::cerr << "Error: " << __FILE__ << "  " << __LINE__ << " - could not open " << fConfName << std::endl;
		exit (EXIT_FAILURE);
	}

	// Read the file
	std::string line;
	Int_t lineNum = 0;
	while (getline(inFile, line)) {
		// std::cout << "Line is: " << line << std::endl;
		// std::cout << "lineNum = " << lineNum << std::endl;
		this->IgnoreComments(line);
		this->ParseLine(line, ++lineNum);
	}
	this->SetFromMap();

	return;
}

/// Erases comments beginning with // or #
void WCSimParameters::IgnoreComments(std::string &str) {
	// std::cout << "IgnoreComments" << std::endl;
	if (str.find("//") != str.npos)
		str.erase(str.find("//"));
	if (str.find("#") != str.npos)
		str.erase(str.find("#"));
	// std::cout << "Without comments " << str << std::endl;
	return;
}

const Bool_t WCSimParameters::IsBlankLine(std::string str) {
	// std::cout << "IsBlankLink" << "   " << str << std::endl;
	if (str.find_first_not_of(' ') == str.npos)
		return true;
	else
		return false;
}

const Bool_t WCSimParameters::IsGoodLine(std::string str) {
	// std::cout << "IsGoodLine: " << str << std::endl;

	Bool_t haveEquals = false;
	Bool_t haveLHS = false;
	Bool_t haveRHS = false;

	// Look for an equals sign
	if (str.find("=") == str.npos) {
		std::cout << "No \"=\" sign found in string: " << std::endl << str << std::endl;
	} else
		haveEquals = true;

	// Look for text on the LHS of = sign:
	std::string tempStr = str;
	tempStr.erase(0, tempStr.find_first_not_of("\t "));
	if (tempStr[0] != '0') {
		haveLHS = true;
	}

	// Look for text on RHS of = sign:
	tempStr = str;
	for (UInt_t rhs = tempStr.find("=") + 1; rhs < tempStr.length(); ++rhs) {
		if (tempStr[rhs] != '\t' && tempStr[rhs] != ' ')
			haveRHS = true;
	}
	// std::cout << "haveEquals = " << haveEquals << std::endl
	//           << "haveLHS = " << haveLHS << std::endl
	//           << "haveRHS = " << haveRHS << std::endl
	//           << "All = " << (haveLHS && haveRHS && haveEquals) << std::endl;

	return (haveEquals && haveLHS && haveRHS);
}

void WCSimParameters::ExtractPair(std::string &lhs, std::string &rhs, std::string str) {
	// std::cout << "ExtractPair" << std::endl;
	UInt_t splitPos = str.find("=");

	// Get left hand side of = sign and strip whitespace
	lhs = str.substr(0, splitPos);
	lhs.erase(std::remove(lhs.begin(), lhs.end(), ' '), lhs.end());
	lhs.erase(std::remove(lhs.begin(), lhs.end(), '\t'), lhs.end());

	// And the other side
	rhs = str.substr(splitPos + 1);
	rhs.erase(std::remove(rhs.begin(), rhs.end(), ' '), rhs.end());
	rhs.erase(std::remove(rhs.begin(), rhs.end(), '\t'), rhs.end());
	// std::cout << "str = " << str << std::endl
	//           << "lhs = " << lhs << std::endl
	//           << "rhs = " << rhs << std::endl;
	return;
}

void WCSimParameters::ParseLine(std::string str, Int_t lineNum) {
	// std::cout << "ParseLine: " << str << std::endl;
	if (!(this->IsGoodLine(str))) {
		std::cerr << "Error: line " << lineNum << "has the wrong format" << std::endl;
		exit (EXIT_FAILURE);
	}

	std::string lhs, rhs;
	this->ExtractPair(lhs, rhs, str);
	this->AddToMap(lhs, rhs);
	return;
}

void WCSimParameters::AddToMap(std::string lhs, std::string rhs) {
	// std::cout << "AddToMap" << std::endl;
	if (fMap.find(lhs) == fMap.end()) {
		fMap[lhs] = rhs;
	} else {
		std::cerr << "Error: map already contains the key: " << lhs << std::endl;
		exit (EXIT_FAILURE);
	}
	std::cout << lhs << "  " << fMap[lhs] << std::endl;
	return;
}

void WCSimParameters::SetFromMap() {
	// std::cout << "SetFromMap" << std::endl;
	std::map<std::string, std::string>::const_iterator itr = fMap.begin();

	// Loop through the map, checking if any of the keys correspond to things we can set
	for (; itr != fMap.end(); ++itr) {
		std::cout << (*itr).first << "   " << (*itr).second << std::endl;
		if ((*itr).first.compare("SlicerClusterDistance") == 0) {
			std::stringstream ss(itr->second);
			float val = 0.0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fSlicerClusterDistance = val;
		} else if ((*itr).first.compare("SlicerMinSize") == 0) {
			std::stringstream ss(itr->second);
			int val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a int" << std::endl;
				exit (EXIT_FAILURE);
			}
			fSlicerMinSize = val;
		} else if ((*itr).first.compare("SlicerChargeCut") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fSlicerChargeCut = val;
		} else if ((*itr).first.compare("SlicerTimeCut") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fSlicerTimeCut = val;
		} else if ((*itr).first.compare("IterateSlicing") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fIterateSlicing = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fIterateSlicing = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("VetoClusterDistance") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fVetoClusterDistance = val;
		} else if ((*itr).first.compare("VetoMinSize") == 0) {
			std::stringstream ss(itr->second);
			int val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a int" << std::endl;
				exit (EXIT_FAILURE);
			}
			fVetoMinSize = val;
		} else if ((*itr).first.compare("VetoMinChargeCut") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fVetoMinChargeCut = val;
		} else if ((*itr).first.compare("VetoMaxChargeCut") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fVetoMaxChargeCut = val;
		} else if ((*itr).first.compare("VetoTimeCut") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fVetoTimeCut = val;
		} else if ((*itr).first.compare("CalculateIntegrals") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fCalculateIntegrals = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fCalculateIntegrals = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("TruncateIntegrals") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fTruncateIntegrals = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fTruncateIntegrals = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("ConstrainExtent") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fConstrainExtent = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fConstrainExtent = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseTransmission") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseTransmission = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseTransmission = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseAngularEfficiency") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseAngularEfficiency = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseAngularEfficiency = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseGlassCathodeReflection") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseGlassCathodeReflection = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseGlassCathodeReflection = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseScatteringTable") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseScatteringTable = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseScatteringTable = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseTime") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseTime = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseTime = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseCharge") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseCharge = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseCharge = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("EqualiseChargeAndTime") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fEqualiseChargeAndTime = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fEqualiseChargeAndTime = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("SaveWCSimRootEvent") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fSaveWCSimRootEvent = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fSaveWCSimRootEvent = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("DigiType") == 0) {
			std::cout << "It's " << ((*itr).second) << std::endl;
			fDigiType = ((*itr).second);
		} else if ((*itr).first.compare("SaveSeedInfo") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fSaveSeedInfo = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fSaveSeedInfo = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("SaveStageInfo") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fSaveStageInfo = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fSaveStageInfo = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("SaveHitComparison") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fSaveHitComparison = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fSaveHitComparison = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseCustomParticleSpeed") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseCustomParticleSpeed = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseCustomParticleSpeed = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseCustomSpeedOfLight") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseCustomSpeedOfLight = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseCustomSpeedOfLight = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseFittedSpeedOfLight") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseFittedSpeedOfLight = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseFittedSpeedOfLight = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("CustomParticleSpeed") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fCustomParticleSpeed = val;
		} else if ((*itr).first.compare("CustomSpeedOfLight") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fCustomSpeedOfLight = val;
		} else if ((*itr).first.compare("FittedSpeedOfLight") == 0) {
			std::stringstream ss(itr->second);
			float val = 0;
			if (!(ss >> val)) {
				std::cerr << "Error: " << (*itr).second << " = " << (*itr).second << " should be a double" << std::endl;
				exit (EXIT_FAILURE);
			}
			fFittedSpeedOfLight = val;
		} else if ((*itr).first.compare("UseSimpleTimeResolution") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseSimpleTimeResolution = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseSimpleTimeResolution = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseSimpleTimeSlew") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseSimpleTimeSlew = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseSimpleTimeSlew = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		} else if ((*itr).first.compare("UseSimpleRefractiveIndex") == 0) {
			if ((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0) {
				fUseSimpleRefractiveIndex = true;
			} else if ((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0) {
				fUseSimpleRefractiveIndex = false;
			} else {
				std::cerr << "Error: " << (*itr).first << " = " << (*itr).second << " should equal true/false or 1/0"
						<< std::endl;
				exit (EXIT_FAILURE);
			}
		}
	}
	fMap.clear();
}
