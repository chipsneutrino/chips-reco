/**
 * \class WCSimDigitizerLikelihood
 * This class is used to model the effect of digitization at
 * the PMT.  The WCSimChargelikelihood class gives the predicted
 * numer of photoelectrons incident on the PMT, and this class models
 * how likely that number is to be converted into the final
 * digitized P.E. recorded by the PMT
 */

#pragma once

#include "Rtypes.h"
#include <string>
#include <map>

class TFile;
class TH2D;
class TH2F;

class WCSimDigitizerLikelihood
{
public:
	/// Enumerate the different models for digitising PMT hits
	enum DigiType_t
	{
		kPoisson, /// A simple Poisson likelihood
		kSK1pe,	  /// A method derived from the default WCSim digitiser that samples the SK1pe distribution
		kPMTSim,  /// A method using the full dynode chain simulation
		kTOT,	  /// Digi simulation that tries to mimic the actual PMTs to be used in CHIPS-10
		kUnknown  /// Error state
	};

	WCSimDigitizerLikelihood();
	virtual ~WCSimDigitizerLikelihood();

	/**
		 * Convert a string to its corresponding WCSimDigitizerLikelihood::DigiType_t enum
		 * @param str Name of digitizer
		 * @return Corresponding DigiType_t enum
		 */
	DigiType_t StringToDigiType(const std::string &str) const;

	/**
		 * Select the digitizer method to use for calculations
		 * @param type Digitzer to use
		 */
	void SetDigiType(WCSimDigitizerLikelihood::DigiType_t type);

	/**
		 * The probability density functions based on the SK 1pe distribution are saved in histograms.
		 * This opens those files and extracts the histograms.
		 */
	void OpenSKPDFs();

	/**
		 * The PDFs based on the dynode chain simulation from WCSim are stored in a 2D histogram
		 * This opens the file containing the histogram and extracts it
		 */
	void OpenPMTSimPDFs();

	/**
		 * The PDFs based on the TOT method for CHIPS-10 are stored in a 2D histogram
		 * This opens the file containing the histogram and extracts it
		 */
	void OpenTOTPDFs();

	/**
		 * Calculate the -2log(likelihood) of the digitizer returning the measured
		 * charge, given the predicted mean number of photoelectrons
		 * @param undigi Predicted mean number of photoelectrons at the PMT
		 * @param digi Digitized P.E. recorded by the PMT
		 * @param PMTName the name of the pmt if different PMTs are required
		 * @return -2 log(likelihood) to measure the digitized charge given the predicted mean number of photoelectrons
		 */
	Double_t GetMinus2LnL(const Double_t &undigi, const Double_t &digi, std::string PMTName);

	/**
		 * Calculate the likelihood of the digitizer returning the measured
		 * charge given the predicted mean number of photoelectrons
		 * @param undigi Predicted mean number of photoelectrons as the PMT
		 * @param digi Digitized P.E. recorded by the PMT
		 * @param PMTName the name of the pmt if different PMTs are required
		 * @return Likelihood to measure the digitized P.E. given the predicted mean number of photoelectrons
		 */
	Double_t GetLikelihood(const Double_t &undigi, const Double_t &digi, std::string PMTName);

	/**
		 * Get the most likely number of digitized P.E. measured by the PMT
		 * given the predicted mean number of photoelectrons hitting it
		 * @param undigi Predicted mean number of photoelectrons hitting the PMT
		 * @return Most likely number of digitized P.E. measured by the PMT
		 */
	Double_t GetExpectation(const Double_t &undigi);

protected:
private:
	/**
		 * Calculate -2log(likelihood) for measured charge, assuming a
		 * Poisson distribution for the digitizer
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @param digi Digitized P.E. recorded by the PMT
		 * @return -2log(likelihood) of getting the recorded P.E.
		 */
	Double_t GetPoissonMinus2LnL(const Double_t &undigi, const Double_t &digi);

	/**
		 * Calculate the likelihood for measured charge, assuming a
		 * Poisson distribution for the digitizer
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @param digi Digitized P.E. recorded by the PMT
		 * @return likelihood of getting the recorded P.E.
		 */
	Double_t GetPoissonLikelihood(const Double_t &undigi, const Double_t &digi);

	/**
		 * Calculate -2log(likelihood) for measured charge, assuming
		 * the sk1pe digitizer
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @param digi Digitized P.E. recorded by the PMT
		 * @return -2log(likelihood) of getting the recorded P.E.
		 */
	Double_t GetSK1peMinus2LnL(const Double_t &undigi, const Double_t &digi);

	/**
		 * Calculate the likelihood for measured charge, assuming
		 * the sk1pe digitizer
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @param digi Digitized P.E. recorded by the PMT
		 * @return likelihood of getting the recorded P.E.
		 */
	Double_t GetSK1peLikelihood(Double_t undigi, const Double_t &digi);

	/**
		 * Mean digitized P.E. returned by the PMT for a given predicted mean
		 * number of photoelectrons, assuming the SK1pe digitizer
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @return Most-likely resulting digitized P.E.
		 */
	Double_t GetSK1peExpectation(const Double_t &undigi);

	/**
		 * Calculate the likelihood for the measured charge using the SK1pe
		 * digitzer, by sampling probability histograms (used for low charges)
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @param digi Digitized P.E. recorded aby the PMT
		 * @return
		 */
	Double_t GetSK1pePickerLikelihood(const Double_t &undigi, const Double_t &digi);

	/**
		 * Calculate the likelihood for the measured charge using the SK1pe
		 * digitzer, by using a Gaussian (x) Exponential (medium P.E.)
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @param digi Digitized P.E. recorded aby the PMT
		 * @return
		 */
	Double_t GetSK1peGausExpoLikelihood(const Double_t &undigi, const Double_t &digi);

	/**
		 * Calculate the likelihood for the measured charge using the SK1pe
		 * digitzer, by using a simple Gaussian (high P.E.)
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @param digi Digitized P.E. recorded aby the PMT
		 * @return
		 */
	Double_t GetSK1peGausLikelihood(const Double_t &undigi, const Double_t &digi);

	/**
		 * Mean digitized P.E. returned by the PMT for a given predicted mean
		 * number of photoelectrons, assuming a Poisson distribution (i.e. same as the
		 * number of photoelectrons in this case)
		 * @param undigi Predicted mean number of photoelectrons at PMT
		 * @return Most-likely resulting digitized P.E.
		 */
	Double_t GetPoissonExpectation(const Double_t &undigi);

	/**
		 * Calculate -2log(likelihood) for measured charge, assuming the full
		 * PMT dynode chain simulation from WCSim
		 * @param undigi Predicted mean number of photoelectrons as the PMT
		 * @param digi Digitised  P.E. recorded by the PMT
		 * @return -2log(likelihood) of getting the recorded P.E.
		 */
	Double_t GetPMTSimMinus2LnL(const Double_t &undigi, const Double_t &digi);

	/**
		 * Calculate the likelihood for measured charge, assuming the full
		 * PMT dynode chain simulation from WCSim
		 * @param undigi Predicted mean number of photoelectrons as the PMT
		 * @param digi Digitised  P.E. recorded by the PMT
		 * @return likelihood of getting the recorded P.E.
		 */
	Double_t GetPMTSimLikelihood(const Double_t &undigi, const Double_t &digi);

	/**
		 * Calculate -2log(likelihood) for measured charge, assuming the full
		 * PMT dynode chain simulation from WCSim
		 * @param undigi Predicted mean number of photoelectrons as the PMT
		 * @param digi Digitised  P.E. recorded by the PMT
		 * @return -2log(likelihood) of getting the recorded P.E.
		 */
	Double_t GetTOTMinus2LnL(const Double_t &undigi, const Double_t &digi, std::string PMTName);

	/**
		 * Calculate the likelihood for measured charge, assuming the
		 * TOT digi method for CHIPS-10
		 * @param undigi Predicted mean number of photoelectrons as the PMT
		 * @param digi Digitised  P.E. recorded by the PMT
		 * @return likelihood of getting the recorded P.E.
		 */
	Double_t GetTOTLikelihood(const Double_t &undigi, const Double_t &digi, std::string PMTName);

	DigiType_t fType; ///< Which digitizer method to use

	std::map<std::string, DigiType_t> fDigiTypeNames; ///< For converting strings to their corresponding DigiType_t

	// WCSim repeatedly samples a 1pe distribution.  For hits < 10pe I've
	// already done this to build a PDF histogram which these variables
	// point to
	TFile *fSK1peFile; ///< File holding probability density histograms for the sub-10pe SK1pe digitizer
	TH2D *fSK1peHist;  ///< Probability density hitogram for the the sub-10pe SK1pe digitizer

	// The WCSim digitizer samples the 1pe distribution repeatedly,
	// then applies a threshold function,
	// then multiplies by an efficiency term
	Double_t fEfficiency; ///< Efficiency value used in the WCSim digitizer

	// WCSim also has a PMT simulation that amplifies the electrons along the full dynode chain and
	// introduces nonlinearity when the current is high.  I parametrised this and filled some histograms
	// so the probability can be looked-up
	TFile *fPMTSimFile; ///< File holding the probability histogram for the dynode chain simulation method
	TH2D *fPMTSimHist;	///< Histogram describing the probability histogram for the dynode chain simulation method

	// WCSim also has a Time-over-threshold digitisation method to mimic the actual PMTs to be used in CHIPS-10
	TFile *fTOTFile;	   ///< File holding the probability histogram for the TOT simulation method
	TH2D *fTOTNikhefHist;  ///< Histogram describing the probability histogram for the TOT simulation method for the Nikhef PMTs
	TH2D *fTOTMadisonHist; ///< Histogram describing the probability histogram for the TOT simulation method for the Madison PMTs

	Double_t fMinimum; ///< Minimum nonzero likelihood from histogram
};
