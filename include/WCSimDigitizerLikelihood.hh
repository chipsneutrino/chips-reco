/**
 * \class WCSimDigitizerLikelihood
 * This class is used to model the effect of digitization at
 * the PMT.  The WCSimChargelikelihood class gives the predicted
 * numer of photoelectrons incident on the PMT, and this class models
 * how likely that number is to be converted into the final
 * digitized P.E. recorded by the PMT
 */

#ifndef WCSIMDIGITIZERLIKELIHOOD_H
#define WCSIMDIGITIZERLIKELIHOOD_H
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
        kSK1pe,  /// A method derived from the default WCSim digitiser that samples the SK 1pe distribution
        kPMTSim, /// A method using the full dynode chain simulation
        kUnknown /// Error state
      };

      WCSimDigitizerLikelihood();
      virtual ~WCSimDigitizerLikelihood();

      /**
       * Convert a string to its corresponding WCSimDigitizerLikelihood::DigiType_t enum
       * @param str Name of digitizer
       * @return Corresponding DigiType_t enum
       */
      DigiType_t StringToDigiType( const std::string &str ) const;

      /**
       * Select the digitizer method to use for calculations
       * @param type Digitzer to use
       */
      void SetDigiType( WCSimDigitizerLikelihood::DigiType_t type );

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
       * Calculate the -2log(likelihood) of the digitizer returning the measured
       * charge, given the predicted mean number of photoelectrons
       * @param undigi Predicted mean number of photoelectrons at the PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return -2 log(likelihood) to measure the digitized charge given the predicted mean number of photoelectrons
       */
      Double_t GetMinus2LnL( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood of the digitizer returning the measured
       * charge given the predicted mean number of photoelectrons
       * @param undigi Predicted mean number of photoelectrons as the PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return Likelihood to measure the digitized P.E. given the predicted mean number of photoelectrons
       */
      Double_t GetLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Get the most likely number of digitized P.E. measured by the PMT
       * given the predicted mean number of photoelectrons hitting it
       * @param undigi Predicted mean number of photoelectrons hitting the PMT
       * @return Most likely number of digitized P.E. measured by the PMT
       */
      Double_t GetExpectation( const Double_t & undigi );

    protected:
    
    private:

      /**
       * Calculate -2log(likelihood) for measured charge, assuming a
       * Poisson distribution for the digitizer
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return -2log(likelihood) of getting the recorded P.E.
       */
      Double_t GetPoissonMinus2LnL( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for measured charge, assuming a
       * Poisson distribution for the digitizer
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return likelihood of getting the recorded P.E.
       */
      Double_t GetPoissonLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Mean digitized P.E. returned by the PMT for a given predicted mean
       * number of photoelectrons, assuming a Poisson distribution (i.e. same as the
       * number of photoelectrons in this case)
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @return Most-likely resulting digitized P.E.
       */
      Double_t GetPoissonExpectation( const Double_t & undigi );
      

      /**
       * Calculate -2log(likelihood) for measured charge, assuming the full
       * PMT dynode chain simulation from WCSim
       * @param undigi Predicted mean number of photoelectrons as the PMT
       * @param digi Digitised  P.E. recorded by the PMT
       * @return -2log(likelihood) of getting the recorded P.E.
       */
      Double_t GetPMTSimMinus2LnL( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for measured charge, assuming the full
       * PMT dynode chain simulation from WCSim
       * @param undigi Predicted mean number of photoelectrons as the PMT
       * @param digi Digitised  P.E. recorded by the PMT
       * @return likelihood of getting the recorded P.E.
       */
      Double_t GetPMTSimLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate -2log(likelihood) for measured charge, assuming
       * the WCSim digitizer
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return -2log(likelihood) of getting the recorded P.E.
       */
      Double_t GetWCSimMinus2LnL( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for measured charge, assuming
       * the WCSim digitizer
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return likelihood of getting the recorded P.E.
       */
      Double_t GetWCSimLikelihood( Double_t undigi, const Double_t &digi );

      /**
       * Mean digitized P.E. returned by the PMT for a given predicted mean
       * number of photoelectrons, assuming the WCSim digitizer
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @return Most-likely resulting digitized P.E.
       */
      Double_t GetWCSimExpectation( const Double_t & undigi );

      /**
       * Calculate the likelihood for the measured charge using the WCSim
       * digitzer, by sampling probability histograms (used for low charges)
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @param digi Digitized P.E. recorded aby the PMT
       * @return
       */
      Double_t GetWCSimPickerLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for the measured charge using the WCSim
       * digitzer, by using a Gaussian (x) Exponential (medium P.E.)
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @param digi Digitized P.E. recorded aby the PMT
       * @return
       */
      Double_t GetWCSimGausExpoLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for the measured charge using the WCSim
       * digitzer, by using a simple Gaussian (high P.E.)
       * @param undigi Predicted mean number of photoelectrons at PMT
       * @param digi Digitized P.E. recorded aby the PMT
       * @return
       */
      Double_t GetWCSimGausLikelihood( const Double_t &undigi, const Double_t &digi );



      DigiType_t fType;  ///< Which digitizer method to use

      std::map<std::string, DigiType_t> fDigiTypeNames; ///< For converting strings to their corresponding DigiType_t
      
      // WCSim repeatedly samples a 1pe distribution.  For hits < 10pe I've
      // already done this to build a PDF histogram which these variables 
      // point to
      TFile * fSK1peFile;  ///< File holding probability density histograms for the sub-10pe WCSim digitizer
      TH2D  * fSK1peHist;  ///< Probabilitiy density hitogram for the the sub-10pe WCSim digitizer

      // WCSim also has a PMT simulation that amplifies the electrons along the full dynode chain and
      // introduces nonlinearity when the current is high.  I parametrised this and filled some histograms
      // so the probability can be looked-up
      TFile * fPMTSimFile; ///< File holding the probability histogram for the dynode chain simulation method
      TH2F  * fPMTSimHist; ///< Histogram describing the probability histogram for the dynode chain simulation method

      // The WCSim digitizer samples the 1pe distribution repeatedly,
      // then applies a threshold function,
      // then multiplies by an efficiency term
      Double_t fEfficiency; ///< Efficiency value used in the WCSim digitizer

      Double_t fMinimum; ///< Minimum nonzero likelihood from histogram
};

#endif // WCSIMDIGITIZERLIKELIHOOD_H
