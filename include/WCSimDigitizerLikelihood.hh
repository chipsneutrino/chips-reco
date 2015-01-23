/**
 * \class WCSimDigitizerLikelihood
 * This class is used to model the effect of digitization at
 * the PMT.  The WCSimChargelikelihood class gives the predicted
 * numer of photons incident on the PMT, and this class models
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


class WCSimDigitizerLikelihood
{
    public:

      /// Enumerate the different models for digitising PMT hits
      enum DigiType_t
      {
        kSimple, /// A simple Poisson likelihood
        kWCSim,  /// A method derived from the default WCSim digitiser
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
       * The WCSim probability density functions are save in histograms.
       * This opens those files and extracts the histograms.
       */
      void OpenPDFs();

      /**
       * Calculate the -2log(likelihood) of the digitizer returning the measured
       * charge, given the predicted mean number of photons
       * @param undigi Predicted mean number of photons at the PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return -2 log(likelihood) to measure the digitized charge given the predicted mean number of photons
       */
      Double_t GetMinus2LnL( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood of the digitizer returning the measured
       * charge given the predicted mean number of photons
       * @param undigi Predicted mean number of photons as the PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return Likelihood to measure the digitized P.E. given the predicted mean number of photons
       */
      Double_t GetLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Get the most likely number of digitized P.E. measured by the PMT
       * given the predicted mean number of photons hitting it
       * @param undigi Predicted mean number of photons hitting the PMT
       * @return Most likely number of digitized P.E. measured by the PMT
       */
      Double_t GetExpectation( const Double_t & undigi );

    protected:
    
    private:

      /**
       * Calculate -2log(likelihood) for measured charge, assuming a
       * Poisson distribution for the digitizer
       * @param undigi Predicted mean number of photons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return -2log(likelihood) of getting the recorded P.E.
       */
      Double_t GetSimpleMinus2LnL( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for measured charge, assuming a
       * Poisson distribution for the digitizer
       * @param undigi Predicted mean number of photons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return likelihood of getting the recorded P.E.
       */
      Double_t GetSimpleLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Mean digitized P.E. returned by the PMT for a given predicted mean
       * number of photons, assuming a Poisson distribution (i.e. same as the
       * number of photons in this case)
       * @param undigi Predicted mean number of photons at PMT
       * @return Most-likely resulting digitized P.E.
       */
      Double_t GetSimpleExpectation( const Double_t & undigi );
      
      /**
       * Calculate -2log(likelihood) for measured charge, assuming
       * the WCSim digitizer
       * @param undigi Predicted mean number of photons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return -2log(likelihood) of getting the recorded P.E.
       */
      Double_t GetWCSimMinus2LnL( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for measured charge, assuming
       * the WCSim digitizer
       * @param undigi Predicted mean number of photons at PMT
       * @param digi Digitized P.E. recorded by the PMT
       * @return likelihood of getting the recorded P.E.
       */
      Double_t GetWCSimLikelihood( Double_t undigi, const Double_t &digi );

      /**
       * Mean digitized P.E. returned by the PMT for a given predicted mean
       * number of photons, assuming the WCSim digitizer
       * @param undigi Predicted mean number of photons at PMT
       * @return Most-likely resulting digitized P.E.
       */
      Double_t GetWCSimExpectation( const Double_t & undigi );

      /**
       * Calculate the likelihood for the measured charge using the WCSim
       * digitzer, by sampling probability histograms (used for low charges)
       * @param undigi Predicted mean number of photons at PMT
       * @param digi Digitized P.E. recorded aby the PMT
       * @return
       */
      Double_t GetWCSimPickerLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for the measured charge using the WCSim
       * digitzer, by using a Gaussian (x) Exponential (medium P.E.)
       * @param undigi Predicted mean number of photons at PMT
       * @param digi Digitized P.E. recorded aby the PMT
       * @return
       */
      Double_t GetWCSimGausExpoLikelihood( const Double_t &undigi, const Double_t &digi );

      /**
       * Calculate the likelihood for the measured charge using the WCSim
       * digitzer, by using a simple Gaussian (high P.E.)
       * @param undigi Predicted mean number of photons at PMT
       * @param digi Digitized P.E. recorded aby the PMT
       * @return
       */
      Double_t GetWCSimGausLikelihood( const Double_t &undigi, const Double_t &digi );



      DigiType_t fType;  ///< Which digitizer method to use

      std::map<std::string, DigiType_t> fDigiTypeNames; ///< For converting strings to their corresponding DigiType_t
      
      // WCSim repeatedly samples a 1pe distribution.  For hits < 10pe I've
      // already done this to build a PDF histogram which these variables 
      // point to
      TFile * fPDFs;    ///< File holding probability density histograms for the sub-10pe WCSim digitizer
      TH2D * fDigiPDF;  ///< Probabilitiy density hitogram for the the sub-10pe WCSim digitizer
          
      // The WCSim digitizer samples the 1pe distribution repeatedly,
      // then applies a threshold function,
      // then multiplies by an efficiency term
      Double_t fEfficiency; ///< Efficiency value used in the WCSim digitizer

      Double_t fMinimum; ///< Minimum nonzero likelihood from histogram
};

#endif // WCSIMDIGITIZERLIKELIHOOD_H
