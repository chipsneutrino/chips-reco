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

      enum DigiType_t
      {
        kSimple,
        kWCSim,
        kUnknown
      };

      WCSimDigitizerLikelihood();
      virtual ~WCSimDigitizerLikelihood();

      DigiType_t StringToDigiType( const std::string &str ) const;

      void SetDigiType( WCSimDigitizerLikelihood::DigiType_t type );
      void OpenPDFs();
      Double_t GetMinus2LnL( const Double_t &undigi, const Double_t &digi );
      Double_t GetLikelihood( const Double_t &undigi, const Double_t &digi );
      Double_t GetExpectation( const Double_t & undigi );
      DigiType_t StrToDigiType(std::string str);
    protected:
    
    private:
      //  Private functions
      ///////////////////////////////////////////////////////////////////////////////////

      Double_t GetSimpleMinus2LnL( const Double_t &undigi, const Double_t &digi );
      Double_t GetSimpleLikelihood( const Double_t &undigi, const Double_t &digi );
      Double_t GetSimpleExpectation( const Double_t & undigi );
      
      Double_t GetWCSimMinus2LnL( const Double_t &undigi, const Double_t &digi );
      Double_t GetWCSimLikelihood( const Double_t &undigi, const Double_t &digi );
      Double_t GetWCSimExpectation( const Double_t & undigi );

      Double_t GetWCSimPickerLikelihood( const Double_t &undigi, const Double_t &digi );
      Double_t GetWCSimGausExpoLikelihood( const Double_t &undigi, const Double_t &digi );
      Double_t GetWCSimGausLikelihood( const Double_t &undigi, const Double_t &digi );


      // Member variables
      ///////////////////////////////////////////////////////////////////////////////////
      
      // Do we use the WCSim-style digitizer, or a simple Poisson
      // Can add something based on real PMT tests later
      DigiType_t fType;

      // For converting strings to their corresponding DigiType_t
      std::map<std::string, DigiType_t> fDigiTypeNames; 
      
      // WCSim repeatedly samples a 1pe distribution.  For hits < 10pe I've
      // already done this to build a PDF histogram which these variables 
      // point to
      TFile * fPDFs;
      TH2D * fDigiPDF;
          
      // The WCSim digitizer samples the 1pe distribution repeatedly,
      // then applies a threshold function,
      // then multiplies by an efficiency term
      Double_t fEfficiency;
};

#endif // WCSIMDIGITIZERLIKELIHOOD_H
