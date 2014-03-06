#ifndef WCSIMDIGITIZERLIKELIHOOD_H
#define WCSIMDIGITIZERLIKELIHOOD_H
#include "Rtypes.h"

class TFile;
class TH2D;

class WCSimDigitizerLikelihood
{
    public:

      enum DigiType
      {
        kSimple,
        kWCSim
      };

      WCSimDigitizerLikelihood( DigiType type = kWCSim);
      virtual ~WCSimDigitizerLikelihood();

      void SetDigiType( WCSimDigitizerLikelihood::DigiType type );
      void OpenPDFs();
      Double_t GetMinus2LnL( const Double_t &undigi, const Double_t &digi );
      Double_t GetLikelihood( const Double_t &undigi, const Double_t &digi );
      Double_t GetExpectation( const Double_t & undigi );

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
      DigiType fType;
	
      
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
