#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "TClonesArray.h"
#include "TCollection.h"
#include "TMath.h"

WCSimLikelihoodFitter::WCSimLikelihoodFitter(WCSimRootEvent * myRootEvent)
{
    fRootEvent = myRootEvent;
    fLikelihoodDigitArray = new WCSimLikelihoodDigitArray(fRootEvent);
//    fRecoDigitList = *(myRecoEvent->GetDigitList());
    std::cout << "Calculated LnL()" << std::endl;
    //ctor
}

WCSimLikelihoodFitter::~WCSimLikelihoodFitter()
{
}

void WCSimLikelihoodFitter::Minimize2LnL()
{
    this->Calc2LnL();
    // Then minimize it
    return ;
}

double WCSimLikelihoodFitter::Calc2LnL()  
{
    return (this->Charge2LnL() + this->Time2LnL());
}

double WCSimLikelihoodFitter::Charge2LnL( )
{
    
    double Charge2LnL = 1.0;

    fTrack = new WCSimLikelihoodTrack(1,1,1,1,1,1,1);
    std::cout << "Calculating charge LnL" << std::endl;
    WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood( fLikelihoodDigitArray, fTrack );
    Charge2LnL = myChargeLikelihood->Calc2LnL();
  //        double myProb = TMath::Poisson( Q, mu );
  //        ChargeLnL += TMath::Log( myProb );
      //}
   // }
    std::cout << "Returning  -2 * charge LnL" << std::endl;
    return Charge2LnL;
}

double WCSimLikelihoodFitter::CalcMuDirect()
{   
    double muDirect = 0.0;
    return muDirect;
}

double WCSimLikelihoodFitter::CalcMuIndirect()
{
    double muIndirect = 1.0;
    return muIndirect;
}

double WCSimLikelihoodFitter::Time2LnL( )
{
    std::cout << "Returning time LnL" << std::endl;
    return 1.0;
}
