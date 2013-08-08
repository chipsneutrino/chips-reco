#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimChargeLikelihood.hh"
#include "WCSimTimeLikelihood.hh"
#include "TClonesArray.h"
#include "TCollection.h"
#include "TMath.h"

WCSimLikelihoodFitter::WCSimLikelihoodFitter(WCSimRootEvent * myRootEvent)
{
    fRootEvent = myRootEvent;
    fLikelihoodDigitArray = new WCSimLikelihoodDigitArray(fRootEvent);
//    fRecoDigitList = *(myRecoEvent->GetDigitList());
    //std::cout << "Calculated LnL()" << std::endl;
    //ctor
}

WCSimLikelihoodFitter::~WCSimLikelihoodFitter()
{
}

void WCSimLikelihoodFitter::Minimize2LnL()
{
    //TODO: should set up some WCSimLikelihoodTrack!
    this->Calc2LnL();
    // Then minimize it
    return ;
}

double WCSimLikelihoodFitter::Calc2LnL()  
{
  //TODO: check if there is a defined track?  //mp
  //  like: if(!fTrack) ShoutOnUser();
  return (this->Charge2LnL() + this->Time2LnL());
}

double WCSimLikelihoodFitter::Calc2LnL(WCSimLikelihoodTrack * myTrack)
{
  return (this->Charge2LnL( myTrack ) + this->Time2LnL());
}

double WCSimLikelihoodFitter::Charge2LnL( WCSimLikelihoodTrack * myTrack )
{
  fTrack = myTrack;
  return (this->Charge2LnL());
}


double WCSimLikelihoodFitter::Charge2LnL( )
{
    if(fTrack == NULL)
    {
      std::cerr << "WCSimLikelihoodFitter::Charge2LnL() - Error: fTrack = NULL" << std::endl;
      exit(EXIT_FAILURE);
    }
    double Charge2LnL = 1.0;

    //std::cout << "Calculating charge LnL" << std::endl;

    WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood( fLikelihoodDigitArray, fTrack );
    Charge2LnL = myChargeLikelihood->Calc2LnL();
    delete myChargeLikelihood;
  //        double myProb = TMath::Poisson( Q, mu );
  //        ChargeLnL += TMath::Log( myProb );
      //}
   // }
    //std::cout << "Returning  -2 * charge LnL" << std::endl;
    std::cout << "-2 * ChargeLnL = " << Charge2LnL << std::endl;
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
    WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood( fLikelihoodDigitArray, fTrack );
    WCSimTimeLikelihood * myTimeLikelihood = new WCSimTimeLikelihood( fLikelihoodDigitArray, fTrack, myChargeLikelihood );
    Double_t time2LnL = myTimeLikelihood->Calc2LnL();
    delete myTimeLikelihood;
    delete myChargeLikelihood;

    //std::cout << "Returning -2* time LnL" << std::endl;
    std::cout << "-2 * Time2LnL = " << time2LnL << std::endl;
    return time2LnL;
}
