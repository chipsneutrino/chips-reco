#include "WCSimDataCleaner.hh"

#include "WCSimRecoDigit.hh"
#include "WCSimRecoCluster.hh"
#include "WCSimRecoClusterDigit.hh"

#include <cmath>
#include <iostream>
#include <cassert>

ClassImp(WCSimDataCleaner)

static WCSimDataCleaner* fgDataCleaner = 0;

WCSimDataCleaner* WCSimDataCleaner::Instance()
{
  if( !fgDataCleaner ){
    fgDataCleaner = new WCSimDataCleaner();
  }

  if( !fgDataCleaner ){
    assert(fgDataCleaner);
  }

  if( fgDataCleaner ){

  }

  return fgDataCleaner;
}
  
void WCSimDataCleaner::Config(Int_t config)
{
  WCSimDataCleaner::Instance()->SetConfig(config);
}

void WCSimDataCleaner::MinPulseHeight(Double_t min)
{
  WCSimDataCleaner::Instance()->SetMinPulseHeight(min);
}
  
void WCSimDataCleaner::NeighbourRadius(Double_t radius)
{
  WCSimDataCleaner::Instance()->SetNeighbourRadius(radius);
}
  
void WCSimDataCleaner::NeighbourDigits(Int_t digits)
{
  WCSimDataCleaner::Instance()->SetNeighbourDigits(digits);
}
  
void WCSimDataCleaner::ClusterRadius(Double_t radius)
{
  WCSimDataCleaner::Instance()->SetClusterRadius(radius);
}
  
void WCSimDataCleaner::ClusterDigits(Int_t digits)
{
  WCSimDataCleaner::Instance()->SetClusterDigits(digits);
}
  
void WCSimDataCleaner::TimeWindow(Double_t window)
{
  WCSimDataCleaner::Instance()->SetTimeWindow(window);
}

void WCSimDataCleaner::PrintParameters()
{
  WCSimDataCleaner::Instance()->RunPrintParameters();
}

WCSimDataCleaner::WCSimDataCleaner()
{
  // cleaning mode
  fConfig = WCSimDataCleaner::kPulseHeightAndClusters;

  // default cleaning parameters
  fMinPulseHeight = 1.0;     // minimum pulse height (PEs)
  fNeighbourRadius = 200.0;  // clustering window (cm)
  fMinNeighbourDigits = 2;   // minimum neighbouring digits
  fClusterRadius = 200.0;    // clustering window (cm)
  fMinClusterDigits = 50;    // minimum clustered digits
  fTimeWindow = 25.0;        // timing window (ns)

  // vector of filtered digits
  fFilterAll = new std::vector<WCSimRecoDigit*>;
  fFilterByPulseHeight = new std::vector<WCSimRecoDigit*>;
  fFilterByNeighbours = new std::vector<WCSimRecoDigit*>;
  fFilterByClusters = new std::vector<WCSimRecoDigit*>;
  
  // vector of clusters
  fClusterList = new std::vector<WCSimRecoCluster*>;
}

WCSimDataCleaner::~WCSimDataCleaner()
{
  delete fFilterAll;
  delete fFilterByPulseHeight;
  delete fFilterByNeighbours;
  delete fFilterByClusters;
  delete fClusterList;
}

void WCSimDataCleaner::RunPrintParameters()
{
  std::cout << " *** WCSimDataCleaner::PrintParameters() *** " << std::endl;
  
  std::cout << "  Data Cleaner Parameters: " << std::endl
            << "   Config = " << fConfig<< std::endl
            << "   MinPulseHeight = " << fMinPulseHeight << std::endl
            << "   NeighbourRadius = " << fNeighbourRadius << std::endl
            << "   MinNeighbourDigits = " << fMinNeighbourDigits << std::endl
            << "   ClusterRadius = " << fClusterRadius << std::endl
            << "   MinClusterDigits = " << fMinClusterDigits << std::endl
            << "   TimeWindow = " << fTimeWindow << std::endl;

  return;
}

void WCSimDataCleaner::Reset()
{

  return;
}

std::vector<WCSimRecoDigit*>* WCSimDataCleaner::Run(std::vector<WCSimRecoDigit*>* myDigitList)
{
  std::cout << " *** WCSimDataCleaner::Run(...) *** " << std::endl;

  // input digit list
  // ================
  std::vector<WCSimRecoDigit*>* myInputList = myDigitList;
  std::vector<WCSimRecoDigit*>* myOutputList = myDigitList;

  // filter all digits
  // =================
  myInputList = ResetDigits(myOutputList);
  myOutputList = (std::vector<WCSimRecoDigit*>*)(this->FilterAll(myInputList));
  myOutputList = FilterDigits(myOutputList);
  if( fConfig==WCSimDataCleaner::kNone ) return myOutputList;

  // filter by pulse height
  // ======================
  myInputList = ResetDigits(myOutputList);
  myOutputList = (std::vector<WCSimRecoDigit*>*)(this->FilterByPulseHeight(myInputList));
  myOutputList = FilterDigits(myOutputList);
  if( fConfig==WCSimDataCleaner::kPulseHeight ) return myOutputList;

  // filter using neighbouring digits
  // ================================
  myInputList = ResetDigits(myOutputList);
  myOutputList = (std::vector<WCSimRecoDigit*>*)(this->FilterByNeighbours(myInputList));
  myOutputList = FilterDigits(myOutputList);
  if( fConfig==WCSimDataCleaner::kPulseHeightAndNeighbours ) return myOutputList;

  // filter using clustered digits
  // =============================
  myInputList = ResetDigits(myOutputList);
  myOutputList = (std::vector<WCSimRecoDigit*>*)(this->FilterByClusters(myInputList));
  myOutputList = FilterDigits(myOutputList);
  if( fConfig==WCSimDataCleaner::kPulseHeightAndClusters ) return myOutputList;

  // return vector of filtered digits
  // ================================
  return myOutputList;
}

std::vector<WCSimRecoDigit*>* WCSimDataCleaner::ResetDigits(std::vector<WCSimRecoDigit*>* myDigitList)
{
  for( UInt_t idigit=0; idigit<myDigitList->size(); idigit++ ){
    WCSimRecoDigit* recoDigit = (WCSimRecoDigit*)(myDigitList->at(idigit));
    recoDigit->ResetFilter();
  }

  return myDigitList;
}

std::vector<WCSimRecoDigit*>* WCSimDataCleaner::FilterDigits(std::vector<WCSimRecoDigit*>* myDigitList)
{
  for( UInt_t idigit=0; idigit<myDigitList->size(); idigit++ ){
    WCSimRecoDigit* recoDigit = (WCSimRecoDigit*)(myDigitList->at(idigit));
    recoDigit->PassFilter();
  }

  return myDigitList;
}

std::vector<WCSimRecoDigit*>* WCSimDataCleaner::FilterAll(std::vector<WCSimRecoDigit*>* myDigitList)
{
  // clear vector of filtered digits
  // ==============================
  fFilterAll->clear();

  // filter all digits
  // =================
  for( UInt_t idigit=0; idigit<myDigitList->size(); idigit++ ){
    WCSimRecoDigit* recoDigit = (WCSimRecoDigit*)(myDigitList->at(idigit));
    fFilterAll->push_back(recoDigit);
  }

  // return vector of filtered digits
  // ================================
  std::cout << "  filter all: " << fFilterAll->size() << std::endl;
  
  return fFilterAll;
}

std::vector<WCSimRecoDigit*>* WCSimDataCleaner::FilterByPulseHeight(std::vector<WCSimRecoDigit*>* myDigitList)
{
  // clear vector of filtered digits
  // ===============================
  fFilterByPulseHeight->clear();

  // filter by pulse height
  // ======================
  for( UInt_t idigit=0; idigit<myDigitList->size(); idigit++ ){
    WCSimRecoDigit* recoDigit = (WCSimRecoDigit*)(myDigitList->at(idigit));
    if( recoDigit->GetQPEs()>fMinPulseHeight ){
      fFilterByPulseHeight->push_back(recoDigit);
    }
  }

  // return vector of filtered digits
  // ================================
  std::cout << "  filter by pulse height: " << fFilterByPulseHeight->size() << std::endl;
  
  return fFilterByPulseHeight;
}

std::vector<WCSimRecoDigit*>* WCSimDataCleaner::FilterByNeighbours(std::vector<WCSimRecoDigit*>* myDigitList)
{
  // clear vector of filtered digits
  // ===============================
  fFilterByNeighbours->clear();

  // create array of neighbours
  // ==========================
  Int_t Ndigits = myDigitList->size();

  if( Ndigits<=0 ){
    return fFilterByNeighbours;
  }

  Int_t* numNeighbours = new Int_t[Ndigits];

  for( Int_t idigit=0; idigit<Ndigits; idigit++ ){
    numNeighbours[idigit] = 0;
  }

  // count number of neighbours
  // ==========================
  for( UInt_t idigit1=0; idigit1<myDigitList->size(); idigit1++ ){
    for( UInt_t idigit2=idigit1+1; idigit2<myDigitList->size(); idigit2++ ){
      WCSimRecoDigit* fdigit1 = (WCSimRecoDigit*)(myDigitList->at(idigit1));
      WCSimRecoDigit* fdigit2 = (WCSimRecoDigit*)(myDigitList->at(idigit2));

      Double_t dx = fdigit1->GetX() - fdigit2->GetX();
      Double_t dy = fdigit1->GetY() - fdigit2->GetY();
      Double_t dz = fdigit1->GetZ() - fdigit2->GetZ();
      Double_t dt = fdigit1->GetTime() - fdigit2->GetTime();
      Double_t drsq = dx*dx + dy*dy + dz*dz;

      if( drsq>0.0
       && drsq<fNeighbourRadius*fNeighbourRadius
       && fabs(dt)<fTimeWindow ){
        numNeighbours[idigit1]++;
        numNeighbours[idigit2]++;
      }
    }
  }

  // filter by number of neighbours
  // ==============================
  for( UInt_t idigit=0; idigit<myDigitList->size(); idigit++ ){
    WCSimRecoDigit* fdigit = (WCSimRecoDigit*)(myDigitList->at(idigit));
    if( numNeighbours[idigit]>=fMinNeighbourDigits ){
      fFilterByNeighbours->push_back(fdigit);
    }
  }

  // delete array of neighbours
  // ==========================
  delete [] numNeighbours;

  // return vector of filtered digits
  // ================================
  std::cout << "  filter by neighbours: " << fFilterByNeighbours->size() << std::endl;
  
  return fFilterByNeighbours;
}

std::vector<WCSimRecoDigit*>* WCSimDataCleaner::FilterByClusters(std::vector<WCSimRecoDigit*>* myDigitList)
{
  // clear vector of filtered digits
  // ===============================
  fFilterByClusters->clear();

  // run clustering algorithm
  // ========================
  std::vector<WCSimRecoCluster*>* myClusterList = (std::vector<WCSimRecoCluster*>*)(this->RecoClusters(myDigitList));

  for( UInt_t icluster=0; icluster<myClusterList->size(); icluster++ ){
    WCSimRecoCluster* myCluster = (WCSimRecoCluster*)(myClusterList->at(icluster));
    for( Int_t idigit=0; idigit<myCluster->GetNDigits(); idigit++ ){
      WCSimRecoDigit* myDigit = (WCSimRecoDigit*)(myCluster->GetDigit(idigit));
      fFilterByClusters->push_back(myDigit);
    }
  }
  
  // return vector of filtered digits
  // ================================
  std::cout << "  filter by clusters: " << fFilterByClusters->size() << std::endl;
  
  return fFilterByClusters;
}

std::vector<WCSimRecoCluster*>* WCSimDataCleaner::RecoClusters(std::vector<WCSimRecoDigit*>* myDigitList)
{  

  // delete cluster digits
  // =====================
  for( UInt_t i=0; i<vClusterDigitList.size(); i++ ){
    delete (WCSimRecoClusterDigit*)(vClusterDigitList.at(i));
  }
  vClusterDigitList.clear();

  // delete clusters
  // ===============
  for( UInt_t i=0; i<vClusterList.size(); i++ ){
    delete (WCSimRecoCluster*)(vClusterList.at(i));
  }
  vClusterList.clear();

  // clear vector clusters
  // =====================
  fClusterList->clear();

  // make cluster digits
  // ===================
  for( UInt_t idigit=0; idigit<myDigitList->size(); idigit++ ){
    WCSimRecoDigit* recoDigit = (WCSimRecoDigit*)(myDigitList->at(idigit));
    WCSimRecoClusterDigit* clusterDigit = new WCSimRecoClusterDigit(recoDigit);
    vClusterDigitList.push_back(clusterDigit);
  }

  // run clustering algorithm
  // ========================
  for( UInt_t idigit1=0; idigit1<vClusterDigitList.size(); idigit1++ ){
    for( UInt_t idigit2=idigit1+1; idigit2<vClusterDigitList.size(); idigit2++ ){

      WCSimRecoClusterDigit* fdigit1 = (WCSimRecoClusterDigit*)(vClusterDigitList.at(idigit1));
      WCSimRecoClusterDigit* fdigit2 = (WCSimRecoClusterDigit*)(vClusterDigitList.at(idigit2));

      Double_t dx = fdigit1->GetX() - fdigit2->GetX();
      Double_t dy = fdigit1->GetY() - fdigit2->GetY();
      Double_t dz = fdigit1->GetZ() - fdigit2->GetZ();
      Double_t dt = fdigit1->GetTime() - fdigit2->GetTime();
      Double_t drsq = dx*dx + dy*dy + dz*dz;

      if( drsq>0.0
       && drsq<fClusterRadius*fClusterRadius
       && fabs(dt)<fTimeWindow ){
        fdigit1->AddClusterDigit(fdigit2);
        fdigit2->AddClusterDigit(fdigit1);
      }
    }
  }

  // collect up clusters
  // ===================
  Bool_t carryon = 0;

  for( UInt_t idigit=0; idigit<vClusterDigitList.size(); idigit++ ){
    WCSimRecoClusterDigit* fdigit = (WCSimRecoClusterDigit*)(vClusterDigitList.at(idigit));

    if( fdigit->IsClustered()==0
     && fdigit->GetNClusterDigits()>0 ){
        
      vClusterDigitCollection.clear();
      vClusterDigitCollection.push_back(fdigit);
      fdigit->SetClustered();

      carryon = 1;
      while( carryon ){
        carryon = 0;

        for( UInt_t jdigit=0; jdigit<vClusterDigitCollection.size(); jdigit++ ){
          WCSimRecoClusterDigit* cdigit = (WCSimRecoClusterDigit*)(vClusterDigitCollection.at(jdigit));

          if( cdigit->IsAllClustered()==0 ){
            for( Int_t kdigit=0; kdigit<cdigit->GetNClusterDigits(); kdigit++ ){
              WCSimRecoClusterDigit* cdigitnew = (WCSimRecoClusterDigit*)(cdigit->GetClusterDigit(kdigit));
                
              if( cdigitnew->IsClustered()==0 ){
                vClusterDigitCollection.push_back(cdigitnew);
                cdigitnew->SetClustered();
                carryon = 1;
              }
            }
          }
        }
      }
    
      if( (Int_t)vClusterDigitCollection.size()>=fMinClusterDigits ){
        WCSimRecoCluster* cluster = new WCSimRecoCluster();
        fClusterList->push_back(cluster);
        vClusterList.push_back(cluster);

        for( UInt_t jdigit=0; jdigit<vClusterDigitCollection.size(); jdigit++ ){
          WCSimRecoClusterDigit* cdigit = (WCSimRecoClusterDigit*)(vClusterDigitCollection.at(jdigit));
          WCSimRecoDigit* recodigit = (WCSimRecoDigit*)(cdigit->GetRecoDigit());
          cluster->AddDigit(recodigit);        
        }
      }
    }
  }

  // return vector of clusters
  // =========================
  return fClusterList;
}


