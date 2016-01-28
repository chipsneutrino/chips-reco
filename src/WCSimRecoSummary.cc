/*
 * WCSimRecoSummary.cc
 *
 *  Created on: 16 Feb 2015
 *      Author: ajperch
 */

#include "WCSimRecoSummary.hh"

#include <iostream>
#include <vector>

#include <TVector3.h>
#include <TLorentzVector.h>

#include "WCSimRecoSummary.hh"

ClassImp(WCSimRecoSummary)

// Standard constructor
WCSimRecoSummary::WCSimRecoSummary() : TObject() {
  this->ResetValues();
}

// Copy constructor
WCSimRecoSummary::WCSimRecoSummary(const WCSimRecoSummary &ts) : TObject(ts) {
  fVertices = ts.GetVertices();
  fVerticesT = ts.GetVerticesT();
  fPrimaryPDGs = ts.GetPrimaryPDGs();
  fPrimaryEnergies = ts.GetPrimaryEnergies();
  fPrimaryDirs = ts.GetPrimaryDirs();
}

// Destructor
WCSimRecoSummary::~WCSimRecoSummary(){

}

void WCSimRecoSummary::ResetValues(){
  fVertices.clear();
  fVerticesT.clear();

  fPrimaryPDGs.clear();
  fPrimaryEnergies.clear();
  fPrimaryDirs.clear();
}

// Get and set the vertex information
void WCSimRecoSummary::AddVertex(TVector3 vtx, double t){
  fVertices.push_back(vtx);
  fVerticesT.push_back(t);
}

void WCSimRecoSummary::AddVertex(double x, double y, double z, double t){
  fVertices.push_back(TVector3(x,y,z));
  fVerticesT.push_back(t);
}

TVector3 WCSimRecoSummary::GetVertex(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fVertices[p];
  }
  else{
    std::cerr << "== Request for vertex " << p << " of [0..." << this->GetNPrimaries()-1 << "]" << std::endl;
    return -999;
  }
}

std::vector<TVector3> WCSimRecoSummary::GetVertices() const{
  return fVertices;
}

double WCSimRecoSummary::GetVertexX(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fVertices[p].X();
  }
  else{
    std::cerr << "== Request for vertex x " << p << " of [0..." << this->GetNPrimaries()-1 << "]" << std::endl;
    return -999;
  }
}

double WCSimRecoSummary::GetVertexY(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fVertices[p].Y();
  }
  else{
    std::cerr << "== Request for vertex y " << p << " of [0..." << this->GetNPrimaries()-1 << "]" << std::endl;
    return -999;
  }
}

double WCSimRecoSummary::GetVertexZ(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fVertices[p].Z();
  }
  else{
    std::cerr << "== Request for vertex z " << p << " of [0..." << this->GetNPrimaries()-1 << "]" << std::endl;
    return -999;
  }
}

double WCSimRecoSummary::GetVertexT(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fVerticesT[p];
  }
  else{
    std::cerr << "== Request for vertex t " << p << " of [0..." << this->GetNPrimaries()-1 << "]" << std::endl;
    return -999;
  }
}

std::vector<double> WCSimRecoSummary::GetVerticesT() const{
  return fVerticesT;
}

// Primary particle functions
void WCSimRecoSummary::AddPrimary(int pdg, double en, TVector3 dir){
  fPrimaryPDGs.push_back(pdg);
  fPrimaryEnergies.push_back(en);
  fPrimaryDirs.push_back(dir);
}

void WCSimRecoSummary::AddPrimary(int pdg, double en, double dx, double dy, double dz){
  fPrimaryPDGs.push_back(pdg);
  fPrimaryEnergies.push_back(en);
  TVector3 tmpDir(dx,dy,dz);
  fPrimaryDirs.push_back(tmpDir);
}

int WCSimRecoSummary::GetPrimaryPDG(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fPrimaryPDGs[p];
  }
  else{
    std::cerr << "== Request for primary particle " << p << " of [0..." << this->GetNPrimaries()-1 << "]" << std::endl;
    return -999;
  }
}

double WCSimRecoSummary::GetPrimaryEnergy(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fPrimaryEnergies[p];
  }
  else{
    std::cerr << "== Request for primary particle " << p << " of [0..." << this->GetNPrimaries()-1 <<  "]" << std::endl;
    return -999.;
  }

}

TVector3 WCSimRecoSummary::GetPrimaryDir(unsigned int p) const{
  if(p < this->GetNPrimaries()){
    return fPrimaryDirs[p];
  }
  else{
    std::cerr << "== Request for primary particle " << p << " of [0..." << this->GetNPrimaries()-1 <<  "]" << std::endl;
    return TVector3(-999.,-999.,-999.);
  }

}

std::vector<int> WCSimRecoSummary::GetPrimaryPDGs() const{
  return fPrimaryPDGs;
}

std::vector<double> WCSimRecoSummary::GetPrimaryEnergies() const{
  return fPrimaryEnergies;
}

std::vector<TVector3> WCSimRecoSummary::GetPrimaryDirs() const{
  return fPrimaryDirs;
}

unsigned int WCSimRecoSummary::GetNPrimaries() const{
  return fPrimaryPDGs.size();
}

bool WCSimRecoSummary::HasCommonVertex() const{
  bool commonVertex = true;
  // Always true for one track. For more, let's find out.
  if(this->GetNPrimaries() > 1){
    for(unsigned int p = 1; p < this->GetNPrimaries(); ++p){
      if(!(fVertices[p] == fVertices[p-1])){
        commonVertex = false;
        break;
      }
    }
  }
  return commonVertex;  
}


