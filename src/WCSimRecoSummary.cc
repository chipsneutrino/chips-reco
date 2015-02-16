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
  fVertex = ts.GetVertex();

  fPrimaryPDGs = ts.GetPrimaryPDGs();
  fPrimaryEnergies = ts.GetPrimaryEnergies();
  fPrimaryDirs = ts.GetPrimaryDirs();
}

// Destructor
WCSimRecoSummary::~WCSimRecoSummary(){

}

void WCSimRecoSummary::ResetValues(){
  fVertex = TVector3(-999.,-999.,-999.);

  fPrimaryPDGs.clear();
  fPrimaryEnergies.clear();
  fPrimaryDirs.clear();
}

// Get and set the vertex information
TVector3 WCSimRecoSummary::GetVertex() const{
  return fVertex;
}

void WCSimRecoSummary::SetVertex(TVector3 vtx){
  fVertex = vtx;
}

void WCSimRecoSummary::SetVertex(double x, double y, double z){
  fVertex = TVector3(x,y,z);
}

// Get the vertex components
double WCSimRecoSummary::GetVertexX() const{
  return fVertex.X();
}

double WCSimRecoSummary::GetVertexY() const{
  return fVertex.Y();
}

double WCSimRecoSummary::GetVertexZ() const{
  return fVertex.Z();
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

