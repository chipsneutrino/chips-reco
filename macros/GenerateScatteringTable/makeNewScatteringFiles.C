// C++ 
#include <algorithm>
#include <string>
#include <sstream>

// Root
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TObject.h>
#include <TVector2.h>
#include <TVector3.h>

// WCSim
#include "/home/lwhitehead/work/CHIPS/recon/test/WCSim/include/WCSimRootEvent.hh"
#include "/home/lwhitehead/work/CHIPS/recon/test/WCSim/include/WCSimRootGeom.hh"
#include "/home/lwhitehead/work/CHIPS/recon/test/WCSim/include/WCSimTruthSummary.hh"

void makeNewScatteringFiles(std::string inputs, std::string output){

  TChain* data = new TChain("wcsimT");
  TChain* geom = new TChain("wcsimGeoT");

  data->Add(inputs.c_str());
  geom->Add(inputs.c_str());

  // Link an event to the chain
  WCSimRootEvent* event = new WCSimRootEvent();
  data->SetBranchAddress("wcsimrootevent",&event);

  // Link a geometry instance to the chain. Get the 0th entry,
  // as this is the only one that exists.
  WCSimRootGeom* geoInst = new WCSimRootGeom();
  geom->SetBranchAddress("wcsimrootgeom",&geoInst);
  geom->GetEntry(0);

  // Set up the output tree.
  // We have six coordinates plus the charge to store
  double vtxR; // Vertex radial position
  double vtxZ; // Vertex Z position
  double v2pZeta; // 2D (disc) angle between vtx and pmt
  double pmtPos; // PMT Z position for barrel, or R position for caps
  double dirTheta; // Source theta direction
  double dirPhi; // Source phi direction
  // Also want the charge
  double charge;
  // Need the pmt location, too.
  int pmtLoc;

  TTree *outTree = new TTree("scatteringNtp","scatteringNtp");
  outTree->Branch("vtxR",&vtxR,"vtxR/D");
  outTree->Branch("vtxZ",&vtxZ,"vtxZ/D");
  outTree->Branch("v2pZeta",&v2pZeta,"v2pZeta/D");
  outTree->Branch("pmtPos",&pmtPos,"pmtPos/D");
  outTree->Branch("dirTheta",&dirTheta,"dirTheta/D");
  outTree->Branch("dirPhi",&dirPhi,"dirPhi/D");
  outTree->Branch("charge",&charge,"charge/D");
  outTree->Branch("pmtLoc",&pmtLoc,"pmtLoc/I");

  // Loop over the events
  for(int i = 0; i < data->GetEntries(); ++i){

    data->GetEntry(i);

    WCSimRootTrigger* trig = event->GetTrigger(0);
    WCSimTruthSummary truthSummary = event->GetTruthSummary();

    // Event level quantities
    TVector3 vtxPos = truthSummary.GetVertex(); // Vertex position
    TVector3 vtxPos2D(vtxPos.X(),vtxPos.Y(),0); // Vertex position in the circular plane
    TVector3 beamDir = truthSummary.GetBeamDir(); // Beam direction

    // Calculate the three vertex coordinates
    vtxR = sqrt(vtxPos.X()*vtxPos.X() + vtxPos.Y()*vtxPos.Y());
    vtxZ = vtxPos.Z();

    // Calculate the two directions
    dirTheta = TMath::ACos(beamDir.Z());
    dirPhi = TMath::ATan2(beamDir.Y(),beamDir.X());
   
    int nHits = trig->GetNcherenkovdigihits();
  
    for(int h = 0; h < nHits; ++h){
      TObject* element = (trig->GetCherenkovDigiHits())->At(h);
      WCSimRootCherenkovDigiHit* hit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);

      // Charge is the easy bit!
      charge = hit->GetQ();

      WCSimRootPMT pmt = geoInst->GetPMT(hit->GetTubeId());
      // PMT position
      TVector3 pmtVec(pmt.GetPosition(0)*10.,pmt.GetPosition(1)*10.,pmt.GetPosition(2)*10.);
      TVector3 pmtVec2D(pmtVec.X(),pmtVec.Y(),0.0);
  
      // Angle in 2D (disc) between vtx and pmt.
      v2pZeta = vtxPos2D.Angle(pmtVec2D);

      // PMT coordinate (R for caps, Z for barrel).
      pmtLoc = pmt.GetCylLoc();
      if(pmtLoc == 0 || pmtLoc == 2){
        pmtPos = sqrt(pmtVec.X()*pmtVec.X() + pmtVec.Y()*pmtVec.Y());
      }
      else{
        pmtPos = pmtVec.Z();
      }

      outTree->Fill();
    }

  }

  TFile outFile(output.c_str(),"RECREATE");
  outTree->Write();
  outFile.Close();
}


