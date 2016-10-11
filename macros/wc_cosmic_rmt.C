#include <iostream>
#include <string>
#include <sstream>
#include <utility>

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>

void MakeNicePlots(TH1D* h);

void wc_cosmic_rmt(std::string fileName, int run){

  gStyle->SetOptStat(11);

  TChain *cTrue = new TChain("fTruthTree");
  TChain *cReco = new TChain("fFitTree");
  
  cTrue->Add(fileName.c_str());
  cReco->Add(fileName.c_str());

	Int_t fTrueEvent;
  Double_t fTrueVtxX;
	Double_t fTrueVtxY;
	Double_t fTrueVtxZ;
	Double_t fTrueVtxT;
	Double_t fTrueDirTheta;
	Double_t fTrueDirPhi;
	Double_t fTrueEnergy;
	Int_t fTruePDG;
	Int_t fTrueEscape;

	Int_t fRecoEvent;
	Double_t fRecoVtxX;
	Double_t fRecoVtxY;
	Double_t fRecoVtxZ;
	Double_t fRecoVtxT;
	Double_t fRecoDirTheta;
	Double_t fRecoDirPhi;
	Double_t fRecoEnergy;
	Int_t fRecoPDG;

  cTrue->SetBranchAddress("event", &fTrueEvent);
	cTrue->SetBranchAddress("vtxX", &fTrueVtxX);
	cTrue->SetBranchAddress("vtxY", &fTrueVtxY);
	cTrue->SetBranchAddress("vtxZ", &fTrueVtxZ);
	cTrue->SetBranchAddress("vtxT", &fTrueVtxT);
	cTrue->SetBranchAddress("theta", &fTrueDirTheta);
	cTrue->SetBranchAddress("phi", &fTrueDirPhi);
	cTrue->SetBranchAddress("energy", &fTrueEnergy);
	cTrue->SetBranchAddress("pdg", &fTruePDG);
	cTrue->SetBranchAddress("escapes", &fTrueEscape);
  
	cReco->SetBranchAddress("event", &fRecoEvent);
	cReco->SetBranchAddress("vtxX", &fRecoVtxX);
	cReco->SetBranchAddress("vtxY", &fRecoVtxY);
	cReco->SetBranchAddress("vtxZ", &fRecoVtxZ);
	cReco->SetBranchAddress("vtxT", &fRecoVtxT);
	cReco->SetBranchAddress("theta", &fRecoDirTheta);
	cReco->SetBranchAddress("phi", &fRecoDirPhi);
	cReco->SetBranchAddress("energy", &fRecoEnergy);
	cReco->SetBranchAddress("pdg", &fRecoPDG);

  // Find the event limits to loop over.
  cTrue->GetEntry(0);
  int firstEvent = fTrueEvent;
  cTrue->GetEntry(cTrue->GetEntries() - 1);
  int finalEvent = fTrueEvent;

  std::cout << std::endl << "Searching for true events from " << firstEvent << " to " << finalEvent << std::endl << std::endl;

  double hMax = 2000;
  double hMin = -hMax;
  double hMaxZoom = 500;
  double hMinZoom = -hMaxZoom;

  // Define some reco minus true histograms
  TH1D* hElecVertexX = new TH1D("hElecVertexX","",50,hMin,hMax);
  TH1D* hElecVertexY = new TH1D("hElecVertexY","",50,hMin,hMax);
  TH1D* hElecVertexZ = new TH1D("hElecVertexZ","",50,hMin,hMax);
  TH1D* hElecDirTheta = new TH1D("hElecDirTheta","",50,-TMath::Pi(),TMath::Pi());
  TH1D* hElecDirPhi = new TH1D("hElecDirPhi","",50,-2*TMath::Pi(),2*TMath::Pi());
  TH1D* hElecVtxDist = new TH1D("hElecVtxDist","",50,0,hMax);
  TH1D* hElecVertexT = new TH1D("hElecVertexT","",50,-20,20);
  TH1D* hElecEnergy = new TH1D("hElecEnergy","",50,-2500,2500);
  TH1D* hElecAngle = new TH1D("hElecAngle","",50,0,45);

  TH1D* hMuonVertexX = new TH1D("hMuonVertexX","",50,hMin,hMax);
  TH1D* hMuonVertexY = new TH1D("hMuonVertexY","",50,hMin,hMax);
  TH1D* hMuonVertexZ = new TH1D("hMuonVertexZ","",50,hMin,hMax);
  TH1D* hMuonDirTheta = new TH1D("hMuonDirTheta","",50,-TMath::Pi(),TMath::Pi());
  TH1D* hMuonDirPhi = new TH1D("hMuonDirPhi","",50,-2*TMath::Pi(),2*TMath::Pi());
  TH1D* hMuonVtxDist = new TH1D("hMuonVtxDist","",50,0,hMax);
  TH1D* hMuonVertexT = new TH1D("hMuonVertexT","",50,-20,20);
  TH1D* hMuonEnergy = new TH1D("hMuonEnergy","",50,-2500,2500);
  TH1D* hMuonAngle = new TH1D("hMuonAngle","",50,0,45);

  // Zoomed plots
  TH1D* hElecVertexXZoom = new TH1D("hElecVertexXZoom","",50,hMinZoom,hMaxZoom);
  TH1D* hElecVertexYZoom = new TH1D("hElecVertexYZoom","",50,hMinZoom,hMaxZoom);
  TH1D* hElecVertexZZoom = new TH1D("hElecVertexZZoom","",50,hMinZoom,hMaxZoom);
  TH1D* hElecDirThetaZoom = new TH1D("hElecDirThetaZoom","",50,-0.25*TMath::Pi(),0.25*TMath::Pi());
  TH1D* hElecDirPhiZoom = new TH1D("hElecDirPhiZoom","",50,-0.5*TMath::Pi(),0.5*TMath::Pi());
  TH1D* hElecVtxDistZoom = new TH1D("hElecVtxDistZoom","",50,0,hMaxZoom);
  TH1D* hElecVertexTZoom = new TH1D("hElecVertexTZoom","",81,-20.125,20.125);
  TH1D* hElecEnergyZoom = new TH1D("hElecEnergyZoom","",51,-2550,2550);
  TH1D* hElecAngleZoom = new TH1D("hElecAngleZoom","",45,0,45);

  TH1D* hMuonVertexXZoom = new TH1D("hMuonVertexXZoom","",50,hMinZoom,hMaxZoom);
  TH1D* hMuonVertexYZoom = new TH1D("hMuonVertexYZoom","",50,hMinZoom,hMaxZoom);
  TH1D* hMuonVertexZZoom = new TH1D("hMuonVertexZZoom","",50,hMinZoom,hMaxZoom);
  TH1D* hMuonDirThetaZoom = new TH1D("hMuonDirThetaZoom","",50,-0.25*TMath::Pi(),0.25*TMath::Pi());
  TH1D* hMuonDirPhiZoom = new TH1D("hMuonDirPhiZoom","",50,-0.5*TMath::Pi(),0.5*TMath::Pi());
  TH1D* hMuonVtxDistZoom = new TH1D("hMuonVtxDistZoom","",50,0,hMaxZoom);
  TH1D* hMuonVertexTZoom = new TH1D("hMuonVertexTZoom","",81,-20.125,20.125);
  TH1D* hMuonEnergyZoom = new TH1D("hMuonEnergyZoom","",51,-2550,2550);
  TH1D* hMuonAngleZoom = new TH1D("hMuonAngleZoom","",45,0,25);

  // Get the first entries
  // Get the first entries
  cTrue->GetEntry(0);
  cReco->GetEntry(0);
  int currentTrueEntry = 0;
  int currentRecoEntry = 0;

  int cTotalMatched = 0;
  int cContainedMuon = 0;
  int cExitingElec = 0;
  int cElecEnergy = 0;
  int cMuonEnergy = 0;

  // Loop over the events and try to match the true and reco objects
  for(int i = firstEvent; i <= finalEvent; ++i){

    int trueElecIndex = -999;
    int trueMuonIndex = -999;
    int recoElecIndex = -999;
    int recoMuonIndex = -999;

    bool trueElecEscape = false;
    bool trueMuonEscape = false;
    double recoElecEnergy = -999.9;
    double recoMuonEnergy = -999.9;

    std::cout << " == Main loop: " << i << " :: " << fTrueEvent << ", " << fRecoEvent << std::endl;

    // Truth info
    while((fTrueEvent == i) && (currentTrueEntry < cTrue->GetEntries())){
//      std::cout << "Checking true entry " << currentTrueEntry << " which is from event " << fTrueEvent
//                << " and has pdg code " << fTruePDG << std::endl;
      if(fTruePDG == 11){
        trueElecIndex = currentTrueEntry;
        trueElecEscape = fTrueEscape;
      }
      else if(fTruePDG == 13){
        trueMuonIndex = currentTrueEntry;
        trueMuonEscape = fTrueEscape;
      }
      ++currentTrueEntry;
      cTrue->GetEntry(currentTrueEntry);
    }
    bool trueSuccess = (trueElecIndex != -999) && (trueMuonIndex != -999);

    if(trueSuccess){
//      std::cout << "True electron = " << trueElecIndex << " and true muon = " << trueMuonIndex << std::endl;
    }
    else{
//      std::cout << "This event didn't have an cosmic overlay." << std::endl;
    }

    // Reco info
    while((fRecoEvent == i) && (currentRecoEntry < cReco->GetEntries())){
//      std::cout << "Checking reco entry " << currentRecoEntry << " which is from event " << fRecoEvent
//                << " and has pdg code " << fRecoPDG << std::endl;
      if(fRecoPDG == 11){
        recoElecIndex = currentRecoEntry;
        recoElecEnergy = fRecoEnergy;
      }
      else if(fRecoPDG == 13){
        recoMuonIndex = currentRecoEntry;
        recoMuonEnergy = fRecoEnergy;
      }
      ++currentRecoEntry;
      cReco->GetEntry(currentRecoEntry);
    }
    bool recoSuccess = (recoElecIndex != -999) && (recoMuonIndex != -999);

    if(recoSuccess){
//      std::cout << "Reco electron = " << recoElecIndex << " and reco muon = " << recoMuonIndex << std::endl;
    }
    else{
//      std::cout << "This event didn't fit two tracks for some reason." << std::endl;
    }

    // If we have found true and reco tracks, time to fill some plots.
    if(trueSuccess && recoSuccess){

      std::cout << "Matched true and reco tracks - filling plots." << std::endl;
      ++cTotalMatched;

      if(!trueMuonEscape){
        std::cout << "Skipping event due to contained muon" << std::endl;
        ++cContainedMuon;
        continue;
      }
      if(trueElecEscape){
        std::cout << "Skipping event due to escaping electron" << std::endl;
        ++cExitingElec;
        continue;
      }
      if(recoElecEnergy < 550 || recoElecEnergy > 4950){
        std::cout << "Skipping event due to railed electron energy" << std::endl;
        ++cElecEnergy;
        continue;
      }
      if(recoMuonEnergy < 550){
        std::cout << "Skipping event due to low muon energy" << std::endl;
        ++cMuonEnergy;
        continue;
      }

      // Electrons first
      cTrue->GetEntry(trueElecIndex);
      cReco->GetEntry(recoElecIndex);

      hElecVertexX->Fill(fRecoVtxX-fTrueVtxX);
      hElecVertexY->Fill(fRecoVtxY-fTrueVtxY);
      hElecVertexZ->Fill(fRecoVtxZ-fTrueVtxZ);
      hElecDirTheta->Fill(fRecoDirTheta-fTrueDirTheta);
      hElecDirPhi->Fill(fRecoDirPhi-fTrueDirPhi);
      hElecVertexT->Fill(fRecoVtxT-fTrueVtxT);
      hElecEnergy->Fill(fRecoEnergy-fTrueEnergy);

      hElecVertexXZoom->Fill(fRecoVtxX-fTrueVtxX);
      hElecVertexYZoom->Fill(fRecoVtxY-fTrueVtxY);
      hElecVertexZZoom->Fill(fRecoVtxZ-fTrueVtxZ);
      hElecDirThetaZoom->Fill(fRecoDirTheta-fTrueDirTheta);
      hElecDirPhiZoom->Fill(fRecoDirPhi-fTrueDirPhi);
      hElecVertexTZoom->Fill(fRecoVtxT-fTrueVtxT);
      hElecEnergyZoom->Fill(fRecoEnergy-fTrueEnergy);

      TVector3 trueVtx(fTrueVtxX,fTrueVtxY,fTrueVtxZ);
      TVector3 recoVtx(fRecoVtxX,fRecoVtxY,fRecoVtxZ);
      hElecVtxDist->Fill((recoVtx-trueVtx).Mag());
      hElecVtxDistZoom->Fill((recoVtx-trueVtx).Mag());

      TVector3 trueDir;
      TVector3 recoDir;
      trueDir.SetMagThetaPhi(1,fTrueDirTheta,fTrueDirPhi);
      recoDir.SetMagThetaPhi(1,fRecoDirTheta,fRecoDirPhi);
      hElecAngle->Fill(recoDir.Angle(trueDir)*180.0/TMath::Pi());
      hElecAngleZoom->Fill(recoDir.Angle(trueDir)*180/TMath::Pi());

      // Now muons
      cTrue->GetEntry(trueMuonIndex);
      cReco->GetEntry(recoMuonIndex);

      hMuonVertexX->Fill(fRecoVtxX-fTrueVtxX);
      hMuonVertexY->Fill(fRecoVtxY-fTrueVtxY);
      hMuonVertexZ->Fill(fRecoVtxZ-fTrueVtxZ);
      hMuonDirTheta->Fill(fRecoDirTheta-fTrueDirTheta);
      hMuonDirPhi->Fill(fRecoDirPhi-fTrueDirPhi);
      hMuonVertexT->Fill(fRecoVtxT-fTrueVtxT);
      hMuonEnergy->Fill(fRecoEnergy-fTrueEnergy);

      hMuonVertexXZoom->Fill(fRecoVtxX-fTrueVtxX);
      hMuonVertexYZoom->Fill(fRecoVtxY-fTrueVtxY);
      hMuonVertexZZoom->Fill(fRecoVtxZ-fTrueVtxZ);
      hMuonDirThetaZoom->Fill(fRecoDirTheta-fTrueDirTheta);
      hMuonDirPhiZoom->Fill(fRecoDirPhi-fTrueDirPhi);
      hMuonVertexTZoom->Fill(fRecoVtxT-fTrueVtxT);
      hMuonEnergyZoom->Fill(fRecoEnergy-fTrueEnergy);

      trueVtx = TVector3(fTrueVtxX,fTrueVtxY,fTrueVtxZ);
      recoVtx = TVector3(fRecoVtxX,fRecoVtxY,fRecoVtxZ);
      hMuonVtxDist->Fill((recoVtx-trueVtx).Mag());
      hMuonVtxDistZoom->Fill((recoVtx-trueVtx).Mag());

      trueDir.SetMagThetaPhi(1,fTrueDirTheta,fTrueDirPhi);
      recoDir.SetMagThetaPhi(1,fRecoDirTheta,fRecoDirPhi);
      hMuonAngle->Fill(recoDir.Angle(trueDir)*180.0/TMath::Pi());
      hMuonAngleZoom->Fill(recoDir.Angle(trueDir)*180.0/TMath::Pi());

      // Reset the chains to the correct entry
      cTrue->GetEntry(currentTrueEntry);
      cReco->GetEntry(currentRecoEntry);

//      std::cout << "Plots filled, event numbers set: " << fTrueEvent << ", " << fRecoEvent << std::endl;
    }
    
  }

  MakeNicePlots(hElecVertexX);
  hElecVertexX->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Position X (cm)");
  MakeNicePlots(hElecVertexY);
  hElecVertexY->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Position Y (cm)");
  MakeNicePlots(hElecVertexZ);
  hElecVertexZ->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Position Z (cm)");
  MakeNicePlots(hElecDirTheta);
  hElecDirTheta->GetXaxis()->SetTitle("(Reco - True) Electron Direction Theta (rad)");
  MakeNicePlots(hElecDirPhi);
  hElecDirPhi->GetXaxis()->SetTitle("(Reco - True) Electron Direction Phi (rad)");
  MakeNicePlots(hElecVtxDist);
  hElecVtxDist->GetXaxis()->SetTitle("Electron |(Reco Vertex - True Vertex)| (cm)");
  MakeNicePlots(hElecVertexT);
  hElecVertexT->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Time (ns)");
  MakeNicePlots(hElecEnergy);
  hElecEnergy->GetXaxis()->SetTitle("(Reco - True) Electron Energy (MeV)");
  MakeNicePlots(hElecAngle);
  hElecAngle->GetXaxis()->SetTitle("Angle between Reco and True Directions (^o)");

  MakeNicePlots(hElecVertexXZoom);
  hElecVertexXZoom->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Position X (cm)");
  MakeNicePlots(hElecVertexYZoom);
  hElecVertexYZoom->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Position Y (cm)");
  MakeNicePlots(hElecVertexZZoom);
  hElecVertexZZoom->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Position Z (cm)");
  MakeNicePlots(hElecDirThetaZoom);
  hElecDirThetaZoom->GetXaxis()->SetTitle("(Reco - True) Electron Direction Theta (rad)");
  MakeNicePlots(hElecDirPhiZoom);
  hElecDirPhiZoom->GetXaxis()->SetTitle("(Reco - True) Electron Direction Phi (rad)");
  MakeNicePlots(hElecVtxDistZoom);
  hElecVtxDistZoom->GetXaxis()->SetTitle("Electron |(Reco Vertex - True Vertex)| (cm)");
  MakeNicePlots(hElecVertexTZoom);
  hElecVertexTZoom->GetXaxis()->SetTitle("(Reco - True) Electron Vertex Time (ns)");
  MakeNicePlots(hElecEnergyZoom);
  hElecEnergyZoom->GetXaxis()->SetTitle("(Reco - True) Electron Energy (MeV)");
  MakeNicePlots(hElecAngleZoom);
  hElecAngleZoom->GetXaxis()->SetTitle("Angle between Reco and True Directions (^o)");

  MakeNicePlots(hMuonVertexX);
  hMuonVertexX->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Position X (cm)");
  MakeNicePlots(hMuonVertexY);
  hMuonVertexY->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Position Y (cm)");
  MakeNicePlots(hMuonVertexZ);
  hMuonVertexZ->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Position Z (cm)");
  MakeNicePlots(hMuonDirTheta);
  hMuonDirTheta->GetXaxis()->SetTitle("(Reco - True) Cosmic Direction Theta (rad)");
  MakeNicePlots(hMuonDirPhi);
  hMuonDirPhi->GetXaxis()->SetTitle("(Reco - True) Cosmic Direction Phi (rad)");
  MakeNicePlots(hMuonVtxDist);
  hMuonVtxDist->GetXaxis()->SetTitle("Cosmic |(Reco Vertex - True Vertex)| (cm)");
  MakeNicePlots(hMuonVertexT);
  hMuonVertexT->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Time (ns)");
  MakeNicePlots(hMuonEnergy);
  hMuonEnergy->GetXaxis()->SetTitle("(Reco - True) Cosmic Energy (MeV)");
  MakeNicePlots(hMuonAngle);
  hMuonAngle->GetXaxis()->SetTitle("Angle between Reco and True Directions (^o)");

  MakeNicePlots(hMuonVertexXZoom);
  hMuonVertexXZoom->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Position X (cm)");
  MakeNicePlots(hMuonVertexYZoom);
  hMuonVertexYZoom->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Position Y (cm)");
  MakeNicePlots(hMuonVertexZZoom);
  hMuonVertexZZoom->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Position Z (cm)");
  MakeNicePlots(hMuonDirThetaZoom);
  hMuonDirThetaZoom->GetXaxis()->SetTitle("(Reco - True) Cosmic Direction Theta (rad)");
  MakeNicePlots(hMuonDirPhiZoom);
  hMuonDirPhiZoom->GetXaxis()->SetTitle("(Reco - True) Cosmic Direction Phi (rad)");
  MakeNicePlots(hMuonVtxDistZoom);
  hMuonVtxDistZoom->GetXaxis()->SetTitle("Cosmic |(Reco Vertex - True Vertex)| (cm)");
  MakeNicePlots(hMuonVertexTZoom);
  hMuonVertexTZoom->GetXaxis()->SetTitle("(Reco - True) Cosmic Vertex Time (ns)");
  MakeNicePlots(hMuonEnergyZoom);
  hMuonEnergyZoom->GetXaxis()->SetTitle("(Reco - True) Cosmic Energy (MeV)");
  MakeNicePlots(hMuonAngleZoom);
  hMuonAngleZoom->GetXaxis()->SetTitle("Angle between Reco and True Directions (^o)");

  std::stringstream outName;
  outName << "~/work/CHIPS/recon/finalCosmic/plots/overlayResolutions_" << run << ".root";
  TFile output(outName.str().c_str(),"RECREATE");
  output.cd();
  hElecVertexX->Write();
  hElecVertexY->Write();
  hElecVertexZ->Write();
  hElecDirTheta->Write();
  hElecDirPhi->Write();
  hElecVtxDist->Write();
  hElecVertexT->Write();
  hElecEnergy->Write();
  hElecAngle->Write();

  hMuonVertexX->Write();
  hMuonVertexY->Write();
  hMuonVertexZ->Write();
  hMuonDirTheta->Write();
  hMuonDirPhi->Write();
  hMuonVtxDist->Write();
  hMuonVertexT->Write();
  hMuonEnergy->Write();
  hMuonAngle->Write();

  hElecVertexXZoom->Write();
  hElecVertexYZoom->Write();
  hElecVertexZZoom->Write();
  hElecDirThetaZoom->Write();
  hElecDirPhiZoom->Write();
  hElecVtxDistZoom->Write();
  hElecVertexTZoom->Write();
  hElecEnergyZoom->Write();
  hElecAngleZoom->Write();

  hMuonVertexXZoom->Write();
  hMuonVertexYZoom->Write();
  hMuonVertexZZoom->Write();
  hMuonDirThetaZoom->Write();
  hMuonDirPhiZoom->Write();
  hMuonVtxDistZoom->Write();
  hMuonVertexTZoom->Write();
  hMuonEnergyZoom->Write();
  hMuonAngleZoom->Write();

  output.Close();

  std::cout << "Total number of matches = " << cTotalMatched << std::endl << std::endl;
  std::cout << "== Skipped events summary ==" << std::endl;
  std::cout << "- Contained muon  : " << cContainedMuon << std::endl;
  std::cout << "- Exiting electron: " << cExitingElec << std::endl;
  std::cout << "- Electron energy : " << cElecEnergy << std::endl;
  std::cout << "- Muon energy     : " << cMuonEnergy << std::endl;
}

void MakeNicePlots(TH1D* h){

  h->GetXaxis()->CenterTitle(1);

}


