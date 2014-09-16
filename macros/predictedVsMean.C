#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"

void predictedVsMean(){
//  if(str == "")	TFile * f = new TFile("../likelihood.root");
//	else
  Double_t zMax=0.0;

  TString str("../likelihoodDoubleBins.root");
  TFile * f = new TFile(str.Data());
  
  TFile *f2 = new TFile("testminim.root");
  TTree *t2 = (TTree*)f2->Get("likelihood");
  
  TTree * t = (TTree*)f->Get("likelihood");
	Double_t X,Y,Z,Qpred, Qdigi, Qundigi;
  Double_t X2, Y2, Z2;
	t->SetBranchAddress("X",&X);
	t->SetBranchAddress("Y",&Y);
	t->SetBranchAddress("Z",&Z);
	t->SetBranchAddress("Qpred",&Qpred);
  t2->SetBranchAddress("Qdigi",&Qdigi);
  t2->SetBranchAddress("Qundigi",&Qundigi);
  t2->SetBranchAddress("X",&X2);
  t2->SetBranchAddress("Y",&Y2);
  t2->SetBranchAddress("Z",&Z2);


	TH2D * hObsUndigiTop = new TH2D("hObsUndigiTop","Observed undigitized charge (top)",58,-2075,2075,58,-2075,2075);
	TH2D * hObsUndigiBottom = new TH2D("hObsUndigiBottom","Observed undigitized charge (bottom)",58,-2075,2075,58,-2075,2075);
	TH2D * hObsUndigiWall = new TH2D("hObsUndigiWall","Observed undigitized charge (wall)",64,0,2.0*TMath::Pi(), 64, -1075., 1075.);
	TH2D * hObsDigiTop = new TH2D("hObsDigiTop","Observed digitized charge (top)",58,-2075,2075,58,-2075,2075);
	TH2D * hObsDigiBottom = new TH2D("hObsDigiBottom","Observed digitized charge (bottom)",58,-2075,2075,58,-2075,2075);
	TH2D * hObsDigiWall = new TH2D("hObsDigiWall","Observed digitized charge (wall)",64,0,2.0*TMath::Pi(), 64, -1075., 1075.);
	TH2D * hPredTop = new TH2D("hPredTop","Expected charge (top)",58,-2075,2075,58,-2075,2075);
	TH2D * hPredBottom = new TH2D("hPredBottom","Expected charge (bottom)",58,-2075,2075,58,-2075,2075);
	TH2D * hPredWall = new TH2D("hPredWall","Expected charge (wall)",64,0.0,2.0*TMath::Pi(), 64, -1075., 1075.);

  TH1D * hRatioTop = new TH1D("hRatioTop","Ratio mean/predicted for top", 120,0,6);
  TH1D * hRatioBottom = new TH1D("hRatioBottom","Ratio mean/predicted for bottom", 150,0,10);
  TH1D * hRatioWall = new TH1D("hRatioWall","Ratio mean/predicted for wall", 150,0,10);

  TH2D * h2DRatioTop = new TH2D("h2DRatioTop","Ratio mean/predicted for top", 58, -2075, 2075, 58, -2075, 2075);
  TH2D * h2DRatioBottom = new TH2D("h2DRatioBottom","Ratio mean/predicted for bottom", 58, -2075, 2075, 58, -2075, 2075);
  TH2D * h2DRatioWall = new TH2D("h2DRatioWall","Ratio mean/predicted for wall", 150, 0.0, 2.0*TMath::Pi(), 150, -2075., 2075.);
  
  std::cout << t->GetEntries() << "   " << t2->GetEntries() << std::endl;

	for(int i =0; i < t->GetEntries(); ++i)
  {
		t->GetEntry(i);
    if(Z > 995)
    { 
			hPredTop->Fill(X,Y,Qpred);
		}
    else if(Z < -995)
    {
      hPredBottom->Fill(X,Y,Qpred);
	  }
    else
    {
      Double_t phi = TMath::ATan2(Y,X);
      if(phi < 0) phi += 2*TMath::Pi();
      hPredWall->Fill(phi, Z, Qpred);
    }
  }    
	
  for(int i = 0; i < t2->GetEntries(); ++i)
  {
  	t2->GetEntry(i);
  	if(Z2 > 995)
    { 
		  hObsDigiTop->Fill(X2,Y2,Qdigi);
		  hObsUndigiTop->Fill(X2,Y2,Qundigi);
		}
    else if(Z2 < -995)
    {
      hObsUndigiBottom->Fill(X2,Y2,Qundigi);
      hObsDigiBottom->Fill(X2,Y2,Qdigi);
	  }
    else
    {
      Double_t phi = TMath::ATan2(Y2,X2);
      if(phi < 0) phi += 2*TMath::Pi();
      hObsUndigiWall->Fill(phi, Z2, Qundigi);
      hObsDigiWall->Fill(phi, Z2, Qdigi);
    }
  }

  for(Int_t i = 1; i < hPredTop->GetNbinsX(); ++i)
  { 
    for(Int_t j=1; j < hPredTop->GetNbinsY(); ++j)
    {
      if(hPredTop->GetBinContent(i,j) > 0.5)
      {
          h2DRatioTop->SetBinContent(i,j, hObsUndigiTop->GetBinContent(i,j)/hPredTop->GetBinContent(i,j));
          hRatioTop->Fill(hObsDigiTop->GetBinContent(i,j)/hPredTop->GetBinContent(i,j));
      }
    }
  }
  for(Int_t i = 1; i < hPredBottom->GetNbinsX(); ++i)
  { 
    for(Int_t j=1; j < hPredBottom->GetNbinsY(); ++j)
    {
      if(hPredBottom->GetBinContent(i,j) > 0.5)
      {
          h2DRatioBottom->SetBinContent(i,j,hObsUndigiBottom->GetBinContent(i,j)/hPredBottom->GetBinContent(i,j));
          hRatioBottom->Fill(hObsUndigiBottom->GetBinContent(i,j)/hPredBottom->GetBinContent(i,j));
      }
    }
  }

  for(Int_t i = 1; i < hPredWall->GetNbinsX(); ++i)
  { 
    for(Int_t j=1; j < hPredWall->GetNbinsY(); ++j)
    {
      if(hPredWall->GetBinContent(i,j) > 0.5)
      {
          hRatioWall->Fill(hObsUndigiWall->GetBinContent(i,j)/hPredWall->GetBinContent(i,j));
          h2DRatioWall->SetBinContent( i,j, hObsUndigiWall->GetBinContent(i,j)/hPredWall->GetBinContent(i,j));

      }
    }
  }

//      if(Qpred>0 && Qobs>0)
//      {
//
//      TVector3 hit(X, Y, Z);
//      TVector3 vtx(0.,0.,0.);
//      TVector3 dir(0.,0.,1.);
//
//      Double_t r = (hit-vtx).Mag();
//      Double_t cosTheta = TMath::Cos( (hit-vtx).Angle(dir));
//      //std::cout << r << "   " << cosTheta << std::endl;
//      hObsPredDistance->Fill( r, Qobs/Qpred);
//      hObsPredAngle->Fill(cosTheta, Qobs,Qpred);
//
//      }
    
//      if(Qpred>0 && Qobs>0)
//      {
//
//      TVector3 hit(X, Y, Z);
//      TVector3 vtx(0.,0.,0.);
//      TVector3 dir(0.,0.,1.);
//
//      Double_t r = (hit-vtx).Mag();
//      Double_t cosTheta = TMath::Cos( (hit-vtx).Angle(dir));
//      //std::cout << r << "   " << cosTheta << std::endl;
//      hObsPredDistance->Fill( r, Qobs/Qpred);
//      hObsPredAngle->Fill(cosTheta, Qobs,Qpred);
//
//      }

    hObsDigiTop->GetXaxis()->SetTitle("hit x (cm)");
    hObsDigiTop->GetYaxis()->SetTitle("hit y (cm)");
    hObsDigiBottom->GetXaxis()->SetTitle("hit x (cm)");
    hObsDigiBottom->GetYaxis()->SetTitle("hit y (cm)");
    hObsDigiWall->GetXaxis()->SetTitle("hit #phi (radians)");
    hObsDigiWall->GetYaxis()->SetTitle("hit z (cm)");

    hObsUndigiTop->GetXaxis()->SetTitle("hit x (cm)");
    hObsUndigiTop->GetYaxis()->SetTitle("hit y (cm)");
    hObsUndigiBottom->GetXaxis()->SetTitle("hit x (cm)");
    hObsUndigiBottom->GetYaxis()->SetTitle("hit y (cm)");
    hObsUndigiWall->GetXaxis()->SetTitle("hit #phi (radians)");
    hObsUndigiWall->GetYaxis()->SetTitle("hit z (cm)");

    hPredTop->GetXaxis()->SetTitle("hit x (cm)");
    hPredTop->GetYaxis()->SetTitle("hit y (cm)");
    hPredBottom->GetXaxis()->SetTitle("hit x (cm)");
    hPredBottom->GetYaxis()->SetTitle("hit y (cm)");
    hPredWall->GetXaxis()->SetTitle("hit #phi (radians)");
    hPredWall->GetYaxis()->SetTitle("hit z (cm)");

    hRatioTop->GetXaxis()->SetTitle("Obs/pred");
    hRatioTop->GetYaxis()->SetTitle("Hits");
    hRatioBottom->GetXaxis()->SetTitle("Obs/pred");
    hRatioBottom->GetYaxis()->SetTitle("Hits");
    hRatioWall->GetXaxis()->SetTitle("Obs/pred");
    hRatioWall->GetYaxis()->SetTitle("Hits");

    h2DRatioTop->GetXaxis()->SetTitle("Obs/pred");
    h2DRatioTop->GetYaxis()->SetTitle("Hits");
    h2DRatioBottom->GetXaxis()->SetTitle("Obs/pred");
    h2DRatioBottom->GetYaxis()->SetTitle("Hits");
    h2DRatioWall->GetXaxis()->SetTitle("Obs/pred");
    h2DRatioWall->GetYaxis()->SetTitle("Hits");
  //  for(int xBin = 0; xBin < hObsTop->GetNbinsX(); ++xBin)
  //  {
  //    for(int yBin = 0; yBin < hObsTop->GetNbinsY(); ++yBin)
  //    {
  //      if( hPredTop->GetBinContent(xBin, yBin) && hObsTop->GetBinContent(xBin, yBin))
  //      {
  //        hObsTopOverPredTop->SetBinContent(xBin, yBin, hObsTop->GetBinContent(xBin, yBin) / hPredTop->GetBinContent(xBin, yBin));
  //        hObsPred->Fill(hObsTopOverPredTop->GetBinContent(xBin, yBin));
  //      }
  //    }
  //  }



  gStyle->SetOptStat(0);
//	TCanvas * c1 = new TCanvas("c1","c1",1200,900);
//	c1->Divide(3,3);
//	c1->cd(1);
//	hObsUndigiTop->Draw("COLZ");
////  hObsDigiTop->Draw("COLZ");
//	c1->cd(2);
//	hObsUndigiBottom->Draw("COLZ");
////  hObsDigiBottom->Draw("COLZ");
//  c1->cd(3);
//  hObsUndigiWall->Draw("COLZ");
////  hObsDigiWall->Draw("COLZ");
//  c1->cd(4);
//  hPredTop->Draw("COLZ");
//  c1->cd(5);
//  hPredBottom->Draw("COLZ");
//  c1->cd(6);
//  hPredWall->Draw("COLZ");
//  c1->cd(7);
//  hRatioTop->Draw();
////  h2DRatioTop->Draw("COLZ");
//  c1->cd(8);
//  hRatioBottom->Draw();
//  c1->cd(9);
//  hRatioWall->Draw();
//  c1->SaveAs(str.Data()+TString(".tabulated.png"));

  zMax = hObsDigiTop->GetMaximum();
  if( hPredTop->GetMaximum() > zMax) zMax = hPredTop->GetMaximum();

	TCanvas * c1 = new TCanvas("c1","c1",640,640);
	c1->cd();
//	hObsUndigiTop->Draw("COLZ");
  hObsDigiTop->GetZaxis()->SetRangeUser(0,zMax);
  hObsDigiTop->GetZaxis()->SetTitle("charge (p.e.)");
  hObsDigiTop->GetZaxis()->SetTitleOffset(0.78);
  hObsDigiTop->GetYaxis()->SetTitleOffset(1.43);
  hObsDigiTop->GetXaxis()->SetLabelSize(0.04);
  hObsDigiTop->GetYaxis()->SetLabelSize(0.04);
  hObsDigiTop->GetZaxis()->SetLabelSize(0.04);
  hObsDigiTop->Draw("COLZ");
  c1->SaveAs("observed.png");

	TCanvas * c2 = new TCanvas("c2","c2",640,640);
	c2->cd();
  hPredTop->GetZaxis()->SetRangeUser(0,zMax);
  hPredTop->GetZaxis()->SetTitleOffset(0.78);
  hPredTop->GetZaxis()->SetTitle("charge (p.e.)");
  hPredTop->GetYaxis()->SetTitleOffset(1.43);
  hPredTop->GetXaxis()->SetLabelSize(0.04);
  hPredTop->GetYaxis()->SetLabelSize(0.04);
  hPredTop->GetZaxis()->SetLabelSize(0.04);
  hPredTop->Draw("COLZ");
  c2->SaveAs("predicted.png");
  c2->SaveAs(str.Data()+TString(".four.C"));

	return;
}
