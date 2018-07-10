void wc_integralLookupTester()
{
  // Load libraries
  // ==============
  gApplication->ProcessLine(".except");
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
  TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
  gSystem->Load(libWCSimRoot.Data());
  gSystem->Load(libWCSimAnalysis.Data());

  float E = 2500;
  TrackType::Type type = TrackType::MuonLike;
  double R0 = 1290.51;
  double cosTheta0 = 0.795948;

  TGraph * gRhoS0Int = new TGraph();
  TGraph * gRhoS1Int = new TGraph();
  TGraph * gRhoS2Int = new TGraph();
  gRhoS0Int->SetName("gRhoS0Int");
  gRhoS1Int->SetName("gRhoS1Int");
  gRhoS2Int->SetName("gRhoS2Int");

  TGraph * gRhoGS0Int = new TGraph();
  TGraph * gRhoGS1Int = new TGraph();
  TGraph * gRhoGS2Int = new TGraph();
  gRhoGS0Int->SetName("gRhoGS0Int");
  gRhoGS1Int->SetName("gRhoGS1Int");
  gRhoGS2Int->SetName("gRhoGS2Int");

  while( E <= 2800 )
  {

    gRhoGS0Int->SetPoint(gRhoGS0Int->GetN(), E, WCSimIntegralLookupReader::Instance()->GetRhoGIntegral(type, E, 0, R0, cosTheta0));
    gRhoGS1Int->SetPoint(gRhoGS1Int->GetN(), E, WCSimIntegralLookupReader::Instance()->GetRhoGSIntegral(type, E, 0, R0, cosTheta0));
    gRhoGS2Int->SetPoint(gRhoGS2Int->GetN(), E, WCSimIntegralLookupReader::Instance()->GetRhoGSSIntegral(type, E, 0, R0, cosTheta0));
    std::cout << "I[0] = " << WCSimIntegralLookupReader::Instance()->GetRhoGIntegral(type, E, 0, R0, cosTheta0) << std::endl;
    std::cout << "I[1] = " << WCSimIntegralLookupReader::Instance()->GetRhoGSIntegral(type, E, 0, R0, cosTheta0) << std::endl;
    std::cout << "I[2] = " << WCSimIntegralLookupReader::Instance()->GetRhoGSSIntegral(type, E, 0, R0, cosTheta0) << std::endl;

    gRhoS0Int->SetPoint(gRhoS0Int->GetN(), E, WCSimIntegralLookupReader::Instance()->GetRhoIntegral(type, E, 0));
    gRhoS1Int->SetPoint(gRhoS1Int->GetN(), E, WCSimIntegralLookupReader::Instance()->GetRhoSIntegral(type, E, 0));
    gRhoS2Int->SetPoint(gRhoS2Int->GetN(), E, WCSimIntegralLookupReader::Instance()->GetRhoSSIntegral(type, E, 0));

    E += 1;
    std::cout << "E = " << E << std::endl;
  }

  gRhoS0Int->SetMarkerColor(kRed);
  gRhoS1Int->SetMarkerColor(kBlue);
  gRhoS2Int->SetMarkerColor(kGreen-2);

  gRhoGS0Int->SetMarkerColor(kRed);
  gRhoGS1Int->SetMarkerColor(kBlue);
  gRhoGS2Int->SetMarkerColor(kGreen-2);

  TCanvas * can = new TCanvas("can","",1200,600);
  TString title = TString::Format("R_{0} = %dcm and cos#theta_{0} = %.05f", R0, cosTheta0);
  can->Divide(2,1);
  can->cd(1);
  gRhoS0Int->Draw("ALP");
  can->cd(2);
  gRhoGS0Int->Draw("ALP");
  gRhoS0Int->SetTitle(title.Data());
  gRhoGS0Int->SetTitle(title.Data());
  gRhoS0Int->GetXaxis()->SetTitle("Energy (MeV)");
  gRhoS0Int->GetYaxis()->SetTitle("Integral of #rho(s)");
  gRhoGS0Int->GetXaxis()->SetTitle("Energy (MeV)");
  gRhoGS0Int->GetYaxis()->SetTitle("Integral of #rho(s)#timesg(s,#theta)");
  can->SaveAs("integrals_S0.png");
  can->SaveAs("integrals_S0.root");

  can->cd(1);
  gRhoS1Int->Draw("ALP");
  can->cd(2);
  gRhoGS1Int->Draw("ALP");
  gRhoS1Int->SetTitle(title.Data());
  gRhoGS1Int->SetTitle(title.Data());
  gRhoS1Int->GetXaxis()->SetTitle("Energy (MeV)");
  gRhoS1Int->GetYaxis()->SetTitle("Integral of s#times#rho(s)");
  gRhoGS1Int->GetXaxis()->SetTitle("Energy (MeV)");
  gRhoGS1Int->GetYaxis()->SetTitle("Integral of s#times#rho(s)#times g(s,#theta)");
  can->SaveAs("integrals_S1.png");
  can->SaveAs("integrals_S1.root");

  can->cd(1);
  gRhoS2Int->Draw("ALP");
  can->cd(2);
  gRhoGS2Int->Draw("ALP");
  gRhoS2Int->SetTitle(title.Data());
  gRhoGS2Int->SetTitle(title.Data());
  gRhoS2Int->GetXaxis()->SetTitle("Energy (MeV)");
  gRhoS2Int->GetYaxis()->SetTitle("Integral of s^{2}#times#rho(s)");
  gRhoGS2Int->GetXaxis()->SetTitle("Energy (MeV)");
  gRhoGS2Int->GetYaxis()->SetTitle("Integral of s^{2}#times#rho(s)#times g(s,#theta)");
  can->SaveAs("integrals_S2.png");
  can->SaveAs("integrals_S2.root");

}
