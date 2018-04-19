void combine(){
    //gSystem->Load("libGeom");
    //gSystem->Load("libEve");
    //gSystem->Load("libMinuit");
    //TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
    //TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
    //gSystem->Load(libWCSimRoot.Data());
    //gSystem->Load(libWCSimAnalysis.Data());

    //TFile * nikhefFile = new TFile("tot/verticalNorm_nikhef.root", "READ");
    //TFile * MadisonFile = new TFile("tot/verticalNorm_madison.root", "READ");
    //TFile * singleFile = new TFile("tot/verticalNorm.root", "READ");

    //TFile * nikhefFile = new TFile("tot/horizontalNorm_nikhef.root", "READ");
    //TFile * MadisonFile = new TFile("tot/horizontalNorm_madison.root", "READ");
    TFile * singleFile = new TFile("tot/horizontalNorm.root", "READ");

    //TH2D * nikhefHist = nikhefFile->Get("digiPDF");
    //TH2D * madisonHist = MadisonFile->Get("digiPDF");

    TH2D * nikhefHist = singleFile->Get("digiPDF");
    //TH2D * madisonHist = singleFile->Get("digiPDF");
    TH2D * madisonHist = nikhefHist->Clone();

    nikhefHist->SetName("digiPDF_nikef");
    nikhefHist->SetTitle("digiPDF_nikef");

    madisonHist->SetName("digiPDF_madison");
    madisonHist->SetTitle("digiPDF_madison");

    std::cout << nikhefHist->GetXaxis()->GetXmin() << std::endl;
    std::cout << nikhefHist->GetYaxis()->GetXmin() << std::endl;

    std::cout << nikhefHist->Interpolate(0, 1) << std::endl;
    std::cout << nikhefHist->Interpolate(1, 0) << std::endl;



    //TFile * mainOutput = new TFile("tot/TOTPDFs_vertical.root","RECREATE");
    //TFile * mainOutput = new TFile("tot/TOTPDFs_vertical_nonPoisson.root","RECREATE");
    //TFile * mainOutput = new TFile("tot/TOTPDFs_horizontal.root","RECREATE");
    //TFile * mainOutput = new TFile("tot/TOTPDFs_horizontal_nonPoisson.root","RECREATE");

    //nikhefHist->Write();
    //madisonHist->Write();

    //mainOutput->Close();

}
