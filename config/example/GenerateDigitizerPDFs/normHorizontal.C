void normHorizontal(){

    //const char* inputName = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/sk1pe/preNorm.root";
    //const char* inputName = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/pmtSim/preNorm.root";
    //const char* inputName = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/tot/preNorm_nikhef.root";
    //const char* inputName = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/tot/preNorm_madison.root";
    const char* inputName = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/tot/preNorm.root";

    TFile * inputFile = new TFile(inputName, "READ");

    //Get pointer to the TTree within the mainFile...
    TH2D *fProbHisto = (TH2D*)inputFile->Get("digiPDF");

    for(int iBinY = 1; iBinY <= fProbHisto->GetNbinsY(); ++iBinY)
    {
    	double totInRow = 0.0;
    	for(int iBinX = 1; iBinX <= fProbHisto->GetNbinsX(); ++iBinX)
    	{
    		totInRow += fProbHisto->GetBinContent(iBinX, iBinY);
    	}

    	for(int iBinX = 1; iBinX <= fProbHisto->GetNbinsX(); ++iBinX)
    	{
    		double tempBinContent = fProbHisto->GetBinContent(iBinX, iBinY);
    		fProbHisto->SetBinContent(iBinX, iBinY, (tempBinContent/totInRow));
    	}
    }

    //TFile * mainOutput = new TFile("/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/sk1pe/horizontalNorm.root","RECREATE");
    //TFile * mainOutput = new TFile("/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/pmtSim/horizontalNorm.root","RECREATE");
    //TFile * mainOutput = new TFile("/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/tot/horizontalNorm_nikhef.root","RECREATE");
    //TFile * mainOutput = new TFile("/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/tot/horizontalNorm_madison.root","RECREATE");
    TFile * mainOutput = new TFile("/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/macros/GenerateDigitizerPDFs/tot/horizontalNorm.root","RECREATE");
    fProbHisto->Write();
    mainOutput->Close();
}
