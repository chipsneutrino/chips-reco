void RunDigitizerPDFMaker(int muBin = 5.0, int type = 0, int nThrows = 1000000, int nBins = 1000, double min = 0.0, double max = 10.0)
{
	gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
    TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
    gSystem->Load(libWCSimRoot.Data());
    gSystem->Load(libWCSimAnalysis.Data());

	WCSimDigitizerPDFMaker maker;
	maker.SetBinning(nBins, min, max);
	maker.SetNumThrows(nThrows);
	maker.SetType(type);
	maker.SetTotPoisson(false);

	double mu = min + (max - min)/(nBins) * muBin;
	//std::cout << "muBin   = " << muBin << std::endl;
	//std::cout << "mu      = " << mu << std::endl;
	//std::cout << "nThrows = " << nThrows << std::endl;
	//std::cout << "nBins   = " << nBins << std::endl;
	//std::cout << "min     = " << min << std::endl;
	//std::cout << "max     = " << max << std::endl;
	maker.SetMu(mu);
	maker.Run();
}
