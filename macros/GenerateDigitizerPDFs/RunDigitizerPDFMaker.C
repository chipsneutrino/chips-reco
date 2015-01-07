void RunDigitizerPDFMaker(int muBin, int nThrows = 1000000, int nBins = 1000, double min = 0.0, double max = 10.0)
{
	// gApplication->ProcessLine(".L WCSimDigitizerPDFMaker.cc+");
  gSystem->Load("WCSimDigitizerPDFMaker_cc.so");
	WCSimDigitizerPDFMaker maker;
	maker.SetBinning(nBins, min, max);
	maker.SetNumThrows(nThrows);

	double mu = min + (max - min)/(nBins) * muBin;
  std::cout << "muBin   = " << muBin << std::endl;
  std::cout << "mu      = " << mu << std::endl;
  std::cout << "nThrows = " << nThrows << std::endl;
  std::cout << "nBins   = " << nBins << std::endl;
  std::cout << "min     = " << min << std::endl;
  std::cout << "max     = " << max << std::endl;
	maker.SetMu(mu);
	maker.Run();
}
