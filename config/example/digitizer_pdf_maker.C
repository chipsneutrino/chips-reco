R__LOAD_LIBRARY(libWCSimRoot.so)
R__LOAD_LIBRARY(libWCSimAnalysisRoot.so)
R__LOAD_LIBRARY(libGeom.so)
R__LOAD_LIBRARY(libEve.so)
R__LOAD_LIBRARY(libMinuit.so)

// 0 = sk1pe, 1 = chipspmt, 2 = tot
void digitizer_pdf_maker(int type = 0, int nThrows = 100000, int nBins = 1000, double min = 0.0, double max = 10.0)
{
	std::cout << "Type: " << type << ", Throws: " << nThrows << ", nBins: "
		<< nBins << ", min: " << min << ", max: " << max << std::endl;

	WCSimDigitizerPDFMaker maker;
	maker.SetBinning(nBins, min, max);
	maker.SetNumThrows(nThrows);
	maker.SetType(type);
	maker.Run();
}
