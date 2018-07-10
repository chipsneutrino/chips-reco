void runMakeNewScatteringFiles(const std::string input, const std::string output){

  gSystem->Load("/home/lwhitehead/work/CHIPS/recon/test/WCSim/libWCSimRoot.so");

  std::string arglist = "(input, output)";
  gROOT->LoadMacro("makeNewScatteringFiles.C++O");
  gROOT->ProcessLineFast(("makeNewScatteringFiles" + arglist).c_str());


}
