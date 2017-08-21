/////////////////////////////////////////////////////////////////////////////////////////
//  Macro to run wc_readTMVA.C  for Neural Network Event PID Classification            //  
//  (NueCC vs NumuCC) and (NueCC vs NC) using  WCSimAnalys output Tree                 //
//                                                                                     //
//  S. Germani - UCL (s.germani@ucl.ac.uk)                                             //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

void run_wcTMVA(const char * inputFileName_mu, const char * inputFileName_el, const char * outputFileName){

  // Load all Libraries and compile wc_readTMVA macro if necessary
  gROOT->ProcessLine(".x LoadLibs.C" );

  //Run macro to get Neural Networ PIDs
  wc_readTMVA( inputFileName_mu, inputFileName_el, outputFileName );
  
}
