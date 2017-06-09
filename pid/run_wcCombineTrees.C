/////////////////////////////////////////////////////////////////////////////////////////
//  Macro to run wc_combineTrees.C to conbine Electron and Muon hypotheses needed      //
//  for Neural Network training used for PID Classification.                           //  
//  (NueCC vs NumuCC) and (NueCC vs NC) using  WCSimAnalys output Tree                 //
//                                                                                     //
//  S. Germani - UCL (s.germani@ucl.ac.uk)                                             //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

void run_wcCombineTrees(const char * inputFileName_mu, const char * inputFileName_el, const char * outputFileName){

  // Load all Libraries and compile wc_readTMVA macro if necessary
  gROOT->ProcessLine(".x LoadLibs.C" );

  //Run macro to combine Trees
  wc_combineTrees( inputFileName_mu, inputFileName_el, outputFileName );
  
}
