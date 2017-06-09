/////////////////////////////////////////////////////////////////////////////////////////
//  Class and main function to combine WCSimAnalys output Tree for Electron and        //
//  Muon Hipitheses in order to train a TMVA Neural Network  for PID classification    //
//  (NueCC vs NumuCC) and (NueCC vs NC).                                               //
//  Adapted from uncommited code by A.J.Perch (UCL)                                    // 
//                                                                                     // 
//  S. Germani - UCL (s.germani@ucl.ac.uk)                                             // 
//                                                                                     //  
/////////////////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TString.h"
#include "TStopwatch.h"


#include "WCSimOutputTree.hh"
#include "WCSimRecoSummary.hh"
#include "WCSimTrueEvent.hh"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <iostream>
#include <cstdlib>


class PIDTreeEntry
{
    public:
  PIDTreeEntry(TTree * tree_mu = NULL, TTree * tree_el = NULL)
          {
	  assert(tree_mu != NULL && tree_el != NULL);
	  assert(tree_mu->GetEntriesFast() == tree_el->GetEntriesFast());
	  SetTreeBranches(tree_mu, tree_el);

	  hitInfo_el   = new HitInfo();  
	  recoInfo_el  = new RecoInfo(); 
	  truthInfo_el = new TruthInfo();
	  
	  hitInfo_mu   = new HitInfo();  
	  recoInfo_mu  = new RecoInfo(); 
	  truthInfo_mu = new TruthInfo();
        }
    
        ~PIDTreeEntry() {}

  Int_t NHits(){

    return  ( hitInfo_mu->fNHits  == hitInfo_el->fNHits  ) ? hitInfo_mu->fNHits  : -999;
  }

  Float_t TotalQ(){
    return  ( hitInfo_mu->fTotalQ  == hitInfo_el->fTotalQ  ) ? hitInfo_mu->fTotalQ  : -999;
  }

  Float_t DeltaCharge2LnL(){
    return  ( recoInfo_el->fCharge2LnL - recoInfo_mu->fCharge2LnL );
  }

  Float_t DeltaTime2LnL(){
    return  ( recoInfo_el->fTime2LnL - recoInfo_mu->fTime2LnL );
  }

  Float_t TotalPredQel(){
    return recoInfo_el->fPredQInRing + recoInfo_el->fPredQOutsideRing + recoInfo_el->fPredQInRingHole;
  }

  Float_t TotalPredQmu(){
    return recoInfo_mu->fPredQInRing + recoInfo_mu->fPredQOutsideRing + recoInfo_mu->fPredQInRingHole;
  }

  Bool_t Veto(){
    return ( hitInfo_mu->fVeto  == hitInfo_el->fVeto  ) ? hitInfo_mu->fVeto  : 1;
  }  
  void GetEntry( int entry){
    
    bhi_mu->GetEntry(entry);
    bri_mu->GetEntry(entry);
    bti_mu->GetEntry(entry);
    bhi_el->GetEntry(entry);
    bri_el->GetEntry(entry);
    bti_el->GetEntry(entry);
    nHits = NHits();
    deltaCharge2LnL = DeltaCharge2LnL();
    deltaTime2LnL   = DeltaTime2LnL();
    totalQ          = TotalQ();
    totalPredQ_mu   = TotalPredQmu();
    totalPredQ_el   = TotalPredQel();
    veto            = Veto();

  }
  
  void SetTreeBranches(TTree * tree_mu, TTree * tree_el)
          {
	  
	  ///// Muon Like  ...   	  
	  bhi_mu =  tree_mu->GetBranch("HitInfo");    bhi_mu->SetAddress(&hitInfo_mu);   
	  bri_mu =  tree_mu->GetBranch("RecoInfo");   bri_mu->SetAddress(&recoInfo_mu);  
	  bti_mu =  tree_mu->GetBranch("TruthInfo");  bti_mu->SetAddress(&truthInfo_mu); 
	  	  

	  ///// Electron Like ...
	  bhi_el =  tree_el->GetBranch("HitInfo");    bhi_el->SetAddress(&hitInfo_el);   
	  bri_el =  tree_el->GetBranch("RecoInfo");   bri_el->SetAddress(&recoInfo_el);  
	  bti_el =  tree_el->GetBranch("TruthInfo");  bti_el->SetAddress(&truthInfo_el); 
	}

  void Clear(){
    
    hitInfo_mu->Clear();
    recoInfo_mu->Clear();
    truthInfo_mu->Clear();

    hitInfo_el->Clear();
    recoInfo_el->Clear();
    truthInfo_el->Clear();
  }

  Int_t          nHits;
  Float_t        deltaCharge2LnL;
  Float_t        deltaTime2LnL;
  Float_t        totalQ;
  Float_t        totalPredQ_mu;
  Float_t        totalPredQ_el;
  Bool_t         veto;

  //Declaration of leaves types
  HitInfo     *hitInfo_el;  
  RecoInfo    *recoInfo_el; 
  TruthInfo   *truthInfo_el;
	                                      
  HitInfo     *hitInfo_mu;  
  RecoInfo    *recoInfo_mu; 
  TruthInfo   *truthInfo_mu;
  
  /// Branches
  TBranch *bhi_el; 
  TBranch *bri_el; 
  TBranch *bti_el; 
  TBranch *bhi_mu; 
  TBranch *bri_mu; 
  TBranch *bti_mu; 

};

TTree* MakeOutputTree(PIDTreeEntry* treeEntry);

void wc_combineTrees(const char * inputFileName_mu, const char * inputFileName_el, const char * outputFileName){

    TMVA::Tools::Instance();

    TFile  * inputFile_mu = new TFile(inputFileName_mu, "READ");
	if (inputFile_mu == 0x0 || inputFile_mu->IsZombie()) {
		std::cout << "ERROR: could not read MuonLike input file: " << inputFileName_mu <<std::endl;
		exit(1);
	}

    TFile  * inputFile_el = new TFile(inputFileName_el, "READ");
	if (inputFile_el == 0x0 || inputFile_el->IsZombie()) {
		std::cout << "ERROR: could not read ElectronLike input file: " << inputFileName_el <<std::endl;
		exit(1);
	}

    
    

    // Load the input variables from the Muon and Electron Like Trees
    TTree * inputTree_mu = static_cast<TTree*>(inputFile_mu->Get("fResultsTree"));
    TTree * inputTree_el = static_cast<TTree*>(inputFile_el->Get("fResultsTree"));

    PIDTreeEntry * treeEntry = new PIDTreeEntry(inputTree_mu, inputTree_el);

    // Make an output tree that copies the input tree and adds the ANN variables
    TFile * outputFile = new TFile(outputFileName, "RECREATE");
    float annNueCCQEvsNumuCCQE = 0.0;
    float annNueCCQEvsNC       = 0.0;
    TTree * outputTree = MakeOutputTree(treeEntry);
    outputTree->Branch("annNueCCQEvsNumuCCQE", &annNueCCQEvsNumuCCQE);
    outputTree->Branch("annNueCCQEvsNC", &annNueCCQEvsNC);

    // Have to do some conversions between doubles and floats in the tree
    float recoEOverQ_el = 0.0;
    float recoEOverQ_mu = 0.0;
    float chargeLikelihoodHitRatio_el = 0.0;
    float chargeLikelihoodHitRatio_mu = 0.0;
    float deltaCharge2LnL = 0.0;
    float deltaTime2LnL = 0.0;
    float fracQOutsideRing_mu = 0.0;
    float fracPredQOutsideRing_mu = 0.0;
    float predictedChargeOverTotalCharge_mu = 0.0;
    float predictedChargeOverTotalCharge_el = 0.0;
    float nHits = 0.0;

    
    // These are functions of variables for the Neural Network, so let's add them to the ouput one explicitly
    outputTree->Branch("recoEOverQ_el", &recoEOverQ_el);
    outputTree->Branch("recoEOverQ_mu", &recoEOverQ_mu);
    outputTree->Branch("chargeLikelihoodHitRatio_el", &chargeLikelihoodHitRatio_el);
    outputTree->Branch("chargeLikelihoodHitRatio_mu", &chargeLikelihoodHitRatio_mu);
    outputTree->Branch("predictedChargeOverTotalCharge_mu", &predictedChargeOverTotalCharge_mu);
    outputTree->Branch("predictedChargeOverTotalCharge_el", &predictedChargeOverTotalCharge_el);
    
    int numEvents_mu = inputTree_mu->GetEntries();
    int numEvents_el = inputTree_el->GetEntries();
    int numEvents    = ( numEvents_mu == numEvents_el )? numEvents_mu : -999;

    cout << " Total Events : " << numEvents << endl;

    if( numEvents <=0 ) { 
      cerr << "ERROR: Empty input file(s) or File Size mismatch!" << endl;
      cerr << "EXIT!"  << endl;
      exit(1);
    }

    int lastPercent = 0;
    TStopwatch sw;
    sw.Start();
    for(Int_t iEvent = 0; iEvent < numEvents; ++iEvent)
    {
        if( (100*iEvent/numEvents) > lastPercent )
        {
            lastPercent = (100*iEvent/numEvents);
            std::cout << "\r ***  - Completed " << lastPercent << "%   [";
            int j = 0;
            while(j < lastPercent/5)
            {
                std::cout << "=";
                j++;
            }
            while(j < 20)
            {
                std::cout << ".";
                j++;
            }
            std::cout << "]" << std::flush;
        }

	treeEntry->GetEntry(iEvent);

        nHits = (float)(treeEntry->NHits());
	//	cout << "nHits " << nHits << endl;
	if(nHits < 50){ continue; }
        deltaCharge2LnL = (float)(treeEntry->DeltaCharge2LnL());
        deltaTime2LnL   = (float)(treeEntry->DeltaTime2LnL());


        fracQOutsideRing_mu     = (float)(treeEntry->recoInfo_mu->fFracTotalQOutsideRing);
        fracPredQOutsideRing_mu = (float)(treeEntry->recoInfo_mu->fFracPredQOutsideRing);


        if(treeEntry->TotalQ() > 0)
        {
 	    recoEOverQ_el = treeEntry->recoInfo_el->fEnergy/(treeEntry->TotalQ());
	    recoEOverQ_mu = treeEntry->recoInfo_mu->fEnergy/(treeEntry->TotalQ());
            predictedChargeOverTotalCharge_mu = treeEntry->TotalPredQmu()/treeEntry->TotalQ();
            predictedChargeOverTotalCharge_el = treeEntry->TotalPredQel()/treeEntry->TotalQ();
        }
        else
        {
            recoEOverQ_el = -999.9;
            recoEOverQ_mu = -999.9;
            predictedChargeOverTotalCharge_mu = -999.9;
            predictedChargeOverTotalCharge_el = -999.9;
        }

        if(nHits)
        {
            chargeLikelihoodHitRatio_el = (treeEntry->recoInfo_el->fCharge2LnL)/(nHits);
            chargeLikelihoodHitRatio_mu = (treeEntry->recoInfo_mu->fCharge2LnL)/(nHits);
        }
        else
        {
            chargeLikelihoodHitRatio_el = -999.9;
            chargeLikelihoodHitRatio_mu = -999.9;
        }

        outputTree->Fill();
	
	treeEntry->Clear();
	
    }
    sw.Stop();
    std::cout << std::endl << " *** Processed " << numEvents << " in " << sw.RealTime() << "s" << std::endl;

    outputTree->Write();
    outputFile->Close();
    std::cout << " *** Combined Tree saved in file: " << outputFileName << std::endl;

};

TTree * MakeOutputTree(PIDTreeEntry * treeEntry)
{
    TTree * outputTree = new TTree("PIDTree_ann","PIDTree_ann");
    //outputTree->Branch("trueType",&(treeEntry->trueType));
    //outputTree->Branch("trueENu",&(treeEntry->trueENu));
    //outputTree->Branch("trueCCEvent",&(treeEntry->trueCCEvent));
    //outputTree->Branch("trueNCEvent",&(treeEntry->trueNCEvent));
    //outputTree->Branch("trueQEEEvent",&(treeEntry->trueQEEEvent));
    //outputTree->Branch("trueResEvent",&(treeEntry->trueResEvent));
    //outputTree->Branch("trueDISEvent",&(treeEntry->trueDISEvent));
    //outputTree->Branch("trueCohEvent",&(treeEntry->trueCohEvent));
    //outputTree->Branch("trueNueElectronElasticEvent",&(treeEntry->trueNueElectronElasticEvent));
    //outputTree->Branch("trueInverseMuonDecayEvent",&(treeEntry->trueInverseMuonDecayEvent));

    outputTree->Branch("charge2LnL_mu",&(treeEntry->recoInfo_mu->fCharge2LnL));
    outputTree->Branch("time2LnL_mu",&(treeEntry->recoInfo_mu->fTime2LnL));
    outputTree->Branch("recoE_mu",&(treeEntry->recoInfo_mu->fEnergy));

    outputTree->Branch("charge2LnL_el",&(treeEntry->recoInfo_el->fCharge2LnL));
    outputTree->Branch("time2LnL_el",&(treeEntry->recoInfo_el->fTime2LnL));
    outputTree->Branch("recoE_el",&(treeEntry->recoInfo_el->fEnergy));

    outputTree->Branch("deltaCharge2LnL",&(treeEntry->deltaCharge2LnL));

    
    outputTree->Branch("deltaTime2LnL",&(treeEntry->deltaTime2LnL));
    outputTree->Branch("totalQ_mu",&(treeEntry->totalQ));

    outputTree->Branch("totalQOutsideRing_mu",&(treeEntry->recoInfo_mu->fTotalQOutsideRing));
    outputTree->Branch("totalQInRing_mu",&(treeEntry->recoInfo_mu->fTotalQInRing));
    outputTree->Branch("totalQInRingHole_mu",&(treeEntry->recoInfo_mu->fTotalQInRingHole));
    outputTree->Branch("fracQOutsideRing_mu",&(treeEntry->recoInfo_mu->fFracTotalQOutsideRing));
    outputTree->Branch("fracQInRing_mu",&(treeEntry->recoInfo_mu->fFracTotalQInRing));
    outputTree->Branch("fracQInRingHole_mu",&(treeEntry->recoInfo_mu->fFracTotalQInRingHole));

    outputTree->Branch("totalQOutsideRing_el",&(treeEntry->recoInfo_el->fTotalQOutsideRing));
    outputTree->Branch("totalQInRing_el",&(treeEntry->recoInfo_el->fTotalQInRing));
    outputTree->Branch("totalQInRingHole_el",&(treeEntry->recoInfo_el->fTotalQInRingHole));
    outputTree->Branch("fracQOutsideRing_el",&(treeEntry->recoInfo_el->fFracTotalQOutsideRing));
    outputTree->Branch("fracQInRing_el",&(treeEntry->recoInfo_el->fFracTotalQInRing));
    outputTree->Branch("fracQInRingHole_el",&(treeEntry->recoInfo_el->fFracTotalQInRingHole));


    outputTree->Branch("nHits",&(treeEntry->nHits));

    outputTree->Branch("nHitsOutsideRing_mu",&(treeEntry->recoInfo_mu->fNHitsOutsideRing));
    outputTree->Branch("nHitsInRing_mu",&(treeEntry->recoInfo_mu->fNHitsInRing));
    outputTree->Branch("nHitsInRingHole_mu",&(treeEntry->recoInfo_mu->fNHitsInRingHole));
    outputTree->Branch("fracHitsOutsideRing_mu",&(treeEntry->recoInfo_mu->fFracNHitsOutsideRing));
    outputTree->Branch("fracHitsInRing_mu",&(treeEntry->recoInfo_mu->fFracNHitsInRing));
    outputTree->Branch("fracHitsInRingHole_mu",&(treeEntry->recoInfo_mu->fFracNHitsInRingHole));
    outputTree->Branch("totalPredQ_mu",&(treeEntry->totalPredQ_mu));
    outputTree->Branch("totalPredQOutsideRing_mu",&(treeEntry->recoInfo_mu->fPredQOutsideRing));
    outputTree->Branch("totalPredQInRing_mu",&(treeEntry->recoInfo_mu->fPredQInRing));
    outputTree->Branch("totalPredQInRingHole_mu",&(treeEntry->recoInfo_mu->fPredQInRingHole));
    outputTree->Branch("fracPredQOutsideRing_mu",&(treeEntry->recoInfo_mu->fFracPredQOutsideRing));
    outputTree->Branch("fracPredQInRing_mu",&(treeEntry->recoInfo_mu->fFracPredQInRing));
    outputTree->Branch("fracPredQInRingHole_mu",&(treeEntry->recoInfo_mu->fFracPredQInRingHole));
    outputTree->Branch("predictedChargeOverTotalCharge_mu",&(treeEntry->recoInfo_mu->fPredictedOverTotalCharge));


    outputTree->Branch("nHitsOutsideRing_el",&(treeEntry->recoInfo_el->fNHitsOutsideRing));
    outputTree->Branch("nHitsInRing_el",&(treeEntry->recoInfo_el->fNHitsInRing));
    outputTree->Branch("nHitsInRingHole_el",&(treeEntry->recoInfo_el->fNHitsInRingHole));
    outputTree->Branch("fracHitsOutsideRing_el",&(treeEntry->recoInfo_el->fFracNHitsOutsideRing));
    outputTree->Branch("fracHitsInRing_el",&(treeEntry->recoInfo_el->fFracNHitsInRing));
    outputTree->Branch("fracHitsInRingHole_el",&(treeEntry->recoInfo_el->fFracNHitsInRingHole));
    outputTree->Branch("totalPredQ_el",&(treeEntry->totalPredQ_el));
    outputTree->Branch("totalPredQOutsideRing_el",&(treeEntry->recoInfo_el->fPredQOutsideRing));
    outputTree->Branch("totalPredQInRing_el",&(treeEntry->recoInfo_el->fPredQInRing));
    outputTree->Branch("totalPredQInRingHole_el",&(treeEntry->recoInfo_el->fPredQInRingHole));
    outputTree->Branch("fracPredQOutsideRing_el",&(treeEntry->recoInfo_el->fFracPredQOutsideRing));
    outputTree->Branch("fracPredQInRing_el",&(treeEntry->recoInfo_el->fFracPredQInRing));
    outputTree->Branch("fracPredQInRingHole_el",&(treeEntry->recoInfo_el->fFracPredQInRingHole));
    outputTree->Branch("predictedChargeOverTotalCharge_el",&(treeEntry->recoInfo_el->fPredictedOverTotalCharge));

    outputTree->Branch("vtxRho_mu",&(treeEntry->recoInfo_mu->fVtxRho));
    outputTree->Branch("vtxZ_mu",&(treeEntry->recoInfo_mu->fVtxZ));
    outputTree->Branch("endRho_mu",&(treeEntry->recoInfo_mu->fEndRho));
    outputTree->Branch("endZ_mu",&(treeEntry->recoInfo_mu->fEndZ));
    outputTree->Branch("dirX_mu",&(treeEntry->recoInfo_mu->fDirX));
    outputTree->Branch("escapes_mu",&(treeEntry->recoInfo_mu->fEscapes));

    outputTree->Branch("vtxRho_el",&(treeEntry->recoInfo_el->fVtxRho));
    outputTree->Branch("vtxZ_el",&(treeEntry->recoInfo_el->fVtxZ));
    outputTree->Branch("endRho_el",&(treeEntry->recoInfo_el->fEndRho));
    outputTree->Branch("endZ_el",&(treeEntry->recoInfo_el->fEndZ));
    outputTree->Branch("dirX_el",&(treeEntry->recoInfo_el->fDirX));
    outputTree->Branch("escapes_el",&(treeEntry->recoInfo_el->fEscapes));


    outputTree->Branch("veto",&(treeEntry->veto));

    //outputTree->Branch("truePDG",&(treeEntry->truePDG));

    return outputTree;
};
