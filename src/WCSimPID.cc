/*
 * WCSimPID.cc
 *
 *  Created on: 5 August 2018
 *      Author: jtingey
 */

#include "WCSimPID.hh"
#include "TList.h"
#include "TFile.h"
#include <TSystem.h>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <iostream>
#include <cstdlib>

#include "WCSimPIDTree.hh"

#ifndef REFLEX_DICTIONARY
ClassImp (WCSimPID);
#endif

WCSimPID::WCSimPID() {
	fElTrainDir = "";
	fMuTrainDir = "";
	fNcTrainDir = "";

	fElAllDir = "";
	fMuAllDir = "";
}

WCSimPID::~WCSimPID() {
	// TODO Auto-generated destructor stub
}

// Public Methods
void WCSimPID::Train() {
	std::cout << " *** WCSimPID::Train() *** " << std::endl;

	// Check the training directories have been set...
	assert(fElTrainDir != "" && fMuTrainDir != "" && fNcTrainDir != "");

	// Define the names of the training files we create and then use...
	const char * elFile = "elTraining.root";
	const char * muFile = "muTraining.root";

	const char * elFile_NCSelect = "elTraining_NCSelect.root";
	const char * ncFile_NCSelect = "ncTraining_NCSelect.root";

	// Make the training files for the Nue vs Numu network
	MakeTrainingFile(fElTrainDir, elFile);
	MakeTrainingFile(fMuTrainDir, muFile);

	// Train the Nue vs Numu network...
	TrainTMVA_electronMuonCCQE(elFile, muFile);

	// Make the training files for the Nue vs NC network
	MakeTrainingFile_NCSelect(fElTrainDir, elFile_NCSelect);
	MakeTrainingFile_NCSelect(fNcTrainDir, ncFile_NCSelect);

	// Train the Nue vs NC network...
	TrainTMVA_electronCCQEvsNC(elFile_NCSelect, ncFile_NCSelect);
}

void WCSimPID::Read() {
	std::cout << " *** WCSimPID::Read() *** " << std::endl;

	// Check the reading directories have been set...
	assert(fElAllDir != "" && fMuAllDir != "");

	// Define the names of the output files...
	const char * elAllFile = "elAllOutput.root";
	const char * muAllFile = "muAllOutput.root";

	// Read to file for the two directories...
	ReadToFile(fElAllDir, elAllFile);
	ReadToFile(fMuAllDir, muAllFile);

	// Find the optimal cuts to maximise eff*pur
	ScanCuts(elAllFile, muAllFile);
}

// Private Methods
void WCSimPID::MakeCombinedFile(std::string trainDir, const char * trainingFileName) {
	std::cout << "Making combined File from -> " << trainDir << std::endl;

	// We want to fill a TChain with all the files in the directory that are worth using...
	TChain * chain = new TChain("fResultsTree");

	char* dir = gSystem->ExpandPathName(trainDir.c_str());
    void* dirp = gSystem->OpenDirectory(dir);
    const char* entry;
    TString str;
	int fileSkipNoTree = 0;
	int fileSkipBadBranches = 0;
	int totalFiles = 0;
    while((entry = (char*)gSystem->GetDirEntry(dirp))){
        str = entry;
        if(str.EndsWith("_tree.root")){
			totalFiles += 1;
			//std::cout << "Processing File -> " << entry << " ..." << std::endl;

			// Open the file and check a couple of things...
			// 1) It contains the TTree "fResultsTree"
			// 2) "fResultsTree" contains 6 branches (EventHeader, TruthInfo, RecoSummary_ElectronLike, PIDInfo_ElectronLike, RecoSummary_MuonLike, PIDInfo_MuonLike)
			TFile * inputFile = new TFile(gSystem->ConcatFileName(dir, entry), "READ");
			if(!inputFile->GetListOfKeys()->Contains("fResultsTree")){ fileSkipNoTree += 1; continue;}
			TTree * resultsTree = (TTree*) (inputFile->Get("fResultsTree"));
			if(resultsTree->GetListOfBranches()->GetEntries() != 6) { fileSkipBadBranches += 1; continue; }

			// If it passes the tests, close the file and add to the TChain so we can use its events...
			inputFile->Close();
			delete inputFile;
			chain->Add(gSystem->ConcatFileName(dir, entry));
		}
	}

	// Now we merge the TChain into a new TFile with the name given by trainingFileName...
	chain->Merge(trainingFileName);

	std::cout << "Total Files -> " << totalFiles << ", Used -> " << (totalFiles-(fileSkipNoTree+fileSkipBadBranches)) << std::endl;
	std::cout << "fileSkipNoTree -> " << fileSkipNoTree << ", fileSkipBadBranches -> " << fileSkipBadBranches << std::endl;
	std::cout << "Chain contains " << chain->GetEntries() << " events" << std::endl;

	delete chain;
}

void WCSimPID::MakeTrainingFile(std::string trainDir, const char * trainingFileName) {

	// Make the combined file first
	MakeCombinedFile(trainDir, trainingFileName);

	// Now we add the PIDTree to the training file which merges various parameters from the resultsTree into one tree...
	TFile * trainingFile = new TFile(trainingFileName,"UPDATE");
	assert(trainingFile->GetListOfKeys()->Contains("fResultsTree"));
    TTree * resultsTree = (TTree*) (trainingFile->Get("fResultsTree"));
	int numEntries = resultsTree->GetEntries();

	WCSimPIDTree * pidTree = new WCSimPIDTree(resultsTree);

	int totalCut = 0;
	// Loop over the events
	for(Int_t iEvent = 0; iEvent < numEntries; ++iEvent) {
		if (!pidTree->GetEntry(iEvent)) { 
			totalCut += 1;
			continue; 
		}
		pidTree->FillTree();
		pidTree->Clear();
	}
	std::cout << "Total Events -> " << numEntries << ", used -> " << (numEntries-totalCut) << std::endl;

	TTree * tree = pidTree->GetOutputTree();
	tree->Write();
	trainingFile->Close();

	delete pidTree;
	delete trainingFile;
}

void WCSimPID::MakeTrainingFile_NCSelect(std::string trainDir, const char * trainingFileName) {

	// Make the combined file first
	MakeCombinedFile(trainDir, trainingFileName);

	// Now we add the PIDTree to the training file which merges various parameters from the resultsTree into one tree...
	TFile * trainingFile = new TFile(trainingFileName,"UPDATE");
	assert(trainingFile->GetListOfKeys()->Contains("fResultsTree"));
    TTree * resultsTree = (TTree*) (trainingFile->Get("fResultsTree"));
	int numEntries = resultsTree->GetEntries();

	WCSimPIDTree * pidTree = new WCSimPIDTree(resultsTree);

	// Define the variables for the TMVA reader
	float annNueCCQEvsNumuCCQE = 0.0;

	float nHits;
	float totalQ;
    float deltaCharge2LnL;
    float deltaTime2LnL;
    float deltaCharge2LnLOverNHits;

	float charge2LnL_el;
	float charge2LnL_mu;
	float time2LnL_el;
	float time2LnL_mu;

	float fracHitsOutsideRing_el;
    float fracHitsOutsideRing_mu;
	float fracHitsInRing_el;
    float fracHitsInRing_mu;
	float fracQOutsideRing_el;
    float fracQOutsideRing_mu;
	float fracQInRing_el;
    float fracQInRing_mu;
	float fracPredQOutsideRing_el;
    float fracPredQOutsideRing_mu;
	float fracPredQInRing_el;
    float fracPredQInRing_mu;

    float predictedChargeOverTotalCharge_mu;
    float predictedChargeOverTotalCharge_el;
    float recoEOverQ_mu;
    float recoEOverQ_el;

	// Reader for the nue CCQE vs. numu CCQE ANN
	TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

	reader->AddVariable("nHits",&nHits);
	reader->AddVariable("totalQ",&nHits);
	reader->AddVariable("deltaCharge2LnL",&deltaCharge2LnL);
	reader->AddVariable("deltaTime2LnL",&deltaTime2LnL);
	reader->AddVariable("deltaCharge2LnLOverNHits",&deltaCharge2LnLOverNHits);

	reader->AddVariable("charge2LnL_el",&charge2LnL_el);
	reader->AddVariable("charge2LnL_mu",&charge2LnL_mu);
	reader->AddVariable("time2LnL_el",&time2LnL_el);
	reader->AddVariable("time2LnL_mu",&time2LnL_mu);

	reader->AddVariable("fracHitsOutsideRing_el",&fracHitsOutsideRing_el);
	reader->AddVariable("fracHitsOutsideRing_mu",&fracHitsOutsideRing_mu);
	reader->AddVariable("fracHitsInRing_el",&fracHitsInRing_el);
	reader->AddVariable("fracHitsInRing_mu",&fracHitsInRing_mu);
	reader->AddVariable("fracQOutsideRing_el",&fracQOutsideRing_el);
	reader->AddVariable("fracQOutsideRing_mu",&fracQOutsideRing_mu);
	reader->AddVariable("fracQInRing_el",&fracQInRing_el);
	reader->AddVariable("fracQInRing_mu",&fracQInRing_mu);
	reader->AddVariable("fracPredQOutsideRing_el",&fracPredQOutsideRing_el);
	reader->AddVariable("fracPredQOutsideRing_mu",&fracPredQOutsideRing_mu);
	reader->AddVariable("fracPredQInRing_el",&fracPredQInRing_el);
	reader->AddVariable("fracPredQInRing_mu",&fracPredQInRing_mu);

	reader->AddVariable("predictedChargeOverTotalCharge_el",&predictedChargeOverTotalCharge_el);
	reader->AddVariable("predictedChargeOverTotalCharge_mu",&predictedChargeOverTotalCharge_mu);
	reader->AddVariable("recoEOverQ_el",&recoEOverQ_el);
	reader->AddVariable("recoEOverQ_mu",&recoEOverQ_mu);
	
	// Set the weight file
	reader->BookMVA("MLP", "weights/tmva_nueCCQE_vs_numuCCQE_MLP.weights.xml");

	int totalCut = 0;
	// Loop over the events
	for(Int_t iEvent = 0; iEvent < numEntries; ++iEvent) {
		if (!pidTree->GetEntry(iEvent)) { 
			totalCut += 1;
			continue; 
		}

		// Fill the needed variables for the readers from the PIDTree
		nHits = (float)(pidTree->GetNHits());
		totalQ = pidTree->GetTotalQ();
		deltaCharge2LnL = pidTree->GetDeltaCharge2LnL();
		deltaTime2LnL = pidTree->GetDeltaTime2LnL();
		deltaCharge2LnLOverNHits = pidTree->GetDeltaCharge2LnLOverNHits();

		charge2LnL_el = pidTree->GetCharge2LnL_el();
		charge2LnL_mu = pidTree->GetCharge2LnL_mu();
		time2LnL_el = pidTree->GetTime2LnL_el();
		time2LnL_mu = pidTree->GetTime2LnL_mu();

		fracHitsOutsideRing_el = pidTree->GetFracHitsOutsideRing_el();
		fracHitsOutsideRing_mu = pidTree->GetFracHitsOutsideRing_mu();
		fracHitsInRing_el = pidTree->GetFracHitsInRing_el();
		fracHitsInRing_mu = pidTree->GetFracHitsInRing_mu();
		fracQOutsideRing_el = pidTree->GetFracQOutsideRing_el();
		fracQOutsideRing_mu = pidTree->GetFracQOutsideRing_mu();
		fracQInRing_el = pidTree->GetFracQInRing_el();
		fracQInRing_mu = pidTree->GetFracQInRing_mu();
		fracPredQOutsideRing_el = pidTree->GetFracPredQOutsideRing_el();
		fracPredQOutsideRing_mu = pidTree->GetFracPredQOutsideRing_mu();
		fracPredQInRing_el = pidTree->GetFracPredQInRing_el();
		fracPredQInRing_mu = pidTree->GetFracPredQInRing_mu();

		predictedChargeOverTotalCharge_mu = pidTree->GetPredictedChargeOverTotalCharge_mu();
		predictedChargeOverTotalCharge_el = pidTree->GetPredictedChargeOverTotalCharge_el();
		recoEOverQ_mu = pidTree->GetRecoEOverQ_mu();
		recoEOverQ_el = pidTree->GetRecoEOverQ_el();

		// Evaluate the TMVA reader and set value in pidTree...
		annNueCCQEvsNumuCCQE = reader->EvaluateMVA("MLP");
		pidTree->SetNueVsNumu(annNueCCQEvsNumuCCQE);

		pidTree->FillTree();
		pidTree->Clear();
	}
	std::cout << "Total Events -> " << numEntries << ", used -> " << (numEntries-totalCut) << std::endl;

	TTree * tree = pidTree->GetOutputTree();
	tree->Write();
	trainingFile->Close();

	delete pidTree;
	delete trainingFile;
}

void WCSimPID::TrainTMVA_electronMuonCCQE(const char * sigFileName, const char * bkgFileName) {
	std::cout << " *** WCSimPID::TrainTMVA_electronMuonCCQE() *** " << std::endl;

    // Fill the signal and background TChains...
    TChain * sig = new TChain("PIDTree_ann"); 
    sig->Add(sigFileName);
	TChain * bkgNumuCC = new TChain("PIDTree_ann");
    bkgNumuCC->Add(bkgFileName);

    // Open the output file and create the TMVA factory object...
    TFile* outputFile = TFile::Open("weights/tmva_nueCCQE_vs_numuCCQE.root", "RECREATE");
    TMVA::Factory* factory = new TMVA::Factory("tmva_nueCCQE_vs_numuCCQE", outputFile, "");

    // global event weights (see below for setting event-wise weights)
    double sigWeight = 1.0; double bkgNumuCCWeight = 1.0;
    factory->AddSignalTree(sig, sigWeight);
    factory->AddBackgroundTree(bkgNumuCC, bkgNumuCCWeight);

	/*
    factory->AddVariable("nHits",'F');
    factory->AddVariable("deltaCharge2LnL",'F');
    factory->AddVariable("deltaTime2LnL",'F');
    factory->AddVariable("deltaCharge2LnLOverNHits",'F');
    factory->AddVariable("fracHitsOutsideRing_mu",'F');
    factory->AddVariable("fracHitsInRing_mu",'F');
    factory->AddVariable("fracHitsOutsideRing_el",'F');
    factory->AddVariable("fracHitsInRing_el",'F');
    factory->AddVariable("fracPredQOutsideRing_mu",'F');
    factory->AddVariable("fracPredQOutsideRing_el",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_mu",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_el",'F');
    factory->AddVariable("recoEOverQ_mu",'F');
    factory->AddVariable("recoEOverQ_el",'F');
	*/

	// Variables for Nue vs Numu ANN...
    factory->AddVariable("nHits",'F');
	factory->AddVariable("totalQ",'F');
    factory->AddVariable("deltaCharge2LnL",'F');
    factory->AddVariable("deltaTime2LnL",'F');
    factory->AddVariable("deltaCharge2LnLOverNHits",'F');	

	factory->AddVariable("charge2LnL_el",'F');
	factory->AddVariable("charge2LnL_mu",'F');
	factory->AddVariable("time2LnL_el",'F');
    factory->AddVariable("time2LnL_mu",'F');

	factory->AddVariable("fracHitsOutsideRing_el",'F');
    factory->AddVariable("fracHitsOutsideRing_mu",'F');
	factory->AddVariable("fracHitsInRing_el",'F');
    factory->AddVariable("fracHitsInRing_mu",'F');
    factory->AddVariable("fracQOutsideRing_el",'F');
    factory->AddVariable("fracQOutsideRing_mu",'F');
	factory->AddVariable("fracQInRing_el",'F');
    factory->AddVariable("fracQInRing_mu",'F');
	factory->AddVariable("fracPredQOutsideRing_el",'F');
    factory->AddVariable("fracPredQOutsideRing_mu",'F');
	factory->AddVariable("fracPredQInRing_el",'F');
    factory->AddVariable("fracPredQInRing_mu",'F');
    
	factory->AddVariable("predictedChargeOverTotalCharge_el",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_mu",'F');
    factory->AddVariable("recoEOverQ_el",'F');
    factory->AddVariable("recoEOverQ_mu",'F');
    
    //factory->AddVariable("nRings",'F');
    //factory->AddVariable("firstRingHeight",'F');
    //factory->AddVariable("lastRingHeight",'F');

    // We apply the preselection...
	//TCut mycut = "preselected == 1 && !veto";
    TCut mycut = "preselected == 1";

    // Use half of the events for training, half for testing
    factory->PrepareTrainingAndTestTree(mycut, mycut, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

    // Book MVA method (see TMVA manual) for more methods...  
    factory->BookMethod(TMVA::Types::kMLP, "MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5" );

    // Train, test and evaluate all methods
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    // Save the output and finish up
    outputFile->Close();
    delete factory;
    return;
}

void WCSimPID::TrainTMVA_electronCCQEvsNC(const char * sigFileName, const char * bkgFileName) {
	std::cout << " *** WCSimPID::TrainTMVA_electronCCQEvsNC() *** " << std::endl;

	// Fill the signal and background TChains...
    TChain * sig = new TChain("PIDTree_ann");
    sig->Add(sigFileName);
    TChain * bkgNumuNC = new TChain("PIDTree_ann");
    bkgNumuNC->Add(bkgFileName);

    // Open the output file and create the TMVA factory object...
    TFile* outputFile = TFile::Open("weights/tmva_nueCCQE_vs_numuNC.root", "RECREATE");
    TMVA::Factory* factory = new TMVA::Factory("tmva_nueCCQE_vs_numuNC", outputFile, "");

    // global event weights (see below for setting event-wise weights)
    double sigWeight = 1.0; double bkgNumuNCWeight = 1.0;
    factory->AddSignalTree(sig, sigWeight);
    factory->AddBackgroundTree(bkgNumuNC, bkgNumuNCWeight);

	/*
    factory->AddVariable("nHits",'F');
    factory->AddVariable("deltaCharge2LnL",'F');
    factory->AddVariable("deltaTime2LnL",'F');
    factory->AddVariable("deltaCharge2LnLOverNHits",'F');
    factory->AddVariable("fracHitsOutsideRing_mu",'F');
    factory->AddVariable("fracHitsInRing_mu",'F');
    factory->AddVariable("fracHitsOutsideRing_el",'F');
    factory->AddVariable("fracHitsInRing_el",'F');
    factory->AddVariable("fracPredQOutsideRing_mu",'F');
    factory->AddVariable("fracPredQOutsideRing_el",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_mu",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_el",'F');
    factory->AddVariable("recoEOverQ_mu",'F');
    factory->AddVariable("recoEOverQ_el",'F');
    factory->AddVariable("dirX_el",'F');
    factory->AddVariable("fracHitsDownstream",'F');
	*/

	// Variables for Nue vs Numu ANN...
    factory->AddVariable("nHits",'F');
	factory->AddVariable("totalQ",'F');
    factory->AddVariable("deltaCharge2LnL",'F');
    factory->AddVariable("deltaTime2LnL",'F');
    factory->AddVariable("deltaCharge2LnLOverNHits",'F');	

	factory->AddVariable("charge2LnL_el",'F');
	factory->AddVariable("charge2LnL_mu",'F');
	factory->AddVariable("time2LnL_el",'F');
    factory->AddVariable("time2LnL_mu",'F');

	factory->AddVariable("fracHitsOutsideRing_el",'F');
    factory->AddVariable("fracHitsOutsideRing_mu",'F');
	factory->AddVariable("fracHitsInRing_el",'F');
    factory->AddVariable("fracHitsInRing_mu",'F');
    factory->AddVariable("fracQOutsideRing_el",'F');
    factory->AddVariable("fracQOutsideRing_mu",'F');
	factory->AddVariable("fracQInRing_el",'F');
    factory->AddVariable("fracQInRing_mu",'F');
	factory->AddVariable("fracPredQOutsideRing_el",'F');
    factory->AddVariable("fracPredQOutsideRing_mu",'F');
	factory->AddVariable("fracPredQInRing_el",'F');
    factory->AddVariable("fracPredQInRing_mu",'F');
    
	factory->AddVariable("predictedChargeOverTotalCharge_el",'F');
    factory->AddVariable("predictedChargeOverTotalCharge_mu",'F');
    factory->AddVariable("recoEOverQ_el",'F');
    factory->AddVariable("recoEOverQ_mu",'F');

	factory->AddVariable("dirX_el",'F');
	factory->AddVariable("dirX_mu",'F');
	factory->AddVariable("dirY_el",'F');
	factory->AddVariable("dirY_mu",'F');
	factory->AddVariable("dirZ_el",'F');
	factory->AddVariable("dirZ_mu",'F');

	factory->AddVariable("fracHitsDownstream",'F');
	factory->AddVariable("fracQDownstream",'F');
	factory->AddVariable("fracHitsBelowMid",'F');
	factory->AddVariable("fracQBelowMid",'F');

    //factory->AddVariable("nRings",'F');
    //factory->AddVariable("firstRingHeight",'F');
    //factory->AddVariable("lastRingHeight",'F');
    
    // We apply the preselection...
	//TCut mycut = "preselected == 1 && annNueCCQEvsNumuCCQE>0.8 && !veto";
    TCut mycut = "preselected == 1 && annNueCCQEvsNumuCCQE>0.8";
	//TCut mycut = "preselected == 1";

    // Use half of the events for training, half for testing
    factory->PrepareTrainingAndTestTree(mycut, mycut, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

    // Book MVA method (see TMVA manual) for more methods...  
    factory->BookMethod(TMVA::Types::kMLP, "MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5" );

    // Train, test and evaluate all methods
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();    

    // Save the output and finish up
    outputFile->Close();
    delete factory;
    return;
}

void WCSimPID::MakePlots(const char * sigFileName, const char * bkgFileName) {
	std::cout << " *** WCSimPID::MakePlots() *** " << std::endl;

	// Now we read the PIDTree's and make plots of sig vs background
	//TFile * sigFile = new TFile(sigFileName,"READ");
	//assert(sigFile->GetListOfKeys()->Contains("PIDTree_ann"));
    //TTree * sigTree = (TTree*) (sigFile->Get("PIDTree_ann"));
	//int sigFileEntries = sigTree->GetEntries();

	//TFile * bkgFile = new TFile(bkgFileName,"READ");
	//assert(bkgFile->GetListOfKeys()->Contains("PIDTree_ann"));
    //TTree * bkgTree = (TTree*) (bkgFile->Get("PIDTree_ann"));
	//int bkgFileEntries = bkgTree->GetEntries();

	//sigFile->Close();
	//bkgFile->Close();
}

void WCSimPID::ReadToFile(std::string readDir, const char * readFileName) {

	// Make the combined file first
	MakeCombinedFile(readDir, readFileName);

	// Now we add the PIDTree to the reading file and Evaluate the two ANNs
	TFile * readingFile = new TFile(readFileName,"UPDATE");
	assert(readingFile->GetListOfKeys()->Contains("fResultsTree"));
    TTree * resultsTree = (TTree*) (readingFile->Get("fResultsTree"));
	int numEntries = resultsTree->GetEntries();

	WCSimPIDTree * pidTree = new WCSimPIDTree(resultsTree);

	// Define the variables for the TMVA reader
	float annNueCCQEvsNumuCCQE = 0.0;
	float annNueCCQEvsNC       = 0.0;

	float nHits = 0.0;
	float totalQ = 0.0;
    float deltaCharge2LnL = 0.0;
    float deltaTime2LnL = 0.0;
    float deltaCharge2LnLOverNHits = 0.0;

	float charge2LnL_el = 0.0;
	float charge2LnL_mu = 0.0;
	float time2LnL_el = 0.0;
	float time2LnL_mu = 0.0;

	float fracHitsOutsideRing_el = 0.0;
    float fracHitsOutsideRing_mu = 0.0;
	float fracHitsInRing_el = 0.0;
    float fracHitsInRing_mu = 0.0;
	float fracQOutsideRing_el = 0.0;
    float fracQOutsideRing_mu = 0.0;
	float fracQInRing_el = 0.0;
    float fracQInRing_mu = 0.0;
	float fracPredQOutsideRing_el = 0.0;
    float fracPredQOutsideRing_mu = 0.0;
	float fracPredQInRing_el = 0.0;
    float fracPredQInRing_mu = 0.0;

    float predictedChargeOverTotalCharge_mu = 0.0;
    float predictedChargeOverTotalCharge_el = 0.0;
    float recoEOverQ_mu = 0.0;
    float recoEOverQ_el = 0.0;

	float dirX_el = 0.0;
	float dirX_mu = 0.0;
	float dirY_el = 0.0;
	float dirY_mu = 0.0;
	float dirZ_el = 0.0;
	float dirZ_mu = 0.0;

	float fracHitsDownstream = 0.0;
	float fracQDownstream = 0.0;
	float fracHitsBelowMid = 0.0;
	float fracQBelowMid = 0.0;

	// Readers for the nue CCQE vs. numu CCQE MVAs and the nue CCQE vs NC trees
	TMVA::Reader* readerElMu = new TMVA::Reader("!Color:!Silent");
	TMVA::Reader* readerElNc = new TMVA::Reader("!Color:!Silent");

	// Add the variables to the readers...
    readerElMu->AddVariable("nHits",&nHits);
	readerElMu->AddVariable("totalQ",&totalQ);
    readerElMu->AddVariable("deltaCharge2LnL",&deltaCharge2LnL);
    readerElMu->AddVariable("deltaTime2LnL",&deltaTime2LnL);
    readerElMu->AddVariable("deltaCharge2LnLOverNHits",&deltaCharge2LnLOverNHits);	
	readerElMu->AddVariable("charge2LnL_el",&charge2LnL_el);
	readerElMu->AddVariable("charge2LnL_mu",&charge2LnL_mu);
	readerElMu->AddVariable("time2LnL_el",&time2LnL_el);
    readerElMu->AddVariable("time2LnL_mu",&time2LnL_mu);
	readerElMu->AddVariable("fracHitsOutsideRing_el",&fracHitsOutsideRing_el);
    readerElMu->AddVariable("fracHitsOutsideRing_mu",&fracHitsOutsideRing_mu);
	readerElMu->AddVariable("fracHitsInRing_el",&fracHitsInRing_el);
    readerElMu->AddVariable("fracHitsInRing_mu",&fracHitsInRing_mu);
    readerElMu->AddVariable("fracQOutsideRing_el",&fracQOutsideRing_el);
    readerElMu->AddVariable("fracQOutsideRing_mu",&fracQOutsideRing_mu);
	readerElMu->AddVariable("fracQInRing_el",&fracQInRing_el);
    readerElMu->AddVariable("fracQInRing_mu",&fracQInRing_mu);
	readerElMu->AddVariable("fracPredQOutsideRing_el",&fracPredQOutsideRing_el);
    readerElMu->AddVariable("fracPredQOutsideRing_mu",&fracPredQOutsideRing_mu);
	readerElMu->AddVariable("fracPredQInRing_el",&fracPredQInRing_el);
    readerElMu->AddVariable("fracPredQInRing_mu",&fracPredQInRing_mu);
	readerElMu->AddVariable("predictedChargeOverTotalCharge_el",&predictedChargeOverTotalCharge_el);
    readerElMu->AddVariable("predictedChargeOverTotalCharge_mu",&predictedChargeOverTotalCharge_mu);
    readerElMu->AddVariable("recoEOverQ_el",&recoEOverQ_el);
    readerElMu->AddVariable("recoEOverQ_mu",&recoEOverQ_mu);

    readerElNc->AddVariable("nHits",&nHits);
	readerElNc->AddVariable("totalQ",&totalQ);
    readerElNc->AddVariable("deltaCharge2LnL",&deltaCharge2LnL);
    readerElNc->AddVariable("deltaTime2LnL",&deltaTime2LnL);
    readerElNc->AddVariable("deltaCharge2LnLOverNHits",&deltaCharge2LnLOverNHits);	
	readerElNc->AddVariable("charge2LnL_el",&charge2LnL_el);
	readerElNc->AddVariable("charge2LnL_mu",&charge2LnL_mu);
	readerElNc->AddVariable("time2LnL_el",&time2LnL_el);
    readerElNc->AddVariable("time2LnL_mu",&time2LnL_mu);
	readerElNc->AddVariable("fracHitsOutsideRing_el",&fracHitsOutsideRing_el);
    readerElNc->AddVariable("fracHitsOutsideRing_mu",&fracHitsOutsideRing_mu);
	readerElNc->AddVariable("fracHitsInRing_el",&fracHitsInRing_el);
    readerElNc->AddVariable("fracHitsInRing_mu",&fracHitsInRing_mu);
    readerElNc->AddVariable("fracQOutsideRing_el",&fracQOutsideRing_el);
    readerElNc->AddVariable("fracQOutsideRing_mu",&fracQOutsideRing_mu);
	readerElNc->AddVariable("fracQInRing_el",&fracQInRing_el);
    readerElNc->AddVariable("fracQInRing_mu",&fracQInRing_mu);
	readerElNc->AddVariable("fracPredQOutsideRing_el",&fracPredQOutsideRing_el);
    readerElNc->AddVariable("fracPredQOutsideRing_mu",&fracPredQOutsideRing_mu);
	readerElNc->AddVariable("fracPredQInRing_el",&fracPredQInRing_el);
    readerElNc->AddVariable("fracPredQInRing_mu",&fracPredQInRing_mu);
	readerElNc->AddVariable("predictedChargeOverTotalCharge_el",&predictedChargeOverTotalCharge_el);
    readerElNc->AddVariable("predictedChargeOverTotalCharge_mu",&predictedChargeOverTotalCharge_mu);
    readerElNc->AddVariable("recoEOverQ_el",&recoEOverQ_el);
    readerElNc->AddVariable("recoEOverQ_mu",&recoEOverQ_mu);
	readerElNc->AddVariable("dirX_el",&dirX_el);
	readerElNc->AddVariable("dirX_mu",&dirX_mu);
	readerElNc->AddVariable("dirY_el",&dirY_el);
	readerElNc->AddVariable("dirY_mu",&dirY_mu);
	readerElNc->AddVariable("dirZ_el",&dirZ_el);
	readerElNc->AddVariable("dirZ_mu",&dirZ_mu);
	readerElNc->AddVariable("fracHitsDownstream",&fracHitsDownstream);
	readerElNc->AddVariable("fracQDownstream",&fracQDownstream);
	readerElNc->AddVariable("fracHitsBelowMid",&fracHitsBelowMid);
	readerElNc->AddVariable("fracQBelowMid",&fracQBelowMid);

	// Book the readers with their weight files...
	readerElMu->BookMVA("MLP", "weights/tmva_nueCCQE_vs_numuCCQE_MLP.weights.xml");
	readerElNc->BookMVA("MLP", "weights/tmva_nueCCQE_vs_numuNC_MLP.weights.xml");

	int totalCut = 0;
	// Loop over the events
	for(Int_t iEvent = 0; iEvent < numEntries; ++iEvent) {
		if (!pidTree->GetEntry(iEvent)) { 
			totalCut += 1;
			continue; 
		}

		// Fill the needed variables for the readers from the PIDTree
		nHits = (float)(pidTree->GetNHits());
		totalQ = pidTree->GetTotalQ();
		deltaCharge2LnL = pidTree->GetDeltaCharge2LnL();
		deltaTime2LnL = pidTree->GetDeltaTime2LnL();
		deltaCharge2LnLOverNHits = pidTree->GetDeltaCharge2LnLOverNHits();

		charge2LnL_el = pidTree->GetCharge2LnL_el();
		charge2LnL_mu = pidTree->GetCharge2LnL_mu();
		time2LnL_el = pidTree->GetTime2LnL_el();
		time2LnL_mu = pidTree->GetTime2LnL_mu();

		fracHitsOutsideRing_el = pidTree->GetFracHitsOutsideRing_el();
		fracHitsOutsideRing_mu = pidTree->GetFracHitsOutsideRing_mu();
		fracHitsInRing_el = pidTree->GetFracHitsInRing_el();
		fracHitsInRing_mu = pidTree->GetFracHitsInRing_mu();
		fracQOutsideRing_el = pidTree->GetFracQOutsideRing_el();
		fracQOutsideRing_mu = pidTree->GetFracQOutsideRing_mu();
		fracQInRing_el = pidTree->GetFracQInRing_el();
		fracQInRing_mu = pidTree->GetFracQInRing_mu();
		fracPredQOutsideRing_el = pidTree->GetFracPredQOutsideRing_el();
		fracPredQOutsideRing_mu = pidTree->GetFracPredQOutsideRing_mu();
		fracPredQInRing_el = pidTree->GetFracPredQInRing_el();
		fracPredQInRing_mu = pidTree->GetFracPredQInRing_mu();

		predictedChargeOverTotalCharge_mu = pidTree->GetPredictedChargeOverTotalCharge_mu();
		predictedChargeOverTotalCharge_el = pidTree->GetPredictedChargeOverTotalCharge_el();
		recoEOverQ_mu = pidTree->GetRecoEOverQ_mu();
		recoEOverQ_el = pidTree->GetRecoEOverQ_el();

		dirX_el = pidTree->GetDirX_el();
		dirX_mu = pidTree->GetDirX_mu();
		dirY_el = pidTree->GetDirY_el();
		dirY_mu = pidTree->GetDirY_mu();
		dirZ_el = pidTree->GetDirZ_el();
		dirZ_mu = pidTree->GetDirZ_mu();

		fracHitsDownstream = pidTree->GetFracHitsDownstream();
		fracQDownstream = pidTree->GetFracQDownstream();
		fracHitsBelowMid = pidTree->GetFracHitsBelowMid();
		fracQBelowMid = pidTree->GetFracQBelowMid();

		// Evaluate the TMVA reader and set values in the pidTree...
		annNueCCQEvsNumuCCQE = readerElMu->EvaluateMVA("MLP");
		annNueCCQEvsNC = readerElNc->EvaluateMVA("MLP");

		pidTree->SetNueVsNumu(annNueCCQEvsNumuCCQE);
		pidTree->SetNueVsNC(annNueCCQEvsNC);

		pidTree->FillTree();
		pidTree->Clear();
	}
	std::cout << "Total Events -> " << numEntries << ", used -> " << (numEntries-totalCut) << std::endl;

	TTree * tree = pidTree->GetOutputTree();
	tree->Write();
	readingFile->Close();

	delete pidTree;
	delete readingFile;	
}

void WCSimPID::ScanCuts(const char * elAllFile, const char * muAllFile) {
	std::cout << " *** WCSimPID::ScanCuts() *** " << std::endl;

	// Settings...
	double normNumEvents = 1000;
  	double scaleNue = 0.05;
  	double scaleNumu = 0.95;

	double min = 0.7;
	double max = 1.0;
	int nBins = 30;
	double width = (max-min)/(double)nBins;

	TChain * PIDTree_ann = new TChain("PIDTree_ann", "PIDTree_ann");
  	PIDTree_ann->Add(elAllFile);
  	PIDTree_ann->Add(muAllFile);

	// Number of pure CC and NC events
    int numNueCCEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueCCEvent");
    int numNueNCEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueNCEvent");
    int numNumuCCEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueCCEvent");
    int numNumuNCEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueNCEvent");

    // Number of each type of CC event 
    int numNueCCQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueCCEvent && trueQEEvent");
    int numNumuCCQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueCCEvent && trueQEEvent");
    int numNueCCnonQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueCCEvent && !trueQEEvent");
    int numNumuCCnonQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueCCEvent && !trueQEEvent");

    double weightNueNCEvents  = normNumEvents * (scaleNue * 0.3)/(double)numNueNCEvents;
    double weightNumuNCEvents = normNumEvents * (scaleNumu * 0.3)/(double)numNumuNCEvents;
    double weightNueCCQEEvents  = normNumEvents * (scaleNue * 0.7 * 0.2)/(double)numNueCCQEEvents;
    double weightNumuCCQEEvents = normNumEvents * (scaleNumu * 0.7 * 0.2)/(double)numNumuCCQEEvents;
    double weightNueCCnonQEEvents  = normNumEvents * (scaleNue * 0.7 * 0.8)/(double)numNueCCnonQEEvents;
    double weightNumuCCnonQEEvents = normNumEvents * (scaleNumu * 0.7 * 0.8)/(double)numNumuCCnonQEEvents;

	double maxCCQEValue = 0.0;
	double xCutCCQEOpt = 0.0;
	double yCutCCQEOpt = 0.0;
	double maxCCValue = 0.0;
	double xCutCCOpt = 0.0;
	double yCutCCOpt = 0.0;

	for(int xBin=1; xBin<=nBins; xBin++) {
		double xCut = min + ((double)xBin*width);
		for (int yBin=0; yBin<=nBins; yBin++) {
			double yCut = min + ((double)yBin*width);

			// Our cut to select each sample is comes from joining together the preselection and ANN 
            // cuts, and then adding the cut for this specific sample
			
			// First the nue CCQEs
			TString nueCCQECutString = Form("(trueBeamPDG==12 && trueCCEvent && trueQEEvent) && (preselected && (annNueCCQEvsNC > %f) && (annNueCCQEvsNumuCCQE > %f))", yCut, xCut);
  			double numNueCCQECut = PIDTree_ann->GetEntries(nueCCQECutString.Data());

			// Then the numu CCQEs
			TString numuCCQECutString = Form("(trueBeamPDG==14 && trueCCEvent && trueQEEvent) && (preselected && (annNueCCQEvsNC > %f) && (annNueCCQEvsNumuCCQE > %f))", yCut, xCut);
  			double numNumuCCQECut = PIDTree_ann->GetEntries(numuCCQECutString.Data());

			// Then the nue NCs
			TString nueNCCutString = Form("(trueBeamPDG==12 && trueNCEvent) && (preselected && (annNueCCQEvsNC > %f) && (annNueCCQEvsNumuCCQE > %f))", yCut, xCut);
  			double numNueNCCut = PIDTree_ann->GetEntries(nueNCCutString.Data());

			// Then the numu NCs
			TString numuNCCutString = Form("(trueBeamPDG==14 && trueNCEvent) && (preselected && (annNueCCQEvsNC > %f) && (annNueCCQEvsNumuCCQE > %f))", yCut, xCut);
  			double numNumuNCCut = PIDTree_ann->GetEntries(numuNCCutString.Data());

			// The the nue CC but non-CCQEs
			TString nueCCnonQECutString = Form("(trueBeamPDG==12 && trueCCEvent && !trueQEEvent) && (preselected && (annNueCCQEvsNC > %f) && (annNueCCQEvsNumuCCQE > %f))", yCut, xCut);
  			double numNueCCnonQECut = PIDTree_ann->GetEntries(nueCCnonQECutString.Data());

			// And finally the numu CC but non-CCQEs
			TString numuCCnonQECutString = Form("(trueBeamPDG==14 && trueCCEvent && !trueQEEvent) && (preselected && (annNueCCQEvsNC > %f) && (annNueCCQEvsNumuCCQE > %f))", yCut, xCut);
  			double numNumuCCnonQECut = PIDTree_ann->GetEntries(numuCCnonQECutString.Data());

			// Now we work out the efficiencies: the numerator is the weighted number of selected events, 
        	// and the denominator is the weighted number of total events in this sample.  The error is the
            // error on the mean of a binomially distributed variable: err = sqrt(p * (1-p) / n)

            // First we'll do the efficiency for selecting nue CCQE events
            //double nueCCQEEff = (numNueCCQECut*weightNueCCQEEvents)/(numNueCCQEEvents*weightNueCCQEEvents);
			// Error example if we wanted to use it...
            //errorNueCCQE = math.sqrt(nueCCQEEff * (1-nueCCQEEff) / (numNueCCQEEvents*weightNueCCQEEvents))

			// Work out totals...
			double nueCCAll   = ((numNueCCQECut*weightNueCCQEEvents) + (numNueCCnonQECut*weightNueCCnonQEEvents));
            double nueCCQEAll = (numNueCCQECut*weightNueCCQEEvents);
            double numuCCAll  = ((numNumuCCQECut*weightNumuCCQEEvents) + (numNumuCCnonQECut*weightNumuCCnonQEEvents));
            double nuNCAll    = ((numNumuNCCut*weightNumuNCEvents) + (numNueNCCut*weightNueNCEvents));

			// Work out Nue CCQE Efficiency and Purity
			double nueCCQEEff = nueCCQEAll/(numNueCCQEEvents*weightNueCCQEEvents);
			double nueCCQEPur = nueCCQEAll/(nueCCAll+numuCCAll+nuNCAll);

			// Work out Nue CC Efficiency and Purity
			double nueCCEff = nueCCAll/((numNueCCQEEvents*weightNueCCQEEvents) + (numNueCCnonQEEvents*weightNueCCnonQEEvents));
			double nueCCPur = nueCCAll/(nueCCAll+numuCCAll+nuNCAll);

            // Work out eff*pur
			double fomCCQE = nueCCQEEff * nueCCQEPur;
			double fomCC = nueCCEff * nueCCPur;

			// Is it the best value so far?
			if (fomCCQE >= maxCCQEValue) {
				maxCCQEValue = fomCCQE;
				xCutCCQEOpt = xCut;
				yCutCCQEOpt = yCut;
			}

			if (fomCC >= maxCCValue) {
				maxCCValue = fomCC;
				xCutCCOpt = xCut;
				yCutCCOpt = yCut;
			}

            //double fom_cc = nueCCAll/sqrt(nueCCAll+numuCCAll+nuNCAll);
            //double fom_ccqe = nueCCQEAll/sqrt(nueCCAll+numuCCAll+nuNCAll);		
		}
	}

	std::cout << "CCQE Opt (" << maxCCQEValue << ") at x=" << xCutCCQEOpt << " y=" << yCutCCQEOpt << std::endl;
	std::cout << "CC Opt (" << maxCCValue << ") at x=" << xCutCCOpt << " y=" << yCutCCOpt << std::endl;
}
