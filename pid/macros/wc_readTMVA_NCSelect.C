/////////////////////////////////////////////////////////////////////////////////////////
//  ROOT Macro for reading the TMVA PID neural network weights for PID classification  //
//  (NueCC vs NumuCC) and combining the Trees so the output files can be               //
//     used for (NueCC vs NC) training                                                 //
//  Adapted from uncommited code by A.J.Perch (UCL)                                    //
//  Adapted from code by S. Germani - UCL (s.germani@ucl.ac.uk)                        //                                                           //
//       J. Tingey UCL (j.tingey.16@ucl.ac.uk)                                         //
/////////////////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TString.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <iostream>
#include <cstdlib>
#include <string>

void wc_readTMVA_NCSelect(const char * inputFileName_mu, const char * inputFileName_el, const char * baseOutputName, const char * PIDDir){

	const int NUM_READERS = 2;

	TMVA::Tools::Instance();

	// Load Root related Libraries
	gSystem->Load("libGeom");
	gSystem->Load("libEve");
	gSystem->Load("libMinuit.so");

	// WCSim And WCSimAnalysis Libraries
	gSystem->Load("$WCSIMHOME/libWCSimRoot.so");
	gSystem->Load("$WCSIMANAHOME/lib/libWCSimAnalysis.so");

	// Put muon fit files from .txt file into vector
    std::ifstream muFile(inputFileName_mu);
    std::string str;
    std::vector<string> muFileVec;
    while(std::getline(muFile, str)){ muFileVec.push_back(str); }
    muFile.close();

    // Put electron fit files from .txt file into vector
    std::ifstream elFile(inputFileName_el);
    std::vector<string> elFileVec;
    while(std::getline(elFile, str)){ elFileVec.push_back(str); }
    elFile.close();

    // Check the vectors and the same size
    assert(muFileVec.size() == elFileVec.size());
    std::cout << "Files to combine-> " << muFileVec.size() << std::endl;

    // Loop over the pairs of files, combine into PID tree and read TMVA
    int numSkipped = 0;
    for(int f=0; f<muFileVec.size(); f++)
    {
    	// Open muon fit file
		TFile  * inputFile_mu = new TFile(muFileVec[f].c_str(), "READ");
		if (inputFile_mu == 0x0 || inputFile_mu->IsZombie()) {
			std::cout << "ERROR: could not read MuonLike input file: " << muFileVec[f].c_str() <<std::endl;
			exit(1);
		}

		// Open electron fit file
		TFile  * inputFile_el = new TFile(elFileVec[f].c_str(), "READ");
		if (inputFile_el == 0x0 || inputFile_el->IsZombie()) {
			std::cout << "ERROR: could not read ElectronLike input file: " << elFileVec[f].c_str() <<std::endl;
			exit(1);
		}

		// Load the input variables from the Muon and Electron Like Trees
		TTree * inputTree_mu = static_cast<TTree*>(inputFile_mu->Get("fResultsTree"));
		TTree * inputTree_el = static_cast<TTree*>(inputFile_el->Get("fResultsTree"));

		// Check the trees are not NULL and have the same number of entries, otherwise skip
		if(inputTree_mu != NULL && inputTree_el != NULL)
		{
			int numEvents_mu = inputTree_mu->GetEntries();
			int numEvents_el = inputTree_el->GetEntries();
			int numEvents = ( numEvents_mu == numEvents_el )? numEvents_mu : -999;

			if(numEvents_mu == numEvents_el)
			{
				// Create a WCSimPIDTree object that deals with reading the WCSimAnalysis output trees and makes the outputTree
			    WCSimPIDTree * treeEntry = new WCSimPIDTree(inputTree_mu, inputTree_el);

				// Make the temporary output name from the baseOutputName
				// NOTE: Use hadd afterwards to combine these together, run_wcPID.py does this for you
				std::string outputName(baseOutputName);
				ostringstream convert;
				convert << f;
				outputName = outputName + convert.str() + "_ann.root";
				TFile * outputFile = new TFile(outputName.c_str(), "RECREATE");

			    // Define the variables for the TMVA reader
			    // Have to do some conversions between doubles and floats in the tree
			    float annNueCCQEvsNumuCCQE = 0.0;

			    float nHits = 0.0;
			    float deltaCharge2LnL = 0.0;
			    float deltaTime2LnL = 0.0;
			    float deltaCharge2LnLOverNHits = 0.0;

			    float fracHitsOutsideRing_mu = 0.0;
			    float fracHitsInRing_mu = 0.0;
			    float fracHitsInRingHole_mu = 0.0;
			    float fracHitsOutsideRing_el = 0.0;
			    float fracHitsInRing_el = 0.0;
			    float fracHitsInRingHole_el = 0.0;

			    float fracPredQOutsideRing_mu = 0.0;
				float fracPredQOutsideRing_el = 0.0;

				float predictedChargeOverTotalCharge_mu = 0.0;
				float predictedChargeOverTotalCharge_el = 0.0;

				float recoEOverQ_mu = 0.0;
				float recoEOverQ_el = 0.0;

				float nRings = 0.0;
				float firstRingHeight = 0.0;
				float lastRingHeight = 0.0;

			    // Reader for the nue CCQE vs. numu CCQE MVAs
			    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

			    // Setup the reader. Add all the tree variables to the MVA and load the weights
			    // nue CCQE vs. numu variables
			    reader->AddVariable("nHits",&nHits);
			    reader->AddVariable("deltaCharge2LnL",&deltaCharge2LnL);
			    reader->AddVariable("deltaTime2LnL",&deltaTime2LnL);
			    reader->AddVariable("deltaCharge2LnLOverNHits",&deltaCharge2LnLOverNHits);
			    reader->AddVariable("fracHitsOutsideRing_mu",&fracHitsOutsideRing_mu);
			    reader->AddVariable("fracHitsInRing_mu",&fracHitsInRing_mu);
			    reader->AddVariable("fracHitsInRingHole_mu",&fracHitsInRingHole_mu);
			    reader->AddVariable("fracHitsOutsideRing_el",&fracHitsOutsideRing_el);
			    reader->AddVariable("fracHitsInRing_el",&fracHitsInRing_el);
			    reader->AddVariable("fracHitsInRingHole_el",&fracHitsInRingHole_el);
			    reader->AddVariable("fracPredQOutsideRing_mu",&fracPredQOutsideRing_mu);
			    reader->AddVariable("fracPredQOutsideRing_el",&fracPredQOutsideRing_el);
			    reader->AddVariable("predictedChargeOverTotalCharge_mu",&predictedChargeOverTotalCharge_mu);
			    reader->AddVariable("predictedChargeOverTotalCharge_el",&predictedChargeOverTotalCharge_el);
			    reader->AddVariable("recoEOverQ_mu",&recoEOverQ_mu);
			    reader->AddVariable("recoEOverQ_el",&recoEOverQ_el);

			    reader->AddVariable("nRings",&nRings);
			    reader->AddVariable("firstRingHeight",&firstRingHeight);
			    reader->AddVariable("lastRingHeight",&lastRingHeight);

			    // Set the weight file
			    std::string weightsDirectory(PIDDir);
			    TString weightPrefix;
			    weightPrefix = TString("/weights/tmva_nueCCQE_vs_numuCCQE_MLP.weights.xml");
			    TString weightFile = TString(weightsDirectory) + weightPrefix;
			    std::cout << "Weights in " << weightFile << std::endl;
			    reader->BookMVA("MLP", weightFile);

			    // Loop over the events
			    for(Int_t iEvent = 0; iEvent < numEvents; ++iEvent)
			    {
			    	// Get the PIDTree entry
					treeEntry->GetEntry(iEvent);

					// Fill the needed variables for the readers from the PIDTree
					nHits = (float)(treeEntry->GetNHits());
					deltaCharge2LnL = treeEntry->GetDeltaCharge2LnL();
					deltaTime2LnL = treeEntry->GetDeltaTime2LnL();
					deltaCharge2LnLOverNHits = treeEntry->GetDeltaCharge2LnLOverNHits();

					fracHitsOutsideRing_mu = treeEntry->GetFracHitsOutsideRing_mu();
					fracHitsInRing_mu = treeEntry->GetFracHitsInRing_mu();
					fracHitsInRingHole_mu = treeEntry->GetFracHitsInRingHole_mu();
					fracHitsOutsideRing_el = treeEntry->GetFracHitsOutsideRing_el();
					fracHitsInRing_el = treeEntry->GetFracHitsInRing_el();
					fracHitsInRingHole_el = treeEntry->GetFracHitsInRingHole_el();

					fracPredQOutsideRing_mu = treeEntry->GetFracPredQOutsideRing_mu();
					fracPredQOutsideRing_el = treeEntry->GetFracPredQOutsideRing_el();

					predictedChargeOverTotalCharge_mu = treeEntry->GetPredictedChargeOverTotalCharge_mu();
					predictedChargeOverTotalCharge_el = treeEntry->GetPredictedChargeOverTotalCharge_el();

					recoEOverQ_mu = treeEntry->GetRecoEOverQ_mu();
					recoEOverQ_el = treeEntry->GetRecoEOverQ_el();

					nRings = (float)(treeEntry->GetNRings());
					firstRingHeight = treeEntry->GetFirstRingHeight();
					lastRingHeight = treeEntry->GetLastRingHeight();

					// Evaluate the TMVA readers
					annNueCCQEvsNumuCCQE = reader->EvaluateMVA("MLP");

					treeEntry->SetNueVsNumu(annNueCCQEvsNumuCCQE);

					treeEntry->FillTree();
					treeEntry->Clear();

			    }

			    //Get the output Tree and save everything to file
			    TTree * outputTree = treeEntry->GetOutputTree();

			    outputTree->Write();
			    outputFile->Close();
			}
			else{ numSkipped += 1; }
		}
		else{ numSkipped += 1; }

		inputFile_mu->Close();
		inputFile_el->Close();
    }

    std::cout << "Read and Combined/Read " << (muFileVec.size() - numSkipped) << " files." << std::endl;
    std::cout << "Skipped " << numSkipped << " files." << std::endl;

};
