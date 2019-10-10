/////////////////////////////////////////////////////////////////////////////////////////
//  ROOT Macro for combining WCSimAnalysis output tree for Electron and Muon           //
//  hypothesis, into a WCSimPIDTree in order to train the TMVA PID Neural Network      //
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

void wc_combineTrees(const char * inputFileName_mu, const char * inputFileName_el, const char * baseOutputName){

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

    // Loop over the pairs of files and combine into PIDTree
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

				// Loop over the events
				for(Int_t iEvent = 0; iEvent < numEvents; ++iEvent)
				{
					treeEntry->GetEntry(iEvent);
					treeEntry->FillTree();
					treeEntry->Clear();
				}

				//Get the output Tree and save to file
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

    std::cout << "Combined " << (muFileVec.size() - numSkipped) << " files." << std::endl;
    std::cout << "Skipped " << numSkipped << " files." << std::endl;

};
