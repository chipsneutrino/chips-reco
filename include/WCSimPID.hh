/*
 * WCSimPID.hh
 *
 *  Created on: 3 March 2018
 *      Author: jtingey
 */

#include "TObject.h"
#include "TTree.h"
#include "TChain.h"

#include <iostream>
#include <cstdlib>

#pragma once

class WCSimPID : public TObject
{
public:
	WCSimPID();
	~WCSimPID();

	void Train();
	void Read();

	void SetTrainingDirs(std::string elTrainDir, std::string muTrainDir, std::string ncTrainDir)
	{
		fElTrainDir = elTrainDir;
		fMuTrainDir = muTrainDir;
		fNcTrainDir = ncTrainDir;
	}
	void PrintTrainingDirs()
	{
		std::cout << "Electron CCQE Training Dir -> " << fElTrainDir << std::endl;
		std::cout << "Muon CCQE Training Dir -> " << fMuTrainDir << std::endl;
		std::cout << "NC Training Dir -> " << fNcTrainDir << std::endl;
	}

	void SetReadingDirs(std::string elAllDir, std::string muAllDir)
	{
		fElAllDir = elAllDir;
		fMuAllDir = muAllDir;
	}
	void PrintReadingDirs()
	{
		std::cout << "Electron All Reading Dir -> " << fElAllDir << std::endl;
		std::cout << "Muon All Reading Dir -> " << fMuAllDir << std::endl;
	}

private:
	// General Methods
	void MakeCombinedFile(std::string trainDir, const char *trainingFileName);

	// Training Methods
	void MakeTrainingFile(std::string trainDir, const char *trainingFileName);
	void MakeTrainingFile_NCSelect(std::string trainDir, const char *trainingFileName);

	void TrainTMVA_electronMuonCCQE(const char *sigFileName, const char *bkgFileName);
	void TrainTMVA_electronCCQEvsNC(const char *sigFileName, const char *bkgFileName);

	void MakePlots(const char *sigFileName, const char *bkgFileName);

	// Reading Methods
	void ReadToFile(std::string readDir, const char *readFileName);
	void ScanCuts(const char *elAllFile, const char *muAllFile, bool makePlots);

	std::string fElTrainDir;
	std::string fMuTrainDir;
	std::string fNcTrainDir;

	std::string fElAllDir;
	std::string fMuAllDir;

	ClassDef(WCSimPID, 1);
};
