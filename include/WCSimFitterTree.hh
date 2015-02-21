/*
 * WCSimFitterTree.hh
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#ifndef INCLUDE_WCSIMFITTERTREE_HH_
#define INCLUDE_WCSIMFITTERTREE_HH_

#include "WCSimLikelihoodTrack.hh"
#include <vector>
class TFile;
class TString;
class TTree;
class WCSimRootGeom;
class WCSimLikelihoodRecoEvent;
class WCSimTrueEvent;
class WCSimRootEvent;
class WCSimRecoSummary;


class WCSimFitterTree {
public:
	WCSimFitterTree();
	virtual ~WCSimFitterTree();

	void MakeTree();
	void SaveTree();
	void SetSaveFileName(TString saveName);
	TString GetSaveFileName() const;
	TTree * GetFitTree();
	TTree * GetTruthTree();
	TTree * GetRecoTree();
	TTree * GetGeoTree();

	void Fill(Int_t iEvent,
			  std::vector<WCSimLikelihoodTrack> bestFitTracks,
			  std::vector<WCSimLikelihoodTrack*> trueTracks,
			  Double_t minimum);
	void FillFitTrack(WCSimLikelihoodTrack track, Double_t twoLnL);
	void FillTrueTrack(WCSimLikelihoodTrack track);

private:

	void MakeRecoSummary(std::vector<WCSimLikelihoodTrack> bestFitTracks);

	TFile * fSaveFile;
	TString fSaveFileName;
	TTree * fTrueTree;
	TTree * fFitTree;
	TTree * fRecoSummaryTree;
	TTree * fWCSimTree;
	TTree * fGeoTree;

	// The best-fit variables
	Int_t fEvent;
	Int_t fFitTrackNum;
	Double_t f2LnL;
	Double_t fFitVtxX;
	Double_t fFitVtxY;
	Double_t fFitVtxZ;
	Double_t fFitVtxT;
	Double_t fFitDirTheta;
	Double_t fFitDirPhi;
	Double_t fFitEnergy;
	Int_t fFitPDG;

	// The truth variables:
	Int_t fTrueTrackNum;
	Double_t fTrueVtxX;
	Double_t fTrueVtxY;
	Double_t fTrueVtxZ;
	Double_t fTrueVtxT;
	Double_t fTrueDirTheta;
	Double_t fTrueDirPhi;
	Double_t fTrueEnergy;
	Int_t fTruePDG;

	WCSimRootGeom * fGeometry;
	WCSimRecoSummary * fRecoSummary;
	WCSimRootEvent * fWCSimRootEvent;


	std::vector<WCSimLikelihoodTrack> fFitTracks;
	std::vector<WCSimLikelihoodTrack> fTruthTracks;

};

#endif /* INCLUDE_WCSIMFITTERTREE_HH_ */
