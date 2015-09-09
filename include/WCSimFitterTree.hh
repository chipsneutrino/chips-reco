/*
 * WCSimFitterTree.hh
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#ifndef INCLUDE_WCSIMFITTERTREE_HH_
#define INCLUDE_WCSIMFITTERTREE_HH_

#include "WCSimLikelihoodTrackBase.hh"
#include <vector>
class TFile;
class TString;
class TTree;
class WCSimRootGeom;
class WCSimLikelihoodRecoEvent;
class WCSimLikelihoodDigitArray;
class WCSimTrueEvent;
class WCSimRootEvent;
class WCSimRecoSummary;

class HitComparison {
  public:
    HitComparison() : pmtID(-999), pmtX(-999.9), pmtY(-999.9), pmtZ(-999.9), trueT(-999.9), predT(-999.9), trueQ(-999.9), predQ(-999.9), correctPredQ(-999.9), minus2LnL(-999.9), correctMinus2LnL(-999.9), charge2LnL(-999.9), time2LnL(-999.9) {};
    HitComparison(int pmt, double x, double y, double z, double truthT, double predictT, double truthQ, double predictQ, double correctPred, double m2LnL, double correctm2LnL, double q2LnL, double t2LnL) :
    pmtID(pmt), pmtX(x), pmtY(y), pmtZ(z), trueT(truthT), predT(predictT), trueQ(truthQ), predQ(predictQ), correctPredQ(correctPred), minus2LnL(m2LnL), correctMinus2LnL(correctm2LnL), charge2LnL(q2LnL), time2LnL(t2LnL) {};
    void Set(int pmt, double x, double y, double z, double truthT, double predictT, double truthQ, double predictQ, double correctPred, double m2LnL, double correctm2LnL, double q2LnL, double t2LnL) {
      pmtID = pmt;
      pmtX = x;
      pmtY = y;
      pmtZ = z;
      trueQ = truthQ;
      predQ = predictQ;
      trueT = truthT;
      predT = predictT;
      correctPredQ = correctPred;
      minus2LnL = m2LnL;
      correctMinus2LnL = correctm2LnL;
      charge2LnL = q2LnL;
      time2LnL = t2LnL;
    };

    virtual ~HitComparison(){};

    int pmtID;
    double pmtX;
    double pmtY;
    double pmtZ;
    double trueT;
    double predT;
    double trueQ;
    double predQ;
    double correctPredQ;
    double minus2LnL;
    double correctMinus2LnL;
    double charge2LnL;
    double time2LnL;
};

class WCSimFitterTree {
public:
	WCSimFitterTree(const TString &saveFileName);
	virtual ~WCSimFitterTree();

	void MakeTree();
	void SaveTree();
	void SetSaveFileName(TString saveName);
	TString GetSaveFileName() const;
	TTree * GetFitTree();
	TTree * GetTruthTree();
	TTree * GetRecoTree();
	TTree * GetGeoTree();
  TTree * GetHitComparisonTree();

	void Fill(Int_t iEvent,
			  std::vector<WCSimLikelihoodTrackBase*> bestFitTracks,
			  std::vector<Bool_t> bestFitEscapes,
			  std::vector<WCSimLikelihoodTrackBase*> trueTracks,
			  std::vector<Bool_t> trueTrackEscapes,
			  Double_t minimum);
	void FillFitTrack(WCSimLikelihoodTrackBase * track, Bool_t escapes, Double_t twoLnL);
	void FillTrueTrack(WCSimLikelihoodTrackBase * track, Bool_t escapes);
	void FillHitComparison(
               const int &event, 
               WCSimLikelihoodDigitArray * digitArray,
						   const std::vector<double> &predictedCharges,
			   	   	   	   const std::vector<double> &correctPredictedCharges,
						   const std::vector<double> &measuredCharges,
               const std::vector<double> &predictedTimes,
						   const std::vector<double> &total2LnLs,
						   const std::vector<double> &correct2LnLs,
               const std::vector<double> &charge2LnLs,
               const std::vector<double> &time2LnLs);
  void FillRecoFailures(const int &event);

private:
  	  void MakeSaveFileName();
	void MakeRecoSummary(std::vector<WCSimLikelihoodTrackBase *> bestFitTracks);

	TFile * fSaveFile;
	TString fSaveFileName;
	TTree * fTrueTree;
	TTree * fFitTree;
	TTree * fRecoSummaryTree;
	TTree * fWCSimTree;
	TTree * fGeoTree;
	TTree * fComparisonTree;
  TTree * fHitComparisonTree;
  TTree * fRecoFailureTree;

  // Failed reco events
  int fFailedEvent;

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
  Double_t fFitConversionDistance;
	Bool_t fFitEscapes;
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
	Bool_t fTrueEscapes;
	Int_t fTruePDG;

	WCSimRootGeom * fGeometry;
	WCSimRecoSummary * fRecoSummary;
	WCSimRootEvent * fWCSimRootEvent;


	std::vector<WCSimLikelihoodTrackBase *> fFitTracks;
	std::vector<WCSimLikelihoodTrackBase *> fTruthTracks;
  HitComparison * fHitComparison;

};

#endif /* INCLUDE_WCSIMFITTERTREE_HH_ */
