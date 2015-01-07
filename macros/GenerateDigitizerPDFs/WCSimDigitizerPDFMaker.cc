/*
 * WCSimDigitizerPDFMaker.cpp
 *
 *  Created on: 5 Jan 2015
 *      Author: andy
 */

#include "WCSimDigitizerPDFMaker.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

#include <algorithm>
#include <iostream>

WCSimDigitizerPDFMaker::WCSimDigitizerPDFMaker() {

	// TODO Auto-generated constructor stub
	fRandom = new TRandom3();
	fRandom->SetSeed(0);

	BuildQPE0();

	fMu = 0;
	fNumThrows = 1000000;
	fEfficiency = 0.985; // Match to "efficiency" in WCSimWCDigitizer.cc

	fProbHisto = 0x0;
	SetBinning();
}

WCSimDigitizerPDFMaker::~WCSimDigitizerPDFMaker() {
	// TODO Auto-generated destructor stub
	if(fProbHisto != 0x0) { delete fProbHisto; }
	delete fRandom;
}

void WCSimDigitizerPDFMaker::SetBinning(int bins, double min, double max) {
	fNumChargeBins = bins;
	fChargeMin = min;
	fChargeMax = max;
}

int WCSimDigitizerPDFMaker::GetBins() const {
	return fNumChargeBins;
}

double WCSimDigitizerPDFMaker::GetMin() const {
	return fChargeMin;
}

double WCSimDigitizerPDFMaker::GetMax() const {
	return fChargeMax;
}

void WCSimDigitizerPDFMaker::SetMu(double mu) {
	fMu = mu;
}

double WCSimDigitizerPDFMaker::GetMu() const {
	return fMu;
}

void WCSimDigitizerPDFMaker::SetNumThrows(int num) {
	fNumThrows = num;
}

int WCSimDigitizerPDFMaker::GetNumThrows() const {
	return fNumThrows;
}

void WCSimDigitizerPDFMaker::Run() {
	MakeHisto();
	LoopDigitize();
	NormHistogram();
	SaveHistogram();
}

int WCSimDigitizerPDFMaker::ThrowPoisson() {

  int myRandom = fRandom->Poisson(fMu);
  fDebug->Fill(myRandom);
	return myRandom;
}

double WCSimDigitizerPDFMaker::Rn1pe() {
	double random = fRandom->Rndm();
	double random2 = fRandom->Rndm();
	int i;
	for(i = 0; i < 501; i++){
		if (random <= fQPE0[i]) break;
	}
	if(i==500)
		random = fRandom->Rndm();

	return (double(i-50) + random2)/22.83;
}

void WCSimDigitizerPDFMaker::BuildQPE0() {

	double tempArr[501] = {
		// 1
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
		0.000000, 0.000129, 0.000754, 0.004060, 0.028471,
		// 2
		0.068449, 0.115679, 0.164646, 0.203466, 0.235631,
		0.262351, 0.282064, 0.303341, 0.320618, 0.338317,
		0.357825, 0.371980, 0.385820, 0.398838, 0.413595,
		0.428590, 0.444387, 0.461685, 0.482383, 0.502369,
		0.520779, 0.540011, 0.559293, 0.579354, 0.599337,
		0.619580, 0.639859, 0.659807, 0.679810, 0.699620,
		0.718792, 0.737382, 0.755309, 0.772042, 0.788232,
		0.803316, 0.817861, 0.831148, 0.844339, 0.855532,
		0.866693, 0.876604, 0.886067, 0.894473, 0.902150,
		0.909515, 0.915983, 0.922050, 0.927418, 0.932492,
		// 3
		0.936951, 0.940941, 0.944660, 0.948004, 0.951090,
		0.953833, 0.956576, 0.958886, 0.961134, 0.963116,
		0.964930, 0.966562, 0.968008, 0.969424, 0.970687,
		0.971783, 0.972867, 0.973903, 0.974906, 0.975784,
		0.976632, 0.977438, 0.978190, 0.978891, 0.979543,
		0.980124, 0.980666, 0.981255, 0.981770, 0.982227,
		0.982701, 0.983146, 0.983566, 0.983975, 0.984357,
		0.984713, 0.985094, 0.985404, 0.985739, 0.986049,
		0.986339, 0.986630, 0.986922, 0.987176, 0.987431,
		0.987655, 0.987922, 0.988173, 0.988414, 0.988639,
		// 4
		0.988856, 0.989065, 0.989273, 0.989475, 0.989662,
		0.989828, 0.990007, 0.990172, 0.990327, 0.990497,
		0.990645, 0.990797, 0.990981, 0.991135, 0.991272,
		0.991413, 0.991550, 0.991673, 0.991805, 0.991928,
		0.992063, 0.992173, 0.992296, 0.992406, 0.992514,
		0.992632, 0.992733, 0.992837, 0.992954, 0.993046,
		0.993148, 0.993246, 0.993354, 0.993458, 0.993549,
		0.993656, 0.993744, 0.993836, 0.993936, 0.994033,
		0.994134, 0.994222, 0.994307, 0.994413, 0.994495,
		0.994572, 0.994659, 0.994739, 0.994816, 0.994886,
		// 5
		0.994970, 0.995032, 0.995110, 0.995178, 0.995250,
		0.995321, 0.995383, 0.995464, 0.995532, 0.995609,
		0.995674, 0.995750, 0.995821, 0.995889, 0.995952,
		0.996010, 0.996071, 0.996153, 0.996218, 0.996283,
		0.996335, 0.996384, 0.996431, 0.996484, 0.996537,
		0.996597, 0.996655, 0.996701, 0.996745, 0.996802,
		0.996860, 0.996917, 0.996962, 0.997014, 0.997079,
		0.997114, 0.997165, 0.997204, 0.997250, 0.997295,
		0.997335, 0.997379, 0.997418, 0.997454, 0.997488,
		0.997530, 0.997573, 0.997606, 0.997648, 0.997685,
		// 6
		0.997725, 0.997762, 0.997795, 0.997835, 0.997866,
		0.997898, 0.997941, 0.997966, 0.997997, 0.998039,
		0.998065, 0.998104, 0.998128, 0.998153, 0.998179,
		0.998205, 0.998223, 0.998254, 0.998293, 0.998319,
		0.998346, 0.998374, 0.998397, 0.998414, 0.998432,
		0.998456, 0.998482, 0.998511, 0.998532, 0.998553,
		0.998571, 0.998594, 0.998614, 0.998638, 0.998669,
		0.998693, 0.998715, 0.998743, 0.998762, 0.998793,
		0.998812, 0.998834, 0.998857, 0.998872, 0.998888,
		0.998904, 0.998926, 0.998946, 0.998963, 0.998983,
		// 7
		0.999007, 0.999027, 0.999044, 0.999064, 0.999079,
		0.999096, 0.999120, 0.999133, 0.999152, 0.999160,
		0.999174, 0.999188, 0.999206, 0.999221, 0.999234,
		0.999248, 0.999263, 0.999276, 0.999286, 0.999300,
		0.999313, 0.999321, 0.999331, 0.999347, 0.999356,
		0.999369, 0.999381, 0.999394, 0.999402, 0.999415,
		0.999427, 0.999433, 0.999446, 0.999458, 0.999472,
		0.999484, 0.999499, 0.999513, 0.999522, 0.999532,
		0.999540, 0.999550, 0.999559, 0.999567, 0.999574,
		0.999588, 0.999599, 0.999613, 0.999618, 0.999627,
		// 8
		0.999635, 0.999639, 0.999652, 0.999662, 0.999667,
		0.999671, 0.999678, 0.999682, 0.999688, 0.999693,
		0.999698, 0.999701, 0.999706, 0.999711, 0.999718,
		0.999722, 0.999727, 0.999732, 0.999737, 0.999740,
		0.999746, 0.999750, 0.999754, 0.999763, 0.999766,
		0.999769, 0.999774, 0.999780, 0.999784, 0.999788,
		0.999796, 0.999803, 0.999807, 0.999809, 0.999815,
		0.999820, 0.999827, 0.999830, 0.999833, 0.999833,
		0.999836, 0.999839, 0.999842, 0.999845, 0.999850,
		0.999853, 0.999857, 0.999860, 0.999865, 0.999870,
		// 9
		0.999873, 0.999877, 0.999880, 0.999882, 0.999883,
		0.999886, 0.999888, 0.999889, 0.999895, 0.999896,
		0.999897, 0.999901, 0.999902, 0.999905, 0.999907,
		0.999907, 0.999909, 0.999911, 0.999911, 0.999912,
		0.999913, 0.999914, 0.999917, 0.999919, 0.999921,
		0.999923, 0.999927, 0.999929, 0.999931, 0.999933,
		0.999936, 0.999942, 0.999942, 0.999944, 0.999947,
		0.999947, 0.999948, 0.999949, 0.999952, 0.999955,
		0.999957, 0.999957, 0.999961, 0.999962, 0.999963,
		0.999963, 0.999963, 0.999964, 0.999965, 0.999965,
		// 10
		0.999965, 0.999965, 0.999966, 0.999968, 0.999969,
		0.999971, 0.999972, 0.999972, 0.999973, 0.999975,
		0.999975, 0.999975, 0.999975, 0.999975, 0.999975,
		0.999975, 0.999979, 0.999979, 0.999980, 0.999982,
		0.999983, 0.999985, 0.999986, 0.999987, 0.999987,
		0.999988, 0.999989, 0.999989, 0.999989, 0.999989,
		0.999990, 0.999990, 0.999992, 0.999993, 0.999994,
		0.999994, 0.999994, 0.999994, 0.999994, 0.999995,
		0.999995, 0.999995, 0.999996, 0.999996, 0.999996,
		0.999996, 0.999998, 0.999999, 1.000000, 1.000000,
		// Dummy element for noticing if the loop reached the end of the array
		0.0
	};
	std::copy(tempArr, tempArr + 501, fQPE0);
	return;
}

void WCSimDigitizerPDFMaker::LoopDigitize() {
	for(int iThrow = 0; iThrow < fNumThrows; ++iThrow)
	{
		Digitize();
	}
}

void WCSimDigitizerPDFMaker::Digitize() {

	// std::cout << "Throwing from a Poisson with mean " << fMu << std::endl;
	int nPhotons = this->ThrowPoisson();
	// std::cout << "Got " << nPhotons << " photons" << std::endl;
	int iflag = 0;
	double peSmeared = 0;
	for(int iPhoton = 0; iPhoton < nPhotons; ++iPhoton)
	{
		// std::cout << "Getting a random 1pe sample" << std::endl;
		peSmeared += Rn1pe();
		// std::cout << "Got it!" << std::endl;
	}
	// std::cout << "Now doing the threshold" << std::endl;
	this->Threshold(peSmeared,iflag);
	// std::cout << "Threshold complete" << std::endl;
	peSmeared *= fEfficiency; // MC tuning correction
  
  if( iflag == 0 )
  {
    double Q = (peSmeared > 0.5) ? peSmeared : 0.5;
    peSmeared = Q;
  }

	fProbHisto->Fill(fMu, peSmeared);
}

void WCSimDigitizerPDFMaker::NormHistogram() {
	double integral = fProbHisto->Integral();
	if(integral > 0)
	{
		fProbHisto->Scale(1.0/integral);
	}
	return;
}

void WCSimDigitizerPDFMaker::SaveHistogram() {
	fProbHisto->SetDirectory(0);
	fDebug->SetDirectory(0);
	TString fileName = TString::Format("digitizerLikelihood_%f.root", fMu);
	TFile * f = new TFile(fileName.Data(), "RECREATE");
	fProbHisto->Write();
  fDebug->Write();
	f->Close();
	delete f;
}

void WCSimDigitizerPDFMaker::Threshold(double& pe, int& iflag) {
	double x = pe+0.1;
	iflag=0;

	double thr;
	double RDUMMY;
	double err;
	if ( x<1.1) {
	  thr = std::min(1.0,
			 -0.06374+x*(3.748+x*(-63.23+x*(452.0+x*(-1449.0+x*(2513.0
									+x*(-2529.+x*(1472.0+x*(-452.2+x*(51.34+x*2.370))))))))));
	} else {
	  thr = 1.0;
	}
	RDUMMY = fRandom->Rndm();
	if (thr < RDUMMY) {
	  pe = 0.0;
	  iflag = 1;
	}
	else {
	  err = fRandom->Gaus(0.0,0.03);
	  /////      call rngaus(0.0, 0.03, err);
	  pe = pe+err;
	}
}

void WCSimDigitizerPDFMaker::MakeHisto()
{
  double binWidth = 0.5 * ( fChargeMax - fChargeMin ) / fNumChargeBins;
  std::cout << "binWidth = " << binWidth << std::endl;

	fProbHisto = new TH2D("fProbHisto","Digitizer PDF",
						   fNumChargeBins, fChargeMin-binWidth, fChargeMax-binWidth,
						   fNumChargeBins, fChargeMin-binWidth, fChargeMax-binWidth);
  fDebug = new TH1D("fDebug","Poisson picks",(int)(ceil(fChargeMax) - floor(fChargeMin)), fChargeMin - binWidth, fChargeMax - binWidth);
  fDebug->GetXaxis()->SetTitle("Q");
  fDebug->GetYaxis()->SetTitle("Events");
	fProbHisto->GetXaxis()->SetTitle("Predicted mean number of photons, #mu");
	fProbHisto->GetYaxis()->SetTitle("Digitized charge in p.e.");
	fProbHisto->GetZaxis()->SetTitle("P(digitized | mean photons)");
	return;
}
