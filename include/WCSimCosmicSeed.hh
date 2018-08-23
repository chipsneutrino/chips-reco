#ifndef WCSIMCOSMICSEED_HH
#define WCSIMCOSMICSEED_HH

#include <vector>

#include <TVector3.h>

class WCSimRecoDigit;
class WCSimRecoClusteringUtil;

class WCSimCosmicSeed {

	public:

		WCSimCosmicSeed();
		WCSimCosmicSeed(std::vector<WCSimRecoDigit*>* fVetoDigitArray);
		~WCSimCosmicSeed();

		void SetVetoDigits(std::vector<WCSimRecoDigit*>* digits) {
			fVetoDigits = digits;
		}

		void SetInnerDigits(std::vector<WCSimRecoDigit*>* digits) {
			fInnerDigits = digits;
		}

		TVector3 GetFittedVtx() const {
			return fFittedVtx;
		}

		TVector3 GetFittedDir() const {
			return fFittedDir;
		}

		double GetFittedVtxT() const {
			return fFittedVtxT;
		}

		bool GetSuccess() const {
			return fSuccess;
		}

		void CalcSeedVtxAndDir();

		bool HasVetoDigits();

		// If we have a seed, check to see if a given position is inside so that we can
		// use the seed to shadow a region of the detector.
		bool IsHitShadowed(double hx, double hy, double hz);

	private:

		void RunSimpleEntrySeed();
		void RunSimpleCrossingSeed();
		void GetSimpleDirectionFromInner();

		void RunClusteringSeed();

		// Generate the ring points
		void GenerateRingPoints();

		WCSimRecoDigit* GetBiggestClusterHit(std::vector<WCSimRecoDigit*> clust);

		std::vector<WCSimRecoDigit*>* fVetoDigits;
		std::vector<WCSimRecoDigit*>* fInnerDigits;

		TVector3 fFittedVtx;
		TVector3 fFittedDir;
		double fFittedVtxT;

		bool fSuccess;

		// Store the 360 points that form the ring outline
		std::vector<double> fRingPosX;
		std::vector<double> fRingPosY;
		std::vector<double> fRingPosZ;

		double fMaxDistBetweenRingPoints;

		WCSimRecoClusteringUtil *fClusteringUtil;

		ClassDef(WCSimCosmicSeed,0)

};

#endif

