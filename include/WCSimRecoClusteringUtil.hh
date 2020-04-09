#ifndef WCSIM_RECOCLUSTERINGUTIL_HH
#define WCSIM_RECOCLUSTERINGUTIL_HH

#include <vector>

class WCSimRecoDigit;

class WCSimRecoClusteringUtil
{
public:
	WCSimRecoClusteringUtil();
	WCSimRecoClusteringUtil(std::vector<WCSimRecoDigit *> digitVec, double qCut, double tCut, double dCut,
							unsigned int minSize);
	~WCSimRecoClusteringUtil();

	void SetDigitVector(std::vector<WCSimRecoDigit *> digitVec)
	{
		fDigitVec = digitVec;
	}

	void SetChargeCut(double qCut)
	{
		fChargeCut = qCut;
	}

	void SetTimeCut(double tCut)
	{
		fTimeCut = tCut;
	}

	void SetDistanceCut(double dCut)
	{
		fDistanceCut = dCut;
	}

	void SetMinSize(unsigned int minSize)
	{
		fMinSliceSize = minSize;
	}

	void PerformClustering();

	std::vector<std::vector<WCSimRecoDigit *>> GetFullSlicedDigits()
	{
		return fFullSlicedDigits;
	}

private:
	void Reset();

	// Zeroth step applies the charge cut.
	void ApplyChargeCut();

	// First step is to order all the digits by time and look for time windows.
	void SliceByTime();

	// Second step is to cluster spatially
	void SliceInSpace();

	// Spatial clustering
	bool AddDigitIfClose(WCSimRecoDigit *, std::vector<WCSimRecoDigit *> &);
	double GetDistanceBetweenDigits(WCSimRecoDigit *, WCSimRecoDigit *);

	unsigned int GetFirstUnclusteredDigit();
	bool AreAllDigitsClustered();

	// Member variables
	double fChargeCut;
	double fTimeCut;
	double fDistanceCut;
	unsigned int fMinSliceSize;

	std::vector<bool> fIsDigitClustered;

	std::vector<WCSimRecoDigit *> fDigitVec;
	std::vector<WCSimRecoDigit *> fChargeCutDigits;
	std::vector<std::vector<WCSimRecoDigit *>> fTimeSlicedDigits;
	std::vector<std::vector<WCSimRecoDigit *>> fFullSlicedDigits;
};

#endif
