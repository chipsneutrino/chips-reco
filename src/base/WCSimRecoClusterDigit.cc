#include "WCSimRecoClusterDigit.hh"

#include "WCSimRecoObjectTable.hh"

ClassImp(WCSimRecoClusterDigit)

WCSimRecoClusterDigit::WCSimRecoClusterDigit(WCSimRecoDigit *recoDigit)
{
	fIsClustered = 0;
	fIsAllClustered = 0;

	fRecoDigit = recoDigit;

	fClusterDigitList = new std::vector<WCSimRecoClusterDigit *>;

	WCSimRecoObjectTable::Instance()->NewClusterDigit();
}

WCSimRecoClusterDigit::~WCSimRecoClusterDigit()
{
	delete fClusterDigitList;

	WCSimRecoObjectTable::Instance()->DeleteClusterDigit();
}

void WCSimRecoClusterDigit::AddClusterDigit(WCSimRecoClusterDigit *clusterDigit)
{
	fClusterDigitList->push_back(clusterDigit);
}

Int_t WCSimRecoClusterDigit::GetNClusterDigits()
{
	return fClusterDigitList->size();
}

WCSimRecoClusterDigit *WCSimRecoClusterDigit::GetClusterDigit(Int_t idigit)
{
	return (WCSimRecoClusterDigit *)(fClusterDigitList->at(idigit));
}

std::vector<WCSimRecoClusterDigit *> *WCSimRecoClusterDigit::GetClusterDigitList()
{
	return fClusterDigitList;
}

Bool_t WCSimRecoClusterDigit::IsAllClustered()
{
	if (fIsAllClustered == 0)
	{
		fIsAllClustered = 1;
		for (UInt_t n = 0; n < fClusterDigitList->size(); n++)
		{
			WCSimRecoClusterDigit *clusterDigit = (WCSimRecoClusterDigit *)(fClusterDigitList->at(n));
			if (clusterDigit->IsClustered() == 0)
				fIsAllClustered = 0;
		}
	}

	return fIsAllClustered;
}
