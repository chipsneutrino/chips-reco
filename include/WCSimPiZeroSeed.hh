/*
 * WCSimPiZeroSeed.hh
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#ifndef WCSIMPIZEROSEED_HH_
#define WCSIMPIZEROSEED_HH_
class WCSimLikelihoodTrackBase;
#include "TObject.h"

class WCSimPiZeroSeed: public TObject {
	public:
		WCSimPiZeroSeed(WCSimLikelihoodTrackBase * track1, WCSimLikelihoodTrackBase* track2, double minus2LnL);
		WCSimPiZeroSeed(const WCSimPiZeroSeed &other);
		WCSimPiZeroSeed& operator =(const WCSimPiZeroSeed &rhs);
		virtual ~WCSimPiZeroSeed();

		WCSimLikelihoodTrackBase* GetTrack1() const;
		WCSimLikelihoodTrackBase* GetTrack2() const;
		double GetMinus2LnL() const;

		void Print();
		friend bool operator <(const WCSimPiZeroSeed &lhs, const WCSimPiZeroSeed &rhs);
		friend bool operator >(const WCSimPiZeroSeed &lhs, const WCSimPiZeroSeed &rhs);
		friend bool operator <=(const WCSimPiZeroSeed &lhs, const WCSimPiZeroSeed &rhs);
		friend bool operator >=(const WCSimPiZeroSeed &lhs, const WCSimPiZeroSeed &rhs);

	private:
		WCSimLikelihoodTrackBase * fTrack1;
		WCSimLikelihoodTrackBase * fTrack2;

		double fMinus2LnL;ClassDef(WCSimPiZeroSeed,0);
};

#endif /* WCSIMPIZEROSEED_HH_ */
