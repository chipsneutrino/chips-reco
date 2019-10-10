/*
 * WCSimImageMaker.hh
 *
 *  Created on: 5 November 2018
 *      Author: jtingey
 */

#include "TObject.h"

#include <iostream>
#include <cstdlib>

#ifndef WCSIMIMAGEMAKER_HH_
#define WCSIMIMAGEMAKER_HH_

class WCSimImageMaker: public TObject {
	public:
		WCSimImageMaker();
		~WCSimImageMaker();

	private:
		ClassDef(WCSimImageMaker,1);
};

#endif /* WCSIMIMAGEMAKER_HH_ */
