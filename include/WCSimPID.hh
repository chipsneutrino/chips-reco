/*
 * WCSimPID.hh
 *
 *  Created on: 3 March 2018
 *      Author: jtingey
 */

#include "TObject.h"

#include <iostream>
#include <cstdlib>

#ifndef WCSIMPID_HH_
#define WCSIMPID_HH_

class WCSimPID : public TObject
{
public:
	WCSimPID();
	~WCSimPID();

	void Echo();

	ClassDef(WCSimPID,1);
};

#endif /* WCSIMPID_HH_ */
