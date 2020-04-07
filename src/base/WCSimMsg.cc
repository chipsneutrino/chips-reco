#include "WCSimMsg.hh"

#include <cassert>

ClassImp (WCSimMsg)

static WCSimMsg* fgMsg = 0;

WCSimMsg* WCSimMsg::Instance() {
	if (!fgMsg) {
		fgMsg = new WCSimMsg();
	}

	if (!fgMsg) {
		assert (fgMsg);
	}

	if (fgMsg) {

	}

	return fgMsg;
}

WCSimMsg::WCSimMsg() {
	fMinLevel = WCSimMsg::kVerbose;
}

WCSimMsg::~WCSimMsg() {

}

void WCSimMsg::SetLevel(WCSimMsg::MsgLevel_t minLevel) {
	return WCSimMsg::Instance()->SetMsgLevel(minLevel);
}

