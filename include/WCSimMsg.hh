#ifndef WCSIMMSG_HH
#define WCSIMMSG_HH

#include "TObject.h"

class WCSimMsg : public TObject
{

public:
	typedef enum EMsgLevel
	{
		kVerbose = 0,
		kDebug = 1,
		kNotice = 2,
		kSynopsis = 3,
		kInfo = 5,
		kWarning = 7,
		kError = 8,
		kFatal = 9
	} MsgLevel_t;

	static WCSimMsg *Instance();

	static void SetLevel(WCSimMsg::MsgLevel_t minLevel);

	void SetMsgLevel(WCSimMsg::MsgLevel_t minLevel)
	{
		fMinLevel = minLevel;
	}

	Bool_t CheckMsgLevel(WCSimMsg::MsgLevel_t minLevel)
	{
		if (minLevel >= fMinLevel)
			return 1;
		else
			return 0;
	}

private:
	WCSimMsg();
	~WCSimMsg();

	WCSimMsg::MsgLevel_t fMinLevel;

	ClassDef(WCSimMsg, 0)
};

#endif
