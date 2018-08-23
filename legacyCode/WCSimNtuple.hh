#ifndef WCSIMNTUPLE_HH
#define WCSIMNTUPLE_HH

#include "TObject.h"
#include "TString.h"

#include "WCSimRecoEvent.hh"
#include "WCSimTrueEvent.hh"

class WCSimNtuple: public TObject {
	public:
		WCSimNtuple();
		virtual ~WCSimNtuple();

		virtual void Fill(WCSimTrueEvent* trueEvent, WCSimRecoEvent* recoEvent) {
		}

		virtual void OpenFile(const char* filename) {
		}
		virtual void CloseFile() {
		}

		virtual void Reset() {
		}

		void SetFileName(const char* filename) {
			fNtpFileName = filename;
		}

	protected:
		TString fNtpFileName;

	private:
		ClassDef(WCSimNtuple,0)

};

#endif
