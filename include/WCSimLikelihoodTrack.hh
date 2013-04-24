#ifndef WCSIMLIKELIHOODTRACK_H
#define WCSIMLIKELIHOODTRACK_H

#include "TObject.h"

class WCSimLikelihoodTrack : public TObject
{
    public:
        WCSimLikelihoodTrack();
        WCSimLikelihoodTrack( double x, double y, double z, double t, double theta, double phi, double E );
        virtual ~WCSimLikelihoodTrack();
		// Setters
		void SetX(double x);
		void SetY(double y);
		void SetZ(double z);
		void SetT(double t);
		void SetTheta(double th);
		void SetPhi(double phi);
		void SetE(double E);
		
		// Getters
		const double GetX();
		const double GetY();
		const double GetZ();
		const double GetT();
		const double GetTheta();
		const double GetPhi();
		const double GetE();
 
    protected:
    private:
			double fVtx[3];	// vertex position
			double fT0;			// vertex time
			double fTheta0;	// polar angle to z axis
			double fPhi0;		// azimuthal angle to z axis
			double fE0;			// kinetic energy




		ClassDef(WCSimLikelihoodTrack,1)
};

#endif

