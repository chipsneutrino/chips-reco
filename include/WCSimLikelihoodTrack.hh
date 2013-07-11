#ifndef WCSIMLIKELIHOODTRACK_H
#define WCSIMLIKELIHOODTRACK_H

#include "TObject.h"
#include "TVector3.h"
#include <string>

class WCSimLikelihoodTrack : public TObject
{
    
    public:
        // To what kind of particle does the track correspond?
        enum TrackType
        {
            Unknown = 0,
            ElectronLike = 1,
            MuonLike = 2
        };
        
        std::string TrackTypeToString( WCSimLikelihoodTrack::TrackType myType );


    
        WCSimLikelihoodTrack();
        WCSimLikelihoodTrack( double x, double y, double z, double t, double theta, double phi, double E, WCSimLikelihoodTrack::TrackType type );
        virtual ~WCSimLikelihoodTrack();
		// Setters
		void SetX(double x);
		void SetY(double y);
		void SetZ(double z);
		void SetT(double t);
		void SetTheta(double th);
		void SetPhi(double phi);
		void SetE(double E);
        void SetType(WCSimLikelihoodTrack::TrackType type);
		
		// Getters
		const double GetX();
		const double GetY();
		const double GetZ();
        const TVector3 GetVtx();
		const double GetT();
		const double GetTheta();
		const double GetPhi();
        const double GetDirX();
        const double GetDirY();
        const double GetDirZ();
        const TVector3 GetDir();
		const double GetE();
        const WCSimLikelihoodTrack::TrackType GetType();
 
    protected:
    private:
			double fVtx[3];	// vertex position
			double fT0;			// vertex time
			double fTheta0;	// polar angle to z axis
			double fPhi0;		// azimuthal angle to z axis
			double fE0;			// kinetic energy
      WCSimLikelihoodTrack::TrackType fType;




		ClassDef(WCSimLikelihoodTrack,1)
};

#endif

