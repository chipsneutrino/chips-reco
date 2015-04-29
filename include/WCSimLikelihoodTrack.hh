/**
 * \class WCSimLikelihoodTrack
 * This class represents a single charged particle
 * track.  It holds its initial position and direction
 * of propagation, as well as its energy and the particle's
 * type.
 */

#ifndef WCSIMLIKELIHOODTRACK_H
#define WCSIMLIKELIHOODTRACK_H

#include "TObject.h"
#include "TVector3.h"
#include "WCSimTrackParameterEnums.hh"
#include <string>
#include <cstdlib>
class WCSimTrueTrack;

class WCSimLikelihoodTrack : public TObject
{
    
    public:
        /**
         * To what kind of particle does the track correspond?
         */
        enum TrackType
        {
            Unknown = 0, //!< Track type is unknown
            ElectronLike = 1, //!< An electron-like track (could be a photon)
            MuonLike = 2 //!< A muon-like track
        };

        /**
         * Get the text name of a give track type
         * @param myType Track type
         * @return Name corresponding to the enum
         */
        static std::string TrackTypeToString( WCSimLikelihoodTrack::TrackType myType );

        static TrackType GetTypeFromPDG(int pdg)
        {
        	if( abs(pdg) == 11 ) { return ElectronLike; }
        	if( abs(pdg) == 13 ) { return MuonLike; }

        	return Unknown;
        }

        static TrackType GetTypeFromName( const char * name )
        {
        	std::string nameStr(name);
        	if( nameStr == "MuonLike" ){ return MuonLike; }
        	else if ( nameStr == "ElectronLike" ) { return ElectronLike; }

        	return Unknown;
        }
    
        WCSimLikelihoodTrack();
        /**
         * Constructo
         * @param x Track vertex x
         * @param y Track vertex y
         * @param z Track vertex z
         * @param t Track start time
         * @param theta Angle to z axis
         * @param phi Azimuthal angle
         * @param E Track energy
         * @param type Particle type
         */
        WCSimLikelihoodTrack( double x, double y, double z, double t,
                              double theta, double phi, double E,
                              WCSimLikelihoodTrack::TrackType type );

        /**
         * Constructor from true WCSim track object
         * @param trueTrack Track to convert into a WCSimLikelihoodTrack
         */
        WCSimLikelihoodTrack( WCSimTrueTrack * trueTrack);

        virtual ~WCSimLikelihoodTrack();
		    void Print();
            
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
		    double GetX() const;
		    double GetY() const;
		    double GetZ() const;
        TVector3 GetVtx() const;
		    double GetT() const;
		    double GetTheta() const;
		    double GetPhi() const;
        double GetDirX() const;
        double GetDirY() const;
        double GetDirZ() const;
        TVector3 GetDir() const;
		    double GetE() const;
        WCSimLikelihoodTrack::TrackType GetType() const;
        double GetTrackParameter(FitterParameterType::Type type) const;

        static bool EnergyGreaterThanOrEqual(const WCSimLikelihoodTrack &a, const WCSimLikelihoodTrack &b);
        static bool EnergyGreaterThanOrEqualPtrs(WCSimLikelihoodTrack *a, WCSimLikelihoodTrack *b);
        /**
         * Function to calculate the position a distance s away from the track vertex in
         * the direction of travel of the track
         * @param s Distance to propagate from the vertex
         * @return The position of the particle after it propagates a distance s from
         *         the vertex
         */
        TVector3 GetPropagatedPos(Double_t s) const;
        Int_t GetPDG() const;
        bool operator == (const WCSimLikelihoodTrack &b) const;
 
    protected:
    private:
      void LoadEmissionProfile(WCSimLikelihoodTrack * myTrack);
			double fVtx[3];	    ///< Vertex position
			double fT0;			///< Vertex time
			double fTheta0;	    ///< Polar angle to z axis
			double fPhi0;		///< Azimuthal angle
			double fE0;			///< Kinetic energy
			WCSimLikelihoodTrack::TrackType fType; ///< Particle type (e, mu, pi...)




		ClassDef(WCSimLikelihoodTrack,1)
};

#endif

