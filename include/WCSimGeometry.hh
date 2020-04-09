#pragma once

#include "TObject.h"
#include "TVector3.h"

#include "WCSimRootGeom.hh"
#include "WCSimTrackParameterEnums.hh"
class WCSimEmissionProfileManager;

class WCSimGeometry : public TObject
{

public:
	typedef enum EGeoConfiguration
	{
		kDUSEL_Unknown = 0,
		kDUSEL_100kton_10inch_12percent = 1001012,
		kDUSEL_100kton_10inch_15percent = 1001015,
		kDUSEL_100kton_10inch_30percent = 1001030,
		kDUSEL_100kton_10inch_40percent = 1001040,
		kDUSEL_200kton_10inch_12percent = 2001012,
		kDUSEL_200kton_12inch_10percent = 2001210,
		kDUSEL_200kton_12inch_14percent = 2001214,
		kDUSEL_150ktonMailbox_10inch_15percent = 101501015,
		kDUSEL_150ktonMailbox_10inch_30percent = 101501030,
		kDUSEL_150ktonMailbox_10inch_40percent = 101501040
	} GeoConfiguration_t;

	typedef enum EGeoType
	{
		kUnknown = -1,
		kCylinder = 0,
		kMailBox = 1
	} GeoType_t;

	typedef enum EGeoRegion
	{
		kTop = 0,
		kSide = 1,
		kBottom = 2,
		kVeto = 3,
		kFront = 10,
		kBack = 12,
		kLeft = 20,
		kRight = 22
	} GeoRegion_t;

	typedef enum ERingRegion
	{
		kOutside = 0,
		kInside = 1,
		kInsideHole = 2
	} RingRegion_t;

	const char *AsString(GeoConfiguration_t config);

	static WCSimGeometry *Instance();
	static void BuildGeometry(WCSimRootGeom *geom);
	static void PrintGeometry();
	static void WriteGeometry(const char *filename = "wcsim.geometry.root");
	static Bool_t TouchGeometry();
	static void Reset();

	void SetGeometry(WCSimRootGeom *geom);
	void WriteToFile(const char *filename = "wcsim.geometry.root");

	// Lookup method for WCSimGeometry
	// ===============================
	WCSimRootGeom *GetWCSimGeometry()
	{
		return fWCSimRootGeom;
	}

	// Lookup method for Configuration
	// ===============================
	Int_t GetGeoConfig()
	{
		return fGeoConfig;
	}

	// Lookup methods for Geometry
	// ===========================
	Int_t GetGeoType()
	{
		return fGeoType;
	}

	Bool_t IsCylinder()
	{
		return fCylinder;
	}
	Double_t GetCylRadius()
	{
		return fCylRadius;
	}
	Double_t GetCylLength()
	{
		return fCylLength;
	}
	Double_t GetCylFiducialRadius()
	{
		return fCylFiducialRadius;
	}
	Double_t GetCylFiducialLength()
	{
		return fCylFiducialLength;
	}

	Bool_t IsMailBox()
	{
		return fMailBox;
	}
	Double_t GetMailBoxX()
	{
		return fMailBoxX;
	}
	Double_t GetMailBoxY()
	{
		return fMailBoxY;
	}
	Double_t GetMailBoxZ()
	{
		return fMailBoxZ;
	}
	Double_t GetMailBoxFiducialX()
	{
		return fMailBoxFiducialX;
	}
	Double_t GetMailBoxFiducialY()
	{
		return fMailBoxFiducialY;
	}
	Double_t GetMailBoxFiducialZ()
	{
		return fMailBoxFiducialZ;
	}

	Double_t GetArea()
	{
		return fDetectorArea;
	}
	Double_t GetVolume()
	{
		return fDetectorVolume;
	}
	Double_t GetFiducialVolume()
	{
		return fDetectorFiducialVolume;
	}

	// Lookup methods for PMTs
	// =======================
	Int_t GetNumPMTs()
	{
		return fPMTs;
	}
	Int_t GetNumVetoPMTs()
	{
		return fVetoPMTs;
	}
	Int_t GetNumInnerPMTs()
	{
		return fPMTs - fVetoPMTs;
	}
	Double_t GetPMTRadius()
	{
		return fPMTRadius;
	}
	Double_t GetPMTCoverage()
	{
		return fPMTCoverage;
	}
	Double_t GetPMTSeparation()
	{
		return fPMTSeparation;
	}

	Int_t GetRegion(Int_t tube);
	Double_t GetX(Int_t tube);
	Double_t GetY(Int_t tube);
	Double_t GetZ(Int_t tube);

	// Calculate Distance to Edge of Detector
	// ======================================
	Bool_t InsideDetector(Double_t x, Double_t y, Double_t z);
	Bool_t InsideFiducialVolume(Double_t x, Double_t y, Double_t z);
	Bool_t InsideDetector(Double_t vx, Double_t vy, Double_t vz, Double_t ex, Double_t ey, Double_t ez);

	Double_t DistanceToEdge(Double_t x, Double_t y, Double_t z);

	void ProjectToNearEdge(Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py, Double_t pz,
						   Double_t &x, Double_t &y, Double_t &z, Int_t &region);
	void ProjectToFarEdge(Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py, Double_t pz, Double_t &x,
						  Double_t &y, Double_t &z, Int_t &region);

	Double_t ForwardProjectionToEdge(Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz);
	Double_t BackwardProjectionToEdge(Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz);
	void ProjectMinMaxPosition(TVector3 vtx, TVector3 dir, Double_t &sMin, Double_t &sMax);

	// Projection to Edge of Detector
	// ==============================
	// Input: A position (x0,y0,z0) and direction (px,py,pz)
	// Output: Projection on edge of detector (x,y,z,region)
	//  [Note: if (x0,y0,z0) is outside detector, then track
	//    will pass through both edges of the detector.
	//    If UseFarEdge=0, project forwards to nearest edge.
	//    If UseFarEdge=1, project forwards to furthest edge.]

	// Second attempt, using vectors to simplify algebra
	void ProjectToEdge(Bool_t useFarEdge, Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py,
					   Double_t pz, Double_t &x, Double_t &y, Double_t &z, Int_t &region);

	// First attempt, containing lengthy algebra...
	void ProjectToEdgeOld(Bool_t UseFarEdge, Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py,
						  Double_t pz, Double_t &x, Double_t &y, Double_t &z, Int_t &region);

	// Coordinates for 'Roll-out' Display
	// ==================================
	// Input: A detector hit (x,y,z,region)
	// Output: Its position on the 'roll-out' display (u,v)
	void XYZtoUV(Int_t region, Double_t x, Double_t y, Double_t z, Double_t &u, Double_t &v);

	// Generic Geometry: Reconstruct a Circle from Three Hits
	// ======================================================
	// Input: Three hits (x,y,z) {0,1,2}
	// Output: Their common circle
	//  Centre=(rx,ry,rz), Normal=(nx,ny,nz), Radius=(r)
	static void FindCircle(Double_t x0, Double_t y0, Double_t z0, Double_t x1, Double_t y1, Double_t z1,
						   Double_t x2, Double_t y2, Double_t z2, Double_t &rx, Double_t &ry, Double_t &rz, Double_t &nx,
						   Double_t &ny, Double_t &nz, Double_t &r);

	// Generic Geometry: Calculate Positions of Cherenkov Hits
	// =======================================================
	// Input: A point along a track
	//  Position of point along track: (x0,y0,z0)
	//  Projected position along track [edge of detector]: (xp,yp,zp)
	//  Cherenkov Angle: angle_degrees
	//  Azimuthal Angle: omega_degrees
	// Output: A point on the Cherenkov Ring
	//  Direction to point on Cherenkov ring: (nx,ny,nz)
	//  Rotated position: (rx,ry,rz); distance moved: r
	//   [The rotated position (rx,ry,rz) is calculated by rotating
	//    the input position (xp,yp,zp) through the two angles.]

	// Second attempt, using ROOT TVector3 for rotations
	static void FindCircle(Double_t xp, Double_t yp, Double_t zp, Double_t x0, Double_t y0, Double_t z0,
						   Double_t angle_degrees, Double_t omega_degrees, Double_t &rx, Double_t &ry, Double_t &rz, Double_t &nx,
						   Double_t &ny, Double_t &nz, Double_t &r);

	// First attempt, calculating all rotations by hand
	static void FindCircleOld(Double_t xp, Double_t yp, Double_t zp, Double_t x0, Double_t y0, Double_t z0,
							  Double_t angle_degrees, Double_t omega_degrees, Double_t &rx, Double_t &ry, Double_t &rz, Double_t &nx,
							  Double_t &ny, Double_t &nz, Double_t &r);

	// Generic Geometry: Reconstruct a Vertex from Four Hits
	// =====================================================
	// Input: Four hits (x,y,z,t) {0,1,2,3}
	// Output: Their common vertex (vx,vy,vz,vt) {m,p}
	//  [Note: this involves solving a quadratic, with two solutions.
	//    In some cases, both solutions appear to be valid.]

	// Use ROOT TMatrix to invert the 4x4 matrix
	//  [Take advantage of ROOT algorithm, as analytic solution is slow.]
	static void FindVertex(Double_t x0, Double_t y0, Double_t z0, Double_t t0, Double_t x1, Double_t y1,
						   Double_t z1, Double_t t1, Double_t x2, Double_t y2, Double_t z2, Double_t t2, Double_t x3, Double_t y3,
						   Double_t z3, Double_t t3, Double_t &vxm, Double_t &vym, Double_t &vzm, Double_t &vtm, Double_t &vxp,
						   Double_t &vyp, Double_t &vzp, Double_t &vtp);

	// Cherenkov Geometry: Project PMT Hit back to Point of Cherenkov Emission
	// ========================================================================
	// Input: PMT Hit      (x0,y0,z0),
	//        Track Vertex (vx,vy,vz),
	//        Track End    (ex,ey,ez)
	// Output: Intersection (x,y,z)
	//         Photon Distance L
	// [Note: use Sine Rule for Cherenkov Geometry]
	static void DistanceToIntersectLine(Double_t x0, Double_t y0, Double_t z0, Double_t vx, Double_t vy,
										Double_t vz, Double_t ex, Double_t ey, Double_t ez, Double_t &x, Double_t &y, Double_t &z, Double_t &L);
	static Double_t DistanceToIntersectLine(Double_t x0, Double_t y0, Double_t z0, Double_t sx, Double_t sy,
											Double_t sz, Double_t ex, Double_t ey, Double_t ez, Double_t &x, Double_t &y, Double_t &z);
	static Double_t DistanceToIntersectLine(Double_t *pos, Double_t *start, Double_t *end, Double_t *intersection);

	bool IsOutsideRing(const TVector3 &vtx, const TVector3 &dir, const TVector3 &pmt, const TrackType::Type &type,
					   const double energy);
	bool IsInRing(const TVector3 &vtx, const TVector3 &dir, const TVector3 &pmt, const TrackType::Type &type,
				  const double energy);
	bool IsInRingHole(const TVector3 &vtx, const TVector3 &dir, const TVector3 &pmt, const TrackType::Type &type,
					  const double energy);
	bool GetEscapes(const TVector3 &vtx, const TVector3 &dir, const TrackType::Type &type, const double energy);
	TVector3 GetEnd(const TVector3 &vtx, const TVector3 &dir, const TrackType::Type &type, const double energy);
	RingRegion_t GetPMTRingRegion(const TVector3 &vtx, const TVector3 &dir, const TVector3 &pmt,
								  const TrackType::Type &type, const double energy);

private:
	WCSimGeometry();
	WCSimGeometry(WCSimRootGeom *geom);
	~WCSimGeometry();
	void TouchEmissionProfileManager();

	Double_t fScale;

	Int_t fPMTs;
	Int_t fVetoPMTs;
	Double_t *fPmtX;
	Double_t *fPmtY;
	Double_t *fPmtZ;
	Int_t *fPmtRegion;

	Double_t fPMTRadius;
	Double_t fPMTCoverage;
	Double_t fPMTSeparation;
	Double_t fPMTSurfaceArea;

	Int_t fGeoType;
	Int_t fGeoTypeInput;
	Int_t fGeoConfig;

	Bool_t fCylinder;
	Double_t fCylRadius;
	Double_t fCylLength;
	Double_t fCylFiducialRadius;
	Double_t fCylFiducialLength;

	Bool_t fMailBox;
	Double_t fMailBoxX;
	Double_t fMailBoxY;
	Double_t fMailBoxZ;
	Double_t fMailBoxFiducialX;
	Double_t fMailBoxFiducialY;
	Double_t fMailBoxFiducialZ;

	Double_t fDetectorArea;
	Double_t fDetectorVolume;
	Double_t fDetectorFiducialVolume;

	Int_t fTube;
	Int_t fRegion;
	Double_t fXpos;
	Double_t fYpos;
	Double_t fZpos;

	WCSimRootGeom *fWCSimRootGeom;
	WCSimEmissionProfileManager *fEmissionProfileManager;

	ClassDef(WCSimGeometry, 0)
};
