#include "WCSimGeometry.hh"
#include "WCSimParameters.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimEmissionProfileManager.hh"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TError.h"

#include <iostream>
#include <cmath>
#include <cassert>

ClassImp(WCSimGeometry)

	static WCSimGeometry *fgGeometryHandle = 0;

Bool_t WCSimGeometry::TouchGeometry()
{
	if (WCSimGeometry::Instance()->GetNumPMTs() > 0)
	{
		return true;
	}
	else
	{
		std::cout << " ***  WCSimGeometry::TouchData() *** " << std::endl;
		std::cout << "  <error> need to build geometry... " << std::endl
				  << "    Call: WCSimGeometry::BuildGeometry(WCSimRootGeom*) " << std::endl;
		return false;
	}
}

WCSimGeometry *WCSimGeometry::Instance()
{
	if (!fgGeometryHandle)
	{
		fgGeometryHandle = new WCSimGeometry();
	}

	return fgGeometryHandle;
}

void WCSimGeometry::BuildGeometry(WCSimRootGeom *geom)
{
	return WCSimGeometry::Instance()->SetGeometry(geom);
}

void WCSimGeometry::WriteGeometry(const char *filename)
{
	WCSimGeometry::Instance()->WriteToFile(filename);
}

void WCSimGeometry::PrintGeometry()
{
	std::cout << " *** WCSimGeometry::PrintGeometry() *** " << std::endl;

	if (WCSimGeometry::TouchGeometry())
	{

		// geometry type
		if (fgGeometryHandle->GetGeoType() == WCSimGeometry::kCylinder)
		{
			std::cout << "  Detector Geometry: CYLINDER " << std::endl
					  << "   radius[m]="
					  << 1.0e-2 * fgGeometryHandle->GetCylRadius() << std::endl
					  << "   length[m]="
					  << 1.0e-2 * fgGeometryHandle->GetCylLength() << std::endl
					  << "   fiducial_radius[m]="
					  << 1.0e-2 * fgGeometryHandle->GetCylFiducialRadius() << std::endl
					  << "   fiducial_length[m]="
					  << 1.0e-2 * fgGeometryHandle->GetCylFiducialLength() << std::endl;
			std::cout << "   area[m^2]=" << 1.0e-4 * fgGeometryHandle->GetArea() << std::endl
					  << "   volume[m^3]="
					  << 1.0e-6 * fgGeometryHandle->GetVolume() << std::endl
					  << "   fiducial_volume[m^3]="
					  << 1.0e-6 * fgGeometryHandle->GetFiducialVolume() << std::endl;
		}
		if (fgGeometryHandle->GetGeoType() == WCSimGeometry::kMailBox)
		{
			std::cout << "  Detector Geometry: MAILBOX " << std::endl
					  << "   dx[m]="
					  << 1.0e-2 * fgGeometryHandle->GetMailBoxX() << std::endl
					  << "   dy[m]="
					  << 1.0e-2 * fgGeometryHandle->GetMailBoxY() << std::endl
					  << "   dz[m]="
					  << 1.0e-2 * fgGeometryHandle->GetMailBoxZ() << std::endl
					  << "   fiducial_dx[m]="
					  << 1.0e-2 * fgGeometryHandle->GetMailBoxFiducialX() << std::endl
					  << "   fiducial_dy[m]="
					  << 1.0e-2 * fgGeometryHandle->GetMailBoxFiducialY() << std::endl
					  << "   fiducial_dz[m]="
					  << 1.0e-2 * fgGeometryHandle->GetMailBoxFiducialZ() << std::endl;
			std::cout << "   area[m^2]=" << 1.0e-4 * fgGeometryHandle->GetArea() << std::endl
					  << "   volume[m^3]="
					  << 1.0e-6 * fgGeometryHandle->GetVolume() << std::endl
					  << "   fiducial_volume[m^3]="
					  << 1.0e-6 * fgGeometryHandle->GetFiducialVolume() << std::endl;
		}
		if (fgGeometryHandle->GetGeoType() == WCSimGeometry::kUnknown)
		{
			std::cout << "  Detector Geometry: UNKNOWN " << std::endl;
			std::cout << "   <warning> geometry unknown! " << std::endl;
		}

		// Dteector PMTs
		std::cout << "  Detector PMTs: " << std::endl;
		std::cout << "   number=" << fgGeometryHandle->GetNumPMTs() << std::endl;
		std::cout << "   radius[cm]=" << fgGeometryHandle->GetPMTRadius() << std::endl;
		std::cout << "   separation[cm]=" << fgGeometryHandle->GetPMTSeparation() << std::endl;
		std::cout << "   coverage=" << fgGeometryHandle->GetPMTCoverage() << std::endl;

		// attempt to figure out configuration...
		std::cout << "  Detector Configuration: " << std::endl;
		GeoConfiguration_t config = (GeoConfiguration_t)(fgGeometryHandle->GetGeoConfig());
		std::cout << "   " << fgGeometryHandle->AsString(config) << std::endl;
	}

	else
	{
		std::cout << " <warning> no geometry loaded " << std::endl;
	}

	return;
}

const char *WCSimGeometry::AsString(GeoConfiguration_t config)
{
	switch (config)
	{
	case WCSimGeometry::kDUSEL_Unknown:
		return "DUSEL_Unknown";
	case WCSimGeometry::kDUSEL_100kton_10inch_12percent:
		return "DUSEL_100kton_10inch_12percent";
	case WCSimGeometry::kDUSEL_100kton_10inch_15percent:
		return "DUSEL_100kton_10inch_15percent";
	case WCSimGeometry::kDUSEL_100kton_10inch_30percent:
		return "DUSEL_100kton_10inch_30percent";
	case WCSimGeometry::kDUSEL_100kton_10inch_40percent:
		return "DUSEL_100kton_10inch_40percent";
	case WCSimGeometry::kDUSEL_200kton_10inch_12percent:
		return "DUSEL_200kton_10inch_12percent";
	case WCSimGeometry::kDUSEL_200kton_12inch_10percent:
		return "DUSEL_200kton_12inch_10percent";
	case WCSimGeometry::kDUSEL_200kton_12inch_14percent:
		return "DUSEL_200kton_12inch_14percent";
	case WCSimGeometry::kDUSEL_150ktonMailbox_10inch_15percent:
		return "DUSEL_150ktonMailbox_10inch_15percent";
	case WCSimGeometry::kDUSEL_150ktonMailbox_10inch_30percent:
		return "DUSEL_150ktonMailbox_10inch_30percent";
	case WCSimGeometry::kDUSEL_150ktonMailbox_10inch_40percent:
		return "DUSEL_150ktonMailbox_10inch_40percent";
	default:
		return "DUSEL_Unknown";
	}
}

WCSimGeometry::WCSimGeometry() : fPMTs(0), fPmtX(0), fPmtY(0), fPmtZ(0), fPmtRegion(0), fWCSimRootGeom(0), fEmissionProfileManager(0)
{
	this->SetGeometry(0);
}

WCSimGeometry::WCSimGeometry(WCSimRootGeom *geom) : fPMTs(0), fPmtX(0), fPmtY(0), fPmtZ(0), fPmtRegion(0), fWCSimRootGeom(0)
{
	this->SetGeometry(geom);
}

WCSimGeometry::~WCSimGeometry()
{
	if (fPmtX)
		delete[] fPmtX;
	if (fPmtY)
		delete[] fPmtY;
	if (fPmtZ)
		delete[] fPmtZ;
	if (fPmtRegion)
		delete[] fPmtRegion;
	if (fEmissionProfileManager)
		delete fEmissionProfileManager;
}

Int_t WCSimGeometry::GetRegion(Int_t tube)
{
	if (tube >= 0 && tube < fPMTs)
	{
		return fPmtRegion[tube];
	}

	return -1;
}

Double_t WCSimGeometry::GetX(Int_t tube)
{
	if (tube >= 0 && tube < fPMTs)
	{
		return fPmtX[tube];
	}

	return -99999.9;
}

Double_t WCSimGeometry::GetY(Int_t tube)
{
	if (tube >= 0 && tube < fPMTs)
	{
		return fPmtY[tube];
	}

	return -99999.9;
}

Double_t WCSimGeometry::GetZ(Int_t tube)
{
	if (tube >= 0 && tube < fPMTs)
	{
		return fPmtZ[tube];
	}

	return -99999.9;
}

void WCSimGeometry::Reset()
{
	if (WCSimGeometry::Instance()->GetNumPMTs() > 0)
	{
		std::cout << " *** WCSimGeometry::Reset() *** " << std::endl;
	}

	WCSimGeometry::Instance()->SetGeometry(0);
}

void WCSimGeometry::SetGeometry(WCSimRootGeom *myGeometry)
{
	// Store WCSim Geometry
	// ====================
	fWCSimRootGeom = myGeometry;

	// Reset Geometry
	// ==============
	fPMTs = 0;

	if (fPmtX)
		delete[] fPmtX;
	fPmtX = 0;
	if (fPmtY)
		delete[] fPmtY;
	fPmtY = 0;
	if (fPmtZ)
		delete[] fPmtZ;
	fPmtZ = 0;
	if (fPmtRegion)
		delete[] fPmtRegion;
	fPmtRegion = 0;

	// Reset Geometry
	// ==============
	fGeoType = WCSimGeometry::kUnknown;
	fGeoConfig = WCSimGeometry::kDUSEL_Unknown;

	fCylinder = 0;
	fCylRadius = 0.0;
	fCylLength = 0.0;
	fCylFiducialRadius = 0.0;
	fCylFiducialLength = 0.0;

	fMailBox = 0;
	fMailBoxX = 0.0;
	fMailBoxY = 0.0;
	fMailBoxZ = 0.0;
	fMailBoxFiducialX = 0.0;
	fMailBoxFiducialY = 0.0;
	fMailBoxFiducialZ = 0.0;

	fDetectorArea = 0.0;
	fDetectorVolume = 0.0;
	fDetectorFiducialVolume = 0.0;

	// Sanity Check
	// ============
	if (myGeometry == 0)
	{
		return;
	}

	// Building New Geometry
	// =====================
	std::cout << " *** WCSimGeometry::BuildGeometry(...) *** " << std::endl;

	// Get Geometry Type
	// =================
	fGeoTypeInput = myGeometry->GetGeo_Type();

	std::cout << "   reading geometry: type=" << fGeoTypeInput << std::endl;

	// --- hack ---
	if (fGeoTypeInput < 0 || fGeoTypeInput > 1)
	{
		std::cout << "   <warning> Detected Unknown Geometry: " << fGeoTypeInput << std::endl;
		std::cout << "             Horrible hack! Assume it's a cylinder..." << std::endl;
		std::cout << "             Setting GeoType = 0 " << std::endl;
		fGeoTypeInput = 0;
	}

	if (fGeoTypeInput == 0)
	{ // geotype=0: cylinder
		std::cout << "   building cylindrical geometry " << std::endl;
		fGeoType = WCSimGeometry::kCylinder;
		fCylinder = 1;
		fCylRadius = myGeometry->GetWCCylRadius();
		fCylLength = myGeometry->GetWCCylLength();
		if (fCylRadius < 3000.0)
		{											 // hack for 100kton -> 200kton transition
			fCylFiducialRadius = fCylRadius - 250.0; // note: fiducial volume formerly defined
		}											 //       as 2.5m from cylinder edge, but
		else
		{											 //       is now reduced to 2.0m from edge.
			fCylFiducialRadius = fCylRadius - 200.0; //
		}											 //
		fCylFiducialLength = fCylLength - 400.0;	 //
		fDetectorArea = 2.0 * TMath::Pi() * fCylRadius * (fCylRadius + fCylLength);
		fDetectorVolume = TMath::Pi() * fCylRadius * fCylRadius * fCylLength;
		fDetectorFiducialVolume = TMath::Pi() * fCylFiducialRadius * fCylFiducialRadius * fCylFiducialLength;
	}

	if (fGeoTypeInput == 1)
	{ // geotype=1: mailbox
		std::cout << "   building mailbox geometry " << std::endl;
		fGeoType = WCSimGeometry::kMailBox;
		fMailBox = 1;
		fMailBoxX = myGeometry->GetMailBox_x();
		fMailBoxY = myGeometry->GetMailBox_y();
		fMailBoxZ = myGeometry->GetMailBox_z();
		fMailBoxFiducialX = fMailBoxX - 400.0; // note: fiducial volume defined
		fMailBoxFiducialY = fMailBoxY - 400.0; //       2.0m from mailbox edge
		fMailBoxFiducialZ = fMailBoxZ - 400.0;
		fDetectorArea = 2.0 * (fMailBoxX * fMailBoxY + fMailBoxY * fMailBoxZ + fMailBoxZ * fMailBoxX);
		fDetectorVolume = fMailBoxX * fMailBoxY * fMailBoxZ;
		fDetectorFiducialVolume = fMailBoxFiducialX * fMailBoxFiducialY * fMailBoxFiducialZ;
	}

	// Look up PMT Locations (NB: +1 to stop overflow)
	// ===============================================
	std::cout << "   reading PMTs: " << myGeometry->GetWCNumPMT() << std::endl;
	fPMTs = myGeometry->GetWCNumPMT();
	fVetoPMTs = myGeometry->GetWCNumVetoPMT();

	fPmtX = new Double_t[fPMTs];
	fPmtY = new Double_t[fPMTs];
	fPmtZ = new Double_t[fPMTs];
	fPmtRegion = new Int_t[fPMTs];

	for (Int_t n = 0; n < fPMTs; n++)
	{
		fPmtX[n] = -99999.9;
		fPmtY[n] = -99999.9;
		fPmtZ[n] = -99999.9;
		fPmtRegion[n] = -1;
	}

	for (Int_t n = 0; n < myGeometry->GetWCNumPMT(); n++)
	{
		// get next PMT
		WCSimRootPMT myPmt = myGeometry->GetPMT(n);

		Int_t tube = myPmt.GetTubeNo();
		Int_t region = -1; // myPmt.GetCylLoc();
		Double_t x = myPmt.GetPosition(0);
		Double_t y = myPmt.GetPosition(1);
		Double_t z = myPmt.GetPosition(2);

		// hack the region (for now)
		Double_t delta = 10.0; // cm

		if (fCylinder)
		{
			if (z >= -0.5 * fCylLength + delta && z <= +0.5 * fCylLength - delta)
				region = WCSimGeometry::kSide;
			if (z < -0.5 * fCylLength + delta)
				region = WCSimGeometry::kBottom;
			if (z > +0.5 * fCylLength - delta)
				region = WCSimGeometry::kTop;
			// Veto hits
			if (myPmt.GetCylLoc() == 3)
			{
				region = WCSimGeometry::kVeto;
			}
		}

		if (fMailBox)
		{
			//if( z>=-0.5*fMailBoxZ+delta
			// && z<=+0.5*fMailBoxZ-delta ) region = WCSimGeometry::kSide;
			if (z < -0.5 * fMailBoxZ + delta)
				region = WCSimGeometry::kBottom;
			if (z > +0.5 * fMailBoxZ - delta)
				region = WCSimGeometry::kTop;
			if (x < -0.5 * fMailBoxX + delta)
				region = WCSimGeometry::kBack;
			if (x > +0.5 * fMailBoxX - delta)
				region = WCSimGeometry::kFront;
			if (y < -0.5 * fMailBoxY + delta)
				region = WCSimGeometry::kRight;
			if (y > +0.5 * fMailBoxY - delta)
				region = WCSimGeometry::kLeft;
		}

		// add to array
		if (tube >= 0 && tube < fPMTs)
		{
			fPmtX[tube] = x;
			fPmtY[tube] = y;
			fPmtZ[tube] = z;
			fPmtRegion[tube] = region;
		}
	}

	// Calculate PMT Coverage
	// ======================

	// horrible hack (WCSim switched units on me!)
	if (myGeometry->GetWCPMTRadius() < 100.0)
	{
		fPMTRadius = 1.0 * myGeometry->GetWCPMTRadius(); // [must be cm already]
	}
	else
	{
		fPMTRadius = 0.1 * myGeometry->GetWCPMTRadius(); // [must be mm, so convert to cm]
	}

	fPMTSurfaceArea = 0.0;
	fPMTCoverage = 0.0;
	fPMTSeparation = 0.0;

	if (fPMTs > 0 && fDetectorArea > 0.0)
	{
		fPMTSurfaceArea = fPMTs * TMath::Pi() * fPMTRadius * fPMTRadius;
		fPMTCoverage = fPMTSurfaceArea / fDetectorArea;
		fPMTSeparation = sqrt(fDetectorArea / fPMTs);
	}

	// Determine Detector Configuration
	// ================================
	std::cout << "   determining detector configuration " << std::endl;

	Double_t kton = 0.5 + 1.0e-9 * fDetectorFiducialVolume;
	Int_t kton_rounded = (Int_t)(10.0 * (Int_t)(0.1 * kton));
	Double_t pmtdiam = 0.5 + (2.0 / 2.54) * fPMTRadius;
	Int_t pmtdiam_rounded = (Int_t)(2.0 * (Int_t)(0.5 * pmtdiam));
	Double_t coverage = 0.5 + 100.0 * fPMTCoverage;
	Int_t coverage_rounded = (Int_t)(coverage);

	std::cout << "    DUSEL_" << kton_rounded << "kton";
	if (fgGeometryHandle->GetGeoType() == WCSimGeometry::kMailBox)
		std::cout << "Mailbox";
	std::cout << "_";
	std::cout << pmtdiam_rounded << "inch_" << coverage_rounded << "perCent " << std::endl;

	// cylindrical geometry
	if (fGeoType == WCSimGeometry::kCylinder)
	{
		if (kton_rounded == 100 && pmtdiam_rounded == 10 && coverage_rounded == 12)
			fGeoConfig = WCSimGeometry::kDUSEL_100kton_10inch_12percent;
		if (kton_rounded == 100 && pmtdiam_rounded == 10 && coverage_rounded == 15)
			fGeoConfig = WCSimGeometry::kDUSEL_100kton_10inch_15percent;
		if (kton_rounded == 100 && pmtdiam_rounded == 10 && coverage_rounded == 30)
			fGeoConfig = WCSimGeometry::kDUSEL_100kton_10inch_30percent;
		if (kton_rounded == 100 && pmtdiam_rounded == 10 && coverage_rounded == 40)
			fGeoConfig = WCSimGeometry::kDUSEL_100kton_10inch_40percent;
		if (kton_rounded == 200 && pmtdiam_rounded == 10 && coverage_rounded == 12)
			fGeoConfig = WCSimGeometry::kDUSEL_200kton_10inch_12percent;
		if (kton_rounded == 200 && pmtdiam_rounded == 12 && coverage_rounded == 10)
			fGeoConfig = WCSimGeometry::kDUSEL_200kton_12inch_10percent;
		if (kton_rounded == 200 && pmtdiam_rounded == 12 && coverage_rounded == 14)
			fGeoConfig = WCSimGeometry::kDUSEL_200kton_12inch_14percent;
	}

	// mailbox geometry
	else if (fGeoType == WCSimGeometry::kMailBox)
	{
		if (kton_rounded == 150 && pmtdiam_rounded == 10 && coverage_rounded == 15)
			fGeoConfig = WCSimGeometry::kDUSEL_150ktonMailbox_10inch_15percent;
		if (kton_rounded == 150 && pmtdiam_rounded == 10 && coverage_rounded == 30)
			fGeoConfig = WCSimGeometry::kDUSEL_150ktonMailbox_10inch_30percent;
		if (kton_rounded == 150 && pmtdiam_rounded == 10 && coverage_rounded == 40)
			fGeoConfig = WCSimGeometry::kDUSEL_150ktonMailbox_10inch_40percent;
	}

	// --- sanity check ---
	else
	{
		std::cout << " <error> unknown geometry" << std::endl;
		assert(0);
	}

	// Print Geometry
	// ==============
	PrintGeometry();

	return;
}

void WCSimGeometry::WriteToFile(const char *filename)
{
	std::cout << " *** WCSimGeometry::WriteGeometry() *** " << std::endl;

	if (fPMTs <= 0)
	{
		std::cout << "  warning: no geometry loaded " << std::endl;
		return;
	}

	std::cout << "  writing geometry to: " << filename << std::endl;

	TDirectory *tmpd = 0;
	TFile *fFileGeometry = 0;
	TTree *fTreeGeometry = 0;

	if (!fFileGeometry)
	{
		tmpd = gDirectory;
		fFileGeometry = new TFile(filename, "recreate");
		fTreeGeometry = new TTree("ntuple", "detector geometry");
		fTreeGeometry->Branch("tube", &fTube, "tube/I");
		fTreeGeometry->Branch("region", &fRegion, "region/I");
		fTreeGeometry->Branch("x", &fXpos, "x/D");
		fTreeGeometry->Branch("y", &fYpos, "y/D");
		fTreeGeometry->Branch("z", &fZpos, "z/D");
		fTreeGeometry->SetAutoSave(100);
		gDirectory = tmpd;
	}

	for (Int_t n = 0; n < fPMTs; n++)
	{
		fTube = n;
		fRegion = fPmtRegion[n];
		fXpos = fPmtX[n];
		fYpos = fPmtY[n];
		fZpos = fPmtZ[n];

		if (fRegion >= 0)
		{
			tmpd = gDirectory;
			fFileGeometry->cd();
			fTreeGeometry->Fill();
			gDirectory = tmpd;
		}
	}

	tmpd = gDirectory;
	fFileGeometry->cd();
	fTreeGeometry->Write();
	fFileGeometry->Close();
	fFileGeometry = 0;
	fTreeGeometry = 0;
	gDirectory = tmpd;

	return;
}

Bool_t WCSimGeometry::InsideDetector(Double_t x, Double_t y, Double_t z)
{
	if (fGeoType == WCSimGeometry::kCylinder)
	{
		if (z >= -0.5 * fCylLength && z <= +0.5 * fCylLength && x * x + y * y <= fCylRadius * fCylRadius)
		{
			return 1;
		}
	}
	if (fGeoType == WCSimGeometry::kMailBox)
	{
		if (x >= -0.5 * fMailBoxX && x <= +0.5 * fMailBoxX && y >= -0.5 * fMailBoxY && y <= +0.5 * fMailBoxY && z >= -0.5 * fMailBoxZ && z <= +0.5 * fMailBoxZ)
		{
			return 1;
		}
	}

	return 0;
}

Bool_t WCSimGeometry::InsideFiducialVolume(Double_t x, Double_t y, Double_t z)
{
	if (fGeoType == WCSimGeometry::kCylinder)
	{
		if (z >= -0.5 * fCylFiducialLength && z <= +0.5 * fCylFiducialLength && x * x + y * y <= fCylFiducialRadius * fCylFiducialRadius)
		{
			return 1;
		}
	}
	if (fGeoType == WCSimGeometry::kMailBox)
	{
		if (x >= -0.5 * fMailBoxFiducialX && x <= +0.5 * fMailBoxFiducialX && y >= -0.5 * fMailBoxFiducialY && y <= +0.5 * fMailBoxFiducialY && z >= -0.5 * fMailBoxFiducialZ && z <= +0.5 * fMailBoxFiducialZ)
		{
			return 1;
		}
	}

	return 0;
}

Bool_t WCSimGeometry::InsideDetector(Double_t vx, Double_t vy, Double_t vz, Double_t ex, Double_t ey, Double_t ez)
{
	Double_t dx = ex - vx;
	Double_t dy = ey - vy;
	Double_t dz = ez - vz;
	Double_t ds = sqrt(dx * dx + dy * dy + dz * dz);

	Double_t px = dx / ds;
	Double_t py = dy / ds;
	Double_t pz = dz / ds;

	if (this->ForwardProjectionToEdge(vx, vy, vz, px, py, pz) >= 0.0 && this->BackwardProjectionToEdge(ex, ey, ez, px, py, pz) >= 0.0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

Double_t WCSimGeometry::DistanceToEdge(Double_t x, Double_t y, Double_t z)
{
	// Cylindrical Geometry
	// ====================
	if (fGeoType == WCSimGeometry::kCylinder)
	{

		// inside detector (convention: +ve dr)
		if (this->InsideDetector(x, y, z))
		{
			Double_t dr = 0.0;
			if (fCylRadius > dr)
				dr = fCylRadius;
			if (0.5 * fCylLength > dr)
				dr = 0.5 * fCylLength;
			if (-sqrt(x * x + y * y) + fCylRadius < dr)
				dr = -sqrt(x * x + y * y) + fCylRadius;
			if (-z + 0.5 * fCylLength < dr)
				dr = -z + 0.5 * fCylLength;
			if (+z + 0.5 * fCylLength < dr)
				dr = +z + 0.5 * fCylLength;
			return dr;
		}

		// outside detector (convention: -ve dr)
		else
		{

			// side region
			if (z >= -0.5 * fCylLength && z <= +0.5 * fCylLength)
			{
				return -sqrt(x * x + y * y) + fCylRadius;
			}

			// top region
			if (z <= -0.5 * fCylLength && x * x + y * y < fCylRadius * fCylRadius)
			{
				return +z + 0.5 * fCylLength;
			}
			if (z >= +0.5 * fCylLength && x * x + y * y < fCylRadius * fCylRadius)
			{
				return -z + 0.5 * fCylLength;
			}

			// corner regions
			if (z >= +0.5 * fCylLength && x * x + y * y >= fCylRadius)
			{
				Double_t dr = sqrt(x * x + y * y) - fCylRadius;
				Double_t dz = -z + 0.5 * fCylLength;
				return -sqrt(dr * dr + dz * dz);
			}
			if (z <= -0.5 * fCylLength && x * x + y * y >= fCylRadius)
			{
				Double_t dr = sqrt(x * x + y * y) - fCylRadius;
				Double_t dz = +z + 0.5 * fCylLength;
				return -sqrt(dr * dr + dz * dz);
			}
		}
	}

	// Mailbox Geometry
	// ================
	if (fGeoType == WCSimGeometry::kMailBox)
	{

		// inside detector (convention: +ve dr)
		if (this->InsideDetector(x, y, z))
		{
			Double_t dr = 0.0;
			if (0.5 * fMailBoxX > dr)
				dr = 0.5 * fMailBoxX;
			if (0.5 * fMailBoxY > dr)
				dr = 0.5 * fMailBoxY;
			if (0.5 * fMailBoxZ > dr)
				dr = 0.5 * fMailBoxZ;
			if (-x + 0.5 * fMailBoxX < dr)
				dr = -x + 0.5 * fMailBoxX;
			if (+x + 0.5 * fMailBoxX < dr)
				dr = +x + 0.5 * fMailBoxX;
			if (-y + 0.5 * fMailBoxY < dr)
				dr = -y + 0.5 * fMailBoxY;
			if (+y + 0.5 * fMailBoxY < dr)
				dr = +y + 0.5 * fMailBoxY;
			if (-z + 0.5 * fMailBoxZ < dr)
				dr = -z + 0.5 * fMailBoxZ;
			if (+z + 0.5 * fMailBoxZ < dr)
				dr = +z + 0.5 * fMailBoxZ;
			return dr;
		}

		// outside detector (convention: -ve dr)
		else
		{

			// six edge regions
			if (x <= -0.5 * fMailBoxX && y >= -0.5 * fMailBoxY && y <= +0.5 * fMailBoxY && z >= -0.5 * fMailBoxZ && z <= +0.5 * fMailBoxZ)
			{
				return +x + 0.5 * fMailBoxX;
			}
			if (x >= +0.5 * fMailBoxX && y >= -0.5 * fMailBoxY && y <= +0.5 * fMailBoxY && z >= -0.5 * fMailBoxZ && z <= +0.5 * fMailBoxZ)
			{
				return -x + 0.5 * fMailBoxX;
			}
			if (x >= -0.5 * fMailBoxX && x <= +0.5 * fMailBoxX && y <= -0.5 * fMailBoxY && z >= -0.5 * fMailBoxZ && z <= +0.5 * fMailBoxZ)
			{
				return +y + 0.5 * fMailBoxY;
			}
			if (x >= -0.5 * fMailBoxX && x <= +0.5 * fMailBoxX && y >= +0.5 * fMailBoxY && z >= -0.5 * fMailBoxZ && z <= +0.5 * fMailBoxZ)
			{
				return -y + 0.5 * fMailBoxY;
			}
			if (x >= -0.5 * fMailBoxX && x <= +0.5 * fMailBoxX && y >= -0.5 * fMailBoxY && y <= +0.5 * fMailBoxY && z <= -0.5 * fMailBoxZ)
			{
				return +z + 0.5 * fMailBoxZ;
			}
			if (x >= -0.5 * fMailBoxX && x <= +0.5 * fMailBoxX && y >= -0.5 * fMailBoxY && y <= +0.5 * fMailBoxY && z >= +0.5 * fMailBoxZ)
			{
				return -z + 0.5 * fMailBoxZ;
			}

			// eight corner regions
			if (x <= -0.5 * fMailBoxX && y <= -0.5 * fMailBoxY && z <= -0.5 * fMailBoxZ)
			{
				Double_t dx = +x + 0.5 * fMailBoxX;
				Double_t dy = +y + 0.5 * fMailBoxY;
				Double_t dz = +z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
			if (x >= +0.5 * fMailBoxX && y <= -0.5 * fMailBoxY && z <= -0.5 * fMailBoxZ)
			{
				Double_t dx = -x + 0.5 * fMailBoxX;
				Double_t dy = +y + 0.5 * fMailBoxY;
				Double_t dz = +z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
			if (x <= -0.5 * fMailBoxX && y >= +0.5 * fMailBoxY && z <= -0.5 * fMailBoxZ)
			{
				Double_t dx = +x + 0.5 * fMailBoxX;
				Double_t dy = -y + 0.5 * fMailBoxY;
				Double_t dz = +z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
			if (x <= -0.5 * fMailBoxX && y <= -0.5 * fMailBoxY && z >= +0.5 * fMailBoxZ)
			{
				Double_t dx = +x + 0.5 * fMailBoxX;
				Double_t dy = +y + 0.5 * fMailBoxY;
				Double_t dz = -z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
			if (x >= +0.5 * fMailBoxX && y >= +0.5 * fMailBoxY && z <= -0.5 * fMailBoxZ)
			{
				Double_t dx = -x + 0.5 * fMailBoxX;
				Double_t dy = -y + 0.5 * fMailBoxY;
				Double_t dz = +z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
			if (x >= +0.5 * fMailBoxX && y <= -0.5 * fMailBoxY && z >= +0.5 * fMailBoxZ)
			{
				Double_t dx = -x + 0.5 * fMailBoxX;
				Double_t dy = +y + 0.5 * fMailBoxY;
				Double_t dz = -z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
			if (x <= -0.5 * fMailBoxX && y >= +0.5 * fMailBoxY && z >= +0.5 * fMailBoxZ)
			{
				Double_t dx = +x + 0.5 * fMailBoxX;
				Double_t dy = -y + 0.5 * fMailBoxY;
				Double_t dz = -z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
			if (x >= +0.5 * fMailBoxX && y >= +0.5 * fMailBoxY && z >= +0.5 * fMailBoxZ)
			{
				Double_t dx = -x + 0.5 * fMailBoxX;
				Double_t dy = -y + 0.5 * fMailBoxY;
				Double_t dz = -z + 0.5 * fMailBoxZ;
				return -sqrt(dx * dx + dy * dy + dz * dz);
			}
		}
	}

	return -99999.9;
}

Double_t WCSimGeometry::ForwardProjectionToEdge(Double_t x, Double_t y, Double_t z, Double_t px, Double_t py,
												Double_t pz)
{
	Double_t xproj = 0.0;
	Double_t yproj = 0.0;
	Double_t zproj = 0.0;
	Int_t regionproj = 0;

	this->ProjectToNearEdge(x, y, z, px, py, pz, xproj, yproj, zproj, regionproj);

	if (regionproj > WCSimGeometry::kUnknown)
	{
		return sqrt((xproj - x) * (xproj - x) + (yproj - y) * (yproj - y) + (zproj - z) * (zproj - z));
	}

	this->ProjectToNearEdge(x, y, z, -px, -py, -pz, xproj, yproj, zproj, regionproj);

	if (regionproj > WCSimGeometry::kUnknown)
	{
		return -sqrt((xproj - x) * (xproj - x) + (yproj - y) * (yproj - y) + (zproj - z) * (zproj - z));
	}

	return -99999.9;
}

Double_t WCSimGeometry::BackwardProjectionToEdge(Double_t x, Double_t y, Double_t z, Double_t px, Double_t py,
												 Double_t pz)
{
	Double_t xproj = 0.0;
	Double_t yproj = 0.0;
	Double_t zproj = 0.0;
	Int_t regionproj = 0;

	this->ProjectToNearEdge(x, y, z, -px, -py, -pz, xproj, yproj, zproj, regionproj);

	if (regionproj > WCSimGeometry::kUnknown)
	{
		return sqrt((xproj - x) * (xproj - x) + (yproj - y) * (yproj - y) + (zproj - z) * (zproj - z));
	}

	this->ProjectToNearEdge(x, y, z, px, py, pz, xproj, yproj, zproj, regionproj);

	if (regionproj > WCSimGeometry::kUnknown)
	{
		return -sqrt((xproj - x) * (xproj - x) + (yproj - y) * (yproj - y) + (zproj - z) * (zproj - z));
	}

	return -99999.9;
}

// Function to find the minimum and maximum distances along the track direction from the vertex to remain
// inside the detector. Used mostly by the fit routine that moves the vertex along the track direction.
void WCSimGeometry::ProjectMinMaxPosition(TVector3 vtx, TVector3 dir, Double_t &sMin, Double_t &sMax)
{

	Double_t xproj1 = 0.0;
	Double_t yproj1 = 0.0;
	Double_t zproj1 = 0.0;
	Int_t regionproj1 = 0;
	this->ProjectToNearEdge(vtx.X(), vtx.Y(), vtx.Z(), dir.X(), dir.Y(), dir.Z(), xproj1, yproj1, zproj1, regionproj1);
	Double_t xproj2 = 0.0;
	Double_t yproj2 = 0.0;
	Double_t zproj2 = 0.0;
	Int_t regionproj2 = 0;
	this->ProjectToNearEdge(vtx.X(), vtx.Y(), vtx.Z(), -dir.X(), -dir.Y(), -dir.Z(), xproj2, yproj2, zproj2,
							regionproj2);

	TVector3 proj1 = TVector3(xproj1, yproj1, zproj1);
	TVector3 proj2 = TVector3(xproj2, yproj2, zproj2);

	std::cout << proj1.X() << "," << proj1.Y() << "," << proj1.Z() << std::endl;
	std::cout << "Projection presumably in direction of track:  " << (proj1 - vtx).Angle(dir) << std::endl;
	std::cout << proj2.X() << "," << proj2.Y() << "," << proj2.Z() << std::endl;
	std::cout << "Projection presumably in the other direction: " << (proj2 - vtx).Angle(dir) << std::endl;

	sMax = (TVector3(xproj1, yproj1, zproj1) - vtx).Mag();
	sMin = -1 * (vtx - TVector3(xproj2, yproj2, zproj2)).Mag();
}

void WCSimGeometry::ProjectToNearEdge(Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py, Double_t pz,
									  Double_t &xproj, Double_t &yproj, Double_t &zproj, Int_t &regionproj)
{
	this->ProjectToEdge(0, x0, y0, z0, px, py, pz, xproj, yproj, zproj, regionproj);

	//
	// this->ProjectToEdgeOld(0,
	//                        x0,y0,z0,
	//                        px,py,pz,
	//                        xproj,yproj,zproj,
	//                        regionproj);
	//

	return;
}

void WCSimGeometry::ProjectToFarEdge(Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py, Double_t pz,
									 Double_t &xproj, Double_t &yproj, Double_t &zproj, Int_t &regionproj)
{
	this->ProjectToEdge(1, x0, y0, z0, px, py, pz, xproj, yproj, zproj, regionproj);

	//
	// this->ProjectToEdgeOld(1,
	//                        x0,y0,z0,
	//                        px,py,pz,
	//                        xproj,yproj,zproj,
	//                        regionproj);
	//

	return;
}

void WCSimGeometry::ProjectToEdgeOld(Bool_t useFarEdge, Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py,
									 Double_t pz, Double_t &xproj, Double_t &yproj, Double_t &zproj, Int_t &regionproj)
{
	// default locations
	// =================
	xproj = -99999.9;
	yproj = -99999.9;
	zproj = -99999.9;
	regionproj = WCSimGeometry::kUnknown;

	// CYLINDRICAL GEOMETRY
	// ====================
	if (fGeoType == WCSimGeometry::kCylinder)
	{

		//
		// std::cout << " --- ProjectToEdge --- " << std::endl
		//           << "  Vtx=(" << x0 << "," << y0 << "," << z0 << ") " << std::endl
		//           << "  Dir=(" << px << "," << py << "," << pz << ") --- " << std::endl;
		//

		// projections
		// ===========
		Double_t r = fCylRadius;
		Double_t L = fCylLength;

		Double_t x = 0.0;
		Double_t y = 0.0;
		Double_t z = 0.0;
		Int_t region = WCSimGeometry::kUnknown;

		Double_t mxy = 0.0;
		Double_t mz = 0.0;
		Double_t xlow = 0.0;
		Double_t ylow = 0.0;
		Double_t zlow = 0.0;
		Double_t slow = 0.0;
		Double_t xhigh = 0.0;
		Double_t yhigh = 0.0;
		Double_t zhigh = 0.0;
		Double_t shigh = 0.0;

		Bool_t foundProjectionXY = 0;
		Bool_t foundProjectionX = 0;
		Bool_t foundProjectionY = 0;
		Bool_t foundProjectionZ = 0;
		Bool_t foundProjection = 0;

		// find radial intersections
		// =========================
		if (px != 0.0)
		{
			mxy = py / px;

			// sanity check: is there an intersection?
			if ((1.0 + mxy * mxy) * r * r - (mxy * x0 - y0) * (mxy * x0 - y0) > 0.0)
			{

				// intersection at low x
				xlow = -(x0 + mxy * y0) - sqrt((1.0 + mxy * mxy) * r * r - (mxy * x0 - y0) * (mxy * x0 - y0));
				xlow /= (1.0 + mxy * mxy);
				ylow = mxy * xlow;

				// intersection at high x
				xhigh = -(x0 + mxy * y0) + sqrt((1.0 + mxy * mxy) * r * r - (mxy * x0 - y0) * (mxy * x0 - y0));
				xhigh /= (1.0 + mxy * mxy);
				yhigh = mxy * xhigh;

				// forward or backward
				if (px >= 0)
				{
					if (xlow >= 0.0)
						slow = +1.0;
					else
						slow = -1.0;
					if (xhigh >= 0.0)
						shigh = +1.0;
					else
						shigh = -1.0;
				}
				else
				{
					if (xlow <= 0.0)
						slow = +1.0;
					else
						slow = -1.0;
					if (xhigh <= 0.0)
						shigh = +1.0;
					else
						shigh = -1.0;
				}

				// intersection z values
				mz = pz / sqrt(px * px + py * py);
				zlow = z0 + mz * slow * sqrt(xlow * xlow + ylow * ylow);
				zhigh = z0 + mz * shigh * sqrt(xhigh * xhigh + yhigh * yhigh);

				// found radial projection
				foundProjectionXY = 1;
			}
		}

		// travelling along y-z plane
		else if (py != 0.0)
		{

			// sanity check: is there an intersection?
			if (r * r - x0 * x0 > 0.0)
			{

				// intersection at low y
				xlow = 0.0;
				ylow = -sqrt(r * r - x0 * x0) - y0;

				// intersection at high y
				xhigh = 0.0;
				yhigh = +sqrt(r * r - x0 * x0) - y0;

				// forward or backward
				if (py >= 0)
				{
					if (ylow >= 0.0)
						slow = +1.0;
					else
						slow = -1.0;
					if (yhigh >= 0.0)
						shigh = +1.0;
					else
						shigh = -1.0;
				}
				else
				{
					if (ylow <= 0.0)
						slow = +1.0;
					else
						slow = -1.0;
					if (yhigh <= 0.0)
						shigh = +1.0;
					else
						shigh = -1.0;
				}

				// intersection z values
				mz = pz / sqrt(px * px + py * py);
				zlow = z0 + mz * slow * sqrt(xlow * xlow + ylow * ylow);
				zhigh = z0 + mz * shigh * sqrt(xhigh * xhigh + yhigh * yhigh);

				// found radial projection
				foundProjectionXY = 1;
			}
		}

		// travelling along z-axis
		else if (pz != 0.0)
		{

			// sanity check: is there an intersection?
			if (x0 * x0 + y0 * y0 < r * r)
			{
				xlow = -x0;
				ylow = -y0;
				zlow = -L / 2.0;

				xhigh = -x0;
				yhigh = -y0;
				zhigh = +L / 2.0;

				// found radial projection
				foundProjectionXY = 1;
			}
		}

		// check projection with cylinder
		// ==============================
		if (foundProjectionXY)
		{

			// project backwards along x-axis
			if (px < 0.0)
			{
				if (xlow < 0 || xhigh <= 0)
				{
					foundProjectionX = 1;
				}
			}

			// project forwards along x-axis
			if (px > 0.0)
			{
				if (xlow >= 0 || xhigh > 0)
				{
					foundProjectionX = 1;
				}
			}

			// project backwards along y-axis
			if (py < 0.0)
			{
				if (ylow < 0 || yhigh <= 0)
				{
					foundProjectionY = 1;
				}
			}

			// project forwards along y-axis
			if (py > 0.0)
			{
				if (ylow >= 0 || yhigh > 0)
				{
					foundProjectionY = 1;
				}
			}

			// project backwards along z-axis
			if (pz < 0.0)
			{
				if (z0 > -L / 2.0 && (zlow <= +L / 2.0 || zhigh <= +L / 2.0))
				{
					foundProjectionZ = 1;
				}
			}

			// project forwards along z-axis
			if (pz > 0)
			{
				if (z0 < +L / 2.0 && (zlow >= -L / 2.0 || zhigh >= -L / 2.0))
				{
					foundProjectionZ = 1;
				}
			}

			// overall projection
			if (foundProjectionX && foundProjectionY && foundProjectionZ)
			{
				foundProjection = 1;
			}
		}

		// apply top/bottom correction
		// ===========================
		if (foundProjection)
		{

			if (zlow > +L / 2.0)
			{
				xlow *= (L / 2.0 - z0) / (zlow - z0);
				ylow *= (L / 2.0 - z0) / (zlow - z0);
				zlow = +L / 2.0;
			}

			if (zlow < -L / 2.0)
			{
				xlow *= (-L / 2.0 - z0) / (zlow - z0);
				ylow *= (-L / 2.0 - z0) / (zlow - z0);
				zlow = -L / 2.0;
			}

			if (zhigh > +L / 2.0)
			{
				xhigh *= (L / 2.0 - z0) / (zhigh - z0);
				yhigh *= (L / 2.0 - z0) / (zhigh - z0);
				zhigh = +L / 2.0;
			}

			if (zhigh < -L / 2.0)
			{
				xhigh *= (-L / 2.0 - z0) / (zhigh - z0);
				yhigh *= (-L / 2.0 - z0) / (zhigh - z0);
				zhigh = -L / 2.0;
			}
		}

		//
		// std::cout << "  FoundProjection = " << foundProjection << std::endl;
		// std::cout << "   low(x,y,z) = (" << xlow+x0 << "," << ylow+y0 << "," << zlow << ") " << std::endl;
		// std::cout << "   high(x,y,z) = (" << xhigh+x0 << "," << yhigh+y0 << "," << zhigh << ") " << std::endl;
		//

		// choose intersection point
		// =========================
		if (foundProjection)
		{

			if (px != 0.0)
			{

				// project backwards along x-axis
				if (px < 0)
				{
					if (useFarEdge)
					{
						if (xlow <= 0.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
					}
					else
					{
						if (xhigh < 0.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
						else if (xlow <= 0.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
					}
				}

				// project forwards along x-axis
				if (px > 0.0)
				{
					if (useFarEdge)
					{
						if (xhigh >= 0.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
					}
					else
					{
						if (xlow > 0.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
						else if (xhigh >= 0.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
					}
				}
			}

			else if (py != 0.0)
			{

				// project backwards along y-axis
				if (py < 0.0)
				{
					if (useFarEdge)
					{
						if (ylow <= 0.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
					}
					else
					{
						if (yhigh < 0.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
						else if (ylow <= 0.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
					}
				}

				// project forwards along y-axis
				if (py > 0.0)
				{
					if (useFarEdge)
					{
						if (yhigh >= 0.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
					}
					else
					{
						if (ylow > 0.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
						else if (yhigh >= 0.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
					}
				}
			}

			else
			{

				// project backwards along z-axis
				if (pz < 0.0)
				{
					if (useFarEdge)
					{
						if (z0 >= -L / 2.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
					}
					else
					{
						if (z0 > +L / 2)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
						else if (z0 >= -L / 2.0)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
					}
				}

				// project forwards along z-axis
				if (pz > 0.0)
				{
					if (useFarEdge)
					{
						if (z0 <= +L / 2.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
					}
					else
					{
						if (z0 < -L / 2)
						{
							x = xlow + x0;
							y = ylow + y0;
							z = zlow;
						}
						else if (z0 <= +L / 2.0)
						{
							x = xhigh + x0;
							y = yhigh + y0;
							z = zhigh;
						}
					}
				}
			}

			// assign region
			if (-L / 2.0 < z && z < +L / 2.0)
			{
				region = WCSimGeometry::kSide; // region 1
			}

			if (z >= +L / 2.0)
			{
				region = WCSimGeometry::kTop; // region 0
			}

			if (z <= -L / 2.0)
			{
				region = WCSimGeometry::kBottom; // region 2
			}
		}

		//
		// std::cout << "  Return Projection: x=" << x << " y=" << y << " z=" << z << " region=" << region << std::endl;
		//

		// return projection
		// =================
		if (region > WCSimGeometry::kUnknown)
		{

			xproj = x;
			yproj = y;
			zproj = z;
			regionproj = region;

			return;
		}
	}

	// MAILBOX GEOMETRY
	// ================
	if (fGeoType == WCSimGeometry::kMailBox)
	{

		Double_t dX = fMailBoxX;
		Double_t dY = fMailBoxY;
		Double_t dZ = fMailBoxZ;

		Double_t x_top = -99999.9;
		Double_t y_top = -99999.9;
		Double_t z_top = -99999.9;

		Double_t x_bottom = -99999.9;
		Double_t y_bottom = -99999.9;
		Double_t z_bottom = -99999.9;

		Double_t x_front = -99999.9;
		Double_t y_front = -99999.9;
		Double_t z_front = -99999.9;

		Double_t x_back = -99999.9;
		Double_t y_back = -99999.9;
		Double_t z_back = -99999.9;

		Double_t x_right = -99999.9;
		Double_t y_right = -99999.9;
		Double_t z_right = -99999.9;

		Double_t x_left = -99999.9;
		Double_t y_left = -99999.9;
		Double_t z_left = -99999.9;

		Double_t x = 0.0;
		Double_t y = 0.0;
		Double_t z = 0.0;
		Int_t region = WCSimGeometry::kUnknown;

		// calculate projections in each plane
		// ===================================
		if (pz != 0.0)
		{
			x_top = x0 + (px / pz) * (+0.5 * dZ - z0);
			y_top = y0 + (py / pz) * (+0.5 * dZ - z0);
			z_top = +0.5 * dZ;

			x_bottom = x0 + (px / pz) * (-0.5 * dZ - z0);
			y_bottom = y0 + (py / pz) * (-0.5 * dZ - z0);
			z_bottom = -0.5 * dZ;
		}

		if (py != 0.0)
		{
			z_left = z0 + (pz / py) * (+0.5 * dY - y0);
			x_left = x0 + (px / py) * (+0.5 * dY - y0);
			y_left = +0.5 * dY;

			z_right = z0 + (pz / py) * (-0.5 * dY - y0);
			x_right = x0 + (px / py) * (-0.5 * dY - y0);
			y_right = -0.5 * dY;
		}

		if (px != 0.0)
		{
			y_front = y0 + (py / px) * (+0.5 * dX - x0);
			z_front = z0 + (pz / px) * (+0.5 * dX - x0);
			x_front = +0.5 * dX;

			y_back = y0 + (py / px) * (-0.5 * dX - x0);
			z_back = z0 + (pz / px) * (-0.5 * dX - x0);
			x_back = -0.5 * dX;
		}

		// calculate intersections with box
		// ================================
		if (x_top >= -0.5 * dX && x_top <= +0.5 * dX && y_top >= -0.5 * dY && y_top <= +0.5 * dY)
		{
			if (pz < 0.0 && z0 > z_top)
			{
				if (!useFarEdge)
				{
					x = x_top;
					y = y_top;
					z = z_top;
					region = WCSimGeometry::kTop;
				}
			}
			if (pz > 0.0 && z0 <= z_top)
			{
				if (z0 > z_bottom || useFarEdge)
				{
					x = x_top;
					y = y_top;
					z = z_top;
					region = WCSimGeometry::kTop;
				}
			}
		}

		if (x_bottom >= -0.5 * dX && x_bottom <= +0.5 * dX && y_bottom >= -0.5 * dY && y_bottom <= +0.5 * dY)
		{
			if (pz > 0.0 && z0 < z_bottom)
			{
				if (!useFarEdge)
				{
					x = x_bottom;
					y = y_bottom;
					z = z_bottom;
					region = WCSimGeometry::kBottom;
				}
			}
			if (pz < 0.0 && z0 >= z_bottom)
			{
				if (z0 < z_top || useFarEdge)
				{
					x = x_bottom;
					y = y_bottom;
					z = z_bottom;
					region = WCSimGeometry::kBottom;
				}
			}
		}

		if (y_front >= -0.5 * dY && y_front <= +0.5 * dY && z_front >= -0.5 * dZ && z_front <= +0.5 * dZ)
		{
			if (px < 0.0 && x0 > x_front)
			{
				if (!useFarEdge)
				{
					x = x_front;
					y = y_front;
					z = z_front;
					region = WCSimGeometry::kFront;
				}
			}
			if (px > 0.0 && x0 <= x_front)
			{
				if (x0 > x_back || useFarEdge)
				{
					x = x_front;
					y = y_front;
					z = z_front;
					region = WCSimGeometry::kFront;
				}
			}
		}

		if (y_back >= -0.5 * dY && y_back <= +0.5 * dY && z_back >= -0.5 * dZ && z_back <= +0.5 * dZ)
		{
			if (px > 0.0 && x0 < x_back)
			{
				if (!useFarEdge)
				{
					x = x_back;
					y = y_back;
					z = z_back;
					region = WCSimGeometry::kBack;
				}
			}
			if (px < 0.0 && x0 >= x_back)
			{
				if (x0 < x_front || useFarEdge)
				{
					x = x_back;
					y = y_back;
					z = z_back;
					region = WCSimGeometry::kBack;
				}
			}
		}

		if (z_left >= -0.5 * dZ && z_left <= +0.5 * dZ && x_left >= -0.5 * dX && x_left <= +0.5 * dX)
		{
			if (py < 0.0 && y0 > y_left)
			{
				if (!useFarEdge)
				{
					x = x_left;
					y = y_left;
					z = z_left;
					region = WCSimGeometry::kLeft;
				}
			}
			if (py > 0.0 && y0 <= y_left)
			{
				if (y0 > y_right || useFarEdge)
				{
					x = x_left;
					y = y_left;
					z = z_left;
					region = WCSimGeometry::kLeft;
				}
			}
		}

		if (z_right >= -0.5 * dZ && z_right <= +0.5 * dZ && x_right >= -0.5 * dX && x_right <= +0.5 * dX)
		{
			if (py > 0.0 && y0 < y_right)
			{
				if (!useFarEdge)
				{
					x = x_right;
					y = y_right;
					z = z_right;
					region = WCSimGeometry::kRight;
				}
			}
			if (py < 0.0 && y0 >= y_right)
			{
				if (y0 < y_left || useFarEdge)
				{
					x = x_right;
					y = y_right;
					z = z_right;
					region = WCSimGeometry::kRight;
				}
			}
		}

		// return projection
		// =================
		if (region > WCSimGeometry::kUnknown)
		{

			xproj = x;
			yproj = y;
			zproj = z;
			regionproj = region;

			return;
		}
	}

	return;
}

void WCSimGeometry::ProjectToEdge(Bool_t useFarEdge, Double_t x0, Double_t y0, Double_t z0, Double_t px, Double_t py,
								  Double_t pz, Double_t &xproj, Double_t &yproj, Double_t &zproj, Int_t &regionproj)
{
	// default locations
	// =================
	xproj = -99999.9;
	yproj = -99999.9;
	zproj = -99999.9;
	regionproj = WCSimGeometry::kUnknown;

	Double_t xNear = -99999.9;
	Double_t yNear = -99999.9;
	Double_t zNear = -99999.9;
	Int_t regionNear = WCSimGeometry::kUnknown;

	Double_t xFar = -99999.9;
	Double_t yFar = -99999.9;
	Double_t zFar = -99999.9;
	Int_t regionFar = WCSimGeometry::kUnknown;

	// CYLINDRICAL GEOMETRY
	// ====================
	if (fGeoType == WCSimGeometry::kCylinder)
	{

		Double_t r = fCylRadius;
		Double_t L = fCylLength;

		Bool_t foundProjectionXY = 0;
		Bool_t foundProjectionZ = 0;

		Double_t t1 = 0.0;
		Double_t x1 = 0.0;
		Double_t y1 = 0.0;
		Double_t z1 = 0.0;
		Int_t region1 = -1;

		Double_t t2 = 0.0;
		Double_t x2 = 0.0;
		Double_t y2 = 0.0;
		Double_t z2 = 0.0;
		Int_t region2 = -1;

		Double_t rSq = r * r;
		Double_t r0r0 = x0 * x0 + y0 * y0;
		Double_t r0p = x0 * px + y0 * py;
		Double_t pSq = px * px + py * py;

		// calculate intersection in XY
		if (pSq > 0.0)
		{
			if (r0p * r0p - pSq * (r0r0 - rSq) > 0.0)
			{
				t1 = (-r0p - sqrt(r0p * r0p - pSq * (r0r0 - rSq))) / pSq;
				t2 = (-r0p + sqrt(r0p * r0p - pSq * (r0r0 - rSq))) / pSq;
				foundProjectionXY = 1;
			}
		}

		// propagation along z-axis
		else if (r0r0 <= rSq)
		{

			if (pz > 0)
			{
				t1 = -L / 2.0 - z0;
				t2 = +L / 2.0 - z0;
			}
			else
			{
				t1 = -L / 2.0 + z0;
				t2 = +L / 2.0 + z0;
			}
			foundProjectionXY = 1;
		}

		// found intersection in XY
		if (foundProjectionXY)
		{

			z1 = z0 + t1 * pz;
			z2 = z0 + t2 * pz;
			if ((z1 >= -L / 2.0 && z2 <= +L / 2.0) || (z2 >= -L / 2.0 && z1 <= +L / 2.0))
			{
				foundProjectionZ = 1;
			}
		}

		// found intersection in Z
		if (foundProjectionZ)
		{

			// first intersection
			if (z1 > -L / 2.0 && z1 < +L / 2.0)
			{
				region1 = WCSimGeometry::kSide;
			}
			if (z1 >= +L / 2.0)
			{
				region1 = WCSimGeometry::kTop;
				if (z1 > +L / 2.0)
				{
					z1 = +L / 2.0;
					t1 = (+L / 2.0 - z0) / pz;
				}
			}
			if (z1 <= -L / 2.0)
			{
				region1 = WCSimGeometry::kBottom;
				if (z1 < -L / 2.0)
				{
					z1 = -L / 2.0;
					t1 = (-L / 2.0 - z0) / pz;
				}
			}

			x1 = x0 + t1 * px;
			y1 = y0 + t1 * py;

			// second intersection
			if (z2 > -L / 2.0 && z2 < +L / 2.0)
			{
				region2 = WCSimGeometry::kSide;
			}
			if (z2 >= +L / 2.0)
			{
				region2 = WCSimGeometry::kTop;
				if (z2 > +L / 2.0)
				{
					z2 = +L / 2.0;
					t2 = (+L / 2.0 - z0) / pz;
				}
			}
			if (z2 <= -L / 2.0)
			{
				region2 = WCSimGeometry::kBottom;
				if (z2 < -L / 2.0)
				{
					z2 = -L / 2.0;
					t2 = (-L / 2.0 - z0) / pz;
				}
			}

			x2 = x0 + t2 * px;
			y2 = y0 + t2 * py;

			// near/far projection
			if (t1 >= 0)
			{
				//      if( fabs(t1)<=fabs(t2)){
				xNear = x1;
				yNear = y1;
				zNear = z1;
				regionNear = region1;

				xFar = x2;
				yFar = y2;
				zFar = z2;
				regionFar = region2;
			}
			else if (t2 > 0)
			{
				//      else{
				xNear = x2;
				yNear = y2;
				zNear = z2;
				regionNear = region2;

				xFar = x1;
				yFar = y1;
				zFar = z1;
				regionFar = region1;
			}
		}
	}

	// MAILBOX GEOMETRY
	// ================
	if (fGeoType == WCSimGeometry::kMailBox)
	{

		Double_t dX = fMailBoxX;
		Double_t dY = fMailBoxY;
		Double_t dZ = fMailBoxZ;

		Bool_t foundProjection = 0;

		Double_t t = 0.0;
		Double_t x = 0.0;
		Double_t y = 0.0;
		Double_t z = 0.0;

		Double_t t1 = 0.0;
		Double_t x1 = 0.0;
		Double_t y1 = 0.0;
		Double_t z1 = 0.0;
		Int_t region1 = -1;

		Double_t t2 = 0.0;
		Double_t x2 = 0.0;
		Double_t y2 = 0.0;
		Double_t z2 = 0.0;
		Int_t region2 = -1;

		if (pz != 0)
		{

			// top face (+dZ/2.0)
			t = (+dZ / 2.0 - z0) / pz;
			x = x0 + t * px;
			y = y0 + t * py;
			if (x >= -dX / 2.0 && x <= +dX / 2.0 && y >= -dY / 2.0 && y <= +dY / 2.0)
			{
				if (pz > 0.0)
				{
					t2 = t;
					region2 = WCSimGeometry::kTop;
					foundProjection = 1;
				}
				else
				{
					t1 = t;
					region1 = WCSimGeometry::kTop;
					foundProjection = 1;
				}
			}

			// bottom face (-dZ/2.0)
			t = (-dZ / 2.0 - z0) / pz;
			x = x0 + t * px;
			y = y0 + t * py;
			if (x >= -dX / 2.0 && x <= +dX / 2.0 && y >= -dY / 2.0 && y <= +dY / 2.0)
			{
				if (pz < 0.0)
				{
					t2 = t;
					region2 = WCSimGeometry::kBottom;
					foundProjection = 1;
				}
				else
				{
					t1 = t;
					region1 = WCSimGeometry::kBottom;
					foundProjection = 1;
				}
			}
		}

		if (px != 0)
		{

			// front face (+dX/2.0)
			t = (+dX / 2.0 - x0) / px;
			y = y0 + t * py;
			z = z0 + t * pz;
			if (y >= -dY / 2.0 && y <= +dY / 2.0 && z >= -dZ / 2.0 && z <= +dZ / 2.0)
			{
				if (px > 0.0)
				{
					t2 = t;
					region2 = WCSimGeometry::kFront;
					foundProjection = 1;
				}
				else
				{
					t1 = t;
					region1 = WCSimGeometry::kFront;
					foundProjection = 1;
				}
			}

			// back face (-dX/2.0)
			t = (-dX / 2.0 - x0) / px;
			y = y0 + t * py;
			z = z0 + t * pz;
			if (y >= -dY / 2.0 && y <= +dY / 2.0 && z >= -dZ / 2.0 && z <= +dZ / 2.0)
			{
				if (px < 0.0)
				{
					t2 = t;
					region2 = WCSimGeometry::kBack;
					foundProjection = 1;
				}
				else
				{
					t1 = t;
					region1 = WCSimGeometry::kBack;
					foundProjection = 1;
				}
			}
		}

		if (py != 0)
		{

			// left face (+dY/2.0)
			t = (+dY / 2.0 - y0) / py;
			z = z0 + t * pz;
			x = x0 + t * px;
			if (z >= -dZ / 2.0 && z <= +dZ / 2.0 && x >= -dX / 2.0 && x <= +dX / 2.0)
			{
				if (py > 0.0)
				{
					t2 = t;
					region2 = WCSimGeometry::kLeft;
					foundProjection = 1;
				}
				else
				{
					t1 = t;
					region1 = WCSimGeometry::kLeft;
					foundProjection = 1;
				}
			}

			// right face (-dY/2.0)
			t = (-dY / 2.0 - y0) / py;
			z = z0 + t * pz;
			x = x0 + t * px;
			if (z >= -dZ / 2.0 && z <= +dZ / 2.0 && x >= -dX / 2.0 && x <= +dX / 2.0)
			{
				if (py < 0.0)
				{
					t2 = t;
					region2 = WCSimGeometry::kRight;
					foundProjection = 1;
				}
				else
				{
					t1 = t;
					region1 = WCSimGeometry::kRight;
					foundProjection = 1;
				}
			}
		}

		// found projection
		if (foundProjection == 1)
		{

			x1 = x0 + t1 * px;
			y1 = y0 + t1 * py;
			z1 = z0 + t1 * pz;

			x2 = x0 + t2 * px;
			y2 = y0 + t2 * py;
			z2 = z0 + t2 * pz;

			// near/far projection
			if (t1 >= 0)
			{
				xNear = x1;
				yNear = y1;
				zNear = z1;
				regionNear = region1;

				xFar = x2;
				yFar = y2;
				zFar = z2;
				regionFar = region2;
			}
			else if (t2 > 0)
			{
				xNear = x2;
				yNear = y2;
				zNear = z2;
				regionNear = region2;

				xFar = x2;
				yFar = y2;
				zFar = z2;
				regionFar = region2;
			}
		}
	}

	// set projections
	// ===============
	if (useFarEdge)
	{
		xproj = xFar;
		yproj = yFar;
		zproj = zFar;
		regionproj = regionFar;
	}
	else
	{
		xproj = xNear;
		yproj = yNear;
		zproj = zNear;
		regionproj = regionNear;
	}

	return;
}

void WCSimGeometry::XYZtoUV(Int_t region, Double_t x, Double_t y, Double_t z, Double_t &u, Double_t &v)
{
	// default locations
	u = -99999.9;
	v = -99999.9;

	// Note:
	// ====
	//  assume that beam points along -ve x direction

	// CYLINDRICAL GEOMETRY
	// ====================
	if (fGeoType == WCSimGeometry::kCylinder)
	{

		// PMTs on top face
		if (region == WCSimGeometry::kTop)
		{
			u = y;
			v = +0.5 * fCylLength + fCylRadius + x;
		}

		// PMTs on side
		if (region == WCSimGeometry::kSide)
		{
			Double_t theta = 0.0;
			if (x != 0.0)
				theta = TMath::ATan(-y / x);
			if (x >= 0)
			{
				if (y > 0.0)
					theta += TMath::Pi();
				if (y < 0.0)
					theta -= TMath::Pi();
			}
			u = fCylRadius * theta;
			v = z;
		}

		// PMTs on bottom face
		if (region == WCSimGeometry::kBottom)
		{
			u = y;
			v = -0.5 * fCylLength - fCylRadius - x;
		}
	}

	// MAILBOX GEOMETRY
	// ================
	//
	// Default Coordinate System:
	//
	//               +ve z(Top)
	//                   ^
	//                   |  / +ve y(Left)
	//                   | /
	//    -ve x(Back) -- --> +ve x(Front)
	//                 / |
	//   -ve y(Right) /  |
	//                  -ve z(Bottom)
	//
	// Rolled out as follows:
	//
	//              -----
	//             |  T  |
	//    --- ----- ----- ----- ---
	//   | F |  R  |  B  |  L  | F |
	//    --- ----- ----- ----- ---
	//             |  B  |
	//              -----
	//

	if (fGeoType == WCSimGeometry::kMailBox)
	{

		// PMTs on top face
		if (region == WCSimGeometry::kTop)
		{
			u = y;
			v = +0.5 * fMailBoxZ + 0.5 * fMailBoxX + x;
		}

		// PMTs on front face
		if (region == WCSimGeometry::kFront)
		{
			if (y > 0.0)
				u = +fMailBoxX + fMailBoxY - y;
			if (y < 0.0)
				u = -fMailBoxX - fMailBoxY - y;
			v = z;
		}

		// PMTs on left face
		if (region == WCSimGeometry::kLeft)
		{
			u = +0.5 * fMailBoxX + 0.5 * fMailBoxY + x;
			v = z;
		}

		// PMTs on right face
		if (region == WCSimGeometry::kRight)
		{
			u = -0.5 * fMailBoxX - 0.5 * fMailBoxY - x;
			v = z;
		}

		// PMTs on back face
		if (region == WCSimGeometry::kBack)
		{
			u = y;
			v = z;
		}

		// PMTs on bottom face
		if (region == WCSimGeometry::kBottom)
		{
			u = y;
			v = -0.5 * fMailBoxZ - 0.5 * fMailBoxX - x;
		}
	}

	return;
}

void WCSimGeometry::FindCircle(Double_t x0, Double_t y0, Double_t z0, Double_t x1, Double_t y1, Double_t z1,
							   Double_t x2, Double_t y2, Double_t z2, Double_t &rx, Double_t &ry, Double_t &rz, Double_t &nx, Double_t &ny,
							   Double_t &nz, Double_t &r)
{
	Double_t centre[3] = {0.0, 0.0, 0.0};
	Double_t normal[3] = {0.0, 0.0, 0.0};
	Double_t radius = 0.0;

	Double_t p01[3] = {x0 - x1, y0 - y1, z0 - z1};
	Double_t p12[3] = {x1 - x2, y1 - y2, z1 - z2};
	Double_t p20[3] = {x2 - x0, y2 - y0, z2 - z0};

	Double_t n[3] = {p01[1] * p12[2] - p01[2] * p12[1], p01[2] * p12[0] - p01[0] * p12[2], p01[0] * p12[1] - p01[1] * p12[0]};

	Double_t D2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];

	if (D2 > 0.0)
	{
		Double_t a0 = -(p12[0] * p12[0] + p12[1] * p12[1] + p12[2] * p12[2]) * (p01[0] * p20[0] + p01[1] * p20[1] + p01[2] * p20[2]) / (2.0 * D2);
		Double_t a1 = -(p20[0] * p20[0] + p20[1] * p20[1] + p20[2] * p20[2]) * (p12[0] * p01[0] + p12[1] * p01[1] + p12[2] * p01[2]) / (2.0 * D2);
		Double_t a2 = -(p01[0] * p01[0] + p01[1] * p01[1] + p01[2] * p01[2]) * (p20[0] * p12[0] + p20[1] * p12[1] + p20[2] * p12[2]) / (2.0 * D2);
		Double_t D = sqrt(D2);

		centre[0] = a0 * x0 + a1 * x1 + a2 * x2;
		centre[1] = a0 * y0 + a1 * y1 + a2 * y2;
		centre[2] = a0 * z0 + a1 * z1 + a2 * z2;

		radius = sqrt(
					 (p01[0] * p01[0] + p01[1] * p01[1] + p01[2] * p01[2]) * (p12[0] * p12[0] + p12[1] * p12[1] + p12[2] * p12[2]) * (p20[0] * p20[0] + p20[1] * p20[1] + p20[2] * p20[2])) /
				 (2.0 * D);

		if (n[0] * centre[0] + n[1] * centre[1] + n[2] * centre[2] > 0.0)
		{
			normal[0] = +n[0] / D;
			normal[1] = +n[1] / D;
			normal[2] = +n[2] / D;
		}
		else
		{
			normal[0] = -n[0] / D;
			normal[1] = -n[1] / D;
			normal[2] = -n[2] / D;
		}
	}

	rx = centre[0];
	ry = centre[1];
	rz = centre[2];

	nx = normal[0];
	ny = normal[1];
	nz = normal[2];

	r = radius;

	return;
}

void WCSimGeometry::FindCircleOld(Double_t xp, Double_t yp, Double_t zp, Double_t x0, Double_t y0, Double_t z0,
								  Double_t angle_degrees, Double_t omega_degrees, Double_t &rx, Double_t &ry, Double_t &rz, Double_t &nx,
								  Double_t &ny, Double_t &nz, Double_t &r)
{
	// inputs
	// ======
	Double_t phi = 0.0;
	Double_t theta = 0.0;
	Double_t angle = (TMath::Pi() / 180.0) * angle_degrees;
	Double_t omega = (TMath::Pi() / 180.0) * omega_degrees;

	Double_t x = xp;
	Double_t y = yp;
	Double_t z = zp;

	Double_t xtemp = 0.0;
	Double_t ytemp = 0.0;
	Double_t ztemp = 0.0;

	Double_t dx = 0.0;
	Double_t dy = 0.0;
	Double_t dz = 0.0;
	Double_t ds = 0.0;

	Double_t radius = 0.0;

	// subtract vertex
	// ===============
	x -= x0;
	y -= y0;
	z -= z0;

	// forward rotation (x,y,z)->(x,0,z)
	// =================================
	phi = 0.0;

	if (x != 0.0)
	{
		phi = atan(y / x);
	}
	if (x <= 0.0)
	{
		if (y > 0.0)
			phi += TMath::Pi();
		if (y < 0.0)
			phi -= TMath::Pi();
	}

	xtemp = x;
	ytemp = y;
	ztemp = z;

	x = +xtemp * cos(phi) + ytemp * sin(phi);
	y = -xtemp * sin(phi) + ytemp * cos(phi);
	z = +ztemp;

	// forward rotation (x,0,z)->(0,0,z)
	// =================================
	theta = 0.0;

	if (z != 0.0)
	{
		theta = atan(x / z);
	}
	if (z <= 0.0)
	{
		if (x > 0.0)
			theta += TMath::Pi();
		if (x < 0.0)
			theta -= TMath::Pi();
	}

	xtemp = x;
	ytemp = y;
	ztemp = z;

	x = +xtemp * cos(theta) - ztemp * sin(theta);
	y = +ytemp;
	z = +xtemp * sin(theta) + ztemp * cos(theta);

	// apply rotation: cone angle + azimuthal angle
	// ============================================
	xtemp = x;
	ytemp = y;
	ztemp = z;

	x = ztemp * sin(angle) * cos(omega);
	y = ztemp * sin(angle) * sin(omega);
	z = ztemp * cos(angle);

	radius = ztemp * sin(angle);

	// backward rotation (0,0,z)->(x,0,z)
	// ==================================
	xtemp = x;
	ytemp = y;
	ztemp = z;

	x = +xtemp * cos(-theta) - ztemp * sin(-theta);
	y = +ytemp;
	z = +xtemp * sin(-theta) + ztemp * cos(-theta);

	// backward rotation (x,0,z)->(x,y,z)
	// =================================
	xtemp = x;
	ytemp = y;
	ztemp = z;

	x = +xtemp * cos(-phi) + ytemp * sin(-phi);
	y = -xtemp * sin(-phi) + ytemp * cos(-phi);
	z = +ztemp;

	// add on vertex
	// =============
	x += x0;
	y += y0;
	z += z0;

	// return coordinates
	// ==================
	rx = x;
	ry = y;
	rz = z;

	// return normal
	// =============
	dx = rx - x0;
	dy = ry - y0;
	dz = rz - z0;
	ds = sqrt(dx * dx + dy * dy + dz * dz);

	nx = dx / ds;
	ny = dy / ds;
	nz = dz / ds;

	// return radius
	// =============
	r = radius;

	return;
}

void WCSimGeometry::FindCircle(Double_t xp, Double_t yp, Double_t zp, Double_t x0, Double_t y0, Double_t z0,
							   Double_t angle_degrees, Double_t omega_degrees, Double_t &rx, Double_t &ry, Double_t &rz, Double_t &nx,
							   Double_t &ny, Double_t &nz, Double_t &r)
{
	// (xp,yp,zp) is position of the pmt digit
	// (x0,y0,z0) is position of the vertex

	// default
	// =======
	rx = xp;
	ry = yp;
	rz = zp;

	nx = 0.0;
	ny = 0.0;
	nz = 0.0;

	r = 0.0;

	// inputs
	// ======
	Double_t angle = (TMath::Pi() / 180.0) * angle_degrees; // radians
	Double_t omega = (TMath::Pi() / 180.0) * omega_degrees; // radians

	// Gives the vector from the vertex to the pmt
	Double_t dx = xp - x0;
	Double_t dy = yp - y0;
	Double_t dz = zp - z0;
	Double_t ds = sqrt(dx * dx + dy * dy + dz * dz); // distance from vertex

	// Becomes unit vector from vertex to digit
	Double_t px = 0.0;
	Double_t py = 0.0;
	Double_t pz = 0.0;

	// Normalized direction vector from vertex to pmt
	if (ds > 0.0)
	{
		px = dx / ds;
		py = dy / ds;
		pz = dz / ds;
	}
	else
	{
		return;
	}

	// do the rotations
	// ================
	TVector3 myVtx(x0, y0, z0);
	TVector3 myDir(px, py, pz);

	TVector3 initDir = myDir;
	TVector3 orthDir = myDir.Orthogonal();
	myDir.Rotate(angle, orthDir);
	myDir.Rotate(omega, initDir);
	// Rotate away from the digit by the cone angle, then around by omega, which runs up to 2pi to draw a cone of given angle around the digit

	// outputs
	// =======
	rx = x0 + ds * myDir.x(); // This is making a cone around the pmt hit starting at the vertex
	ry = y0 + ds * myDir.y(); // whose cone angle is "angle"
	rz = z0 + ds * myDir.z(); // We increment omega to step round the edge of the cone
							  // For each step we return the (x,y,z) coordinate of the point on the cone
							  // rim which gets converted to (phi, costheta) elsewehere, and the cone radius
							  // (which gives the Cherenkov angle)

	nx = myDir.x();
	ny = myDir.y(); // n is the unit vector from the vertex to the rim of the cone
	nz = myDir.z();

	r = ds * sin(angle); // the radius from the PMT hit to the edge of the cone

	return;
}

void WCSimGeometry::FindVertex(Double_t x0, Double_t y0, Double_t z0, Double_t t0, Double_t x1, Double_t y1,
							   Double_t z1, Double_t t1, Double_t x2, Double_t y2, Double_t z2, Double_t t2, Double_t x3, Double_t y3,
							   Double_t z3, Double_t t3, Double_t &vxm, Double_t &vym, Double_t &vzm, Double_t &vtm, Double_t &vxp,
							   Double_t &vyp, Double_t &vzp, Double_t &vtp)
{
	// switch off error messages
	// =========================
	// suppress messages like:
	//  Error in <TDecompLU::DecomposeLUCrout>: matrix is singular
	Int_t myIgnoreLevel = gErrorIgnoreLevel;
	gErrorIgnoreLevel = kFatal;

	// default vertex
	// ==============
	vxm = -99999.9;
	vym = -99999.9;
	vzm = -99999.9;
	vtm = -99999.9;

	vxp = -99999.9;
	vyp = -99999.9;
	vzp = -99999.9;
	vtp = -99999.9;

	// speed of light in water
	// =======================
	Double_t c = 29.98 / 1.33; // cm/ns [water]

	// causality checks
	// ================
	if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0) >= c * c * (t1 - t0) * (t1 - t0) && (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1) >= c * c * (t2 - t1) * (t2 - t1) && (x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2) >= c * c * (t3 - t2) * (t3 - t2) && (x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0) + (z2 - z0) * (z2 - z0) >= c * c * (t2 - t0) * (t2 - t0) && (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1) + (z3 - z1) * (z3 - z1) >= c * c * (t3 - t1) * (t3 - t1) && (x3 - x0) * (x3 - x0) + (y3 - y0) * (y3 - y0) + (z3 - z0) * (z3 - z0) >= c * c * (t3 - t0) * (t3 - t0))
	{

		// [Note: for causality, require that |x_{i}-x_{j}| >= c*|t_{i}-t_{j}|
		//        for each pair of points]

		Double_t dx1 = x1 - x0;
		Double_t dy1 = y1 - y0;
		Double_t dz1 = z1 - z0;
		Double_t dt1 = c * (t1 - t0);
		Double_t dx2 = x2 - x0;
		Double_t dy2 = y2 - y0;
		Double_t dz2 = z2 - z0;
		Double_t dt2 = c * (t2 - t0);
		Double_t dx3 = x3 - x0;
		Double_t dy3 = y3 - y0;
		Double_t dz3 = z3 - z0;
		Double_t dt3 = c * (t3 - t0);

		Double_t epsilon = 1.0e-7;

		// check that points don't all lie in a plane
		if (!(fabs(dx1) < epsilon && fabs(dx2) < epsilon && fabs(dx3) < epsilon) && !(fabs(dy1) < epsilon && fabs(dy2) < epsilon && fabs(dy3) < epsilon) && !(fabs(dz1) < epsilon && fabs(dz2) < epsilon && fabs(dz3) < epsilon) && !(fabs(dx1) < epsilon && fabs(dy1) < epsilon && fabs(dz1) < epsilon) && !(fabs(dx2) < epsilon && fabs(dy2) < epsilon && fabs(dz2) < epsilon) && !(fabs(dx3) < epsilon && fabs(dy3) < epsilon && fabs(dz3) < epsilon))
		{

			// [Note: this is a problem for detectors with flat faces!]

			Double_t Mdata[9] = {dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3};

			Double_t Qdata[3] = {0.5 * (dx1 * dx1 + dy1 * dy1 + dz1 * dz1 - dt1 * dt1), 0.5 * (dx2 * dx2 + dy2 * dy2 + dz2 * dz2 - dt2 * dt2), 0.5 * (dx3 * dx3 + dy3 * dy3 + dz3 * dz3 - dt3 * dt3)};

			Double_t Tdata[3] = {dt1, dt2, dt3};

			TMatrixD M(3, 3, Mdata);
			TMatrixD Q(3, 1, Qdata);
			TMatrixD T(3, 1, Tdata);

			if (M.Determinant() != 0.0)
			{

				TMatrixD A(3, 1);
				TMatrixD B(3, 1);

				M.Invert();
				A.Mult(M, T);
				B.Mult(M, Q);

				Double_t ax = A(0, 0);
				Double_t ay = A(1, 0);
				Double_t az = A(2, 0);

				Double_t bx = B(0, 0);
				Double_t by = B(1, 0);
				Double_t bz = B(2, 0);

				Double_t ab = ax * bx + ay * by + az * bz;
				Double_t a2 = ax * ax + ay * ay + az * az;
				Double_t b2 = bx * bx + by * by + bz * bz;

				Double_t qa = a2 - 1.0;
				Double_t qb = 2.0 * ab;
				Double_t qc = b2;

				// check for solutions
				if (qb * qb - 4.0 * qa * qc > 0.0)
				{

					// The common vertex is given by a quadratic equation, which has two solutions.
					// Typically, one solution corresponds to photons travelling forwards in time,
					// and the other solution corresponds to photons travelling backwards in time.
					// However, sometimes there appear to be two valid solutions.

					Double_t ctm = (-qb - sqrt(qb * qb - 4.0 * qa * qc)) / (2.0 * qa);
					Double_t ctp = (-qb + sqrt(qb * qb - 4.0 * qa * qc)) / (2.0 * qa);

					Double_t tm = t0 + ctm / c;
					Double_t xm = x0 + ctm * ax + bx;
					Double_t ym = y0 + ctm * ay + by;
					Double_t zm = z0 + ctm * az + bz;
					Bool_t foundVertexM = 0;

					if (tm < t0 && tm < t1 && tm < t2 && tm < t3)
					{
						vxm = xm;
						vym = ym;
						vzm = zm;
						vtm = tm;
						foundVertexM = 1;
					}

					Double_t tp = t0 + ctp / c;
					Double_t xp = x0 + ctp * ax + bx;
					Double_t yp = y0 + ctp * ay + by;
					Double_t zp = z0 + ctp * az + bz;
					Bool_t foundVertexP = 0;

					if (tp < t0 && tp < t1 && tp < t2 && tp < t3)
					{
						vxp = xp;
						vyp = yp;
						vzp = zp;
						vtp = tp;
						foundVertexP = 1;
					}

					// std::cout << "  Vertex: [a=" << qa << ",b=" << qb << ",c=" << qc << "] [ctm=" << ctm << ",ctp=" << ctp << "] " << std::endl
					//           << "   Low T:  (x,y,z,t)=(" << xm << "," << ym << "," << zm << "," << tm << ") [foundVtx=" << foundVertexM << "] " << std::endl
					//           << "   High T: (x,y,z,t)=(" << xp << "," << yp << "," << zp << "," << tp << ") [foundVtx=" << foundVertexP << "] " << std::endl;
				}
			}
		}
	}

	// switch on error messages
	// ========================
	gErrorIgnoreLevel = myIgnoreLevel;

	return;
}

void WCSimGeometry::DistanceToIntersectLine(Double_t x0, Double_t y0, Double_t z0, Double_t vx, Double_t vy,
											Double_t vz, Double_t ex, Double_t ey, Double_t ez, Double_t &x, Double_t &y, Double_t &z, Double_t &L)
{
	Double_t dx = 0.0;
	Double_t dy = 0.0;
	Double_t dz = 0.0;
	Double_t ds = 0.0;

	Double_t px = 0.0;
	Double_t py = 0.0;
	Double_t pz = 0.0;
	Double_t dsTrack = 0.0;

	Double_t qx = 0.0;
	Double_t qy = 0.0;
	Double_t qz = 0.0;
	Double_t dsPmt = 0.0;

	// vector from vertex to end
	dx = ex - vx;
	dy = ey - vy;
	dz = ez - vz;
	ds = sqrt(dx * dx + dy * dy + dz * dz);

	if (ds > 0.0)
	{
		px = dx / ds;
		py = dy / ds;
		pz = dz / ds;
		dsTrack = ds;
	}

	// vector from vertex to pmt
	dx = x0 - vx;
	dy = y0 - vy;
	dz = z0 - vz;
	ds = sqrt(dx * dx + dy * dy + dz * dz);

	if (ds > 0.0)
	{
		qx = dx / ds;
		qy = dy / ds;
		qz = dz / ds;
		dsPmt = ds;
	}

	// Cherenkov geometry
	// ==================
	Double_t cosphi = px * qx + py * qy + pz * qz;
	Double_t phi = acos(cosphi);
	Double_t theta = WCSimParameters::CherenkovAngle() * TMath::Pi() / 180.0;

	Double_t Ltrack = 0.0;
	Double_t Lphoton = 0.0;

	if (phi < theta)
	{
		Ltrack = +dsPmt * sin(theta - phi) / sin(theta);
		Lphoton = +dsPmt * sin(phi) / sin(theta);
	}
	else
	{
		Ltrack = -dsPmt * sin(phi - theta) / sin(phi);
		Lphoton = -dsPmt * sin(phi) / sin(theta);
	}

	// return intersection
	// ===================
	x = vx + px * Ltrack;
	y = vy + py * Ltrack;
	z = vz + pz * Ltrack;
	L = Lphoton;

	return;
}

Double_t WCSimGeometry::DistanceToIntersectLine(Double_t x0, Double_t y0, Double_t z0, Double_t sx, Double_t sy,
												Double_t sz, Double_t ex, Double_t ey, Double_t ez, Double_t &x, Double_t &y, Double_t &z)
{
	Double_t L = -1.0;

	WCSimGeometry::DistanceToIntersectLine(x0, y0, z0, sx, sy, sz, ex, ey, ez, x, y, z, L);
	return L;
}

Double_t WCSimGeometry::DistanceToIntersectLine(Double_t *pos, Double_t *start, Double_t *end, Double_t *intersection)
{
	Double_t x = -999.9;
	Double_t y = -999.9;
	Double_t z = -999.9;
	Double_t L = -1.0;

	WCSimGeometry::DistanceToIntersectLine(pos[0], pos[1], pos[2], start[0], start[1], start[2], end[0], end[1], end[2],
										   x, y, z, L);

	intersection[0] = x;
	intersection[1] = y;
	intersection[2] = z;

	return L;
}

WCSimGeometry::RingRegion_t WCSimGeometry::GetPMTRingRegion(const TVector3 &vtx, const TVector3 &dir,
															const TVector3 &pmt, const TrackType::Type &type, const double energy)
{
	TouchEmissionProfileManager();

	double length = fEmissionProfileManager->GetTrackLength(type, energy);
	TVector3 stop = vtx + length * dir;

	double startAngle = (pmt - vtx).Angle(dir) * TMath::RadToDeg();
	double stopAngle = (pmt - stop).Angle(dir) * TMath::RadToDeg();
	//std::cout << "Start at (" << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << ") and " << TrackType(type) << " with E = " << energy << " travels " << length << "cm" << std::endl;
	//std::cout << "Stops at (" << stop.X() << ", " << stop.Y() << ", " << stop.Z() << ") so start angle = " << startAngle << " and stopAngle = " << stopAngle << std::endl;

	if (startAngle < 42 && stopAngle > 42)
	{
		return WCSimGeometry::kInside;
	} // It's in the body of the ring
	if (startAngle < 42 && stopAngle < 42)
	{
		return WCSimGeometry::kInsideHole;
	} // It's in the hole in the middle of the ring
	//std::cout << "Outside" << std::endl << std::endl;
	return kOutside; // It isn't in the ring at all
}

bool WCSimGeometry::IsOutsideRing(const TVector3 &vtx, const TVector3 &dir, const TVector3 &pmt,
								  const TrackType::Type &type, const double energy)
{
	int region = (GetPMTRingRegion(vtx, dir, pmt, type, energy));
	return (region == WCSimGeometry::kOutside);
}

bool WCSimGeometry::IsInRing(const TVector3 &vtx, const TVector3 &dir, const TVector3 &pmt, const TrackType::Type &type,
							 const double energy)
{
	return (GetPMTRingRegion(vtx, dir, pmt, type, energy) == WCSimGeometry::kInside);
}

bool WCSimGeometry::IsInRingHole(const TVector3 &vtx, const TVector3 &dir, const TVector3 &pmt,
								 const TrackType::Type &type, const double energy)
{
	return (GetPMTRingRegion(vtx, dir, pmt, type, energy) == WCSimGeometry::kInsideHole);
}

void WCSimGeometry::TouchEmissionProfileManager()
{
	if (fEmissionProfileManager == NULL)
	{
		fEmissionProfileManager = new WCSimEmissionProfileManager();
	}
}

bool WCSimGeometry::GetEscapes(const TVector3 &vtx, const TVector3 &dir, const TrackType::Type &type,
							   const double energy)
{
	bool startIn = InsideDetector(vtx.X(), vtx.Y(), vtx.Z());
	TVector3 end = GetEnd(vtx, dir, type, energy);
	bool endIn = InsideDetector(end.X(), end.Y(), end.Z());
	return !(startIn && endIn);
}

TVector3 WCSimGeometry::GetEnd(const TVector3 &vtx, const TVector3 &dir, const TrackType::Type &type,
							   const double energy)
{
	TouchEmissionProfileManager();
	double length = fEmissionProfileManager->GetTrackLength(type, energy);
	TVector3 end = vtx + length * dir;
	return end;
}
