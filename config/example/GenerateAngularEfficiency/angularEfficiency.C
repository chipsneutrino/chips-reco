///////////////////////////////////////////////////////////////////////////////
//  Author: S. Germani  (s.germani@ucl.ac.uk) and J. Tingey (j.tingey.16@ucl.ac.uk) UCL
//  First  version             19/05/2017
//  Second version             28/11/2017 --> Change rotation scheme: rotate starting plane(position)
//  Third  version			   17/07/2018
// ----------------------------------------------------------------------------
//  This macro makes use of the ROOT Geometry package  to define a PMT and
//  a Light Collector (Winstone Cone) volume.
//  Light rays are generated from a source and tracked, according to the
//  geometrcal optics laws, to the PMT in order to evaluate the overall(PMT+LC)
//  angular efficiency to be used in the Reconstruction.
//  The output is a ROOT file containing histograms with information about the
//  generated and detected photons.
//
// ----------------------  PMT Geometry ---------------------------------------
//
//  The PMT geometry is a spherical shell defined by two parameters as in the
//  Geant4 (WCSim) simulation:
//  - PMT_D             which defines the PMT Diameter
//  - PMT_EXP_HEIGHT    which defines the PMT Exposure Height
//
// ----------------------  LC  Geometry ---------------------------------------
//
//  The Collector geometry is loaded from a configuration file which need to
//  be produced using other tools.
//
//  The Light Collector optical properties are define by the variables:
//  - LC_REFLECTIVITY   which the deifine the probability for a photon to be
//                      reflected (not absorbed) with range 0-1.
//  - LC_ROUGHNESS      which define the sigma for the gaussian smearing of
//                      the components of the normal to the surface to add
//                      a diffusing effect related to the surface roughness.
//                      This is a first attempt to emulate roughness and may
//                      need to be improved.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>

#include "TVector3.h"

#include "TGeometry.h"
#include "TGeoManager.h"

//####################################################################################
//####################  Functions  (defined at the bottom) ###########################
//####################################################################################

void read_lcfile(string fname, Double_t *lcz, Double_t *lcr);

//####################################################################################

//####################################################################################
//####################   Main Function    ############################################
//####################################################################################

void angularEfficiency(int iOpt = -1) {

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////      Load Libraries                         ///////////////////
	////////////////////////////////////////////////////////////////////////////////////

	gSystem->Load("libGeom");

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////      Dealing with Options                   ///////////////////
	////////////////////////////////////////////////////////////////////////////////////

	//////  Check what the chosen Option means

	enum runOptions {
		kDrawOnly, kAngleScan, kNumRunOptions
	};
	string runOptName[] = { "DrawOnly", "AngleScan" };

	int runOption = 0;
	int kVerbose = 0;
	int makeGifs = 0;

	int use88mm = 0;
	int useR6091 = 1;
	int useLC = 1;

	//////  Print out Option meaning
	string geoName;
	string runOptPrint[] = { "Draw Only", "Position Scan - Pencil Beam", "Position Scan - Light Cone", "Angle Scan",
			"Diffuse Posistion", "Diffuse Position and Angle" };

	cout << "\n=========== OPTION SUMMARY ===================" << endl;
	if (useR6091) {
		cout << "    R6091 PMT            " << endl;
		geoName += "3inchPMT";
	}
	if (use88mm) {
		cout << "    88-mm  PMT            " << endl;
		geoName += "88mmPMT";
	}
	if (useLC) {
		cout << "      + Light Collector   " << endl;
		geoName += "_LC";
	}
	cout << "    " << runOptPrint[runOption] << endl;
	if (kVerbose)
		cout << "     Verbose Mode " << endl;
	cout << "==============================================\n" << endl;

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////      Light Collector File(s) and Properties ///////////////////
	////////////////////////////////////////////////////////////////////////////////////

	//////  Generic Properties
	///// Light Collector  Optical Properties
	Double_t LC_REFLECTIVITY = 0.95;  // Using 90% fot now
	Double_t LC_ROUGHNESS = 0.000000001;  // Sigma applied to surface normal versor components

	///////////  PMT  Collection Efficiency ///////////
	Double_t collection_angle[10] = { 0, 10, 20, 30, 40, 50, 60, 70, 74, 90 };
	Double_t collection_eff[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.2 };

	//////  Units
	Double_t mm = 1. / 1000;

	////// Thickness of wrappings
	Double_t alu_thickness;
	Double_t pla_thickness;

	//////  Specific Properties
	//////  Just one file for the time being ...
	string lc_R6091_file = "R6091_LC_v1.txt";
	string lc_88mm_file = "88mm_LC_v2.txt";

	//////  Arrays to store Light Collector dimensions
	Double_t LCz[100];
	Double_t LCrmin[100];

	//////  PMT geometry variables
	Double_t PMT_D;					// PMT Diameter
	Double_t PMT_EXP_HEIGHT;		// PMT exposure height
	Double_t R_SPHERE;				// Sphere Radius
	Double_t UNSENSITIVE_HEIGHT;  	// PMT Unsensitive height
	Double_t LCR;					// Light Collector Maximum Radius

	if (useR6091) {
		std::cout << "Using R6091 PMT..." << std::endl;
		read_lcfile(lc_R6091_file, LCz, LCrmin);
	} else if (use88mm) {
		std::cout << "Using 88mm PMT..." << std::endl;

		read_lcfile(lc_88mm_file, LCz, LCrmin);

		PMT_D = 88.0;
		PMT_EXP_HEIGHT = 25.5;
		R_SPHERE = 52.4;
		UNSENSITIVE_HEIGHT = 2.8 * mm;
		LCR = LCrmin[99];

		alu_thickness = 0.02 * mm;
		pla_thickness = 2 * mm;

	} else {
		std::cout << "Unknown geometry!" << std::endl;

		read_lcfile(lc_R6091_file, LCz, LCrmin);

		PMT_D = 65.0;
		PMT_EXP_HEIGHT = 8.0;
		R_SPHERE = 0.0;
		UNSENSITIVE_HEIGHT = 0.0 * mm;
		LCR = LCrmin[99];

		alu_thickness = 0.02 * mm;
		pla_thickness = 2 * mm;
		return;
	}

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////      Define Geometry and Run Constants      ///////////////////
	////////////////////////////////////////////////////////////////////////////////////

	int Npos = 0;
	int Nphotons = 500000;

	//Double_t AMAX = 90;  // (deg)  Angle scan maximum angle wrt normal
	//Double_t DELTAA = 5;  // (deg)  Angle scan angular step

	Double_t AMAX = 90;  // (deg)  Angle scan maximum angle wrt normal
	Double_t DELTAA = 10;  // (deg)  Angle scan angular step

	if (runOption == kAngleScan) {
		Npos = int(AMAX / DELTAA) + 1;
	}
	if (Npos == 1) {
		Nphotons *= 100; // Random position is inefficient, need more photons
	}

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////     Output File and Histograms              ///////////////////
	////////////////////////////////////////////////////////////////////////////////////

	//////  Output File
	TFile *fout;
	if (runOption > 0)
		fout = TFile::Open(Form("lightCollectorPhotons_%s_%s.root", runOptName[runOption].c_str(), geoName.c_str()),
				"RECREATE");

	// Histograms for every angle tested...
	TH2F *h_sxypos[100];      		// all photons start x-y position

	TH2F *hDet_xypos[100];  		// Detection x-y positions
	TH2F *hDet_first_xypos[100];  	// direct detected photons x-y position
	TH2F *hDet_multiple_xypos[100];  	// reflected detected photons x-y position
	TH2F *hNon_xypos[100];  		// all non detected photons x-y position
	TH2F *hNon_reflect_xypos[100];  	// reflected non detected photons x-y position

	// Histograms across the whole scan...
	TH1F *hndet;             		// Number of detected photons (useful for scans and efficiency comparisons)
	TH1F *hnhit;             		// Number of photons hitting the LC+PMT

	TH1F *hnTop;             		// Number of photons hitting the Top volume
	TH1F *hnPmt;             		// Number of photons hitting the PMT volume
	TH1F *hnSup;             		// Number of photons hitting the Support volume
	TH1F *hnRef;             		// Number of photons hitting the Reflecting volume
	TH1F *hnAbs;             		// Number of photons hitting the Absorbing volume

	if (runOption == kAngleScan) {
		hndet = new TH1F("hndet", ";Beam Angle [deg];Detected Photons", Npos, -DELTAA / 2, AMAX + DELTAA / 2);
		hnhit = new TH1F("hnhit", ";Beam Angle [deg];Photons Hitting LC+PMT", Npos, -DELTAA / 2, AMAX + DELTAA / 2);

		hnTop = new TH1F("hnTop", ";Beam Angle [deg];Photons Hitting Top", Npos, -DELTAA / 2, AMAX + DELTAA / 2);
		hnPmt = new TH1F("hnPmt", ";Beam Angle [deg];Photons Hitting PMT", Npos, -DELTAA / 2, AMAX + DELTAA / 2);
		hnSup = new TH1F("hnSup", ";Beam Angle [deg];Photons Hitting Support", Npos, -DELTAA / 2, AMAX + DELTAA / 2);
		hnRef = new TH1F("hnRef", ";Beam Angle [deg];Photons Hitting Reflective", Npos, -DELTAA / 2, AMAX + DELTAA / 2);
		hnAbs = new TH1F("hnAbs", ";Beam Angle [deg];Photons Hitting Absorbing", Npos, -DELTAA / 2, AMAX + DELTAA / 2);
	}

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////    Define Geometry Objects                  ///////////////////
	////////////////////////////////////////////////////////////////////////////////////

	//////  Material Colours
	const Color_t PlasticColor = 1;
	const Color_t AluColor = 9;
	const Color_t GlassColor = 3;
	const Color_t AcrylicColor = kYellow;

	//////  Geometry Manager
	TGeoManager *geom = new TGeoManager("Assemblies", "Geometry using assemblies");

	//////  Define Materials
	TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98, 13, 2.7);
	TGeoMaterial *matAcrylic = new TGeoMaterial("Acrylic", 13, 6, 1.);

	//////  Define Media
	TGeoMedium *Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
	TGeoMedium *Al = new TGeoMedium("Aluminium", 2, matAl);
	TGeoMedium *Acrylic = new TGeoMedium("Acrylic", 3, matAcrylic);

	//////  Make the top container volume
	TGeoVolume *top = geom->MakeBox("TOP", Vacuum, 500000., 500000., 500000.);
	geom->SetTopVolume(top);

	//////
	//////  Define PMT surface as a Sphere Cap (it allows to simulate both flat and spherical PMTs)
	//////
	double dR = PMT_EXP_HEIGHT * mm;
	double pmtR = PMT_D / 2 * mm;

	double Rsphere;
	double beta;
	if (use88mm) {
		double Rsphere = R_SPHERE * mm;
		double beta = TMath::ASin(2 * (pmtR) * dR / (pow(pmtR, 2) + pow(dR, 2)));
	} else if (useR6091) {
		double beta = 2*TMath::ATan2(dR, pmtR);
		double Rsphere = 0.5*PMT_D/sin(beta);
	}

	//std::cout << R_SPHERE * mm << std::endl;
	//std::cout << TMath::ASin(2 * (pmtR) * dR / (pow(pmtR, 2) + pow(dR, 2))) << std::endl;
	//std::cout << R_SPHERE * mm << std::endl;
	//std::cout << 2*TMath::ATan2(dR, pmtR) << std::endl;

	//   Double_t angle = TMath::ASin(pmtR / Rsphere * 180 / TMath::Pi());
	Double_t angle = beta * 180 / TMath::Pi();

	cout << "PMT size " << dR << " " << pmtR << " " << Rsphere << "  |  " << LCR << endl;

	TGeoVolume *pmt = geom->MakeSphere("PMT", Al, Rsphere - 0.5 * mm, Rsphere, 0, angle, 0, 360);
	pmt->SetLineColor(GlassColor);

	//////  Move PMT to the correct position
	TGeoTranslation *tpmt = new TGeoTranslation("tpmt", 0.0, 0.0, dR - Rsphere + UNSENSITIVE_HEIGHT);
	tpmt->RegisterYourself();

	//////  Add PMT to top container volume
	top->AddNode(pmt, 1, tpmt);

	/////// Define PMT support base ring ...
	TGeoVolume *supvol = gGeoManager->MakePcon("SUPPORT", Al, 0, 360, 2);
	TGeoPcon *suppcon = (TGeoPcon*) (supvol->GetShape());
	suppcon->DefineSection(0, LCz[0] - 5 * mm, 0, LCrmin[0]);
	suppcon->DefineSection(1, LCz[0] - 0.5 * mm, 0, LCrmin[0]);

	//////  Define Light Collector Reflective and Absorbing (Structure) Volume as Poliycones
	TGeoVolume *refvol = gGeoManager->MakePcon("REFPCON", Al, 0, 360, 100);
	TGeoVolume *absvol = gGeoManager->MakePcon("ABSPCON", Al, 0, 360, 100);

	TGeoPcon *refpcon = (TGeoPcon*) (refvol->GetShape());
	TGeoPcon *abspcon = (TGeoPcon*) (absvol->GetShape());

	/////  Assign Polycones surface points
	for (int ip = 0; ip < 100; ip++) {
		refpcon->DefineSection(ip, LCz[ip], LCrmin[ip], LCrmin[ip] + alu_thickness);
		abspcon->DefineSection(ip, LCz[ip], LCrmin[ip] + alu_thickness, LCrmin[ip] + alu_thickness + pla_thickness);
	}  // for ip

	//////  Assign Polycone colours and add to top container volume
	refvol->SetLineColor(AluColor);
	absvol->SetLineColor(PlasticColor);
	if (useLC) {
		top->AddNode(supvol, 1);
		top->AddNode(refvol, 1);
		top->AddNode(absvol, 1);
	}

	////// Close the Geometry  ---------------------------------------------
	geom->CloseGeometry();

	geom->SetVisLevel(4);
	geom->SetVisOption(0);
	top->Draw("ogl");

	if (runOption > 0) {
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////     Generate and Track Photons              ///////////////////
		////////////////////////////////////////////////////////////////////////////////////

		//////  Random Generator
		TRandom3 *tr = new TRandom3();

		//////  Define starting Photon position components
		Double_t px = 0, py = 0, pz = 200 * mm;

		//////  Define starting Photon direction components
		Double_t dx = 0, dy = 0, dz = -1;

		/// default direction is vertical
		Double_t theta = TMath::Pi();
		Double_t phi = 0;

		TCanvas * c1 = new TCanvas("canvas", "canvas", 600, 600);

		double scale = 0;
		double zeroNum = 0;

		//////  Loop over PMT Angular Positions  ---------------------------------
		for (int ipos = 0; ipos < Npos; ipos++) {
			cout << " Running Step " << ipos << " ..." << endl;

			// Start position plots
			h_sxypos[ipos] = new TH2F(Form("h_sxypos%02d", ipos), ";Start X [m];Start Y [m]", 500, -5 * LCR, 5 * LCR,
					100, -5 * LCR, 5 * LCR);

			// End position plots
			hDet_xypos[ipos] = new TH2F(Form("hDet_xypos%02d", ipos), ";End X [m];End Y [m]", 100, -LCR, LCR,
					100, -LCR, LCR);
			hDet_first_xypos[ipos] = new TH2F(Form("hDet_first_xypos%02d", ipos), ";End X [m];End Y [m]", 100, -LCR, LCR,
					100, -LCR, LCR);
			hDet_multiple_xypos[ipos] = new TH2F(Form("hDet_multiple_xypos%02d", ipos), ";End X [m];End Y [m]", 100, -LCR, LCR,
					100, -LCR, LCR);
			hNon_xypos[ipos] = new TH2F(Form("hNon_xypos%02d", ipos), ";End X [m];End Y [m]", 100, -LCR, LCR,
					100, -LCR, LCR);
			hNon_reflect_xypos[ipos] = new TH2F(Form("hNon_reflect_xypos%02d", ipos), ";End X [m];End Y [m]", 100, -LCR, LCR,
					100, -LCR, LCR);

			if (runOption == kAngleScan) {
				////// Define starting  angle
				theta = TMath::Pi() - ipos * DELTAA * TMath::Pi() / 180;
				cout << " --> Photon angle   : " << ipos * DELTAA << " deg " << endl;
			} // if runOption

			//// Get position for filling whole scan histograms
			Double_t nDetHistoX = hndet->GetXaxis()->GetBinCenter(ipos + 1);

			////// Initialise photon counters
			int nDet = 0;
			int nHit = 0;

			int nTop = 0;
			int nPmt = 0;
			int nSup = 0;
			int nRef = 0;
			int nAbs = 0;

			std::vector < TVector3 > refl_points;
			vector<TVector3>::iterator vit;

			//////  Loop over emitted photons     ------------------
			for (int iph = 0; iph < Nphotons; iph++) {

				// Clear the reflections vector
				refl_points.clear();

				//// Random Position if needed
				if (runOption == kAngleScan) {
					// Just have them uniform over large area...
					px = tr->Uniform(-2 * LCR, 2 * LCR);
					py = tr->Uniform(-2 * LCR, 2 * LCR);
					pz = 200 * mm;

					TVector3 *startV = new TVector3(px,py,pz);
					startV->RotateY((-ipos)*DELTAA*TMath::Pi()/180);

					px = startV->X();
					py = startV->Y();
					pz = startV->Z();
				}

				dx = sin(theta) * cos(phi);
				dy = sin(theta) * sin(phi);
				dz = cos(theta);

				geom->ResetState();
				geom->InitTrack(px, py, pz, dx, dy, dz);

				h_sxypos[ipos]->Fill(px, py);

				// ---- Start Tracking -------------

				if (kVerbose)
					cout << " ============    New Track " << iph << " " << geom->GetCurrentVolume()->GetName()
							<< " =========== " << endl;

				//////  Track Photon ....
				bool Stop = false;
				bool detected = false;
				int reflections = 0;

				while (!Stop) {
					Double_t *ppos = (Double_t*) geom->GetCurrentPoint();
					if (kVerbose) cout << " Position  " << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
					if (kVerbose) cout << "---- Tracking to ------->  " << geom->GetCurrentVolume()->GetName() << endl;

					////// Look for next surface
					geom->FindNextBoundary();
					Double_t *dir = (Double_t*) geom->GetCurrentDirection();

					///// Step over to next surface
					geom->Step();

					///// Get new  position
					Double_t *pos = (Double_t*) geom->GetCurrentPoint();
					if (kVerbose) cout << " Direction " << dir[0] << " " << dir[1] << " " << dir[2] << endl;
					if (kVerbose) cout << " Position  " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
					if (kVerbose) cout << "---- Tracking to ------->  " << geom->GetCurrentVolume()->GetName() << endl;

					// If photon is hitting something we should increment the counter and fill the appropriate end histogram...
					// If it hits the top, abs or support volumes we need to stop the photon...

					//std::cout << geom->GetCurrentVolume()->GetName() << std::endl;

					// 1) If the photon hits the "top" volume we can stop it...
					if ((geom->GetCurrentVolume()->GetName()) == top->GetName()) {
						if (reflections == 0) {
							nTop += 1;
							hnTop->Fill(nDetHistoX);
						}
						Stop = true;
					} // if top

					// 2) If the photon hits the "support" volume we can stop it...
					if ((geom->GetCurrentVolume()->GetName()) == supvol->GetName()) {
						if (reflections == 0) {
							nSup += 1;
							hnSup->Fill(nDetHistoX);
						}
						Stop = true;
					} // if sup

					// 3) If the photon hits the "absvol" volume we can stop it...
					if ((geom->GetCurrentVolume()->GetName()) == absvol->GetName()) {
						if (reflections == 0) {
							nAbs += 1;
							hnAbs->Fill(nDetHistoX);
						}
						Stop = true;
					} // if abs

					// 4) If the photon hits any of the PMT+LC record it's first hit...
					if ((geom->GetCurrentVolume()->GetName()) != top->GetName() && reflections == 0) {
						nHit += 1;
						hnhit->Fill(nDetHistoX);
					}

					// 5) If the photon hits the reflective surface it can be reflected...
					if ((geom->GetCurrentVolume()->GetName()) == refvol->GetName()) {
						// emulate non perfect reflectivity
						if (tr->Uniform() > LC_REFLECTIVITY) {
							Stop = true;
							continue;
						}

						if (reflections == 0) {
							nRef += 1;
							hnRef->Fill(nDetHistoX);
						}

						reflections += 1;

						Double_t norm[3], newdir[3];
						refpcon->ComputeNormal(pos, dir, norm);
						for (int i = 0; i < 3; i++)
							norm[i] *= -1;

						////// ADD smearing to norm (emulate non perfectly smooth surface)
						//////  Reasonable smearing need to be quantitatively evaluated or derived from literature !!!!!

						Double_t smear_rms = LC_ROUGHNESS;

						// Add Gaussian smearing to norm components
						norm[0] += tr->Gaus(0, smear_rms);
						norm[1] += tr->Gaus(0, smear_rms);
						norm[2] += tr->Gaus(0, smear_rms);

						// Normalise new vector back to 1
						Double_t nnmod = sqrt(pow(norm[0], 2) + pow(norm[1], 2) + pow(norm[2], 2));
						norm[0] = norm[0] / nnmod;
						norm[1] = norm[1] / nnmod;
						norm[2] = norm[2] / nnmod;

						////// Perform Reflection computing new photon direction
						if(kVerbose)  cout <<  " Normal    " << norm[0] << " " << norm[1] << " " << norm[2] << endl;

						///// Direction and normal scalar product (dir * norm)
						Double_t d_dot_n = dir[0] * norm[0] + dir[1] * norm[1] + dir[2] * norm[2];

						///// Compute reflected direction components
						for (int i = 0; i < 3; i++)
							newdir[i] = dir[i] - 2 * d_dot_n * norm[i];

						if (kVerbose) cout << " Reflection " << newdir[0] << " " << newdir[1] << " " << newdir[2] << endl;

						////// Assign new direction to photon
						geom->SetCurrentPoint(pos);
						geom->SetCurrentDirection(newdir);
						geom->FindNextBoundary();
						geom->Step();

						refl_points.push_back(TVector3(pos[0], pos[1], pos[2]));
					} // if reflector

					// 6) If the photon reaches the PMT we can see if it detected...
					if ((geom->GetCurrentVolume()->GetName()) == pmt->GetName()) {
						if (reflections == 0) {
							nPmt += 1;
							hnPmt->Fill(nDetHistoX);
						}

						double CathAngle = TMath::ACos(
								fabs(pos[2])
										/ sqrt(
												TMath::Power(pos[0], 2) + TMath::Power(pos[1], 2)
														+ TMath::Power(pos[2], 2))) * 180 / TMath::Pi();
						double CollEfficiency = 1;

						for (int iEntry = 0; iEntry < 9; ++iEntry) {
							if (CathAngle >= collection_angle[iEntry] && CathAngle < collection_angle[iEntry + 1]) {
								CollEfficiency = collection_eff[iEntry]
										+ (CathAngle - collection_angle[iEntry])
												/ (collection_angle[iEntry + 1] - collection_angle[iEntry])
												* (collection_eff[iEntry + 1] - collection_eff[iEntry]);
							}
						}

						if(kVerbose)  cout <<  " CathAngle -> " << CathAngle << ", CollEfficiency -> " << CollEfficiency << endl;

						if (CollEfficiency > tr->Uniform()) {
							if (kVerbose) cout << "Photon Detected! " << nDet << endl;
							detected = true;

							hndet->Fill(nDetHistoX);
							nDet += 1;

							// Fill 2D position histograms...
							hDet_xypos[ipos]->Fill(pos[0], pos[1]);

							if (reflections == 0) {
								hDet_first_xypos[ipos]->Fill(pos[0], pos[1]);
							} else if (reflections == 1) {
								hDet_first_xypos[ipos]->Fill(refl_points[0].X(), refl_points[0].Y());
							} else {
								for (vit = refl_points.begin(); vit != refl_points.end(); ++vit) {
									hDet_multiple_xypos[ipos]->Fill(vit->X(), vit->Y());
								}
							}
						}
						Stop = true;
					} // if pmt

					// 7) If the photon has been stopped but not detected add to non detected plots...
					if (Stop && !detected && (geom->GetCurrentVolume()->GetName()) != top->GetName()) {
						hNon_xypos[ipos]->Fill(pos[0], pos[1]);
						if (reflections > 0) {
							for (vit = refl_points.begin(); vit != refl_points.end(); ++vit) {
								hNon_reflect_xypos[ipos]->Fill(vit->X(), vit->Y());
							}
						}
					}
				} // while !Stop
				//std::cout << iph << std::endl;
			} // for iph
			std::cout << "nDet = " << nDet << ", nHit = " << nHit << ", nTop = " << nTop << ", nPmt = " << nPmt << ", nSup = " << nSup << ", nRef = " << nRef << ", nAbs = " << nAbs << std::endl;

			if (makeGifs) {
				hDet_xypos[ipos]->Draw();
				c1->Print("gifs/hDet_xypos.gif+50");  // print canvas to GIF file with 500ms delays
				hDet_first_xypos[ipos]->Draw();
				c1->Print("gifs/hDet_first_xypos.gif+50");  // print canvas to GIF file with 500ms delays
				hDet_multiple_xypos[ipos]->Draw();
				c1->Print("gifs/hDet_multiple_xypos.gif+50");  // print canvas to GIF file with 500ms delays
				hNon_xypos[ipos]->Draw();
				c1->Print("gifs/hNon_xypos.gif+50");  // print canvas to GIF file with 500ms delays
				hNon_reflect_xypos[ipos]->Draw();
				c1->Print("gifs/hNon_reflect_xypos.gif+50");  // print canvas to GIF file with 500ms delays
			}

			if (ipos == 0) {
				zeroNum = double(nDet);
				scale = double(nDet) / double(nHit);
				std::cout << "Scale = " << scale << std::endl;
			}

		} // for pos

		TH1F *heff = (TH1F*) hndet->Clone();
		heff->SetName("heff");
		heff->SetTitle("heff");
		heff->GetYaxis()->SetTitle("Efficiency");
		heff->Scale(1/zeroNum);
		heff->Scale(scale);

		std::cout << std::endl << "******Efficiency histogram printout******" << std::endl;
		std::cout << "*****************************************" << std::endl;
		heff->Print("all");
		std::cout << "*****************************************" << std::endl;

		////// Write Histograms to output file and close it
		fout->Write(0);
		fout->Close();

	} // if runOpt > 0

} // rootGeo

///////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_lcfile(string fname, Double_t *lcz, Double_t *lcr) {

	ifstream lcin;
	lcin.open(fname.c_str());

	int index;
	Double_t z, r;

	string line;
	while (!lcin.eof()) {

		getline(lcin, line);
		if (line[0] != '#') {
			std::stringstream ss(line);
			ss >> index >> z >> r;
			//      std::cout << index << " " << z << " " << r << std::endl;
			lcz[index] = z / 1000.;
			lcr[index] = r / 1000.;

		} // if line
	} // while

} // read_lcfile

