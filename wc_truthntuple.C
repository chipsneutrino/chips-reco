#include <algorithm>
#include <cstdio>
#include <iostream>
#include <vector>

void wc_truthntuple(std::string fPath = "/unix/fnu/ajperch/sampleEvents/WCSimOutput/tauTest.root", std::string oPath = "../recoNtp/nutau_test_TRUTH.root", int nEvt = 2000){
// ../mcSamples/numi_numu_NC_100kT_10pC_1000events_PEperMeV.root"

	// Load the WC Library
	gSystem->Load("~/MINOS/waterCherenkov/WCSim/libWCSimRoot.so");

	TChain *d = new TChain("wcsimT");
	d->Add(fPath.c_str());

	WCSimRootEvent *evt = new WCSimRootEvent();
	d->SetBranchAddress("wcsimrootevent",&evt);

	TChain *g = new TChain("wcsimGeoT");
	g->Add(fPath.c_str());
	WCSimRootGeom *geo = new WCSimRootGeom();
	g->SetBranchAddress("wcsimrootgeom",&geo);
	g->GetEntry(0);

	int numEvt = d->GetEntries();

	std::cout << "= Loaded file " << fPath << " with " << numEvt << " entries." << std::endl;

	if(nEvt!=0 && nEvt < numEvt) numEvt = nEvt;

	// Create the output file and output trees.
	TFile output(oPath.c_str(),"RECREATE");
	TTree *truthTree = new TTree("truthNtp","truthNtp"); // Event level tree
	truthTree->SetAutoFlush(0);	

	

	// Variables that we want.
	int eventID;
	int mode;	// interaction mode
	int neutrinoPDG;
	double neutrinoEnergy;
	double trueVtxX, trueVtxY, trueVtxZ;
	double trueDirX, trueDirY, trueDirZ;
	int nDaughters;
	int trueDaughterPDG[2500] = {0};
	int trueDaughterID[2500] = {0};
	double trueDaughterE[2500] = {0};
	double trueDaughterVtxX[2500] = {0.0};
	double trueDaughterVtxY[2500] = {0.0};
	double trueDaughterVtxZ[2500] = {0.0};
	double trueDaughterDirX[2500] = {0.0};
	double trueDaughterDirY[2500] = {0.0};
	double trueDaughterDirZ[2500] = {0.0};
	int trueDaughterParentPDG[2500] = {0};

	int nHits;
	int nHitsUp; // NHits in -z
	int nHitsDown; // NHits in +z
	double nHitsUDRatio; // Ratio of Upstream to Downstream hits
	double chargeTotal; // Total charge
	double chargeUp; // Charge in -z
	double chargeDown; // Charge in +z
	double chargeUDRatio; // Ratio of Upstream to Downstream charge
	double chargePerPmtUp;
	double chargePerPmtDown;
	double chargePerPmtUDRatio;

	// Temporary variables that we use but don't keep
	int hitTubeID[50000];
	double hitQ[50000];
	double hitX[50000];
	double hitY[50000];
	double hitZ[50000];

	// Event level arrays
	double hitU, hitV;


	// Set all the tree branches
	truthTree->Branch("eventID",&eventID,"eventID/I");
	truthTree->Branch("mode",&mode,"mode/I");
	truthTree->Branch("neutrinoEnergy",&neutrinoEnergy,"neutrinoEnergy/D");
	truthTree->Branch("neutrinoPDG",&neutrinoPDG,"neutrinoPDG/I");
	truthTree->Branch("trueVtxX",&trueVtxX,"trueVtxX/D");
	truthTree->Branch("trueVtxY",&trueVtxY,"trueVtxY/D");
	truthTree->Branch("trueVtxZ",&trueVtxZ,"trueVtxZ/D");
	truthTree->Branch("trueDirX",&trueDirX,"trueDirX/D");
	truthTree->Branch("trueDirY",&trueDirY,"trueDirY/D");
	truthTree->Branch("trueDirZ",&trueDirZ,"trueDirZ/D");
  
	truthTree->Branch("nDaughters",&nDaughters,"nDaughters/I");
	truthTree->Branch("trueDaughterPDG",trueDaughterPDG,"trueDaughterPDG[250]/I");
	truthTree->Branch("trueDaughterID",trueDaughterID,"trueDaughterID[250]/I");
	truthTree->Branch("trueDaughterE",trueDaughterE,"trueDaughterE[250]/D");
	truthTree->Branch("trueDaughterVtxX",trueDaughterVtxX,"trueDaughterVtxX[250]/D");
	truthTree->Branch("trueDaughterVtxY",trueDaughterVtxY,"trueDaughterVtxY[250]/D");
	truthTree->Branch("trueDaughterVtxZ",trueDaughterVtxZ,"trueDaughterVtxZ[250]/D");
	truthTree->Branch("trueDaughterDirX",trueDaughterDirX,"trueDaughterDirX[250]/D");
	truthTree->Branch("trueDaughterDirY",trueDaughterDirY,"trueDaughterDirY[250]/D");
	truthTree->Branch("trueDaughterDirZ",trueDaughterDirZ,"trueDaughterDirZ[250]/D");
	truthTree->Branch("trueDaughterParentPDG",trueDaughterParentPDG,"trueDaughterParentPDG[250]/I");

	truthTree->Branch("nHits",&nHits,"nHits/I");
	truthTree->Branch("nHitsUpstream",&nHitsUp,"nHitsUp/I");
	truthTree->Branch("nHitsDownstream",&nHitsDown,"nHitsDown/I");
	truthTree->Branch("nHitsUDRatio",&nHitsUDRatio,"nHitsUDRatio/D");
	truthTree->Branch("chargeTotal",&chargeTotal,"chargeTotal/D");
	truthTree->Branch("chargeUpstream",&chargeUp,"chargeUp/D");
	truthTree->Branch("chargeDownstream",&chargeDown,"chargeDown/D");
  truthTree->Branch("chargeUDRatio",&chargeUDRatio,"chargeUDRatio/D");
	truthTree->Branch("chargePerPmtUp",&chargePerPmtUp,"chargePerPmtUp/D");
	truthTree->Branch("chargePerPmtDown",&chargePerPmtDown,"chargePerPmtDown/D");
	truthTree->Branch("chargePerPmtUDRatio",&chargePerPmtUDRatio,"chargePerPmtUDRatio/D");
	



	for(int i = 0; i < numEvt; ++i){

	
		d->GetEntry(i);
	
		nHitsUp = 0; nHitsDown = 0;
		chargeUp = 0; chargeDown = 0;
		chargeTotal = 0;



	// Load the trigger. This contains the event information.  Trigger 0 is the initial event
		WCSimRootTrigger *trig = evt->GetTrigger(0);


		// Track 0 is the incoming neutrino - use this for its energy and the primary vertex
		WCSimRootTrack * track = (WCSimRootTrack *)(trig->GetTracks()->At(0));
	

		// Get the number of hits
		nHits = trig->GetNcherenkovdigihits();
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\\	
		for(int c = 0; c < nHits; ++c){
			// Get the hit
			WCSimRootCherenkovDigiHit *hit = ( WCSimRootCherenkovDigiHit *)(trig->GetCherenkovDigiHits()->At(c));
			WCSimRootPMT pmt = geo->GetPMT(hit->GetTubeId());
			output.cd();
			hitTubeID[c] = hit->GetTubeId();
			hitQ[c] = hit->GetQ();
			hitX[c] = pmt.GetPosition(0);
			hitY[c] = pmt.GetPosition(1);
			hitZ[c] = pmt.GetPosition(2);
//			if(i == 9 && c < 10) std::cout << hitQ[c] << "   " << hitX[c] << "   " << hitY[c] << "   " << hitZ[c] << "   " << hitTubeID[c] << std::endl;

			CalcFlatCoords(hitX[c],hitY[c],hitZ[c],hitU,hitV);
//			if(hitQ < 0) continue;
//			if(hitQ > 1e8) continue;
			// Get the hit position from the PMT
	
			// Get rid of the noisy PMT... quite why this is necessary I have no idea.
			if(hitU>-2421&&hitU<-2410&&hitV>-6922&&hitV<-6921){
			 hitQ[c] = 0; hitX[c] = 0; hitY[c] = 0; hitZ[c] = 0; hitTubeID[c] = 0;
			 continue;
			}
		
		// Fill up/down stream variables
			if(hitZ[c] > trig->GetVtx(2)){
				++nHitsDown;
				chargeDown += hitQ[c];
			}
			else{
				++nHitsUp;
				chargeUp += hitQ[c];
			}
			chargeTotal += hitQ[c];
		

		}
		for(int c = nHits; c < 50000; c++){
			hitX[c] = 0;
			hitY[c] = 0;
			hitZ[c] = 0;
			hitQ[c] = 0;
			hitTubeID[c] = 0;
		}
		
			
	
		nHitsUDRatio = nHitsUp / (double)nHitsDown;
		chargeUDRatio = chargeUp / chargeDown;
		chargePerPmtUp = chargeUp/(double)nHitsUp;
		chargePerPmtDown = chargeDown/(double)nHitsDown;
		chargePerPmtUDRatio = chargePerPmtUp/chargePerPmtDown;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





		eventID = i;
		mode = trig->GetMode();
		neutrinoPDG = track->GetIpnu();
		neutrinoEnergy = track->GetE();
		trueVtxX = trig->GetVtx(0);
		trueVtxY = trig->GetVtx(1);
		trueVtxZ = trig->GetVtx(2);
		trueDirX = track->GetDir(0);
		trueDirY = track->GetDir(1);
		trueDirZ = track->GetDir(2);

		int nDaughters = trig->GetNtrack();
		std::cout << std::endl << std::endl << "Event " << i << "   " << nDaughters << " = nDaughters" << std::endl;
		for( int j = 0; j < nDaughters; j++ )
		{
//			std::cout << j << "   " << trig ;
			// trig = (WCSimRootTrigger *) evt->GetTrigger(j+1); // Trigger 0 is the neutrino
//			std::cout << "     "  << trig << "   done!   " << std::endl;
			track = (WCSimRootTrack *) (trig->GetTracks()->At(j));
			if(track->GetIpnu()==13 || track->GetIpnu() == 15)
            {
                if(track->GetIpnu()==13) std::cout << "It's a muon" << std::endl;
                else if(track->GetIpnu()==15) std::cout << "It's a tau" << std::endl;
                std::cout << track->GetFlag() << " = flag" << std::endl;
			    std::cout << track->GetParenttype() << " = parenttype" << std::endl;
			    std::cout << track->GetId() << " = ID" << std::endl;
			    std::cout << track->GetIpnu() << " = Ipnu" << std::endl;
			    std::cout << track->GetE() << " = E" << std::endl;
			    printf("Starts at (%f,%f,%f) and stops at (%f,%f,%f)\n\n", track->GetStart(0), track->GetStart(1), track->GetStart(2), track->GetStop(0), track->GetStop(1), track->GetStop(2));
            }
            trueDaughterPDG[j] = track->GetIpnu();
			trueDaughterParentPDG[j] = track->GetParenttype();
			trueDaughterID[j] = track->GetId();
			trueDaughterE[j] = track->GetE();
			trueDaughterVtxX[j] = trig->GetVtx(0);
			trueDaughterVtxY[j] = trig->GetVtx(1);
			trueDaughterVtxZ[j] = trig->GetVtx(2);
			trueDaughterDirX[j] = track->GetDir(0);
			trueDaughterDirY[j] = track->GetDir(1);
			trueDaughterDirZ[j] = track->GetDir(2);
	
		}
		truthTree->Fill();
	}
		
	truthTree->Write();
	output.Close();
	std::cout << "Wrote tree to " << oPath << std::endl;
	std::cout << "Finished successfully." << std::endl;
}




// Function to calculate the flattened coordinates U and V
void CalcFlatCoords(double fEndX, double fEndY, double fEndZ, double &fEndU, double &fEndV){
	double fHalfWidthX = 4533.3 / 2.0;
	double fHalfWidthY = 4000.0 / 2.0;
	double fHalfWidthZ = 5500.0 / 2.0;

	// Downstream End
	if(fEndZ >= fHalfWidthZ){
		fEndU = fEndX;
		fEndV = fEndY + fHalfWidthZ + fHalfWidthY;
	}
	// Upstream End
	else if(fEndZ <= -fHalfWidthZ){
		fEndU = fEndX;
		fEndV = fEndY - fHalfWidthZ - fHalfWidthY;
	}
	// Top (split this in half)
	else if(fEndY >= fHalfWidthY){
		if(fEndX>0) fEndU = fHalfWidthX + 2*fHalfWidthY + fEndX;
		else fEndU =  -fHalfWidthX - 2*fHalfWidthY + fEndX;
		fEndV = fEndZ;
	}
	// Bottom
	else if(fEndY <= -fHalfWidthY){
		fEndU = fEndX;
		fEndV = fEndZ;
	}
	// Left Side	
	else if(fEndX <= -fHalfWidthX){
		fEndU = -fHalfWidthX - fHalfWidthY + fEndY;
		fEndV = fEndZ;
	}
	// Right Side
	else if(fEndX >= fHalfWidthX){
		fEndU = +fHalfWidthX + fHalfWidthY + fEndY;
		fEndV = fEndZ;
	}
}


