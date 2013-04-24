{
  gStyle->SetOptStat("0"); 
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  gSystem->Load("../WCSim/libWCSimRoot.so");
  gSystem->Load("lib/libWCSimAnalysis.so");


  // Load Data
  // =========
//  WCSimInterface::LoadData("/r01/lbne/water/wcsim_root_files/DUSEL_200kton_12inch_10perCent/muon_minus_001600MeV_200kton.0001.root");
//  WCSimInterface::LoadData("../nuFlux/numu/sim/numi_numu_sim_1.root");
//  WCSimInterface::LoadData("../mcSamples/numi_numu_1_100kT_10pC.root");
  WCSimInterface::LoadData("../mcSamples/piZero.root");

  // Initialize Reconstruction
  // =========================
  WCSimReco* myReco = WCSimRecoFactory::Instance()->MakeReco();

  // Reconstruct Event
  // =================

 // get first entry	 (nb. first entry is 0, but it's got no ring in my data so I've set this to the second for testing - AJP)
  WCSimInterface::LoadEvent(3);

  WCSimRecoEvent* myRecoEvent = WCSimInterface::RecoEvent();
  WCSimTrueEvent* myTrueEvent = WCSimInterface::TrueEvent();


  // apply filter
  myReco->RunFilter(myRecoEvent);


  // True Vertex
  // ===========
  WCSimRecoVertex* myTrueVertex = new WCSimRecoVertex();
  myTrueVertex->SetVertex(myTrueEvent->GetVtxX(),
                          myTrueEvent->GetVtxY(),
                          myTrueEvent->GetVtxZ());

  // Vertex Finder
  // =============
  WCSimVertexFinder* myVertexFinder = WCSimVertexFinder::Instance();
  WCSimRecoVertex* myRecoVertex = myVertexFinder->Run(myRecoEvent);




  // Ring Finder
  // ===========
  WCSimRingFinder* myRingFinder = WCSimRingFinder::Instance();

  // Perform Hough Transform
  // ========================
  WCSimHoughTransformArray* myTransformArray = myRingFinder->HoughTransformArray(myRecoEvent,myRecoVertex);  

  // Find Hough Peak
  // ===============    
  Double_t phi = 0.0;
  Double_t costheta = 0.0;
  Double_t hx = 0.0;
  Double_t hy = 0.0;
  Double_t hz = 0.0;
  Double_t angle = 0.0;
  Double_t height = 0.0;
  Int_t bin = 0;
  
	// MAKE TH3D OF ALL (angle,phi,cos(theta)) 
  TH3D * h3d = new TH3D("h3d","h3d",30,24.5,54.5,120,-180,180,60,-1,1);
  for( int coneAngleBin = 0; coneAngleBin < myTransformArray->GetBins(); coneAngleBin++ ){
		WCSimHoughTransform* myHoughTransform = myTransformArray->GetHoughTransform(coneAngleBin);
		TH2D * h2d = myHoughTransform->GetTH2D("h2d");
		for( int phiBin=0; phiBin < h2d->GetNbinsX(); phiBin++ ){
			for(int thetaBin = 0; thetaBin < h2d->GetNbinsY(); thetaBin++){
				h3d->SetBinContent(coneAngleBin,phiBin,thetaBin,h2d->GetBinContent(phiBin,thetaBin));
			}
		}
		delete h2d;
  }
  TFile * hough3d = new TFile("hough3d.root","RECREATE");
	h3d->Write();

	WCSimHoughTransform* myHoughTransform = myTransformArray->GetHoughTransform( 0 );


	// MAKE A LIST OF THE MAXIMA OF EACH PHI, THETA BIN, INDEXED BY CONE ANGLE
  vector<vector<double>> pointsMax;
  
  for( int phiBin = 0; phiBin < myHoughTransform->GetPhiBins(); phiBin++){
		for( int thetaBin = 0; thetaBin < myHoughTransform->GetThetaBins(); thetaBin++){
			vector<double> point;
			point.push_back(0);
      double max = 0;
			for( int coneBin = 0; coneBin < myTransformArray->GetBins(); coneBin++){
				myHoughTransform = myTransformArray->GetHoughTransform( coneBin );
				if(myHoughTransform->GetEntry(phiBin, thetaBin) > max){
					max = myHoughTransform->GetEntry(phiBin, thetaBin);
					point[0] = coneBin;
				}
			}
     	point.push_back(phiBin);
			point.push_back(thetaBin);
			point.push_back(max);
			pointsMax.push_back(point);
			point.clear();
		}
	}	

  TH3D * hMax3d = new TH3D("hMax3d","hMax3d",120,-180,180,60,-1,1,30,24.5,54.5);
  TH2D * hMax2d = new TH2D("hMax2d","hMax2d",120,-180,180,60,-1,1);
	for(int x = 0; x < pointsMax.size(); x++){
			hMax3d->SetBinContent(pointsMax[x][1],pointsMax[x][2],pointsMax[x][0],pointsMax[x][3]);
			hMax2d->SetBinContent(pointsMax[x][1],pointsMax[x][2],pointsMax[x][3]);
	}

	hMax2d->Write();
	hMax3d->Write();
	hough3d->Close();


  myTransformArray->FindPeak(hx,hy,hz,angle,height);

  Double_t deltaphi = 180.0;
  Int_t bin = myTransformArray->FindBin(42);
  
	// Added this "if" (and another one when drawing histTransform) as when bin == 0
	// this was crashing - AJP 15/1/13
  bool binExists = (bin>0);
  if( binExists ) {
		WCSimHoughTransform* myTransform = myTransformArray->GetHoughTransform(bin);
  	TH2D* histTransform = myTransform->GetRotatedTH2D("hello",deltaphi);
	}

  // Plot Truth
  // ==========
  Double_t px = myTrueEvent->GetDirX();
  Double_t py = myTrueEvent->GetDirY();
  Double_t pz = myTrueEvent->GetDirZ();

  if( px!=0.0 ){
    phi = atan(py/px);
  }
  if( px<=0.0 ){
    if( py>0.0 ) phi += TMath::Pi();
    if( py<0.0 ) phi -= TMath::Pi();
  }

  phi *= 180.0/TMath::Pi();
  phi += deltaphi; if( phi>180.0 ) phi -= 360.0;
  costheta = pz;

  marker = new TMarker(phi,costheta,29);
  marker->SetMarkerColor(2);
  marker->SetMarkerSize(2.5);

  // Plot Reco
  // =========
  if( hx!=0.0 ){
    phi = atan(hy/hx);
  }
  if( hx<=0.0 ){
    if( hy>0.0 ) phi += TMath::Pi();
    if( hy<0.0 ) phi -= TMath::Pi();
  }

  phi *= 180.0/TMath::Pi();
  phi += deltaphi; if( phi>180.0 ) phi -= 360.0;
  costheta = hz;

  markerpeak = new TMarker(phi,costheta,29);
  markerpeak->SetMarkerColor(1);
  markerpeak->SetMarkerSize(2.5);

  // Display Hough Transform
  // =======================
  TCanvas* fHoughCanvas = new TCanvas("wcsim_houghtransform_canvas","WCSim HoughTransform Viewer");
  if(binExists) histTransform->Draw("col");
  
   marker->Draw();      // MC Truth
   markerpeak->Draw();  // Hough Peak

}
