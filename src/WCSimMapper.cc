/*
 * WCSimMapper.cc
 *
 *  Created on: 26 Nov 2019
 *      Author: Josh
 */

// 1) Need to look at Hough maps, which ones shall we use
// 2) Do we want to slice events?
// 3) Do we want to change the filtering, do we want to use it
// 4) Do I just want to save the info for the highest ring?

#include "WCSimMapper.hh"
#include "WCSimInterface.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoSeed.hh"
#include "WCSimRecoRing.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimHoughTransform.hh"
#include "WCSimHoughTransformArray.hh"

ClassImp (WCSimMapper)

WCSimMapper::WCSimMapper(const char* in_file, const char* out_file, int max_events, int category, int pdg_code)
    :input_file_(in_file)
    ,max_events_(max_events)
	,pdg_code_(pdg_code)
	,true_category_(category)
{
	resetVariables(); // Reset all TTree variables

	// Open the output file with a TTree
	map_f_ = new TFile(out_file, "RECREATE", "HitMapFile");

	// Truth info tree
	true_t_ = new TTree("true","true");
	true_t_->Branch("true_category",&true_category_,"true_category_/I");
	true_t_->Branch("true_vtx_x",&true_vtx_x_,"true_vtx_x_/F");
	true_t_->Branch("true_vtx_y",&true_vtx_y_,"true_vtx_y_/F");
	true_t_->Branch("true_vtx_z",&true_vtx_z_,"true_vtx_z_/F");
	true_t_->Branch("true_dir_costheta",&true_dir_costheta_,"true_dir_costheta_/F");
	true_t_->Branch("true_dir_phi",&true_dir_phi_,"true_dir_phi_/F");
	true_t_->Branch("true_nu_energy",&true_nu_energy_,"true_nu_energy_/F");
	true_t_->Branch("true_lep_energy",&true_lep_energy_,"true_lep_energy_/F");

	// Event and reco tree
	reco_t_ = new TTree("reco","reco");
	reco_t_->Branch("raw_num_hits",&raw_num_hits_,"raw_num_hits_/I");
	reco_t_->Branch("raw_total_digi_q",&raw_total_digi_q_,"raw_total_digi_q_/F");
	reco_t_->Branch("filtered_num_hits",&filtered_num_hits_,"filtered_num_hits_/I");
	reco_t_->Branch("filtered_total_digi_q",&filtered_total_digi_q_,"filtered_total_digi_q_/F");
	reco_t_->Branch("num_hough_rings",&num_hough_rings_,"num_hough_rings_/I");
	reco_t_->Branch("first_ring_height",&first_ring_height_,"first_ring_height_/F");
	reco_t_->Branch("last_ring_height",&last_ring_height_,"last_ring_height_/F");
	reco_t_->Branch("vtx_x",&vtx_x_,"vtx_x_/F");
	reco_t_->Branch("vtx_y",&vtx_y_,"vtx_y_/F");
	reco_t_->Branch("vtx_z",&vtx_z_,"vtx_z_/F");
	reco_t_->Branch("vtx_t",&vtx_t_,"vtx_t_/F");
	reco_t_->Branch("dir_theta",&dir_costheta_,"dir_theta_/F");
	reco_t_->Branch("dir_phi",&dir_phi_,"dir_phi_/F");
	reco_t_->Branch("raw_hit_map_origin",&raw_hit_map_origin_,"raw_hit_map_origin_[64][64]/F");
	reco_t_->Branch("raw_charge_map_origin",&raw_charge_map_origin_,"raw_charge_map_origin_[64][64]/F");
	reco_t_->Branch("raw_time_map_origin",&raw_time_map_origin_,"raw_time_map_origin_[64][64]/F");
	reco_t_->Branch("filtered_hit_map_origin",&filtered_hit_map_origin_,"filtered_hit_map_origin_[64][64]/F");
	reco_t_->Branch("filtered_charge_map_origin",&filtered_charge_map_origin_,"filtered_charge_map_origin_[64][64]/F");
	reco_t_->Branch("filtered_time_map_origin",&filtered_time_map_origin_,"filtered_time_map_origin_[64][64]/F");
	reco_t_->Branch("raw_hit_map_vtx",&raw_hit_map_vtx_,"raw_hit_map_vtx_[64][64]/F");
	reco_t_->Branch("raw_charge_map_vtx",&raw_charge_map_vtx_,"raw_charge_map_vtx_[64][64]/F");
	reco_t_->Branch("raw_time_map_vtx",&raw_time_map_vtx_,"raw_time_map_vtx_[64][64]/F");
	reco_t_->Branch("filtered_hit_map_vtx",&filtered_hit_map_vtx_,"filtered_hit_map_vtx_[64][64]/F");
	reco_t_->Branch("filtered_charge_map_vtx",&filtered_charge_map_vtx_,"filtered_charge_map_vtx_[64][64]/F");
	reco_t_->Branch("filtered_time_map_vtx",&filtered_time_map_vtx_,"filtered_time_map_vtx_[64][64]/F");
	//reco_t_->Branch("raw_hit_hough_map_vtx",&raw_hit_hough_map_vtx_,"raw_hit_hough_map_vtx_[64][64]/F");
	reco_t_->Branch("filtered_hit_hough_map_vtx",&filtered_hit_hough_map_vtx_,"filtered_hit_hough_map_vtx_[64][64]/F");
	//reco_t_->Branch("raw_charge_hough_map_vtx",&raw_charge_hough_map_vtx_,"raw_charge_hough_map_vtx_[64][64]/F");
	//reco_t_->Branch("filtered_charge_hough_map_vtx",&filtered_charge_hough_map_vtx_,"filtered_charge_hough_map_vtx_[64][64]/F");
}

void WCSimMapper::run()
{
	WCSimInterface::LoadData(input_file_); // Load the input file

	WCSimRecoSeed* reco = new WCSimRecoSeed(); // Create the reconstruction seeder
	reco->SetNumberOfTracks(1); 

	// Loop through all events and generate hit maps and fill output tree
	int ev = 0;
	int valid_events = 0;
	for(ev; ev<WCSimInterface::GetNumEvents(); ev++)
	{	
		if(valid_events==max_events_) break;

		WCSimInterface::Instance()->BuildEvent(ev);
		if(WCSimInterface::Instance()->CheckEvent())
		{
			reco_event_ = WCSimInterface::RecoEvent();
			reco->Run(reco_event_);
			if(reco_event_->GetNFilterDigits() == 0) continue;

			valid_events ++;

			fillTrueTree();
			fillRecoTree();

			resetVariables();
		}
	}

	// Write the output trees to file
	writeFile(); 
	map_f_->Close();

	std::cout << "Total = " << (ev-1) << ", Valid = " << valid_events << std::endl;

    delete reco;
}

void WCSimMapper::resetVariables()
{
	// True TTree variables (Used as labels in the CVN)
	true_vtx_x_ = -999;
	true_vtx_y_ = -999;
	true_vtx_z_ = -999;
	true_vtx_t_ = -999;
	true_dir_costheta_ = -999;
	true_dir_phi_ = -999;
	true_nu_energy_ = -999;
	true_lep_energy_ = -999;

	// Reco TTree variables (Used as examples in the CVN)
	raw_num_hits_ = -999;
	raw_total_digi_q_ = -999;
	filtered_num_hits_ = -999;
	filtered_total_digi_q_ = -999;
	num_hough_rings_ = -999;
	first_ring_height_ = -999;
	last_ring_height_ = -999;
	vtx_x_ = -999;
	vtx_y_ = -999;
	vtx_z_ = -999;
	vtx_t_ = -999;
	dir_costheta_ = -999;
	dir_phi_ = -999;

	memset(raw_hit_map_origin_, 0, sizeof(raw_hit_map_origin_));
	memset(raw_charge_map_origin_, 0, sizeof(raw_charge_map_origin_));
	memset(raw_time_map_origin_, 0, sizeof(raw_time_map_origin_));
	memset(filtered_hit_map_origin_, 0, sizeof(filtered_hit_map_origin_));
	memset(filtered_charge_map_origin_, 0, sizeof(filtered_charge_map_origin_));
	memset(filtered_time_map_origin_, 0, sizeof(filtered_time_map_origin_));
	memset(raw_hit_map_vtx_, 0, sizeof(raw_hit_map_vtx_));
	memset(raw_charge_map_vtx_, 0, sizeof(raw_charge_map_vtx_));
	memset(raw_time_map_vtx_, 0, sizeof(raw_time_map_vtx_));
	memset(filtered_hit_map_vtx_, 0, sizeof(filtered_hit_map_vtx_));
	memset(filtered_charge_map_vtx_, 0, sizeof(filtered_charge_map_vtx_));
	memset(filtered_time_map_vtx_, 0, sizeof(filtered_time_map_vtx_));
	//memset(raw_hit_hough_map_vtx_, 0, sizeof(raw_hit_hough_map_vtx_));
	memset(filtered_hit_hough_map_vtx_, 0, sizeof(filtered_hit_hough_map_vtx_));
	//memset(raw_charge_hough_map_vtx_, 0, sizeof(raw_charge_hough_map_vtx_));
	//memset(filtered_charge_hough_map_vtx_, 0, sizeof(filtered_charge_hough_map_vtx_));
}

void WCSimMapper::fillTrueTree()
{
	WCSimTruthSummary truth_sum = WCSimInterface::Instance()->GetTruthSummary();
	true_vtx_x_ = truth_sum.GetVertexX()/10;
	true_vtx_y_ = truth_sum.GetVertexY()/10;
	true_vtx_z_ = truth_sum.GetVertexZ()/10;
	true_vtx_t_ = truth_sum.GetVertexT();
	true_nu_energy_ = truth_sum.GetBeamEnergy();
	
	// Find the leading lepton
	for (int p=0; p<(int)truth_sum.GetNPrimaries(); p++) {
		int code = truth_sum.GetPrimaryPDG(p);
		if (code == pdg_code_) {
			TVector3 dir = truth_sum.GetPrimaryDir(p);
			float dirX = dir.X();
			float dirY = dir.Y();
			float dirZ = dir.Z();

			true_dir_phi_ = TMath::ATan2(dirY, dirX) * (180.0/TMath::Pi());
			true_dir_costheta_ = dirZ;

			true_lep_energy_ = truth_sum.GetPrimaryEnergy(p);
		}
	}
	true_t_->Fill();
}

float WCSimMapper::GetCosTheta(float x, float y, float z)
{
	TVector3 vec(x, y, z);
	TVector3 norm = vec.Unit();
	return norm.Z();
}

void WCSimMapper::fillRecoTree()
{
	// First get the Hough transform space, we use the leading ring for all variables
	num_hough_rings_ = reco_event_->GetNRings();
	first_ring_height_ = reco_event_->GetRing(0)->GetHeight();
	last_ring_height_ = reco_event_->GetRing(num_hough_rings_-1)->GetHeight();

	vtx_x_ = reco_event_->GetRing(0)->GetVtxX();
	vtx_y_ = reco_event_->GetRing(0)->GetVtxY();
	vtx_z_ = reco_event_->GetRing(0)->GetVtxZ();
	vtx_t_ = reco_event_->GetVtxTime();

	float dirX = reco_event_->GetRing(0)->GetDirX();
	float dirY = reco_event_->GetRing(0)->GetDirY();
	dir_phi_ = TMath::ATan2(dirY, dirX) * (180.0/TMath::Pi());
	dir_costheta_ = reco_event_->GetRing(0)->GetDirZ();

	WCSimHoughTransformArray* hough_array = reco_event_->GetHoughArray();
	WCSimHoughTransform* hough = hough_array->GetHoughTransform(0);
	TH2D* filtered_hit_hough_map_vtx_h = hough->GetTH2D("h");

	// We can then generate the hit and time maps in both raw and "vertex view" forms
	TH2F* raw_hit_map_origin_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* raw_charge_map_origin_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* raw_time_map_origin_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* filtered_hit_map_origin_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* filtered_charge_map_origin_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* filtered_time_map_origin_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* raw_hit_map_vtx_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* raw_charge_map_vtx_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* raw_time_map_vtx_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* filtered_hit_map_vtx_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* filtered_charge_map_vtx_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);
	TH2F* filtered_time_map_vtx_h = new TH2F("h", "h", 64, -180, 180, 64, -1, 1);

	// Loop through the digi hits
	raw_num_hits_ = reco_event_->GetNDigits();
	raw_total_digi_q_ = 0;
	float first_hit_time = 1000000;
	for(int h=0; h<raw_num_hits_; h++)
	{
		// Get the hit position and deposited charge
		float hit_x = reco_event_->GetDigit(h)->GetX();
		float hit_y = reco_event_->GetDigit(h)->GetY();
		float hit_z = reco_event_->GetDigit(h)->GetZ();
		float hit_q = reco_event_->GetDigit(h)->GetRawQPEs();
		raw_total_digi_q_ += hit_q;

		// Is the hit time the earliest in the event?
		float digi_hit_time = reco_event_->GetDigit(h)->GetRawTime();
		if (digi_hit_time < first_hit_time) 
		{ 
			first_hit_time = digi_hit_time; 
		}

		// Calculate origin phi and theta and fill hists
		float hit_phi_origin = TMath::ATan2(hit_y, hit_x) * (180.0/TMath::Pi());
		float hit_theta_origin = GetCosTheta(hit_x, hit_y, hit_z);

		raw_hit_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, 1);
		raw_charge_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, hit_q);
		int bin_num = raw_time_map_origin_h->FindBin(hit_phi_origin, hit_theta_origin);
		float currentBinTime = raw_time_map_origin_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			raw_time_map_origin_h->SetBinContent(bin_num, digi_hit_time);
		}

		// Calculate the "vertex view" phi and theta and fill hists 
		float hit_phi_vtx = TMath::ATan2((hit_y-vtx_y_), (hit_x-vtx_x_)) * (180.0/TMath::Pi());
		float hit_theta_vtx = GetCosTheta((hit_x-vtx_x_), (hit_y-vtx_y_), (hit_z-vtx_z_));

		raw_hit_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, 1);
		raw_charge_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, hit_q);
		bin_num = raw_time_map_vtx_h->FindBin(hit_phi_vtx, hit_theta_vtx);
		currentBinTime = raw_time_map_vtx_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			raw_time_map_vtx_h->SetBinContent(bin_num, digi_hit_time);
		}
	} 

	// We now need to move all timestamps to be from zero
	for(int i=0; i<(raw_time_map_origin_h->GetNbinsX()+1); i++)
	{
		for(int j=0; j<(raw_time_map_origin_h->GetNbinsY()+1); j++)
		{
			float raw_content = raw_time_map_origin_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				raw_time_map_origin_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
			raw_content = raw_time_map_vtx_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				raw_time_map_vtx_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
		}	
	} 

	// Loop through filtered digi hits
	filtered_num_hits_ = reco_event_->GetNFilterDigits();
	filtered_total_digi_q_ = 0;
	first_hit_time = 1000000;
	for(int h=0; h<filtered_num_hits_; h++)
	{
		// Get the hit position and deposited charge
		float hit_x = reco_event_->GetFilterDigit(h)->GetX();
		float hit_y = reco_event_->GetFilterDigit(h)->GetY();
		float hit_z = reco_event_->GetFilterDigit(h)->GetZ();
		float hit_q = reco_event_->GetFilterDigit(h)->GetRawQPEs();
		filtered_total_digi_q_ += hit_q;

		// Is the hit time the earliest in the event?
		float digi_hit_time = reco_event_->GetFilterDigit(h)->GetRawTime();
		if (digi_hit_time < first_hit_time) 
		{ 
			first_hit_time = digi_hit_time; 
		}

		// Calculate origin phi and theta and fill hists
		float hit_phi_origin = TMath::ATan2(hit_y, hit_x) * (180.0/TMath::Pi());
		float hit_theta_origin = GetCosTheta(hit_x, hit_y, hit_z);

		filtered_hit_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, 1);
		filtered_charge_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, hit_q);
		int bin_num = filtered_time_map_origin_h->FindBin(hit_phi_origin, hit_theta_origin);
		float currentBinTime = filtered_time_map_origin_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			filtered_time_map_origin_h->SetBinContent(bin_num, digi_hit_time);
		}

		// Calculate the "vertex view" phi and theta and fill hists 
		float hit_phi_vtx = TMath::ATan2((hit_y-vtx_y_), (hit_x-vtx_x_)) * (180.0/TMath::Pi());
		float hit_theta_vtx = GetCosTheta((hit_x-vtx_x_), (hit_y-vtx_y_), (hit_z-vtx_z_));

		filtered_hit_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, 1);
		filtered_charge_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, hit_q);
		bin_num = filtered_time_map_vtx_h->FindBin(hit_phi_vtx, hit_theta_vtx);
		currentBinTime = filtered_time_map_vtx_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			filtered_time_map_vtx_h->SetBinContent(bin_num, digi_hit_time);
		}
	}

	// We now need to move all timestamps to be from zero
	for(int i=0; i<(filtered_time_map_origin_h->GetNbinsX()+1); i++)
	{
		for(int j=0; j<(filtered_time_map_origin_h->GetNbinsY()+1); j++)
		{
			float raw_content = filtered_time_map_origin_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				filtered_time_map_origin_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
			raw_content = filtered_time_map_vtx_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				filtered_time_map_vtx_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
		}	
	}

	// We now fill the arrays we save to the reco tree
	for(int i=0; i<64; i++)
	{
		for(int j=0; j<64; j++)
		{
			raw_hit_map_origin_[j][i] = raw_hit_map_origin_h->GetBinContent(i,j);
			raw_charge_map_origin_[j][i] = raw_charge_map_origin_h->GetBinContent(i,j);
			raw_time_map_origin_[j][i] = raw_time_map_origin_h->GetBinContent(i,j);
			filtered_hit_map_origin_[j][i] = filtered_hit_map_origin_h->GetBinContent(i,j);
			filtered_charge_map_origin_[j][i] = filtered_charge_map_origin_h->GetBinContent(i,j);
			filtered_time_map_origin_[j][i] = filtered_time_map_origin_h->GetBinContent(i,j);
			raw_hit_map_vtx_[j][i] = raw_hit_map_vtx_h->GetBinContent(i,j);
			raw_charge_map_vtx_[j][i] = raw_charge_map_vtx_h->GetBinContent(i,j);
			raw_time_map_vtx_[j][i] = raw_time_map_vtx_h->GetBinContent(i,j);
			filtered_hit_map_vtx_[j][i] = filtered_hit_map_vtx_h->GetBinContent(i,j);
			filtered_charge_map_vtx_[j][i] = filtered_charge_map_vtx_h->GetBinContent(i,j);
			filtered_time_map_vtx_[j][i] = filtered_time_map_vtx_h->GetBinContent(i,j);
			
			filtered_hit_hough_map_vtx_[j][i] = filtered_hit_hough_map_vtx_h->GetBinContent(i,j);
		}	
	}

	delete raw_hit_map_origin_h;
	delete raw_charge_map_origin_h;
	delete raw_time_map_origin_h;
	delete filtered_hit_map_origin_h;
	delete filtered_charge_map_origin_h;
	delete filtered_time_map_origin_h;
	delete raw_hit_map_vtx_h;
	delete raw_charge_map_vtx_h;
	delete raw_time_map_vtx_h;
	delete filtered_hit_map_vtx_h;
	delete filtered_charge_map_vtx_h;
	delete filtered_time_map_vtx_h;

	reco_t_->Fill();
}

void WCSimMapper::writeFile()
{
	true_t_->Write();
	reco_t_->Write();
}