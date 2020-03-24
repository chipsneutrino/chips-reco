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

WCSimMapper::WCSimMapper(const char* in_file, const char* out_file, int max_events, bool save_extra)
    :input_file_(in_file)
	,output_file_(out_file)
    ,max_events_(max_events)
	,save_extra_(save_extra)
{
	resetVariables(); // Reset all TTree variables

	// Truth info tree
	true_t_ = new TTree("true","true");
	true_t_->Branch("t_nu",&t_nu_,"t_nu_/I");
	true_t_->Branch("t_code",&t_code_,"t_code_/I");
	true_t_->Branch("t_vtxX",&t_vtxX_,"t_vtxX_/F");
	true_t_->Branch("t_vtxY",&t_vtxY_,"t_vtxY_/F");
	true_t_->Branch("t_vtxZ",&t_vtxZ_,"t_vtxZ_/F");
	true_t_->Branch("t_vtxT",&t_vtxT_,"t_vtxT_/F");
	true_t_->Branch("t_nuEnergy",&t_nuEnergy_,"t_nuEnergy_/F");
	true_t_->Branch("t_p_energies",&t_p_energies_,"t_p_energies_[10]/F");
	true_t_->Branch("t_p_pdgs",&t_p_pdgs_,"t_p_pdgs_[10]/I");
	true_t_->Branch("t_p_dirTheta",&t_p_dirTheta_,"t_p_dirTheta_[10]/F");
	true_t_->Branch("t_p_dirPhi",&t_p_dirPhi_,"t_p_dirPhi_[10]/F");

	// Event and reco tree
	reco_t_ = new TTree("reco","reco");
	reco_t_->Branch("r_raw_num_hits",&r_raw_num_hits_,"r_raw_num_hits_/I");
	reco_t_->Branch("r_raw_total_digi_q",&r_raw_total_digi_q_,"r_raw_total_digi_q_/F");
	reco_t_->Branch("r_filtered_num_hits",&r_filtered_num_hits_,"r_filtered_num_hits_/I");
	reco_t_->Branch("r_filtered_total_digi_q",&r_filtered_total_digi_q_,"r_filtered_total_digi_q_/F");
	reco_t_->Branch("r_num_hough_rings",&r_num_hough_rings_,"r_num_hough_rings_/I");
	reco_t_->Branch("r_first_ring_height",&r_first_ring_height_,"r_first_ring_height_/F");
	reco_t_->Branch("r_last_ring_height",&r_last_ring_height_,"r_last_ring_height_/F");
	reco_t_->Branch("r_vtxX",&r_vtxX_,"r_vtxX_/F");
	reco_t_->Branch("r_vtxY",&r_vtxY_,"r_vtxY_/F");
	reco_t_->Branch("r_vtxZ",&r_vtxZ_,"r_vtxZ_/F");
	reco_t_->Branch("r_vtxT",&r_vtxT_,"r_vtxT_/F");
	reco_t_->Branch("r_dirTheta",&r_dirTheta_,"r_dirTheta_/F");
	reco_t_->Branch("r_dirPhi",&r_dirPhi_,"r_dirPhi_/F");
	if(save_extra_)
	{
		reco_t_->Branch("r_raw_hit_map_origin",&r_raw_hit_map_origin_,"r_raw_hit_map_origin_[64][64]/b");
		reco_t_->Branch("r_raw_charge_map_origin",&r_raw_charge_map_origin_,"r_raw_charge_map_origin_[64][64]/b");
		reco_t_->Branch("r_raw_time_map_origin",&r_raw_time_map_origin_,"r_raw_time_map_origin_[64][64]/b");
		reco_t_->Branch("r_filtered_hit_map_origin",&r_filtered_hit_map_origin_,"r_filtered_hit_map_origin_[64][64]/b");
		reco_t_->Branch("r_filtered_charge_map_origin",&r_filtered_charge_map_origin_,"r_filtered_charge_map_origin_[64][64]/b");
		reco_t_->Branch("r_filtered_time_map_origin",&r_filtered_time_map_origin_,"r_filtered_time_map_origin_[64][64]/b");
		reco_t_->Branch("r_raw_hit_map_vtx",&r_raw_hit_map_vtx_,"r_raw_hit_map_vtx_[64][64]/b");
		reco_t_->Branch("r_filtered_hit_map_vtx",&r_filtered_hit_map_vtx_,"r_filtered_hit_map_vtx_[64][64]/b");
		reco_t_->Branch("r_filtered_charge_map_vtx",&r_filtered_charge_map_vtx_,"r_filtered_charge_map_vtx_[64][64]/b");
		reco_t_->Branch("r_filtered_time_map_vtx",&r_filtered_time_map_vtx_,"r_filtered_time_map_vtx_[64][64]/b");
	}
	reco_t_->Branch("r_raw_charge_map_vtx",&r_raw_charge_map_vtx_,"r_raw_charge_map_vtx_[64][64]/b");
	reco_t_->Branch("r_raw_time_map_vtx",&r_raw_time_map_vtx_,"r_raw_time_map_vtx_[64][64]/b");
	reco_t_->Branch("r_filtered_hit_hough_map_vtx",&r_filtered_hit_hough_map_vtx_,"r_filtered_hit_hough_map_vtx_[64][64]/b");
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
	map_f_ = new TFile(output_file_, "RECREATE", "HitMapFile");
	true_t_->Write();
	reco_t_->Write();
	map_f_->Close();

	std::cout << "Total = " << (ev-1) << ", Valid = " << valid_events << std::endl;

    delete reco;
}

void WCSimMapper::resetVariables()
{
	// True TTree variables (Used as labels in the CVN)
	t_nu_ = -999;
	t_code_ = -999;
	t_vtxX_ = -999;
	t_vtxY_ = -999;
	t_vtxZ_ = -999;
	t_vtxT_ = -999;
	t_nuEnergy_ = -999;
	memset(t_p_energies_, 0, sizeof(t_p_energies_));
	memset(t_p_pdgs_, 0, sizeof(t_p_pdgs_));
	memset(t_p_dirTheta_, 0, sizeof(t_p_dirTheta_));
	memset(t_p_dirPhi_, 0, sizeof(t_p_dirPhi_));

	// Reco TTree variables (Used as examples in the CVN)
	r_raw_num_hits_ = -999;
	r_raw_total_digi_q_ = -999;
	r_filtered_num_hits_ = -999;
	r_filtered_total_digi_q_ = -999;
	r_num_hough_rings_ = -999;
	r_first_ring_height_ = -999;
	r_last_ring_height_ = -999;
	r_vtxX_ = -999;
	r_vtxY_ = -999;
	r_vtxZ_ = -999;
	r_vtxT_ = -999;
	r_dirTheta_ = -999;
	r_dirPhi_ = -999;

	memset(r_raw_hit_map_origin_, 0, sizeof(r_raw_hit_map_origin_));
	memset(r_raw_charge_map_origin_, 0, sizeof(r_raw_charge_map_origin_));
	memset(r_raw_time_map_origin_, 0, sizeof(r_raw_time_map_origin_));
	memset(r_filtered_hit_map_origin_, 0, sizeof(r_filtered_hit_map_origin_));
	memset(r_filtered_charge_map_origin_, 0, sizeof(r_filtered_charge_map_origin_));
	memset(r_filtered_time_map_origin_, 0, sizeof(r_filtered_time_map_origin_));
	memset(r_raw_hit_map_vtx_, 0, sizeof(r_raw_hit_map_vtx_));
	memset(r_raw_charge_map_vtx_, 0, sizeof(r_raw_charge_map_vtx_));
	memset(r_raw_time_map_vtx_, 0, sizeof(r_raw_time_map_vtx_));
	memset(r_filtered_hit_map_vtx_, 0, sizeof(r_filtered_hit_map_vtx_));
	memset(r_filtered_charge_map_vtx_, 0, sizeof(r_filtered_charge_map_vtx_));
	memset(r_filtered_time_map_vtx_, 0, sizeof(r_filtered_time_map_vtx_));
	memset(r_filtered_hit_hough_map_vtx_, 0, sizeof(r_filtered_hit_hough_map_vtx_));
}

void WCSimMapper::fillTrueTree()
{
	WCSimTruthSummary truth_sum = WCSimInterface::Instance()->GetTruthSummary();
	t_nu_ = truth_sum.GetBeamPDG();
	t_code_ = truth_sum.GetInteractionMode();
	t_vtxX_ = truth_sum.GetVertexX()/10;
	t_vtxY_ = truth_sum.GetVertexY()/10;
	t_vtxZ_ = truth_sum.GetVertexZ()/10;
	t_vtxT_ = truth_sum.GetVertexT();
	t_nuEnergy_ = truth_sum.GetBeamEnergy();

	// First we get all the particle info
	std::vector<int> pdgs = truth_sum.GetPrimaryPDGs();
	std::vector<double> energies = truth_sum.GetPrimaryEnergies();
	std::vector<TVector3> dirs = truth_sum.GetPrimaryDirs();

	// Create a vector of tuples with all the primary particle info we need
	std::vector<std::tuple<double, int, float, float>> primaries; 
	for(int p=0; p<dirs.size(); p++)
	{
		float dir_theta = dirs[p].Z();
		float dir_phi = TMath::ATan2(dirs[p].Y(), dirs[p].X())*(180.0/TMath::Pi());
		primaries.push_back(std::make_tuple(energies[p], pdgs[p], dir_theta, dir_phi));
	}

	std::sort(primaries.begin(), primaries.end()); // Sort the primaries vector by energy
	std::reverse(primaries.begin(), primaries.end());  // Get them in descending order of energy

	// Store the top 10 primaries
	for(int p=0; p<10; p++) {
		if(p<primaries.size()) {
			t_p_energies_[p] = std::get<0>(primaries[p]);
			t_p_pdgs_[p] = std::get<1>(primaries[p]);
			t_p_dirTheta_[p] = std::get<2>(primaries[p]);
			t_p_dirPhi_[p] = std::get<3>(primaries[p]);
		} else {
			t_p_energies_[p] = -999;
			t_p_pdgs_[p] = -999;
			t_p_dirTheta_[p] = -999;
			t_p_dirPhi_[p] = -999;
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
	r_num_hough_rings_ = reco_event_->GetNRings();
	r_first_ring_height_ = reco_event_->GetRing(0)->GetHeight();
	r_last_ring_height_ = reco_event_->GetRing(r_num_hough_rings_-1)->GetHeight();

	r_vtxX_ = reco_event_->GetRing(0)->GetVtxX();
	r_vtxY_ = reco_event_->GetRing(0)->GetVtxY();
	r_vtxZ_ = reco_event_->GetRing(0)->GetVtxZ();
	r_vtxT_ = reco_event_->GetVtxTime();

	float dirX = reco_event_->GetRing(0)->GetDirX();
	float dirY = reco_event_->GetRing(0)->GetDirY();
	r_dirPhi_ = TMath::ATan2(dirY, dirX) * (180.0/TMath::Pi());
	r_dirTheta_ = reco_event_->GetRing(0)->GetDirZ();

	WCSimHoughTransformArray* hough_array = reco_event_->GetHoughArray();
	WCSimHoughTransform* hough = hough_array->GetHoughTransform(0);
	TH2D* r_filtered_hit_hough_map_vtx_h = hough->GetTH2D("h");

	// We can then generate the hit and time maps in both raw and "vertex view" forms
	TH2D* r_raw_hit_map_origin_h = new TH2D("r_raw_hit_map_origin_h", "r_raw_hit_map_origin_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_raw_charge_map_origin_h = new TH2D("r_raw_charge_map_origin_h", "r_raw_charge_map_origin_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_raw_time_map_origin_h = new TH2D("r_raw_time_map_origin_h", "r_raw_time_map_origin_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_filtered_hit_map_origin_h = new TH2D("r_filtered_hit_map_origin_h", "r_filtered_hit_map_origin_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_filtered_charge_map_origin_h = new TH2D("r_filtered_charge_map_origin_h", "r_filtered_charge_map_origin_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_filtered_time_map_origin_h = new TH2D("r_filtered_time_map_origin_h", "r_filtered_time_map_origin_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_raw_hit_map_vtx_h = new TH2D("r_raw_hit_map_vtx_h", "r_raw_hit_map_vtx_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_raw_charge_map_vtx_h = new TH2D("r_raw_charge_map_vtx_h", "r_raw_charge_map_vtx_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_raw_time_map_vtx_h = new TH2D("r_raw_time_map_vtx_h", "r_raw_time_map_vtx_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_filtered_hit_map_vtx_h = new TH2D("r_filtered_hit_map_vtx_h", "r_filtered_hit_map_vtx_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_filtered_charge_map_vtx_h = new TH2D("r_filtered_charge_map_vtx_h", "r_filtered_charge_map_vtx_h", 64, -180, 180, 64, -1, 1);
	TH2D* r_filtered_time_map_vtx_h = new TH2D("r_filtered_time_map_vtx_h", "r_filtered_time_map_vtx_h", 64, -180, 180, 64, -1, 1);

	// Loop through the digi hits
	r_raw_num_hits_ = reco_event_->GetNDigits();
	r_raw_total_digi_q_ = 0;
	float first_hit_time = 1000000;
	for(int h=0; h<r_raw_num_hits_; h++)
	{
		// Get the hit position and deposited charge
		float hit_x = reco_event_->GetDigit(h)->GetX();
		float hit_y = reco_event_->GetDigit(h)->GetY();
		float hit_z = reco_event_->GetDigit(h)->GetZ();
		float hit_q = reco_event_->GetDigit(h)->GetRawQPEs();
		r_raw_total_digi_q_ += hit_q;

		// Is the hit time the earliest in the event?
		float digi_hit_time = reco_event_->GetDigit(h)->GetRawTime();
		if (digi_hit_time < first_hit_time) 
		{ 
			first_hit_time = digi_hit_time; 
		}

		// Calculate origin phi and theta and fill hists
		float hit_phi_origin = TMath::ATan2(hit_y, hit_x) * (180.0/TMath::Pi());
		float hit_theta_origin = GetCosTheta(hit_x, hit_y, hit_z);

		r_raw_hit_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, 1);
		r_raw_charge_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, hit_q);
		int bin_num = r_raw_time_map_origin_h->FindBin(hit_phi_origin, hit_theta_origin);
		float currentBinTime = r_raw_time_map_origin_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			r_raw_time_map_origin_h->SetBinContent(bin_num, digi_hit_time);
		}

		// Calculate the "vertex view" phi and theta and fill hists 
		float hit_phi_vtx = TMath::ATan2((hit_y-r_vtxY_), (hit_x-r_vtxX_)) * (180.0/TMath::Pi());
		float hit_theta_vtx = GetCosTheta((hit_x-r_vtxX_), (hit_y-r_vtxY_), (hit_z-r_vtxZ_));

		r_raw_hit_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, 1);
		r_raw_charge_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, hit_q);
		bin_num = r_raw_time_map_vtx_h->FindBin(hit_phi_vtx, hit_theta_vtx);
		currentBinTime = r_raw_time_map_vtx_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			r_raw_time_map_vtx_h->SetBinContent(bin_num, digi_hit_time);
		}
	} 

	// We now need to move all timestamps to be from zero
	for(int i=0; i<(r_raw_time_map_origin_h->GetNbinsX()+1); i++)
	{
		for(int j=0; j<(r_raw_time_map_origin_h->GetNbinsY()+1); j++)
		{
			float raw_content = r_raw_time_map_origin_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				r_raw_time_map_origin_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
			raw_content = r_raw_time_map_vtx_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				r_raw_time_map_vtx_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
		}	
	} 

	// Loop through filtered digi hits
	r_filtered_num_hits_ = reco_event_->GetNFilterDigits();
	r_filtered_total_digi_q_ = 0;
	first_hit_time = 1000000;
	for(int h=0; h<r_filtered_num_hits_; h++)
	{
		// Get the hit position and deposited charge
		float hit_x = reco_event_->GetFilterDigit(h)->GetX();
		float hit_y = reco_event_->GetFilterDigit(h)->GetY();
		float hit_z = reco_event_->GetFilterDigit(h)->GetZ();
		float hit_q = reco_event_->GetFilterDigit(h)->GetRawQPEs();
		r_filtered_total_digi_q_ += hit_q;

		// Is the hit time the earliest in the event?
		float digi_hit_time = reco_event_->GetFilterDigit(h)->GetRawTime();
		if (digi_hit_time < first_hit_time) 
		{ 
			first_hit_time = digi_hit_time; 
		}

		// Calculate origin phi and theta and fill hists
		float hit_phi_origin = TMath::ATan2(hit_y, hit_x) * (180.0/TMath::Pi());
		float hit_theta_origin = GetCosTheta(hit_x, hit_y, hit_z);

		r_filtered_hit_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, 1);
		r_filtered_charge_map_origin_h->Fill(hit_phi_origin, hit_theta_origin, hit_q);
		int bin_num = r_filtered_time_map_origin_h->FindBin(hit_phi_origin, hit_theta_origin);
		float currentBinTime = r_filtered_time_map_origin_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			r_filtered_time_map_origin_h->SetBinContent(bin_num, digi_hit_time);
		}

		// Calculate the "vertex view" phi and theta and fill hists 
		float hit_phi_vtx = TMath::ATan2((hit_y-r_vtxY_), (hit_x-r_vtxX_)) * (180.0/TMath::Pi());
		float hit_theta_vtx = GetCosTheta((hit_x-r_vtxX_), (hit_y-r_vtxY_), (hit_z-r_vtxZ_));

		r_filtered_hit_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, 1);
		r_filtered_charge_map_vtx_h->Fill(hit_phi_vtx, hit_theta_vtx, hit_q);
		bin_num = r_filtered_time_map_vtx_h->FindBin(hit_phi_vtx, hit_theta_vtx);
		currentBinTime = r_filtered_time_map_vtx_h->GetBinContent(bin_num);
		if ((currentBinTime == 0) || (currentBinTime != 0 && digi_hit_time<currentBinTime)) 
		{
			r_filtered_time_map_vtx_h->SetBinContent(bin_num, digi_hit_time);
		}
	}

	// We now need to move all timestamps to be from zero
	for(int i=0; i<(r_filtered_time_map_origin_h->GetNbinsX()+1); i++)
	{
		for(int j=0; j<(r_filtered_time_map_origin_h->GetNbinsY()+1); j++)
		{
			float raw_content = r_filtered_time_map_origin_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				r_filtered_time_map_origin_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
			raw_content = r_filtered_time_map_vtx_h->GetBinContent(i,j);
			if(raw_content != 0.0)
			{
				r_filtered_time_map_vtx_h->SetBinContent(i, j, (raw_content-first_hit_time));
			}
		}	
	}

	std::vector<TH2D*> hists;
	std::vector<float> range_low;
	std::vector<float> range_high;
	
	hists.push_back(r_raw_hit_map_origin_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_raw_charge_map_origin_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_raw_time_map_origin_h);
	range_low.push_back(0.0);
	range_high.push_back(80.0);
	hists.push_back(r_filtered_hit_map_origin_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_filtered_charge_map_origin_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_filtered_time_map_origin_h);
	range_low.push_back(0.0);
	range_high.push_back(80.0);
	hists.push_back(r_raw_hit_map_vtx_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_raw_charge_map_vtx_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_raw_time_map_vtx_h);
	range_low.push_back(0.0);
	range_high.push_back(80.0);
	hists.push_back(r_filtered_hit_map_vtx_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_filtered_charge_map_vtx_h);
	range_low.push_back(0.0);
	range_high.push_back(15.0);
	hists.push_back(r_filtered_time_map_vtx_h);
	range_low.push_back(0.0);
	range_high.push_back(80.0);
	hists.push_back(r_filtered_hit_hough_map_vtx_h);
	range_low.push_back(0.0);
	range_high.push_back(1500.0);

	// First we clip all the histograms so the values lie in the correct ranges
	// Then we normalise to the range [0, 255]
	for(int hist=0; hist<13; hist++)
	{
		for(int i=0; i<64; i++)
		{
			for(int j=0; j<64; j++)
			{
				float content = hists[hist]->GetBinContent(i,j);
				//if(hist == 1) std::cout << content << " -> ";
				if(content > range_high[hist]) content = range_high[hist];
				if(content < range_low[hist]) content = range_low[hist];
				content = content - range_low[hist];
				content = content / (range_high[hist]-range_low[hist]);
				content = content * 255.0;
				content = round(content);
				//UChar_t c_content = round(content);
				//int i_content = c_content;
				//if(hist == 1) std::cout << content << " (" << i_content << ")" << std::endl;
				hists[hist]->SetBinContent(i,j,content);
			}
		}
	}

	// We now fill the arrays we save to the reco tree
	for(int i=0; i<64; i++)
	{
		for(int j=0; j<64; j++)
		{
			r_raw_hit_map_origin_[j][i] = (UChar_t)hists[0]->GetBinContent(i,j);
			r_raw_charge_map_origin_[j][i] = (UChar_t)hists[1]->GetBinContent(i,j);
			r_raw_time_map_origin_[j][i] = (UChar_t)hists[2]->GetBinContent(i,j);
			r_filtered_hit_map_origin_[j][i] = (UChar_t)hists[3]->GetBinContent(i,j);
			r_filtered_charge_map_origin_[j][i] = (UChar_t)hists[4]->GetBinContent(i,j);
			r_filtered_time_map_origin_[j][i] = (UChar_t)hists[5]->GetBinContent(i,j);
			r_raw_hit_map_vtx_[j][i] = (UChar_t)hists[6]->GetBinContent(i,j);
			r_raw_charge_map_vtx_[j][i] = (UChar_t)hists[7]->GetBinContent(i,j);
			r_raw_time_map_vtx_[j][i] = (UChar_t)hists[8]->GetBinContent(i,j);
			r_filtered_hit_map_vtx_[j][i] = (UChar_t)hists[9]->GetBinContent(i,j);
			r_filtered_charge_map_vtx_[j][i] = (UChar_t)hists[10]->GetBinContent(i,j);
			r_filtered_time_map_vtx_[j][i] = (UChar_t)hists[11]->GetBinContent(i,j);
			r_filtered_hit_hough_map_vtx_[j][i] = (UChar_t)hists[12]->GetBinContent(i,j);
		}	
	}

	delete r_raw_hit_map_origin_h;
	delete r_raw_charge_map_origin_h;
	delete r_raw_time_map_origin_h;
	delete r_filtered_hit_map_origin_h;
	delete r_filtered_charge_map_origin_h;
	delete r_filtered_time_map_origin_h;
	delete r_raw_hit_map_vtx_h;
	delete r_raw_charge_map_vtx_h;
	delete r_raw_time_map_vtx_h;
	delete r_filtered_hit_map_vtx_h;
	delete r_filtered_charge_map_vtx_h;
	delete r_filtered_time_map_vtx_h;

	reco_t_->Fill();
}