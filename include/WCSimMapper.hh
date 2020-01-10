/*
 * WCSimMapper.hh
 *
 *  Created on: 26 Nov 2019
 *      Author: Josh
 */

#ifndef WCSIMMAPPER_HH_
#define WCSIMMAPPER_HH_

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"

#include "WCSimRecoEvent.hh"

class WCSimRecoDigit;
class WCSimRecoRing;
class WCSimHoughTransform;
class WCSimHoughTransformArray;

class WCSimMapper 
{
	public:
		WCSimMapper(const char* in_file, const char* out_file, int max_files, int pdg_code);
		virtual ~WCSimMapper() {};
		
        void run();

	private:
		void resetVariables();
		void fillTrueTree();
		void fillRecoTree();
		void writeFile();

		float GetCosTheta(float x, float y, float z);

		// General variables
        const char* input_file_;		///< Input file to process
        int max_events_;				///< Max number of events to process
		int pdg_code_;					///< PDG code to select correct primary track
		TFile *map_f_;					///< Output File
		TTree *true_t_;					///< Truth TTree for CVN labelling
		TTree *reco_t_;					///< Reco TTree for CVN input

		WCSimRecoEvent* reco_event_;	///< Current reco event

		// True TTree variables (Used as labels in the CVN)
		int true_type_;					///< Interaction type (truth)
		float true_vtx_x_;				///< Vertex x-position (truth)
		float true_vtx_y_;				///< Vertex y-position (truth)
		float true_vtx_z_;				///< Vertex z-position (truth)
		float true_vtx_t_;				///< Vertex time (truth)
		float true_dir_costheta_;		///< Lepton track costheta-direction (truth)
		float true_dir_phi_;			///< Lepton track phi-direction (truth)
		float true_nu_energy_;			///< Neutrino energy (truth)
		float true_lep_energy_;			///< Leading lepton energy (truth)

		// Reco TTree variables (Used as examples in the CVN)
		int raw_num_hits_;				///< Raw number of hits in the event
		float raw_total_digi_q_;		///< Raw charge in the event
		int filtered_num_hits_;			///< Filtered number of hits in the event
		float filtered_total_digi_q_;	///< Filtered charge in the event
		int num_hough_rings_;			///< Number of seed rings in the event
		float first_ring_height_;		///< First seed ring hough height
		float last_ring_height_;		///< Lats seed ring hough height
		float vtx_x_;					///< Vertex x-position
		float vtx_y_;					///< Vertex y-position
		float vtx_z_;					///< Vertex z-position
		float vtx_t_;					///< Vertex time
		float dir_costheta_;			///< Lepton track costheta-direction
		float dir_phi_;					///< Lepton track phi-direction

		float raw_hit_map_origin_[64][64];				///< Raw hit map from origin
		float raw_charge_map_origin_[64][64];			///< Raw charge map from origin
		float raw_time_map_origin_[64][64];				///< Raw time map from origin
		float filtered_hit_map_origin_[64][64];			///< Filtered hit number map from origin
		float filtered_charge_map_origin_[64][64];		///< Filtered charge map from origin
		float filtered_time_map_origin_[64][64];			///< Filtered time map from origin

		float raw_hit_map_vtx_[64][64];					///< Raw hit map from reco vertex
		float raw_charge_map_vtx_[64][64];				///< Raw charge map from reco vertex
		float raw_time_map_vtx_[64][64];					///< Raw time map from reco vertex
		float filtered_hit_map_vtx_[64][64];				///< Filtered hit number map from reco vertex
		float filtered_charge_map_vtx_[64][64];			///< Filtered charge map from reco vertex 
		float filtered_time_map_vtx_[64][64];			///< Filtered time map from reco vertex

		//float raw_hit_hough_map_vtx_[64][64];			///< Raw hit Hough map from reco vertex 
		float filtered_hit_hough_map_vtx_[64][64];		///< Filtered hit Hough map from reco vertex 
		//float raw_charge_hough_map_vtx_[64][64];		///< Raw charge Hough map from reco vertex 
		//float filtered_charge_hough_map_vtx_[64][64];	///< Filtered charge Hough map from reco vertex 

		ClassDef(WCSimMapper,0)
};

#endif /* WCSIMMAPPER_HH_ */
