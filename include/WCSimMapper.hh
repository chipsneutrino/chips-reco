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
		WCSimMapper(const char* in_file, const char* out_file, int max_files, int category, int pdg_code);
		virtual ~WCSimMapper() {};
		
        void run();

	private:
		void resetVariables();
		void fillTrueTree();
		void fillRecoTree();
		void writeFile();

		float GetCosTheta(float x, float y, float z);

		// General variables
        const char* input_file_;	///< Input file to process
        int max_events_;			///< Max number of events to process
		int pdg_code_;				///< PDG code to select correct primary track
		TFile *map_f_;				///< Output File
		TTree *true_t_;				///< Truth TTree for CVN labelling
		TTree *reco_t_;				///< Reco TTree for CVN input

		WCSimRecoEvent* reco_event_;///< Current reco event

		// True TTree variables (Used as labels in the CVN)
		int true_n_hits_;			///< Number Cherenkov hits
		int true_category_;			///< PID category (truth)
		float true_vtx_x_;			///< Vertex x-position (truth)
		float true_vtx_y_;			///< Vertex y-position (truth)
		float true_vtx_z_;			///< Vertex z-position (truth)
		float true_vtx_t_;			///< Vertex time (truth)
		float true_dir_costheta_;	///< Lepton track costheta-direction (truth)
		float true_dir_phi_;		///< Lepton track phi-direction (truth)
		float true_nu_energy_;		///< Neutrino energy (truth)
		float true_lep_energy_;		///< Leading lepton energy (truth)

		// Reco TTree variables (Used as examples in the CVN)
		int num_hits_;				///< Number of hits in the event
		int num_filtered_hits_;		///< Number of hits in the event
		float total_digi_q_;		///< Total charge in the event
		int num_hough_rings_;		///< Number of seed rings in the event
		float first_ring_height_;	///< First seed ring hough height
		float last_ring_height_;	///< Lats seed ring hough height
		float vtx_x_;				///< Vertex x-position
		float vtx_y_;				///< Vertex y-position
		float vtx_z_;				///< Vertex z-position
		float vtx_t_;				///< Vertex time
		float dir_costheta_;		///< Lepton track costheta-direction
		float dir_phi_;				///< Lepton track phi-direction
		float hit_map_raw_[64][64];		///< Hit map for the event
		float time_map_raw_[64][64];	///< Hit time map for the event
		float hit_map_vtx_[64][64];		///< Hit map for the event
		float time_map_vtx_[64][64];	///< Hit time map for the event
		float hough_map_vtx_[64][64];	///< Hough map for the event 

		ClassDef(WCSimMapper,0)
};

#endif /* WCSIMMAPPER_HH_ */
