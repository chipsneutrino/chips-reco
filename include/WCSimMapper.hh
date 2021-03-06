/*
 * WCSimMapper.hh
 *
 *  Created on: 26 Nov 2019
 *      Author: Josh
 */

#pragma once

#include <math.h>

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
	WCSimMapper(const char *in_file, const char *out_file, int max_files, bool save_extra,
				float detector_height, float detector_radius);
	virtual ~WCSimMapper(){};

	void run();

private:
	void resetVariables();
	void fillTrees();

	float GetCosTheta(float x, float y, float z);
	void GetXValues(float *x_plus, float *x_minus, float rho, float z, float phi);

	// General variables
	const char *input_file_;  ///< Input file to process
	const char *output_file_; ///< Output file name
	int max_events_;		  ///< Max number of events to process
	bool save_extra_;		  ///< Should we save all the extra maps that we don't need?
	TFile *map_f_;			  ///< Output File
	TTree *true_t_;			  ///< Truth TTree for CVN labelling
	TTree *reco_t_;			  ///< Reco TTree for CVN input
	float height_;			  ///< Detector height
	float radius_;			  ///< Detector radius

	WCSimRecoEvent *reco_event_; ///< Current reco event

	// True TTree variables (Used as labels in the CVN)
	int t_nu_;				 ///< Beam PDG (truth)
	int t_code_;			 ///< Interaction type (truth)
	float t_vtxX_;			 ///< Vertex x-position (truth)
	float t_vtxY_;			 ///< Vertex y-position (truth)
	float t_vtxZ_;			 ///< Vertex z-position (truth)
	float t_vtxT_;			 ///< Vertex time (truth)
	float t_nuEnergy_;		 ///< Neutrino energy (truth)
	float t_p_energies_[10]; ///< Event primaries energies
	int t_p_pdgs_[10];		 ///< Event primaries PDG codes
	float t_p_dirTheta_[10]; ///< Event primaries theta direction
	float t_p_dirPhi_[10];	 ///< Event primaries phi direction

	// Reco TTree variables (Used as examples in the CVN)
	int r_num_hits_;			///< Raw number of hits in the event
	float r_total_digi_q_;		///< Raw charge in the event
	int r_num_hough_rings_;			///< Number of seed rings in the event
	float r_first_ring_height_;		///< First seed ring hough height
	float r_last_ring_height_;		///< Lats seed ring hough height
	float r_vtxX_;					///< Vertex x-position
	float r_vtxY_;					///< Vertex y-position
	float r_vtxZ_;					///< Vertex z-position
	float r_vtxT_;					///< Vertex time
	float r_dirTheta_;				///< Lepton track costheta-direction
	float r_dirPhi_;				///< Lepton track phi-direction

	UChar_t r_charge_map_origin_[64][64];	///< Raw charge map from origin
	UChar_t r_time_map_origin_[64][64];		///< Raw time map from origin
	UChar_t r_charge_map_iso_[64][64];		///< Raw charge map in iso space
	UChar_t r_time_map_iso_[64][64];		///< Raw time map in iso space
	UChar_t r_charge_map_vtx_[64][64];		///< Raw charge map from reco vertex
	UChar_t r_time_map_vtx_[64][64];		///< Raw time map from reco vertex
	UChar_t r_hough_map_vtx_[64][64]; 	///< Raw hit Hough map from reco vertex

	ClassDef(WCSimMapper, 0)
};
