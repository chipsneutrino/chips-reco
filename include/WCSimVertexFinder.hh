#ifndef WCSIMVERTEXFINDER_HH
#define WCSIMVERTEXFINDER_HH

#include "TObject.h"
#include "TMinuit.h"

#include "WCSimRecoVertex.hh"

class WCSimRecoDigit;
class WCSimRecoEvent;

#include <vector>

class WCSimVertexFinder: public TObject {
	public:
		static WCSimVertexFinder* Instance();

		static void PointFitOnly(Bool_t yesno = 1);
		static void SimpleVertexOnly(Bool_t yesno = 1);
		static void UseTrueVertex(Bool_t yesno = 1);
		static void UseTruePosition(Bool_t yesno = 1);
		static void UseTrueDirection(Bool_t yesno = 1);
		static void SeedWithTrueVertex(Bool_t yesno = 1);
		static void SeedWithSimpleVertex(Bool_t yesno = 1);
		static void SeedWithQuadruples(Bool_t yesno = 1);
		static void NumSeeds(Int_t nseeds);

		static void FitWeights(Double_t tw, Double_t cw);

		static void FitTimeParams();
		static void FixTimeParams(Double_t param0);

		static void FitConeParams();
		static void FixConeParams(Double_t param0, Double_t param1, Double_t param2);

		static void FixVertexBias(Double_t bias);

		static void PrintParameters();
		void RunPrintParameters();

		void SetPointFitOnly(Bool_t yesno) {
			fPointFitOnly = yesno;
		}
		void SetSimpleVertexOnly(Bool_t yesno) {
			fSimpleVertexOnly = yesno;
		}
		void UsingTrueVertex(Bool_t yesno) {
			fUseTrueVertex = yesno;
		}
		void UsingTruePosition(Bool_t yesno) {
			fUseTruePosition = yesno;
		}
		void UsingTrueDirection(Bool_t yesno) {
			fUseTrueDirection = yesno;
		}
		void SeedingWithTrueVertex(Bool_t yesno) {
			fSeedWithTrueVertex = yesno;
		}
		void SeedingWithSimpleVertex(Bool_t yesno) {
			fSeedWithSimpleVertex = yesno;
		}
		void SeedingWithQuadruples(Bool_t yesno) {
			fSeedWithQuadruples = yesno;
		}
		void SetNumSeeds(Int_t nseeds) {
			if (nseeds > 0)
				fNumSeeds = nseeds;
		}

		void SetFitTimeParams() {
			fFitTimeParams = 1;
		}
		void SetFixTimeParams(Double_t param0) {
			fFitTimeParams = 0;
			fFixTimeParam0 = param0;
		}

		void SetFitConeParams() {
			fFitConeParams = 1;
		}
		void SetFixConeParams(Double_t param0, Double_t param1, Double_t param2) {
			fFitConeParams = 0;
			fFixConeParam0 = param0;
			fFixConeParam1 = param1;
			fFixConeParam2 = param2;
		}

		void SetFitWeights(Double_t tw, Double_t cw) {
			fTimeFitWeight = tw;
			fConeFitWeight = cw;
		}

		void SetVertexBias(Double_t bias) {
			if (bias > 0.0) {
				fFixVertexBias = 1;
				fVertexBias = bias;
			} else {
				fFixVertexBias = 0;
				fVertexBias = bias;
			}
		}

		// Run Vertex Fitter
		WCSimRecoVertex* Run(WCSimRecoEvent* event);
		WCSimRecoVertex* RunSimpleFit(WCSimRecoEvent* event);
		WCSimRecoVertex* RunPointFit(WCSimRecoEvent* event);
		WCSimRecoVertex* RunExtendedFit(WCSimRecoEvent* event);

		// Run Vertex Fitter (using True Vertex)
		WCSimRecoVertex* RunPointFitFromTruth(WCSimRecoEvent* event);
		WCSimRecoVertex* RunExtendedFitFromTruth(WCSimRecoEvent* event);

		// Run Vertex Fitter (using Given Vertex)
		WCSimRecoVertex* RunPointFit(WCSimRecoEvent* event, Double_t vx, Double_t vy, Double_t vz);
		WCSimRecoVertex* RunPointFit(WCSimRecoEvent* event, Double_t vx, Double_t vy, Double_t vz, Double_t px,
				Double_t py, Double_t pz);
		WCSimRecoVertex* RunExtendedFit(WCSimRecoEvent* event, Double_t vx, Double_t vy, Double_t vz, Double_t px,
				Double_t py, Double_t pz);
		WCSimRecoVertex* Run(WCSimRecoEvent* event, Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py,
				Double_t pz);

		// External Methods: Access to Vertices
		WCSimRecoVertex* GetSimplePosition() {
			return fSimplePosition;
		}
		WCSimRecoVertex* GetSimpleDirection() {
			return fSimpleDirection;
		}
		WCSimRecoVertex* GetPointPosition() {
			return fPointPosition;
		}
		WCSimRecoVertex* GetPointDirection() {
			return fPointDirection;
		}
		WCSimRecoVertex* GetPointVertex() {
			return fPointVertex;
		}
		WCSimRecoVertex* GetExtendedVertex() {
			return fExtendedVertex;
		}

	private:
		// Internal Methods: Build True Vertex
		WCSimRecoVertex* BuildTrueVertex();

		// Internal Methods: Build Dummy Vertex
		WCSimRecoVertex* BuildDummyVertex();

		// Internal Methods: Fit the Vertex and Direction
		WCSimRecoVertex* FindSimplePosition();
		WCSimRecoVertex* FindSeedPosition();
		WCSimRecoVertex* FindSimpleDirection(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FindSeedDirection(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FitPointPosition(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FitPointDirection(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FitPointVertex(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FitExtendedVertex(WCSimRecoVertex* vertex);

		// Internal Methods: Correct Vertex Bias
		WCSimRecoVertex* CorrectExtendedVertex(WCSimRecoVertex* vertex);

		// Internal Methods: Fit the Vertex and Direction (using Given Vertex)
		WCSimRecoVertex* FixSimplePosition(Double_t vx, Double_t vy, Double_t vz);
		WCSimRecoVertex* FixSimpleDirection(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py,
				Double_t pz);
		WCSimRecoVertex* FixPointPosition(Double_t vx, Double_t vy, Double_t vz);
		WCSimRecoVertex* FixPointDirection(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py,
				Double_t pz);
		WCSimRecoVertex* FixPointVertex(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py, Double_t pz);
		WCSimRecoVertex* FixExtendedVertex(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py,
				Double_t pz);

		// Internal Methods: Run Fits with Minuit
		WCSimRecoVertex* FitPointPositionWithMinuit(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FitPointDirectionWithMinuit(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FitPointVertexWithMinuit(WCSimRecoVertex* vertex);
		WCSimRecoVertex* FitExtendedVertexWithMinuit(WCSimRecoVertex* vertex);

		// Internal Methods: Calculate Residuals
		void PointResiduals(Double_t vx, Double_t vy, Double_t vz);
		void PointResiduals(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py, Double_t pz);
		void ExtendedResiduals(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py, Double_t pz);

		// Simple method to just seed the time as the average hit time.
		void GetTimeSeed(WCSimRecoEvent *evt);

	public:

		// External Methods: Point Position Figure of Merit
		void PointPositionChi2(Double_t vx, Double_t vy, Double_t vz, Double_t& vtime, Double_t& chi2);
		void PointPositionChi2(Double_t& vtime, Double_t& chi2);

		// External Methods: Point Direction Figure of Merit
		void PointDirectionChi2(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py, Double_t pz,
				Double_t& vangle, Double_t& chi2);
		void PointDirectionChi2(Double_t& vangle, Double_t& chi2);

		// External Methods: Point Vertex Figure of Merit
		void PointVertexChi2(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py, Double_t pz,
				Double_t& vangle, Double_t& vtime, Double_t& chi2);
		void PointVertexChi2(Double_t& vangle, Double_t& vtime, Double_t& chi2);

		// External Methods: Extended Vertex Figure of Merit
		void ExtendedVertexChi2(Double_t vx, Double_t vy, Double_t vz, Double_t px, Double_t py, Double_t pz,
				Double_t& vangle, Double_t& vtime, Double_t& chi2);
		void ExtendedVertexChi2(Double_t& vangle, Double_t& vtime, Double_t& chi2);

		// External Methods: Penalty Term
		void PenaltyChi2(Double_t vx, Double_t vy, Double_t vz, Double_t& chi2);

		// External Methods: Penalty Term, to Fix Vertex or Direction
		void FixPositionChi2(Double_t vx, Double_t vy, Double_t vz, Double_t& chi2);
		void FixDirectionChi2(Double_t px, Double_t py, Double_t pz, Double_t& chi2);

		//  External Methods: Fit Profile in Angle
		void FitConePropertiesFoM(Double_t& vangle, Double_t& chi2);
		void ConePropertiesFoM(Double_t& chi2);

		void FitPointConePropertiesLnL(Double_t& vangle, Double_t& chi2);
		void FitExtendedConePropertiesLnL(Double_t& vangle, Double_t& chi2);
		void ConePropertiesLnL(Double_t vparam0, Double_t vparam1, Double_t vparam2, Double_t& vangle, Double_t& chi2);

		//  External Methods: Fit Profile in Time
		void FindSimpleTimeProperties(Double_t& vtime);
		void FitTimePropertiesFoM(Double_t& vtime, Double_t& chi2);
		void TimePropertiesFoM(Double_t vtime, Double_t& chi2);

		void FitPointTimePropertiesLnL(Double_t& vtime, Double_t& chi2);
		void FitExtendedTimePropertiesLnL(Double_t& vtime, Double_t& chi2);
		void TimePropertiesLnL(Double_t vtime, Double_t vparam, Double_t& chi2);

		// Count Iterations
		Int_t GetIterations() {
			return fItr;
		}

		void time_fit_itr() {
			fTimeFitItr++;
		}
		void cone_fit_itr() {
			fConeFitItr++;
		}
		void point_position_itr() {
			fPointPosItr++;
		}
		void point_direction_itr() {
			fPointDirItr++;
		}
		void point_vertex_itr() {
			fPointVtxItr++;
		}
		void extended_vertex_itr() {
			fExtendedVtxItr++;
		}

		void time_fit_reset_itr() {
			fTimeFitItr = 0;
		}
		void cone_fit_reset_itr() {
			fConeFitItr = 0;
		}
		void point_position_reset_itr() {
			fPointPosItr = 0;
		}
		void point_direction_reset_itr() {
			fPointDirItr = 0;
		}
		void point_vertex_reset_itr() {
			fPointVtxItr = 0;
		}
		void extended_vertex_reset_itr() {
			fExtendedVtxItr = 0;
		}

		Int_t time_fit_iterations() {
			return fTimeFitItr;
		}
		Int_t cone_fit_iterations() {
			return fConeFitItr;
		}
		Int_t point_position_iterations() {
			return fPointPosItr;
		}
		Int_t point_direction_iterations() {
			return fPointDirItr;
		}
		Int_t point_vertex_iterations() {
			return fPointVtxItr;
		}
		Int_t extended_vertex_iterations() {
			return fExtendedVtxItr;
		}

		// Reset Everything
		void Reset();

	private:

		WCSimVertexFinder();
		~WCSimVertexFinder();

		void Clear();

		Double_t fBaseFOM;
		Bool_t fPointFitOnly;
		Bool_t fSimpleVertexOnly;
		Bool_t fUseTrueVertex;
		Bool_t fUseTruePosition;
		Bool_t fUseTrueDirection;
		Bool_t fSeedWithTrueVertex;
		Bool_t fSeedWithSimpleVertex;
		Bool_t fSeedWithQuadruples;
		Int_t fNumSeeds;

		Double_t fTimeFitWeight;
		Double_t fConeFitWeight;

		Bool_t fFixVertexBias;
		Double_t fVertexBias;

		Double_t fInputVtxX;
		Double_t fInputVtxY;
		Double_t fInputVtxZ;
		Double_t fInputVtxTime;

		Double_t fInputDirX;
		Double_t fInputDirY;
		Double_t fInputDirZ;

		Double_t fInputVtxFOM;

		WCSimRecoVertex* fSimplePosition;
		WCSimRecoVertex* fSimpleDirection;
		WCSimRecoVertex* fPointPosition;
		WCSimRecoVertex* fPointDirection;
		WCSimRecoVertex* fPointVertex;
		WCSimRecoVertex* fExtendedVertex;

		TMinuit* fMinuitPointPosition;
		TMinuit* fMinuitPointDirection;
		TMinuit* fMinuitPointVertex;
		TMinuit* fMinuitExtendedVertex;

		TMinuit* fMinuitTimeFit;
		TMinuit* fMinuitConeFit;

		Bool_t fPass;
		Int_t fItr;

		Int_t fTimeFitItr;
		Int_t fConeFitItr;
		Int_t fPointPosItr;
		Int_t fPointDirItr;
		Int_t fPointVtxItr;
		Int_t fExtendedVtxItr;

		Double_t fSconeA;
		Double_t fSconeB;
		Double_t fSmu;
		Double_t fSel;
		Bool_t fIntegralsDone;

		Double_t fTimeSeed;

	public:

		Double_t fVtxX;
		Double_t fVtxY;
		Double_t fVtxZ;
		Double_t fVtxTime;

		Double_t fDirX;
		Double_t fDirY;
		Double_t fDirZ;

		Double_t fVtxFOM;

		// --- fitting ---
		Bool_t fFitTimeParams;
		Double_t fFixTimeParam0;

		Bool_t fFitConeParams;
		Double_t fFixConeParam0;
		Double_t fFixConeParam1;
		Double_t fFixConeParam2;

		// --- debug ---
		Double_t fTimeParam0;
		Double_t fConeParam0;
		Double_t fConeParam1;
		Double_t fConeParam2;

		ClassDef(WCSimVertexFinder,0)

};

#endif

