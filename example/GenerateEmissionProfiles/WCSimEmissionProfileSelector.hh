//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 25 22:54:22 2014 by ROOT version 5.26/00
// from TTree ntuple/my analysis ntuple
// found on file: /unix/fnu/ajperch/numu_5mrad_LE_1000_photons.root
//////////////////////////////////////////////////////////

#ifndef WCSimEmissionProfileSelector_h
#define WCSimEmissionProfileSelector_h

#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TVector3.h>

class WCSimEmissionProfileSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           eventID;
   Int_t           pdgCode;
   Int_t           trackID;
   Int_t           parentID;
   Int_t           processID;
   Double_t        energy;
   Double_t        wavelength;
   Int_t           optical;
   Int_t           scattered;
   Double_t        vtxX;
   Double_t        vtxY;
   Double_t        vtxZ;
   Double_t        endX;
   Double_t        endY;
   Double_t        endZ;
   Double_t        vtxdirX;
   Double_t        vtxdirY;
   Double_t        vtxdirZ;

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_pdgCode;   //!
   TBranch        *b_trackID;   //!
   TBranch        *b_parentID;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_wavelength;   //!
   TBranch        *b_optical;   //!
   TBranch        *b_scattered;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_endX;   //!
   TBranch        *b_endY;   //!
   TBranch        *b_endZ;   //!
   TBranch        *b_vtxdirX;   //!
   TBranch        *b_vtxdirY;   //!
   TBranch        *b_vtxdirZ;   //!

   WCSimEmissionProfileSelector(TTree * /*tree*/ =0) { }
   virtual ~WCSimEmissionProfileSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree * tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   /** We calculate cos(theta) relative to the vector from the vertex to the
    * centre of momentum of photon production.  This adds a photon vertex to
    * the vectors keeping track of the weighted sum
    */
   void FillComVtx( const Double_t &x, const Double_t &y, const Double_t &z, const Double_t &energy);

   /**
    * Here we keep track of the raw direction of each particle, so we can
    * subtract off the CoM angle later
    */
   void FillRawDir( const Double_t &x, const Double_t &y, const Double_t &z);

   /**
    * Fill the vector of z positions for each photon
    */
   void FillRawPos( const Double_t &z );

   /// Calculate the centre of mass of photon emission (energy-weighted)
   void CalculateCoM();

   /// Fill the 1D and 2D emission profile histograms
   void FillProfiles();

   /// Correctly normalise the 1D and 2D emission profiles
   void NormaliseProfiles();

   /// Save emission profiles into a Root file for WCSimAnalysis to read
   void SaveProfiles();

   /** Calculate the angle theta to the z axis given 3 components of a vector
    * @param x x-coordinate of direction vector
    * @param y y-coordinate of direction vector
    * @param z z-coordinate of direction vector
    * @return the angle theta to the z-axis (0,0,1)
    */
   Double_t Theta( const Double_t& x, const Double_t& y, const Double_t& z );

   /**
    * Get the quantum efficiency of a WCSim PMT at wavelength lambda
    * @param lambda Photon wavelength (nm)
    * @return Quantum efficiency
    */
   Double_t GetQE( const Double_t& lambda );

   /**
    * Process the first entry in the tree to set vector etc.
    */
   void ProcessFirst();

   // Andy's variables

   Int_t    fNBinsS;  /// Number of s bins for emission profile
   Double_t fSMin;    /// Minimum value of s
   Double_t fSMax;    /// Maximum value of s
   Int_t    fNBinsTh; /// Number of cosTheta bins for emission profile
   Double_t fThMin;   /// Minimum value of cosTheta
   Double_t fThMax;   /// Maximum value of cosTheta
   TString fSaveName; /// File name to save the emission profile as

   std::vector<Double_t> fComVtxX;
   std::vector<Double_t> fComVtxY;
   std::vector<Double_t> fComVtxZ;
   std::vector<TVector3> fComDir;
   std::vector<Int_t>    fNEvents;

   std::vector<std::vector<Double_t> > fRawPos;
   std::vector<std::vector<TVector3> > fRawDir;

   Float_t MM_TO_CM;
   Int_t   fEvents;
   Int_t   fMaxTreeEvents;
   Int_t   fOverallEventNum;
  

   TH1D f_hNEvents;
   TH1D f_hS;
   TH2D f_hSTh;
   TH2D f_hG;

   ClassDef(WCSimEmissionProfileSelector,0);
};

#endif

#ifdef WCSimEmissionProfileSelector_cxx
void WCSimEmissionProfileSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("pdgCode", &pdgCode, &b_pdgCode);
   fChain->SetBranchAddress("trackID", &trackID, &b_trackID);
   fChain->SetBranchAddress("parentID", &parentID, &b_parentID);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("wavelength", &wavelength, &b_wavelength);
   fChain->SetBranchAddress("optical", &optical, &b_optical);
   fChain->SetBranchAddress("scattered", &scattered, &b_scattered);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("endX", &endX, &b_endX);
   fChain->SetBranchAddress("endY", &endY, &b_endY);
   fChain->SetBranchAddress("endZ", &endZ, &b_endZ);
   fChain->SetBranchAddress("vtxdirX", &vtxdirX, &b_vtxdirX);
   fChain->SetBranchAddress("vtxdirY", &vtxdirY, &b_vtxdirY);
   fChain->SetBranchAddress("vtxdirZ", &vtxdirZ, &b_vtxdirZ);
}

Bool_t WCSimEmissionProfileSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}





#endif // #ifdef WCSimEmissionProfileSelector_cxx


