void combineMuonFitFiles(){
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
    TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
    gSystem->Load(libWCSimRoot.Data());
    gSystem->Load(libWCSimAnalysis.Data());

    const char* inputDir = "/unix/chips/jtingey/CHIPS/data/cosmicEvents/numu_cr_flux/NuMI/MuonLike";
    const char* outputName = "numu_cr_combined.root";

    //const char* inputDir = "/unix/chips/jtingey/CHIPS/data/cosmicEvents/numu_cc_flux/NuMI/MuonLike";
    //const char* outputName = "numu_cc_combined.root";

    gStyle->SetOptStat(1101);
    //gStyle->SetOptStat(0);

    char* dir = gSystem->ExpandPathName(inputDir);
    void* dirp = gSystem->OpenDirectory(dir);
    const char* entry;
    TString str;
    int n=0;
    int totalEvents = 0;
    int totalEscaped = 0;

    int fNVetoHits;
    int fNHits; //
    double fFracHitsUpstream; //
    double fFracHitsDownstream; //
    //double fFracHitsInBottom; //
    //double fFracHitsInTop; //
    double fFracHitsAboveMid; //
    double fFracHitsBelowMid; //
    double fFracNHitsInRing; //
    double fFracNHitsOutsideRing; //
    double fFracNHitsInRingHole; //
    double fEnergy; //
    double fVtxX; //
    double fVtxY; //
    double fVtxZ; //
    double fVtxRho; //
    double fDirX; //
    double fDirY; //
    double fDirZ; //
    bool fEscapes;

    TTree *PIDOutputTree = new TTree("PIDTree_ann","PIDTree_ann");
    PIDOutputTree->Branch("fNVetoHits",&fNVetoHits);
    PIDOutputTree->Branch("fNHits",&fNHits);
    PIDOutputTree->Branch("fFracHitsUpstream",&fFracHitsUpstream);
    PIDOutputTree->Branch("fFracHitsDownstream",&fFracHitsDownstream);
    //PIDOutputTree->Branch("fFracHitsInBottom",&fFracHitsInBottom);
    //PIDOutputTree->Branch("fFracHitsInTop",&fFracHitsInTop);
    PIDOutputTree->Branch("fFracHitsAboveMid",&fFracHitsAboveMid);
    PIDOutputTree->Branch("fFracHitsBelowMid",&fFracHitsBelowMid);
    PIDOutputTree->Branch("fFracNHitsInRing",&fFracNHitsInRing);
    PIDOutputTree->Branch("fFracNHitsOutsideRing",&fFracNHitsOutsideRing);
    PIDOutputTree->Branch("fFracNHitsInRingHole",&fFracNHitsInRingHole);
    PIDOutputTree->Branch("fEnergy",&fEnergy);
    PIDOutputTree->Branch("fVtxX",&fVtxX);
    PIDOutputTree->Branch("fVtxY",&fVtxY);
    PIDOutputTree->Branch("fVtxZ",&fVtxZ);
    PIDOutputTree->Branch("fVtxRho",&fVtxRho);
    PIDOutputTree->Branch("fDirX",&fDirX);
    PIDOutputTree->Branch("fDirY",&fDirY);
    PIDOutputTree->Branch("fDirZ",&fDirZ);
    PIDOutputTree->Branch("fEscapes",&fEscapes);

    while((entry = (char*)gSystem->GetDirEntry(dirp)) && n<2500){
        str = entry;
        if(str.EndsWith("_tree.root")){
            n++;
            std::cout << "Processing File-> " << entry << std::endl;
            TFile * inputFile = new TFile(gSystem->ConcatFileName(dir, entry), "READ");

            try{

            std::cout << "Checking the file..." << std::endl;
            if(!inputFile->GetListOfKeys()->Contains("fResultsTree")){ std::cout << "Skipping File" << std::endl; continue;}
            TTree * tree;
            tree = (TTree*) (inputFile->Get("fResultsTree"));

            //Set up the PidInfo...
            PidInfo *pidInfo = new PidInfo();
            TBranch *b_pi =  tree->GetBranch("PidInfo");
            b_pi->SetAddress(&pidInfo);

            //std::cout << tree->GetEntries() << std::endl;
            for(int evt=0; evt<tree->GetEntries(); evt++){
                b_pi->GetEntry(evt);

            	totalEvents++;

                fNVetoHits = pidInfo->GetNVetoHits();
                fNHits = pidInfo->GetNHits();
                fFracHitsUpstream = pidInfo->GetFracHitsUpstream();
                fFracHitsDownstream = pidInfo->GetFracHitsDownstream();
                //fFracHitsInBottom = pidInfo->GetFracHitsInBottom();
                //fFracHitsInTop = pidInfo->GetFracHitsInTop();
                fFracHitsAboveMid = pidInfo->GetFracHitsAboveMid();
                fFracHitsBelowMid = pidInfo->GetFracHitsBelowMid();
                fFracNHitsInRing = pidInfo->GetFracNHitsInRing();
                fFracNHitsOutsideRing = pidInfo->GetFracNHitsOutsideRing();
                fFracNHitsInRingHole = pidInfo->GetFracNHitsInRingHole();
                fEnergy = pidInfo->GetEnergy();
                fVtxX = pidInfo->GetVtxX();
                fVtxY = pidInfo->GetVtxY();
                fVtxZ = pidInfo->GetVtxZ();
                fVtxRho = pidInfo->GetVtxRho();
                fDirX = pidInfo->GetDirX();
                fDirY = pidInfo->GetDirY();
                fDirZ = pidInfo->GetDirZ();
                fEscapes = pidInfo->Escapes();

                PIDOutputTree->Fill();

                pidInfo->Clear();
            }

            }
            catch (...) { std::cout << "default exception"; }

            inputFile->Close();
        }
    }

    TFile * mainOutput = new TFile(outputName,"RECREATE");

    PIDOutputTree->Write();

    mainOutput->Close();

    std::cout << "Finished! Total Events-> " << totalEvents << ", numEscaped-> " << totalEscaped << std::endl;
}
