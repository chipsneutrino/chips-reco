void wc_train_read_PID(){
	std::cout << "#### START ####" << std::endl;

	// Load the libraries we need...
	std::cout << "Loading libraries..." << std::endl;
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("libTMVA");
    TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
    TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
    gSystem->Load(libWCSimRoot.Data());
    gSystem->Load(libWCSimAnalysis.Data());

    // Define the training input directories to use full of WCSimAnalysis output files...
    std::string elTrainDir = "/unix/chips/jtingey/CHIPS/data/testing/CHIPS_1200_veto_LC/sk1pe/nue_cc_flux/NuMI/wc_Chips_1200_el_mu/";
    std::string muTrainDir = "/unix/chips/jtingey/CHIPS/data/testing/CHIPS_1200_veto_LC/sk1pe/numu_cc_flux/NuMI/wc_Chips_1200_el_mu/";
    std::string ncTrainDir = "/unix/chips/jtingey/CHIPS/data/testing/CHIPS_1200_veto_LC/sk1pe/numu_nc_flux/NuMI/wc_Chips_1200_el_mu/";

    // Define the reading(sample) input directories to use full of WCSimAnalysis output files...
    std::string elAllDir = "/unix/chips/jtingey/CHIPS/data/testing/CHIPS_1200_veto_LC/sk1pe/nue_all_flux/NuMI/wc_Chips_1200_el_mu/";
    std::string muAllDir = "/unix/chips/jtingey/CHIPS/data/testing/CHIPS_1200_veto_LC/sk1pe/numu_all_flux/NuMI/wc_Chips_1200_el_mu/";

    // Make a new WCSimPID and set the directories
    WCSimPID *pid = new WCSimPID();
    pid->SetTrainingDirs(elTrainDir, muTrainDir, ncTrainDir);
    pid->SetReadingDirs(elAllDir, muAllDir);

    // Train the two ANNs
    pid->Train();

    // Test the two ANNs with the sample
    pid->Read();

	std::cout << "#### END ####" << std::endl;
}
