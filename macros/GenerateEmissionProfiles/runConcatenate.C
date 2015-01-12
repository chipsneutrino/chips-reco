{
	gROOT->ProcessLine(".L WCSimConcatenateEmissionProfiles.cc+");
	WCSimConcatenateEmissionProfiles conc;
	conc.AddFile("ep_nu_mu_1250.root",1250,14);
	conc.AddFile("ep_nu_mu_1500.root",1500,14);
	conc.AddFile("ep_nu_mu_1750.root",1750,14);
	conc.Run();


}
