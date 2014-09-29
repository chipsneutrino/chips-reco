{
	gROOT->ProcessLine(".L WCSimConcatenateEmissionProfiles.cc+");
	WCSimConcatenateEmissionProfiles conc;
	conc.AddFile("emissionProfile.root",1500,14);
	conc.Run();


}
