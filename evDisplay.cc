#include <iostream>
#include <string>
#include <cstring>
#include "TEnv.h"
#include "TApplication.h"
#include "WCSimRecoEvDisplay.hh"

void PrintHelp();

int main(int argc,char* argv[])
{
	// Have to deal with this here since TApplication "usefully" strips out
	// root files from the argument list!
	std::string filename = "";
	bool useTTF = false;

	if(argc > 1){
		for(int i = 1; i < argc; ++i){
			// File name switch 
			if(std::strcmp(argv[i],"-f") == 0){
				if(argc >= i+1){
					filename = argv[i+1];
					std::cout << "== Input file = " << argv[i+1] << std::endl;
					++i;
				}
				else{
						std::cout << "== -f requires a file name!" << std::cout;
					return 0;
				}
			}
			// TrueType fonts switch
			else if(std::strcmp(argv[i],"-t") == 0){
				std::cout << "== Using TTF Fonts" << std::endl;
				useTTF = true;
			}
			// Help switch
			else if (std::strcmp(argv[i],"-h") == 0){
				PrintHelp();
				return 0;
			}
			else{
				std::cout << "Option " << argv[i] << " not recognised, showing usage information." << std::endl;
				PrintHelp();
				return 0;
			}
		}
	}

	if(useTTF){
		gEnv->SetValue("Unix.*.Root.UseTTFonts","true");
	}
	else{
		gEnv->SetValue("Unix.*.Root.UseTTFonts","false");
	}
	TApplication theApp("App", &argc, argv);

  // Popup the GUI...
	WCSimRecoEvDisplay *ev = new WCSimRecoEvDisplay(gClient->GetRoot(),1000,800);
	if(filename != ""){
		ev->SetInputFile(filename);
	}
	theApp.Run();
	return 0;

}

void PrintHelp(){
	std::cout << "Usage instructions for evDisplay" << std::endl;
	std::cout << "\t-h Displays the help message" << std::endl;
	std::cout << "\t-f <filename> Supply an input file" << std::endl;
	std::cout << "\t-t Enable useage of TTF fonts (slow over ssh)" << std::endl; 
}
