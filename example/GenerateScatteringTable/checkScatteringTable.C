#include <iostream>
#include <string>
#include <cmath>

#include <TFile.h>
#include <TAxis.h>
#include <THnSparse.h>

void checkTable(THnSparseF *hTable, std::string id);

int main(int argc, char** argv){

  if(argc != 2){
    std::cerr << "Expect only a single argument: checkScatteringTable <filename>" << std::endl;
    return 1;
  }

  TFile *fInput = new TFile(argv[1],"read"); 
  THnSparseF *hTopTable = (THnSparseF*)fInput->Get("hScatteringTableTop");
  THnSparseF *hBottomTable = (THnSparseF*)fInput->Get("hScatteringTableBarrel");
  THnSparseF *hBarrelTable = (THnSparseF*)fInput->Get("hScatteringTableBottom");

  checkTable(hTopTable,"   Top");
  checkTable(hBottomTable,"Bottom");
  checkTable(hBarrelTable,"Barrel");

  return 0;
}

void checkTable(THnSparseF *hTable, std::string id){

  // Same number of bins in all dimensions
  int nBins = hTable->GetAxis(0)->GetNbins();

  unsigned int nEmptyBins = 0;
  unsigned int nNegBins = 0;
 
  double total = 0.0;
  double min = 999.;
  double max = -999.;
 
  for(int a = 0; a < nBins; ++a){
    for(int b = 0; b < nBins; ++b){
      for(int c = 0; c < nBins; ++c){
        for(int d = 0; d < nBins; ++d){
          for(int e = 0; e < nBins; ++e){
            for(int f = 0; f < nBins; ++f){
    
              int currBin[6] = {a+1,b+1,c+1,d+1,e+1,f+1};

              float content = hTable->GetBinContent(currBin);
              total += content;
              if(content < min) min = content;
              if(content > max) max = content;

              if(content < 0.){
                ++nNegBins;
              }
              else if(content < 1e-4){
                ++nEmptyBins;
              }
            }
          }
        }
      }
    }
  }

  std::cout << id << " Scattering Table: Found " << nEmptyBins << " entries with values < 1e-4 of which " << nNegBins << " were negative, of a total of " << pow(nBins,6) << " bins." << std::endl;
  std::cout << id << " Scattering Table: Average value = " << total / (double)pow(nBins,6) << std::endl;
  std::cout << id << " Scattering Table: Minimum value = " << min << std::endl;
  std::cout << id << " Scattering Table: Maximum value = " << max << std::endl;

}

