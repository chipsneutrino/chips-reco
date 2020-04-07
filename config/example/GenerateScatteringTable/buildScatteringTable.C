#include <iostream>
#include <string>

#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <THnSparse.h>

// Number of bins in each dimension
static const int nBins = 4;

void buildNewScatteringTable(std::string inputAll, std::string inputDirect, std::string output);
void LoopAndFill6DArray(TChain* chain, float qArray[][nBins][nBins][nBins][nBins][nBins], std::vector<std::vector<float> > &binArray, bool add4D, float q4DArray[][nBins][nBins][nBins], std::vector<std::string> &varNames, int location);
int GetBin(double val, std::vector<float> binArray);
void AddArrayToVector(std::vector<std::vector<float> > &binVector, float *array);
void CalculateBinVector(float min, float max, std::vector<float> &binVec);
void GetVarNames(std::vector<std::string> &vecName, bool type);

int main(int argc, char** argv){

  if(argc != 4){
    std::cout << "Usage: buildScatteringTable <normal_files> <direct_files> <output_file>" << std::endl;
    std::cout << "Options: <type> = new or old" << std::endl;
    return 1;
  }
  else{
    buildNewScatteringTable(argv[1],argv[2],argv[3]);
  }
  return 0;
}


void buildNewScatteringTable(std::string inputAll, std::string inputDirect, std::string output){

  // Get the data containing all light
  TChain *dAll = new TChain("scatteringNtp");
  dAll->Add(inputAll.c_str());

  // Get the data containing just the direct light
  TChain *dDirect = new TChain("scatteringNtp");
  dDirect->Add(inputDirect.c_str());

  // Need to fill 3 arrays
  // Bin definitions: vtxR, vtxZ, pmtZ, cosZeta, cosTheta, phi
  // Charge in each of the 6D phase-space from all light
  float chargeAllTop[nBins][nBins][nBins][nBins][nBins][nBins];
  float chargeAllBot[nBins][nBins][nBins][nBins][nBins][nBins];
  float chargeAllBar[nBins][nBins][nBins][nBins][nBins][nBins];
  // Charge in each of the 6D phase-space from direct light
  float chargeDirectTop[nBins][nBins][nBins][nBins][nBins][nBins];
  float chargeDirectBot[nBins][nBins][nBins][nBins][nBins][nBins];
  float chargeDirectBar[nBins][nBins][nBins][nBins][nBins][nBins];
  // Normalisation, where we integrate over cosTheta and phi
  float chargeDirectNormTop[nBins][nBins][nBins][nBins];
  float chargeDirectNormBot[nBins][nBins][nBins][nBins];
  float chargeDirectNormBar[nBins][nBins][nBins][nBins];
  // Initialise arrays
  for(int a = 0; a < nBins; ++a){
    for(int b = 0; b < nBins; ++b){
      for(int c = 0; c < nBins; ++c){
        for(int d = 0; d < nBins; ++d){
          chargeDirectNormTop[a][b][c][d] = 0.0;
          chargeDirectNormBot[a][b][c][d] = 0.0;
          chargeDirectNormBar[a][b][c][d] = 0.0;
          for(int e = 0; e < nBins; ++e){
            for(int f = 0; f < nBins; ++f){
              chargeAllTop[a][b][c][d][e][f] = 0.0;
              chargeAllBot[a][b][c][d][e][f] = 0.0;
              chargeAllBar[a][b][c][d][e][f] = 0.0;
              chargeDirectTop[a][b][c][d][e][f] = 0.0;
              chargeDirectBot[a][b][c][d][e][f] = 0.0;
              chargeDirectBar[a][b][c][d][e][f] = 0.0;
            }
          }
        }
      }
    }
  }

  // Create a 6-dimensional histogram.
  Int_t binsArray[6] = {nBins,nBins,nBins,nBins,nBins,nBins};
  // Index order: vtxR, vtxZ, v2pZeta, pmtPos, dirTheta, dirPhi
  Double_t minArrayCap[6] = {0., -10500., 0., 0., 0., -1*TMath::Pi()};
  Double_t maxArrayCap[6] = {12200., 10500., TMath::Pi(), 12200., TMath::Pi(), TMath::Pi()};
  Double_t minArrayBar[6] = {0., -10500., 0., -10500., 0., -1*TMath::Pi()};
  Double_t maxArrayBar[6] = {12200., 10500., TMath::Pi(), 10500., TMath::Pi(), TMath::Pi()};

  // Store all of the bins in a vector of vectors
  std::vector<std::vector<float> > binVectorCap;
  std::vector<std::vector<float> > binVectorBar;
  // The individual variable bins are hence just vectors.
  std::vector<float> vtxRBins;
  std::vector<float> vtxZBins;
  std::vector<float> v2pZetaBins;
  std::vector<float> pmtPosCapBins; // PMT R
  std::vector<float> pmtPosBarBins; // PMT Z
  std::vector<float> dirThetaBins;
  std::vector<float> dirPhiBins;
  // Calculate the bins
  CalculateBinVector(minArrayCap[0],maxArrayCap[0],vtxRBins);
  CalculateBinVector(minArrayCap[1],maxArrayCap[1],vtxZBins);
  CalculateBinVector(minArrayCap[2],maxArrayCap[2],v2pZetaBins);
  CalculateBinVector(minArrayCap[3],maxArrayCap[3],pmtPosCapBins); // PMT R
  CalculateBinVector(minArrayBar[3],maxArrayBar[3],pmtPosBarBins); // PMT Z
  CalculateBinVector(minArrayCap[4],maxArrayCap[4],dirThetaBins);
  CalculateBinVector(minArrayCap[5],maxArrayCap[5],dirPhiBins);
  // Cap bin vector
  binVectorCap.push_back(vtxRBins);
  binVectorCap.push_back(vtxZBins);
  binVectorCap.push_back(v2pZetaBins);
  binVectorCap.push_back(pmtPosCapBins);
  binVectorCap.push_back(dirThetaBins);
  binVectorCap.push_back(dirPhiBins);
  // Barrel bin vector
  binVectorBar.push_back(vtxRBins);
  binVectorBar.push_back(vtxZBins);
  binVectorBar.push_back(v2pZetaBins);
  binVectorBar.push_back(pmtPosBarBins);
  binVectorBar.push_back(dirThetaBins);
  binVectorBar.push_back(dirPhiBins);

  // The three individual scattering tables
  THnSparseF *hScatteringTableTop = new THnSparseF("hScatteringTableTop","",6,binsArray,minArrayCap,maxArrayCap);
  THnSparseF *hScatteringTableBot = new THnSparseF("hScatteringTableBottom","",6,binsArray,minArrayCap,maxArrayCap);
  THnSparseF *hScatteringTableBar = new THnSparseF("hScatteringTableBarrel","",6,binsArray,minArrayBar,maxArrayBar);

  // Make the variables names
  std::vector<std::string> varNames;
  GetVarNames(varNames,true);

  // Loop over all of the entries in the "all light" case
  LoopAndFill6DArray(dAll,chargeAllTop,binVectorCap,false,chargeDirectNormTop,varNames,0);
  LoopAndFill6DArray(dAll,chargeAllBot,binVectorCap,false,chargeDirectNormBot,varNames,2); 
  LoopAndFill6DArray(dAll,chargeAllBar,binVectorBar,false,chargeDirectNormBar,varNames,1);

  // Loop over all of the entries in the "all light" case
  LoopAndFill6DArray(dDirect,chargeDirectTop,binVectorCap,true,chargeDirectNormTop,varNames,0);
  LoopAndFill6DArray(dDirect,chargeDirectBot,binVectorCap,true,chargeDirectNormBot,varNames,2);
  LoopAndFill6DArray(dDirect,chargeDirectBar,binVectorBar,true,chargeDirectNormBar,varNames,1);

  // Fill the scattering table
  int totalBins = pow(nBins,6);
  double totalTop = 0.;
  double totalBot = 0.;
  double totalBar = 0.;
  for(int a = 0; a < nBins; ++a){
    for(int b = 0; b < nBins; ++b){
      for(int c = 0; c < nBins; ++c){
        for(int d = 0; d < nBins; ++d){
          for(int e = 0; e < nBins; ++e){
            for(int f = 0; f < nBins; ++f){
              // Which bin?
              int binToFill[6] = {a+1,b+1,c+1,d+1,e+1,f+1};
              float scatteringValue;
              // Top
              if(chargeDirectNormTop[a][b][c][d] > 1e-6){
                scatteringValue = (chargeAllTop[a][b][c][d][e][f] - chargeDirectTop[a][b][c][d][e][f]) / chargeDirectNormTop[a][b][c][d];
              }
              else{
                scatteringValue = 0;
                std::cout << "Empty bin in \"isotropic\" direct charge: " << a << ", " << b << ", " << c << ", " << d << std::endl;
              }
              hScatteringTableTop->SetBinContent(binToFill,scatteringValue);
              totalTop += scatteringValue;
              if(scatteringValue < 0){
//                std::cerr << "Negative scattering value (top): " << a << "," << b << "," << c << "," << d << "," << e << "," << f << " :: " << scatteringValue << std::endl;
              }
              // Bottom
              if(chargeDirectNormBot[a][b][c][d] > 1e-6){
                scatteringValue = (chargeAllBot[a][b][c][d][e][f] - chargeDirectBot[a][b][c][d][e][f]) / chargeDirectNormBot[a][b][c][d];
              }
              else{
                scatteringValue = 0;
                std::cout << "Empty bin in \"isotropic\" direct charge: " << a << ", " << b << ", " << c << ", " << d << std::endl;
              }
              hScatteringTableBot->SetBinContent(binToFill,scatteringValue);
              totalBot += scatteringValue;
              if(scatteringValue < 0){
//                std::cerr << "Negative scattering value (bot): " << a << "," << b << "," << c << "," << d << "," << e << "," << f << " :: " << scatteringValue << std::endl;
              }
              // Barrel
              if(chargeDirectNormBar[a][b][c][d] > 1e-6){
                scatteringValue = (chargeAllBar[a][b][c][d][e][f] - chargeDirectBar[a][b][c][d][e][f]) / chargeDirectNormBar[a][b][c][d];
              }
              else{
                scatteringValue = 0;
                std::cout << "Empty bin in \"isotropic\" direct charge: " << a << ", " << b << ", " << c << ", " << d << std::endl;
              }
              hScatteringTableBar->SetBinContent(binToFill,scatteringValue);
              totalBar += scatteringValue;
              if(scatteringValue < 0){
//                std::cerr << "Negative scattering value (bar): " << a << "," << b << "," << c << "," << d << "," << e << "," << f << " :: " << scatteringValue << std::endl;
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Average top    = " << totalTop / (double)totalBins << std::endl;
  std::cout << "Average bottom = " << totalBot / (double)totalBins << std::endl;
  std::cout << "Average barrel = " << totalBar / (double)totalBins << std::endl;
  std::cout << "Table has " << totalBins << " bins" << std::endl;
  std::cout << "Had " << dAll->GetEntries() << " and " << dDirect->GetEntries() << " in each tree" << std::endl;

  TFile *fTable = new TFile(output.c_str(),"recreate");
  hScatteringTableTop->Write();
  hScatteringTableBot->Write();
  hScatteringTableBar->Write();
  fTable->Close();
  delete fTable;
  delete hScatteringTableTop;
  delete hScatteringTableBot;
  delete hScatteringTableBar;
  delete dAll;
  delete dDirect;

}

void LoopAndFill6DArray(TChain* chain, float qArray[][nBins][nBins][nBins][nBins][nBins], std::vector<std::vector<float> > &binArray, bool add4D, float q4DArray[][nBins][nBins][nBins], std::vector<std::string> &varNames, int location){

  // Attach the 7 variables to the tree.
  double vtxR, vtxZ, pmtZ;
  double cosZeta, cosTheta, phi;
  double charge;

  double var1; // This is always vtxR
  double var2; // This is always vtxZ
  double var3; // New files: v2pZeta.  Old files: pmtZ
  double var4; // New files: pmtPos.   Old files: cosZeta
  double var5; // New files: dirTheta. Old files: cosTheta
  double var6; // New files: dirPhi.   Old files: phi
  double var7; // This is always charge
  int var8;    // This is always pmtLoc for the new files. Not used for the old files.

  chain->SetBranchAddress(varNames[0].c_str(),&var1);
  chain->SetBranchAddress(varNames[1].c_str(),&var2);
  chain->SetBranchAddress(varNames[2].c_str(),&var3);
  chain->SetBranchAddress(varNames[3].c_str(),&var4);
  chain->SetBranchAddress(varNames[4].c_str(),&var5);
  chain->SetBranchAddress(varNames[5].c_str(),&var6);
  chain->SetBranchAddress(varNames[6].c_str(),&var7);
  if(varNames.size() == 8){
    chain->SetBranchAddress(varNames[7].c_str(),&var8);
  }

  // Count entries in each bin
  int nEntries[nBins][nBins][nBins][nBins][nBins][nBins];
  for(int a = 0; a < nBins; ++a){
    for(int b = 0; b < nBins; ++b){
      for(int c = 0; c < nBins; ++c){
        for(int d = 0; d < nBins; ++d){
          for(int e = 0; e < nBins; ++e){
            for(int f = 0; f < nBins; ++f){
              nEntries[a][b][c][d][e][f] = 0;
            }
          }
        }
      }
    }
  }

  // Loop
  for(int i = 0; i < chain->GetEntries(); ++i){

    chain->GetEntry(i);

    // If we are using new files, need to worry about pmt location
    if(varNames.size() == 8){
      if(var8 != location) continue;
    }

    int bin0 = GetBin(var1,binArray[0]);
    int bin1 = GetBin(var2,binArray[1]);
    int bin2 = GetBin(var3,binArray[2]);
    int bin3 = GetBin(var4,binArray[3]);
    int bin4 = GetBin(var5,binArray[4]);
    int bin5 = GetBin(var6,binArray[5]);

//    std::cout << "(" << bin0 << "," << bin1 << "," << bin2 << "," << bin3 << "," << bin4 << "," << bin5 << ") ";
//    std::cout << "(" << vtxR << "," << vtxZ << "," << pmtZ << "," << cosZeta << "," << cosTheta << "," << phi << ") " << std::endl;

    qArray[bin0][bin1][bin2][bin3][bin4][bin5] += var7;
    ++nEntries[bin0][bin1][bin2][bin3][bin4][bin5];
    if(add4D){
      q4DArray[bin0][bin1][bin2][bin3] += var7;
    } 
  }

  for(int a = 0; a < nBins; ++a){
    for(int b = 0; b < nBins; ++b){
      for(int c = 0; c < nBins; ++c){
        for(int d = 0; d < nBins; ++d){
          for(int e = 0; e < nBins; ++e){
            for(int f = 0; f < nBins; ++f){
              int entries = nEntries[a][b][c][d][e][f];
              if(entries < 1000){
//                std::cout << "(" << a << "," << b << "," << c << "," << d << "," << e << "," << f << ") " << entries << std::endl;
              }
            }
          }
        }
      }
    }
  }

  // Make sure we get rid of the attached variables before they go out of scope
  chain->ResetBranchAddresses();

}

int GetBin(double val, std::vector<float> binArray){

  int bin = -1;

  for(unsigned int b = 1; b < binArray.size(); ++b){
    if(val < binArray[b]){
      bin = b - 1;
      break;
    }
  }

  return bin;

}

void AddArrayToVector(std::vector<std::vector<float> > &binVector, float *array){
  // Remember we need nBins+1 values to define nBins.
  std::vector<float> temp;
  for(int i = 0; i <= nBins; ++i){
    temp.push_back(array[i]);
  }
  binVector.push_back(temp);
}

void CalculateBinVector(float min, float max, std::vector<float> &binVec){

  // Start with a fresh vector.
  binVec.clear();
  for(int i = 0; i <= nBins; ++i){
    binVec.push_back(min + (i*(max-min))/(float)nBins);
    std::cout << binVec[i] << std::endl;
  }

}

void GetVarNames(std::vector<std::string> &vecName, bool type){

  if(type){
    vecName.push_back("vtxR");
    vecName.push_back("vtxZ");
    vecName.push_back("v2pZeta");
    vecName.push_back("pmtPos");
    vecName.push_back("dirTheta");
    vecName.push_back("dirPhi");
    vecName.push_back("charge");
    vecName.push_back("pmtLoc");
  }

}


