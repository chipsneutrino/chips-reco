#---Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(chips-reco CXX)

#---Build Configuration
set(CMAKE_CXX_STANDARD 11)
set(BUILD_SHARED_LIBS YES CACHE BOOL "")
set(BUILD_STATIC_LIBS NO CACHE BOOL "")
set(CMAKE_VERBOSE_MAKEFILE off)

#---Default to a 'release' build with O2 optimisations
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

#---Find chips-sim
list(APPEND CMAKE_PREFIX_PATH $ENV{CHIPSSIM})
find_package(chips-sim REQUIRED)
if(NOT chips-sim_FOUND)
  message(STATUS "chips-sim package not found.")  
  return()
endif()

#---Find ROOT and define it's cmake macros
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT CONFIG REQUIRED COMPONENTS TMVA)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT package not found.")  
  return()
endif()
include(${ROOT_USE_FILE})

#---Locate sources and headers for this project
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${CHIPSSIM_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})

#---WCSimAnalysisRootDict.cxx and WCSimAnalysisRootDict.pcm generation by rootcling
ROOT_GENERATE_DICTIONARY(WCSimAnalysisRootDict 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimDigitizerPDFMaker.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPID.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPIDTree.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimOutputTree.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimHitComparison.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimTrackParameterEnums.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimTimeLikelihood3.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimFastMath.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimTransmissionFunctionLookup.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPiZeroSeed.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPiZeroSeedGenerator.hh
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPiZeroSeeder.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPiZeroHoughSeeder.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPiZeroSingleElectronSeeder.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPiZeroElectronAdjuster.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimPiZeroFitter.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimDetectorParameters.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimTimePredictor.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimLikelihoodTrackFactory.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimLikelihoodTrackBase.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimFitterTrackParMap.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimIntegralLookupMaker3D.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimIntegralLookup3D.hh
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimIntegralLookupReader.hh  
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoEvDisplay.hh
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoSummary.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimEmissionProfileManager.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimEmissionProfiles.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimFitterConfig.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimFitterInterface.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimMapper.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimFitterParameters.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimFitterPlots.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimDigitizerLikelihood.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimTotalLikelihood.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimLikelihoodTrack.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimLikelihoodDigit.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimLikelihoodDigitArray.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimLikelihoodTuner.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimLikelihoodFitter.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimGeometry.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimInterface.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimParameters.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoObjectTable.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoDigit.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoCluster.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoClusterDigit.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoRing.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoVertex.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoEvent.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimTrueEvent.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimTrueTrack.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimHoughTransform.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimHoughTransformArray.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimDataCleaner.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimVertexFinder.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimVertexGeometry.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRingFinder.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimMsg.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoSeed.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoSlicer.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimScatteringTableManager.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimCosmicSeed.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimCosmicFitter.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimRecoClusteringUtil.hh 
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimChargePredictor.hh
                         LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/WCSimAnalysisRootLinkDef.hh)

#---Make shared object library, crucial for reading ROOT classes
add_library(WCSimAnalysisRoot SHARED 
            ./src/base/WCSimDigitizerPDFMaker.cc 
            ./src/base/WCSimPID.cc 
            ./src/base/WCSimPIDTree.cc 
            ./src/base/WCSimOutputTree.cc 
            ./src/base/WCSimHitComparison.cc 
            ./src/base/WCSimTrackParameterEnums.cc 
            ./src/base/WCSimTimeLikelihood3.cc 
            ./src/base/WCSimFastMath.cc 
            ./src/base/WCSimTransmissionFunctionLookup.cc 
            ./src/base/WCSimPiZeroSeed.cc 
            ./src/base/WCSimPiZeroSeedGenerator.cc
            ./src/base/WCSimPiZeroSeeder.cc 
            ./src/base/WCSimPiZeroHoughSeeder.cc 
            ./src/base/WCSimPiZeroSingleElectronSeeder.cc 
            ./src/base/WCSimPiZeroElectronAdjuster.cc 
            ./src/base/WCSimPiZeroFitter.cc 
            ./src/base/WCSimDetectorParameters.cc 
            ./src/base/WCSimTimePredictor.cc 
            ./src/base/WCSimLikelihoodTrackFactory.cc 
            ./src/base/WCSimLikelihoodTrackBase.cc 
            ./src/base/WCSimFitterTrackParMap.cc 
            ./src/base/WCSimIntegralLookupMaker3D.cc 
            ./src/base/WCSimIntegralLookup3D.cc 
            ./src/base/WCSimIntegralLookupReader.cc
            ./src/base/WCSimRecoEvDisplay.cc
            ./src/base/WCSimRecoSummary.cc 
            ./src/base/WCSimEmissionProfileManager.cc 
            ./src/base/WCSimEmissionProfiles.cc 
            ./src/base/WCSimFitterConfig.cc 
            ./src/base/WCSimFitterInterface.cc 
            ./src/base/WCSimMapper.cc 
            ./src/base/WCSimFitterParameters.cc 
            ./src/base/WCSimFitterPlots.cc 
            ./src/base/WCSimDigitizerLikelihood.cc 
            ./src/base/WCSimTotalLikelihood.cc 
            ./src/base/WCSimLikelihoodTrack.cc 
            ./src/base/WCSimLikelihoodDigit.cc 
            ./src/base/WCSimLikelihoodDigitArray.cc 
            ./src/base/WCSimLikelihoodTuner.cc 
            ./src/base/WCSimLikelihoodFitter.cc 
            ./src/base/WCSimGeometry.cc 
            ./src/base/WCSimInterface.cc 
            ./src/base/WCSimParameters.cc 
            ./src/base/WCSimRecoObjectTable.cc 
            ./src/base/WCSimRecoDigit.cc 
            ./src/base/WCSimRecoCluster.cc 
            ./src/base/WCSimRecoClusterDigit.cc 
            ./src/base/WCSimRecoRing.cc 
            ./src/base/WCSimRecoVertex.cc 
            ./src/base/WCSimRecoEvent.cc 
            ./src/base/WCSimTrueEvent.cc 
            ./src/base/WCSimTrueTrack.cc 
            ./src/base/WCSimHoughTransform.cc 
            ./src/base/WCSimHoughTransformArray.cc 
            ./src/base/WCSimDataCleaner.cc 
            ./src/base/WCSimVertexFinder.cc 
            ./src/base/WCSimVertexGeometry.cc 
            ./src/base/WCSimRingFinder.cc 
            ./src/base/WCSimMsg.cc 
            ./src/base/WCSimRecoSeed.cc 
            ./src/base/WCSimRecoSlicer.cc 
            ./src/base/WCSimScatteringTableManager.cc 
            ./src/base/WCSimCosmicSeed.cc 
            ./src/base/WCSimCosmicFitter.cc 
            ./src/base/WCSimRecoClusteringUtil.cc 
            ./src/base/WCSimChargePredictor.cc
            WCSimAnalysisRootDict.cxx)
target_link_libraries(WCSimAnalysisRoot ${ROOT_LIBRARIES} Spectrum)

#---Add the geometryHelper executable and link to all the required libraries
add_executable(evDisplay ./src/apps/evDisplay.cc ${sources} ${headers})
target_link_libraries(evDisplay ${ROOT_LIBRARIES} ${CHIPSSIM_LIBRARY} Gui EG Spectrum WCSimAnalysisRoot Tree)

#---Download large data files to the config directory
message([<mode>] "message text" ...)

if(EXISTS "./config/electronIntegralsSmall.root")
  message(STATUS "Already have electronIntegralsSmall.root")
else()
  message(STATUS "Downloading electronIntegralsSmall.root")
  file(DOWNLOAD http://www.hep.ucl.ac.uk/~jtingey/chips-reco/electronIntegralsSmall.root ./config/electronIntegralsSmall.root)
endif()

if(EXISTS "./config/emissionProfilesElectrons.root")
  message(STATUS "Already have emissionProfilesElectrons.root")
else()
  message(STATUS "Downloading emissionProfilesElectrons.root")
  file(DOWNLOAD http://www.hep.ucl.ac.uk/~jtingey/chips-reco/emissionProfilesElectrons.root ./config/emissionProfilesElectrons.root)
endif()

if(EXISTS "./config/emissionProfilesMuons.root")
  message(STATUS "Already have emissionProfilesMuons.root")
else()
  message(STATUS "Downloading emissionProfilesMuons.root")
  file(DOWNLOAD http://www.hep.ucl.ac.uk/~jtingey/chips-reco/emissionProfilesMuons.root ./config/emissionProfilesMuons.root)
endif()

if(EXISTS "./config/muonIntegralsSmall.root")
  message(STATUS "Already have muonIntegralsSmall.root")
else()
  message(STATUS "Downloading muonIntegralsSmall.root")
  file(DOWNLOAD http://www.hep.ucl.ac.uk/~jtingey/chips-reco/muonIntegralsSmall.root ./config/muonIntegralsSmall.root)
endif()