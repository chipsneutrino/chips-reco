# --------------------------------------------------------------
# GNUmakefile for WCSimAnalysis
# --------------------------------------------------------------


CXX			  = g++
CXXDEPEND     = -MM
CXXFLAGS      = -g -Wall -fPIC -O3 -std=c++11
LD            = g++ 
LDFLAGS       = -g -O3

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBS     := $(shell root-config --libs) -lSpectrum
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

INCDIR = ./include
SRCDIR = ./src
TMPDIR = ./tmp
LIBDIR = ./lib

INCLUDES = -I$(INCDIR)

WCSIM_INCDIR = ${WCSIMHOME}/include
WCSIM_LIBDIR = ${WCSIMHOME}/geant4/tmp/Linux-g++/WCSim

WCSIM_INCLUDES = -I$(WCSIM_INCDIR)
WCSIM_LDFLAGS += -L$(WCSIM_LIBDIR)
WCSIM_LIBS += -lWCSim

.PHONY: 
	all

all: rootcint shared libWCSimAnalysis.a

ROOTSO := libWCSimAnalysis.so

ROOTDICT := $(SRCDIR)/WCSimAnalysisRootDict.cc

# This is the list of all the ROOT-based classes we need to worry about.
# Assumes that each class has src/*.cc, include/*.hh and tmp/*.o files.
ROOTCLASS := WCSimDigitizerPDFMaker WCSimPID WCSimPIDTree WCSimOutputTree WCSimHitComparison WCSimTrackParameterEnums WCSimTimeLikelihood3 WCSimFastMath WCSimTransmissionFunctionLookup WCSimPiZeroSeed WCSimPiZeroSeedGenerator WCSimPiZeroSeeder WCSimPiZeroHoughSeeder WCSimPiZeroSingleElectronSeeder WCSimPiZeroElectronAdjuster WCSimPiZeroFitter WCSimDetectorParameters WCSimTimePredictor WCSimLikelihoodTrackFactory WCSimLikelihoodTrackBase WCSimFitterTrackParMap WCSimIntegralLookupMaker3D WCSimIntegralLookup3D WCSimIntegralLookupReader WCSimRecoEvDisplay WCSimRecoSummary WCSimEmissionProfileManager WCSimEmissionProfiles WCSimFitterConfig WCSimFitterInterface WCSimMapper WCSimFitterParameters WCSimFitterPlots WCSimDigitizerLikelihood WCSimTotalLikelihood WCSimLikelihoodTrack WCSimLikelihoodDigit WCSimLikelihoodDigitArray WCSimLikelihoodTuner WCSimLikelihoodFitter WCSimGeometry WCSimInterface WCSimParameters WCSimRecoObjectTable WCSimRecoDigit WCSimRecoCluster WCSimRecoClusterDigit WCSimRecoRing WCSimRecoVertex WCSimRecoEvent WCSimTrueEvent WCSimTrueTrack WCSimHoughTransform WCSimHoughTransformArray WCSimDataCleaner WCSimVertexFinder WCSimVertexGeometry WCSimRingFinder WCSimMsg WCSimChargePredictor WCSimRecoSeed WCSimRecoSlicer WCSimScatteringTableManager WCSimCosmicSeed WCSimCosmicFitter WCSimRecoClusteringUtil

# Create the ROOTINC list from the class list, remembering to also add the LinkDef file
ROOTINC = $(ROOTCLASS:%=$(INCDIR)/%.hh)
ROOTINC += $(INCDIR)/WCSimAnalysisRootLinkDef.hh
# For dictionary generation, want a version without the include directory.
ROOTINCNODIR = $(ROOTCLASS:%=%.hh)
ROOTINCNODIR += WCSimAnalysisRootLinkDef.hh
# Now for the ROOTSRC list
ROOTSRC = $(ROOTCLASS:%=$(SRCDIR)/%.cc)
# Finally, the ROOTOBJ list, remembering to add on the RootDict file
ROOTOBJS = $(ROOTCLASS:%=$(TMPDIR)/%.o)
ROOTOBJS += $(TMPDIR)/WCSimAnalysisRootDict.o

ROOTEXTOBJS := $(WCSIM_LIBDIR)/WCSimRootEvent.o $(WCSIM_LIBDIR)/WCSimRootGeom.o $(WCSIM_LIBDIR)/WCSimRootDict.o

$(TMPDIR)/%.o : $(SRCDIR)/%.cc
	@echo "<**Compiling $@**>"
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) -c $< -o $@ 

$(TMPDIR)/%.d: $(SRCDIR)/%.cc
	@echo "<**Depend $@**>"
	@mkdir -p $(@D)
	set -e; $(CXX) $(CXXDEPEND) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) $< \
	| sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

$(ROOTDICT) : $(ROOTSRC) $(ROOTINC)
	rootcint -f $(ROOTDICT) -c -I$(INCDIR) -I$(WCSIM_INCDIR) -I$(shell root-config --incdir) $(ROOTINCNODIR)

rootcint: $(ROOTDICT)

shared: $(ROOTDICT) $(ROOTSRC) $(ROOTINC) $(ROOTOBJS)
	@mkdir -p lib
	$(CXX) -shared $(ROOTLIBS) $(ROOTGLIBS) -O $(ROOTOBJS) -o $(LIBDIR)/$(ROOTSO)

libWCSimAnalysis.a : $(ROOTOBJS)
	$(RM) $@
	ar clq $@ $(ROOTOBJS)
	mv $@ $(LIBDIR)
	
evDisp:	
	$(CXX) `root-config --cflags --glibs --libs --evelibs` -I./include ${WCSIM_INCLUDES} -L${WCSIMANAHOME} -L${WCSIMHOME} -o evDisplay evDisplay.cc ${WCSIMHOME}/src/WCSimRootDict.cc ${WCSIMANAHOME}/src/WCSimAnalysisRootDict.cc -lWCSim -lWCSimAnalysis -lEG -lSpectrum -lMinuit -lTMVA -lMLP -lXMLIO

clean :
	@echo "<**Clean**>"
	rm -f $(SRCDIR)/*~ $(INCDIR)/*~ $(TMPDIR)/*.o $(TMPDIR)/*.d $(TMPDIR)/*.a $(SRCDIR)/WCSimAnalysisRootDict.*

DEPS = $(ROOTOBJS:$(TMPDIR)/%.o=$(TMPDIR)/%.d)

ifeq ($(MAKECMDGOALS),all)
 include $(DEPS)
endif

ifeq ($(MAKECMDGOALS),shared)
 include $(DEPS)
endif
