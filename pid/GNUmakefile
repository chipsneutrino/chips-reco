PROGNAME      = wc_make_all_TMVA
SOURCEFILES   = wc_make_all_TMVA.cc wc_trainTMVA_electronMuonCCQE.cc wc_trainTMVA_electronCCQEvsNC.cc
INCLFILES     = wc_trainTMVA_electronMuonCCQE.hh wc_trainTMVA_electronCCQEvsNC.hh
INCLDIR       = .
INCLUDES      = $(INCLDIR)/$(INCLFILES)
CFLAGS       += -I$(INCLDIR)
OBJS          = $(patsubst %.cc, %.o, $(SOURCEFILES))
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
LDFLAGS       = -O
LIBS         += $(ROOTLIBS)
CFLAGS       += $(ROOTCFLAGS)

#  Not sure why Minuit, TMV aren't being included -- put in by hand
#
LIBS         += -lMinuit
LIBS         += -lTMVA
#LIBS          += -L/nfs/scratch0/cowan/TMVA/lib -lTMVA

%.o: %.cc 
	g++ ${CFLAGS} -c  -g -o $@ $<

$(PROGNAME):    $(OBJS) $(INCLUDES)
	g++ -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

test:
	@echo $(CFLAGS)

clean:	
	-rm -f ${PROGNAME} ${OBJS}
