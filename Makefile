#!/bin/make
# GNUmakefile using implicit rules and standard definitions
SHELL=/bin/sh

# Uncomment the following line for root functionality
# USEROOT = 1

# Uncomment this line for a more verbose scan
# CXXFLAGS += -DVERBOSE

# Undefine to make a "online" version
ONLINE = 1 

# Define to use Gamma-Gamma gates in GeProcessor
# This turns on Gamma-Gamma angular distribution
# and Gamma-Gamma-Gamma gates
#GGATES = 1

# Use gfortran
HHIRF_GFORTRAN = 1

#These will set the analysis used on the waveforms
#Uncomment this line to use the Pulse Fitting routine
PULSEFIT = 1
#Uncomment this line to use the cfd
#DCFD = 1


#------- instruct make to search through these
#------- directories to find files
vpath %.f scan/ 
vpath %.hpp include/
vpath %.h include/
vpath %.icc include/
vpath %.cpp src/
vpath %.o obj/

#DIRA2=/usr/hhirf-intel64  #//at JINR
#DIRB= /usr/acq2/lib
HHIRF_DIR = /usr/hhirf-intel64
#ACQ2_LIBDIR = /usr/hhirf-intel64
ACQ2_LIBDIR = /usr/acq2n/lib

LIBS = $(HHIRF_DIR)/scanorlib.a $(HHIRF_DIR)/orphlib.a\
       $(ACQ2_LIBDIR)/acqlib.a  $(ACQ2_LIBDIR)/ipclib.a

OutPutOpt     = -o # keep whitespace after "-o"
ObjSuf        = o

#------- define file suffixes
fSrcSuf   = f
cSrcSuf   = c
c++SrcSuf = cpp
cxxSrcSuf = cxx

#------- define compilers
#define to compile with gfortran (>=4.2) if required for the hhirf libs
ifeq ($(HHIRF_GFORTRAN), )
FC        = g77
else
FC        = gfortran
endif

GCC       = gcc 
CXX       = g++
LINK.o    = $(FC) $(LDFLAGS)

# -Dnewreadout is needed to account for a change to pixie16 readout
# structure change on 03/20/08.  Remove for backwards compatability
#
# for debug and profiling add options -g -pg
# and remove -O
#------- define basic compiler flags, no warnings on code that is not my own
FFLAGS   += -O3
GCCFLAGS += -fPIC $(CINCLUDEDIRS) -Dnewreadout
#CXXFLAGS += -Wall -fPIC $(CINCLUDEDIRS) -Dnewreadout 
CXXFLAGS += -fPIC $(CINCLUDEDIRS) -Dnewreadout

ifdef ONLINE
CXXFLAGS += -DONLINE
endif

#------- include directories for the pixie c files
CINCLUDEDIRS  = -Iinclude

#------- basic linking instructions
LDLIBS   += -lm -lstdc++
LDLIBS   += -lgsl -lgslcblas
CXXFLAGS += -Dpulsefit
CXXFLAGS += -Ddcfd
#CXXFLAGS += -Dpixel

ifeq ($(FC),gfortran)
FFLAGS	 += -fsecond-underscore
LDLIBS	 += -lgfortran
GCCFLAGS += -O3
CXXFLAGS += -O3 -DLINK_GFORTRAN
else
LDFLAGS += -g77libs
LDLIBS	+= -lg2c
endif

#-------- define file variables -----------------------

# objects from fortran
SET2CCO          = set2cc.$(ObjSuf)
MESSLOGO         = messlog.$(ObjSuf)
MILDATIMO        = mildatim.$(ObjSuf)
SCANORO          = scanor.$(ObjSuf)

#XML parser
PUGIXMLO = pugixml.$(ObjSuf)

# objects from cpp
PIXIEO           = PixieStd.$(ObjSuf)

# ReadBufData
READBUFFDATADFO    = ReadBuffData.RevD.$(ObjSuf)
READBUFFDATAAO    = ReadBuffData.RevA.$(ObjSuf)

BEAMLOGICPROCESSORO  = BeamLogicProcessor.$(ObjSuf)
BETASCINTPROCESSORO = BetaScintProcessor.$(ObjSuf)
BETA4HEN3PROCESSORO = Beta4Hen3Processor.$(ObjSuf)
CALIBRATORO      = Calibrator.$(ObjSuf)
CFDANALYZERO     = CfdAnalyzer.$(ObjSuf)
CHANEVENTO       = ChanEvent.$(ObjSuf)
CHANIDENTIFIERO  = ChanIdentifier.$(ObjSuf)
CORRELATORO      = Correlator.$(ObjSuf)
DETECTORDRIVERO  = DetectorDriver.$(ObjSuf)
DETECTORLIBRARYO = DetectorLibrary.$(ObjSuf)
DETECTORSUMMARYO = DetectorSummary.$(ObjSuf)
DOUBLETRACEO     = DoubleTraceAnalyzer.$(ObjSuf)
DSSDPROCESSORO   = DssdProcessor.$(ObjSuf)
DSSD4SHEPROCESSORO = Dssd4SHEProcessor.$(ObjSuf)
DSSD4JAEAPROCESSORO = Dssd4JAEAProcessor.$(ObjSuf)
NAIPROCESSORO    = NaIProcessor.$(ObjSuf)
PINPROCESSORO    = PINProcessor.$(ObjSuf)
EVENTPROCESSORO  = EventProcessor.$(ObjSuf)
FITTINGANALYZERO = FittingAnalyzer.$(ObjSuf)
GEPROCESSORO     = GeProcessor.$(ObjSuf)
GE4HEN3PROCESSORO= Ge4Hen3Processor.$(ObjSuf)
GECALIBPROCESSORO= GeCalibProcessor.$(ObjSuf)
GLOBALSO         = Globals.$(ObjSuf)
HEN3PROCESSORO   = Hen3Processor.$(ObjSuf)
ISSDPROCESSORO   = ImplantSsdProcessor.$(ObjSuf)
INITIALIZEO      = Initialize.$(ObjSuf)
IONCHAMBERPROCESSORO = IonChamberProcessor.$(ObjSuf)
LIQUIDSCINTPROCESSORO = LiquidScintProcessor.$(ObjSuf)
LOGICPROCESSORO  = LogicProcessor.$(ObjSuf)
MESSENGERO       = Messenger.$(ObjSuf)
MCPPROCESSORO    = McpProcessor.$(ObjSuf)
MTCPROCESSORO    = MtcProcessor.$(ObjSuf)
NEUTRONSCINTPROCESSORO  = NeutronScintProcessor.$(ObjSuf)
NOTEBOOKO		 = Notebook.$(ObjSuf)
PLOTSO           = Plots.$(ObjSuf)
PLOTSREGISTERO   = PlotsRegister.$(ObjSuf)
POSITIONPROCESSORO = PositionProcessor.$(ObjSuf)
RANDOMPOOLO      = RandomPool.$(ObjSuf)
RAWEVENTO        = RawEvent.$(ObjSuf)
SHECORRELATORO   = SheCorrelator.$(ObjSuf)
JAEACORRELATORO  = JAEACorrelator.$(ObjSuf)
ROOTPROCESSORO   = RootProcessor.$(ObjSuf)
PLACEBUILDERO    = PlaceBuilder.$(ObjSuf)
PLACESO          = Places.$(ObjSuf)
PULSERPROCESSORO = PulserProcessor.$(ObjSuf)
SSDPROCESSORO    = SsdProcessor.$(ObjSuf)
STATSDATAO       = StatsData.$(ObjSuf)
TAUANALYZERO     = TauAnalyzer.$(ObjSuf)
TIMINGINFOO      = TimingInformation.$(ObjSuf)
TRIGGERLOGICPROCESSORO = TriggerLogicProcessor.$(ObjSuf)
TRACEO           = Trace.$(ObjSuf)
TRACEEXTRACTERO  = TraceExtracter.$(ObjSuf)
TRACEFILTERO     = TraceFilterer.$(ObjSuf)
TRACESUBO        = TraceAnalyzer.$(ObjSuf)
TREECORRELATORO  = TreeCorrelator.$(ObjSuf)
VANDLEPROCESSORO = VandleProcessor.$(ObjSuf)
VANDLEROOTO      = VandleROOT.$(ObjSuf)
WALKCORRECTORO   = WalkCorrector.$(ObjSuf)
WAVEFORMSUBO     = WaveformAnalyzer.$(ObjSuf)
WAVEFORMSUBO     = WaveformAnalyzer.$(ObjSuf)

ifdef USEROOT
PIXIE = pixie_ldf_c_root$(ExeSuf)
else
ifdef ONLINE
PIXIE = pixie_ldf_c_online$(ExeSuf)
else
PIXIE = pixie_ldf_c$(ExeSuf)
endif
endif

#----- list of objects
# Fortran objects
FORT_OBJS   = \
$(SET2CCO)\
$(MESSLOGO)\
$(MILDATIMO)\
$(SCANORO)
# Important to compile READBUFFDATA first
CXX_OBJS += $(READBUFFDATAAO)
CXX_OBJS += $(READBUFFDATADFO)
# other C++ objects
CXX_OBJS += \
$(PUGIXMLO)\
$(PIXIEO)\
$(BEAMLOGICPROCESSORO)\
$(BETASCINTPROCESSORO)\
$(BETA4HEN3PROCESSORO)\
$(CALIBRATORO)\
$(CORRELATORO)\
$(CHANEVENTO)\
$(CHANIDENTIFIERO)\
$(HISTOGRAMMERO)\
$(DETECTORDRIVERO)\
$(DETECTORLIBRARYO)\
$(DETECTORSUMMARYO)\
$(DOUBLETRACEO)\
$(DSSDPROCESSORO)\
$(DSSD4SHEPROCESSORO)\
$(DSSD4JAEAPROCESSORO)\
$(NAIPROCESSORO)\
$(PINPROCESSORO)\
$(EVENTPROCESSORO)\
$(GEPROCESSORO)\
$(GE4HEN3PROCESSORO)\
$(GECALIBPROCESSORO)\
$(GLOBALSO)\
$(HEN3PROCESSORO)\
$(ISSDPROCESSORO)\
$(INITIALIZEO)\
$(IONCHAMBERPROCESSORO)\
$(LIQUIDSCINTPROCESSORO)\
$(LOGICPROCESSORO)\
$(MESSENGERO)\
$(MCPPROCESSORO)\
$(MTCPROCESSORO)\
$(NEUTRONSCINTPROCESSORO)\
$(NOTEBOOKO)\
$(PLOTSO)\
$(PLOTSREGISTERO)\
$(POSITIONPROCESSORO)\
$(RANDOMPOOLO)\
$(RAWEVENTO)\
$(SHECORRELATORO)\
$(JAEACORRELATORO)\
$(PLACEBUILDERO)\
$(PLACESO)\
$(PULSERPROCESSORO)\
$(ACCUMULATORO)\
$(SSDPROCESSORO)\
$(STATSDATAO)\
$(TAUANALYZERO)\
$(TIMINGINFOO)\
$(TRIGGERLOGICPROCESSORO)\
$(TRACEO)\
$(TRACEEXTRACTERO)\
$(TRACEFILTERO)\
$(TRACESUBO)\
$(TREECORRELATORO)\
$(VANDLEPROCESSORO)\
$(WALKCORRECTORO)\
$(WAVEFORMSUBO)\
$(WAVEFORMSUBO) \
$(VANLDEPROCESSORO)\
$(FITTINGANALYZERO) \
$(CFDANALYZERO)  \

ifdef USEROOT
CXX_OBJS  += $(ROOTPROCESSORO) $(VANDLEROOTO) $(SCINTROOTO)
endif

#PROGRAMS = $(PIXIE)

#------------ adjust compilation if ROOT capability is desired -------
ifdef USEROOT
ROOTCONFIG   := root-config
#no uncomment ROOTCLFAGS   := $(filter-out pthread,$(ROOTCFLAGS))
CXXFLAGS     += $(shell $(ROOTCONFIG) --cflags) -Duseroot
LDFLAGS      += $(shell $(ROOTCONFIG) --ldflags)
LDLIBS       := $(shell $(ROOTCONFIG) --libs)
endif
#---------- Update some information about the object files 
FORT_OBJDIR = obj/fortran
FORT_OBJS_W_DIR = $(addprefix $(FORT_OBJDIR)/,$(FORT_OBJS))
CXX_OBJDIR = obj/c++
CXX_OBJS_W_DIR = $(addprefix $(CXX_OBJDIR)/,$(CXX_OBJS))


#------------ Compile with Gamma-Gamma gates support in GeProcessor
ifdef GGATES
CXXFLAGS	+= -DGGATES
endif

#--------- Add to list of known file suffixes
.SUFFIXES: .$(cxxSrcSuf) .$(fSrcSuf) .$(c++SrcSuf) .$(cSrcSuf)

.phony: all clean
all:     $(FORT_OBJS_W_DIR) $(CXX_OBJS_W_DIR) $(PIXIE)

$(FORT_OBJS_W_DIR): | $(FORT_OBJDIR)

$(FORT_OBJDIR):
	mkdir -p $(FORT_OBJDIR)

$(FORT_OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

$(CXX_OBJS_W_DIR): | $(CXX_OBJDIR)

$(CXX_OBJDIR):
	mkdir -p $(CXX_OBJDIR)

$(CXX_OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



#----------- link all created objects together
#----------- to create pixie_ldf_c program
$(PIXIE): $(FORT_OBJS_W_DIR) $(CXX_OBJS_W_DIR) $(LIBS)
	$(LINK.o) $^ -o $@ $(LDLIBS)
#----------- remove all objects, core and .so file
clean:
	@echo "Cleaning up..."
#	@rm -f $(CXX_OBJS_W_DIR) $(PIXIE) core *~ src/*~ include/*~ scan/*~ config/*~
	@rm -f $(FORT_OBJS_W_DIR) $(CXX_OBJS_W_DIR) $(PIXIE) core *~ src/*~ include/*~ scan/*~ config/*~

tidy:
	@echo "Tidying up..."
	@rm -f $(CXX_OBJS_W_DIR) core *~ src/*~ include/*~ scan/*~


