all : testMcbFit CombMcbFit libFit.so select

LIBBOOST = -lboost_system -lboost_program_options -lboost_filesystem
LIBROOT = `root-config --libs` -lMinuit   -lRooFit -lMinuit -lRooFitCore
LIBS = $(LIBROOT) $(LIBBOOST)

ROOTFLAGS = `root-config --cflags`
CXXFLAGS =   -std=c++1y -fPIC $(ROOTFLAGS)

BINDIR=$(HOME)/work/bin

testMcbFit : testMcbFit.cpp RooMcbPdf.h   RooMcbPdf.o
		$(CXX)  -o $@   $< RooMcbPdf.o  $(CXXFLAGS) $(LIBS) 

CombMcbFit : CombMcbFit.cpp RooMcbPdf.h   RooMcbPdf.o
		$(CXX)  -o $@   $< RooMcbPdf.o  $(CXXFLAGS) $(LIBS)

ROOT_ADAPTORS = RootEvent.o RootMC.o RootMdc.o RootEmc.o RootTof.o RootDedx.o RootMuc.o RootNeutral.o mctopo/RootMCTopo.o

libFit.so : RooMcbPdf.o libFitDict.o  $(ROOT_ADAPTORS)
		$(CXX) $^  -shared -o $@ -fPIC   $(LIBS)

RooMcbPdf.o : RooMcbPdf.cpp
		$(CXX) -fPIC -o RooMcbPdf.o $(CXXFLAGS) -c RooMcbPdf.cpp $(LIBS)


select : select.cpp $(ROOT_ADAPTORS) libFitDict.o RooMcbPdf.o
		$(CXX) $< -o $(BINDIR)/kksel  $(CXXFLAGS) $(ROOT_ADAPTORS)      $(LIBS)
	
libFitDict.o :  libFitDict.cxx
		$(CXX) -o $@  $(CXXFLAGS) -c $^

libFitDict.cxx : libFit.h libFitLinkDef.h
		rootcint -f $@ -c $^

.C.o :
		g++ -c $< -o $@ $(CXXFLAGS) $(LIBS)

clean :
			rm -f *.o *.so *~ *.pcm testMcbFit CombMcbFit libFitDict.o libFitDict.cxx libFit.so  $(ROOT_ADAPTORS)
