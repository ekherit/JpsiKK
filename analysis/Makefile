all : testMcbFit CombMcbFit libFit.so select

#LIBS = `root-config --libs` -lMinuit  -lboost_program_options -lRooFit -lMinuit -lRooFitCore -lgsl -lgslcblas

LIBBOOST = -lboost_system -lboost_program_options -lboost_filesystem
LIBROOT = `root-config --libs` -lMinuit   -lRooFit -lMinuit -lRooFitCore
LIBS = $(LIBROOT) $(LIBBOOST)

ROOTFLAGS = `root-config --cflags`
CXXFLAGS =   -std=c++1y -fPIC $(ROOTFLAGS)

BINDIR=$(HOME)/work/bin

testMcbFit : testMcbFit.cpp RooMcbPdf.h   RooMcbPdf.o
		g++ -o testMcbFit   RooMcbPdf.o testMcbFit.cpp $(CXXFLAGS) $(LIBS) 

CombMcbFit : CombMcbFit.cpp RooMcbPdf.h   RooMcbPdf.o
		g++  -o CombMcbFit   RooMcbPdf.o CombMcbFit.cpp $(CXXFLAGS) $(LIBS)

libFit.so : RooMcbPdf.o libFitDict.o   RootEvent.C RootMC.C
		g++ -fPIC -shared  RooMcbPdf.o -o libFit.so libFitDict.o $(LIBS)

RooMcbPdf.o : RooMcbPdf.cpp
		g++ -fPIC -o RooMcbPdf.o $(CXXFLAGS) -c RooMcbPdf.cpp $(LIBS)

ROOT_ADAPTORS = RootEvent.o RootMC.o RootMdc.o RootEmc.o RootTof.o RootDedx.o RootMuc.o RootNeutral.o RootMCTopo.o

select : select.cpp $(ROOT_ADAPTORS)
		g++ -std=c++1y  $(ROOT_ADAPTORS) libFitDict.o RooMcbPdf.o -o select select.cpp $(CXXFLAGS) $(LIBS)
	
libFitDict.o :  libFitDict.cxx
		g++ -o libFitDict.o $(CXXFLAGS) -c libFitDict.cxx

libFitDict.cxx : libFit.h libFitLinkDef.h
		rootcint -f libFitDict.cxx -c libFit.h libFitLinkDef.h


RootEvent.o : 
	  g++ -o RootEvent.o $(CXXFLAGS) $(LIBS) -c RootEvent.C

RootMC.o : 
	  g++ -o RootMC.o $(CXXFLAGS) $(LIBS)  -c RootMC.C

RootMdc.o : 
	  g++ -o RootMdc.o $(CXXFLAGS) $(LIBS) -c RootMdc.C

RootEmc.o : 
	  g++ -o RootEmc.o $(CXXFLAGS) $(LIBS) -c RootEmc.C

RootTof.o : 
	  g++ -o RootTof.o $(CXXFLAGS) $(LIBS) -c RootTof.C

RootDedx.o : 
	  g++ -o RootDedx.o $(CXXFLAGS) $(LIBS) -c RootDedx.C

RootMCTopo.o : 
	  g++ -o RootMCTopo.o $(CXXFLAGS) $(LIBS) -c RootMCTopo.C

RootMuc.o : 
	  g++ -o RootMuc.o $(CXXFLAGS) $(LIBS) -c RootMuc.C

RootNeutral.o : 
	  g++ -o RootNeutral.o $(CXXFLAGS) $(LIBS) -c RootNeutral.C



clean :
			rm -f *.o *.so *~ *.pcm testMcbFit CombMcbFit libFitDict.o libFitDict.cxx libFit.so  $(ROOT_ADAPTORS)
