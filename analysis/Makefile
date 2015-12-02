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

libFit.so : RooMcbPdf.o libFitDict.o  Selector.C RootEvent.C RootMC.C
		g++ -fPIC -shared  RooMcbPdf.o -o libFit.so libFitDict.o $(LIBS)

RooMcbPdf.o : RooMcbPdf.cpp
		g++ -fPIC -o RooMcbPdf.o $(CXXFLAGS) -c RooMcbPdf.cpp $(LIBS)


select : select.cpp RootEvent.o Selector.o 
		g++ -std=c++1y Selector.o RootEvent.o RootMC.o libFitDict.o RooMcbPdf.o -o select select.cpp $(CXXFLAGS) $(LIBS)
	
libFitDict.o :  libFitDict.cxx
		g++ -o libFitDict.o $(CXXFLAGS) -c libFitDict.cxx

libFitDict.cxx : libFit.h libFitLinkDef.h
		rootcint -f libFitDict.cxx -c libFit.h libFitLinkDef.h




RootEvent.o : 
	  g++ -o RootEvent.o $(CXXFLAGS) $(LIBS) -c RootEvent.C

Selector.o : 
	  g++ -o Selector.o $(CXXFLAGS) $(LIBS)  -c Selector.C

RootMC.o : 
	  g++ -o RootMC.o $(CXXFLAGS) $(LIBS)  -c RootMC.C




clean :
			rm -f *.o *.so *~ *.pcm testMcbFit CombMcbFit libFitDict.o libFitDict.cxx libFit.so
