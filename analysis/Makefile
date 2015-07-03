all : testMcbFit

LIBS = `root-config --libs` -lMinuit  -lboost_program_options -lRooFit -lMinuit -lRooFitCore
CXXFLAGS = `root-config --cflags`  -std=c++1y -fPIC

BINDIR=$(HOME)/work/bin

testMcbFit : testMcbFit.cpp RooMcbPdf.h   RooMcbPdf.o
		g++ $(LIBS) -o testMcbFit   RooMcbPdf.o testMcbFit.cpp $(CXXFLAGS)

libFit.so : RooMcbPdf.o libFitDict.o 
		g++ -fPIC -shared  RooMcbPdf.o -o libFit.so libFitDict.o $(LIBS)

RooMcbPdf.o : RooMcbPdf.cpp
		g++ -fPIC -o RooMcbPdf.o $(CXXFLAGS) -c RooMcbPdf.cpp

libFit.cxx : libFit.h libFitLinkDef.h
		rootcint -f libFitDict.cxx -c libFit.h libFitLinkDef.h
	
libFitDict.o :  libFitDict.cxx
		g++ -o libFitDict.o $(CXXFLAGS) -c libFitDict.cxx




clean :
			rm *.o *.so *~ *.pcm
