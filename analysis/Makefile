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

libFit.so : RooMcbPdf.o libFitDict.o   RootEvent.C RootMC.C
		$(CXX) -fPIC -shared -o $@ RooMcbPdf.o  libFitDict.o $(LIBS)

RooMcbPdf.o : RooMcbPdf.cpp
		$(CXX) -fPIC -o RooMcbPdf.o $(CXXFLAGS) -c RooMcbPdf.cpp $(LIBS)

ROOT_ADAPTORS = RootEvent.o RootMC.o RootMdc.o RootEmc.o RootTof.o RootDedx.o RootMuc.o RootNeutral.o mctopo/RootMCTopo.o

select : select.cpp $(ROOT_ADAPTORS)
		$(CXX) -std=c++1y  $(ROOT_ADAPTORS) libFitDict.o RooMcbPdf.o  -o $@ $< $(CXXFLAGS) $(LIBS)
	
libFitDict.o :  libFitDict.cxx
		$(CXX) -o $@  $(CXXFLAGS) -c $^

libFitDict.cxx : libFit.h libFitLinkDef.h
		rootcint -f $@ -c $^

.C.o :
		g++ -o $@ $(CXXFLAGS) $(LIBS) -c $<

#RootEvent.o : 
#	  g++ -o RootEvent.o $(CXXFLAGS) $(LIBS) -c RootEvent.C
#
#RootMC.o : 
#	  g++ -o RootMC.o $(CXXFLAGS) $(LIBS)  -c RootMC.C
#
#RootMdc.o : 
#	  g++ -o RootMdc.o $(CXXFLAGS) $(LIBS) -c RootMdc.C
#
#RootEmc.o : 
#	  g++ -o RootEmc.o $(CXXFLAGS) $(LIBS) -c RootEmc.C
#
#RootTof.o : 
#	  g++ -o RootTof.o $(CXXFLAGS) $(LIBS) -c RootTof.C
#
#RootDedx.o : 
#	  g++ -o RootDedx.o $(CXXFLAGS) $(LIBS) -c RootDedx.C
#
#RootMCTopo.o : 
#	  g++ -o RootMCTopo.o $(CXXFLAGS) $(LIBS) -c RootMCTopo.C
#
#RootMuc.o : 
#	  g++ -o RootMuc.o $(CXXFLAGS) $(LIBS) -c RootMuc.C
#
#RootNeutral.o : 
#	  g++ -o RootNeutral.o $(CXXFLAGS) $(LIBS) -c RootNeutral.C



clean :
			rm -f *.o *.so *~ *.pcm testMcbFit CombMcbFit libFitDict.o libFitDict.cxx libFit.so  $(ROOT_ADAPTORS)
