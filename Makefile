CXX             = g++
#CXXFLAGS        = -fPIC -Iheader
CXXFLAGS        = -g -Iheader
LDFLAGS        = -L. -lm -lRNA

all: strDiversity

strDiversity: measures.o misc.o  structuralDiversity.o
	${CXX} ${CXXFLAGS} -o strDiversity structuralDiversity.o measures.o misc.o ${LDFLAGS}
	
structuralDiversity.o: structuralDiversity.cpp 
	${CXX} ${CXXFLAGS} -c structuralDiversity.cpp
	
measures.o: measures.cpp
	${CXX} ${CXXFLAGS} -c measures.cpp 

misc.o : misc.cpp
	${CXX} ${CXXFLAGS} -c misc.cpp

clean:
	rm -fr *.o strDiversity 	
