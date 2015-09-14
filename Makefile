CXX             = g++
CXXFLAGS        = -fPIC
LDFLAGS        = -L. -lm -lRNA

all: strDiversity

strDiversity: structuralDiversity.o measures.o misc.o
	${CXX} ${CXXFLAGS} -o strDiversity structuralDiversity.o measures.o misc.o ${LDFLAGS}
	
structuralDiversity.o: structuralDiversity.cpp 
	${CXX} ${CXXFLAGS} -c structuralDiversity.cpp
	
measures.o: measures.cpp
	${CXX} ${CXXFLAGS} -c measures.cpp 

misc.o : misc.cpp
	${CXX} ${CXXFLAGS} -c misc.cpp

clean:
	rm -fr *.o strDiversity 	
