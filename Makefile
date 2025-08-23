CC=g++
CFLAGS=-Wall -Werror -fopenmp -lm -std=c++11
DEBUGFLAGS=-g
RELFLAGS=-O3 -mavx2

HEADERS = EulerSolve/block.h EulerSolve/matrix.h EulerSolve/Structures.h EulerSolve/Tools.h
SOURCES = EulerSolve/block.cpp EulerSolve/matrix.cpp EulerSolve/Structures.cpp EulerSolve/Tools.cpp EulerSolve/FiniteVolume.cpp CudaCFD/CudaCFD/main.cpp 

debug: $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) $(CFLAGS) $(DEBUGFLAGS) -o cudaCFD.exe
	
release: $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) $(CFLAGS) $(RELFLAGS) -o cudaCFD.exe
