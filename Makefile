CC=g++
CFLAGS=-Wall -Werror -fopenmp -lm -std=c++11
DEBUGFLAGS=-g
RELFLAGS=-O3 -mavx2

HEADERS = src/block.h src/matrix.h src/Structures.h src/Tools.h
SOURCES = src/block.cpp src/matrix.cpp src/Structures.cpp src/Tools.cpp src/FiniteVolume.cpp src/main.cpp 

debug: $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) $(CFLAGS) $(DEBUGFLAGS) -o cudaCFD.exe
	
release: $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) $(CFLAGS) $(RELFLAGS) -o cudaCFD.exe
	
clean:
	rm cudaCFD.exe -rf
