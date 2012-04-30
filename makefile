#CPP = g++ -g -fopenmp
CPP = /usr/local/stow/openmpi-1.2.6/bin/mpic++ -g

all: haploi  

haploi: main.o GenoHaploDB.o MUTDP.o HapAssembler.o  
	$(CPP) -O main.o GenoHaploDB.o MUTDP.o HapAssembler.o -o Haploi 

GenoHaploDB.o: GenoHaploDB.cpp
	$(CPP) -c GenoHaploDB.cpp

MUTDP.o: MUTDP.cpp GenoHaploDB.cpp util.h
	$(CPP) -c MUTDP.cpp

HapAssembler.o: MUTDP.cpp GenoHaploDB.cpp util.h
	$(CPP) -c HapAssembler.cpp

main.o: main.cpp MUTDP.cpp GenoHaploDB.cpp HapAssembler.cpp 
	$(CPP) -c main.cpp

clean:
	rm -rf *.o Haploi 

