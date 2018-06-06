CC=g++

CPPFLAGS = -g -Wall -O3

SLR: main.o contig.o aligningFromBam.o scaffoldgraph.o scaffolding.o
	$(CC) -o $@ $^ ./lp/liblpsolve55.a -I $(BAMTOOLS_HOME)/include/bamtools/ $(BAMTOOLS_HOME)/lib/libbamtools.a -lm -ldl -lz
	
SORT-Contig: contig.o sortContigSet.o 
	$(CC) -o $@ $^
	
main.o: main.cpp contig.h aligningFromBam.h scaffoldgraph.h scaffolding.h
	$(CC) -c main.cpp -I $(BAMTOOLS_HOME)/include/bamtools/

contig.o: contig.cpp contig.h
	$(CC) -c contig.cpp
	
sortContigSet.o: sortContigSet.cpp contig.h
	$(CC) -c sortContigSet.cpp
	
aligningFromBam.o: aligningFromBam.cpp contig.h aligningFromBam.h
	$(CC) -c aligningFromBam.cpp -I $(BAMTOOLS_HOME)/include/bamtools/

scaffoldgraph.o: scaffoldgraph.cpp contig.h scaffoldgraph.h
	$(CC) -c scaffoldgraph.cpp -I $(BAMTOOLS_HOME)/include/bamtools/
	
scaffolding.o: scaffolding.cpp contig.h scaffoldgraph.h scaffolding.h
	$(CC) -c scaffolding.cpp -I $(BAMTOOLS_HOME)/include/bamtools/
	

all: SLR SORT-Contig
	rm -f *.o

clean:
	rm -f *.o
	rm SLR
	rm SORT-Contig
