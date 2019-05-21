CC=g++

CPPFLAGS = -g -Wall -O3

SLR: main.o contig.o aligningFromBam.o scaffoldgraph.o scaffolding.o
	$(CC) -o $@ $^ ./lp/liblpsolve55.a -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lm -ldl -lz -no-pie
	
SORT-Contig: contig.o sortContigSet.o 
	$(CC) -o $@ $^ -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lz
	
SLR-unique-ambiguous: main-unique.o contig.o aligningFromBam.o aligningFromBamOfScaffold.o scaffoldgraph.o scaffolding.o
	$(CC) -o $@ $^ ./lp/liblpsolve55.a -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lm -ldl -lz -no-pie

main-unique.o: main-unique.cpp contig.h aligningFromBam.h scaffoldgraph.h scaffolding.h aligningFromBamOfScaffold.h 
	$(CC) -c main-unique.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

main.o: main.cpp contig.h aligningFromBam.h scaffoldgraph.h scaffolding.h
	$(CC) -c main.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

contig.o: contig.cpp contig.h
	$(CC) -c contig.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
sortContigSet.o: sortContigSet.cpp contig.h
	$(CC) -c sortContigSet.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
aligningFromBam.o: aligningFromBam.cpp contig.h aligningFromBam.h
	$(CC) -c aligningFromBam.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

scaffoldgraph.o: scaffoldgraph.cpp contig.h scaffoldgraph.h
	$(CC) -c scaffoldgraph.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
scaffolding.o: scaffolding.cpp contig.h scaffoldgraph.h scaffolding.h
	$(CC) -c scaffolding.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
aligningFromBamOfScaffold.o: aligningFromBamOfScaffold.cpp contig.h scaffolding.h aligningFromBamOfScaffold.h
	$(CC) -c aligningFromBamOfScaffold.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	

all: SLR SORT-Contig SLR-unique-ambiguous
	rm -f *.o

clean:
	rm -f *.o
	rm SLR
	rm SORT-Contig
	rm SLR-unique-ambiguous
