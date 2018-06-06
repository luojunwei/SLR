#ifndef CONTIG_H_INCLUDED 
#define CONTIG_H_INCLUDED 


#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

using namespace std;

typedef struct ArrayTwoElementList{
	int readPosition;
	int contigPosition;
	ArrayTwoElementList * next;
}ArrayTwoElementList;


typedef struct ReadAligning{
	int readIndex;
	int readStartPosition;
	int readEndPosition;
	int contigStartPosition;
	int contigEndPosition;
	ArrayTwoElementList * kmerAligning;
	ReadAligning * next;
}ReadAligning;

typedef struct Contig{
	char * contig;
	int contigLength;
	bool shortContig;
	int realContigIndex;
	ReadAligning * readAligning;
	Contig(){
		contig = NULL;
		contigLength = 0;
		shortContig = true;
		readAligning = NULL;
		realContigIndex = -1;
	}
}Contig;

typedef struct ContigSetHead{
	Contig * contigSet;
	long int contigCount;
	long int allContigLength;
	bool * repeatContigIndex;
	Contig * minContigSet;
	long int minContigCount;
	long int minAllContigLength;
	bool * visited;
	ContigSetHead(){
		contigSet = NULL;
		contigCount = 0;
		allContigLength = 0;
	}
	
}ContigSetHead;




ContigSetHead * GetContigSet(char * contigSetFile, int contigLengthThreshold);

void SortContigSet(char * contigSetFile, char * sortContigSetFile);

void GetContigIDandPosition(ContigSetHead * contigSetHead, int tempPosition, int * P);

void DeleteTailOfContigSet(ContigSetHead * contigSetHead, int tailLength);

void DeleteArrayTwoElementList(ArrayTwoElementList * arrayTwoElementList);

char * ReverseComplement(char * temp);







































#endif