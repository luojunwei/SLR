#ifndef READ_H_INCLUDED 
#define READ_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "contig.h"
#include "kmer.h"

using namespace std;

typedef struct Read{
	char * read;
	int readLength;
	int allocateReadLength;
}Read;

typedef struct ReadSetHead{
	Read * readSet;
	int readCount;
	int realReadCount;
	int kmerLength;
	int step;
}ReadSetHead;


void GetReadToKmer(ContigSetHead * contigSetHead, KmerHashTableHead * kmerHashTableHead, Read * read, int kmerLength, int step, char * kmer, char * tempKmer, FILE * fp, int readStart);

void GetReadSetToKmer(Read * readSet, int readCount, ContigSetHead * contigSetHead, KmerHashTableHead * kmerHashTableHead, int kmerLength, int step, int readLengthCutOff, char * kmer, char * tempKmer, int readStart, char * file, char * line);

void GetReadSet(char * readLibraryFileName, int startIndex, ReadSetHead * readSetHead, char * read, int maxSize);

char * ReverseComplement(char * temp);






















#endif