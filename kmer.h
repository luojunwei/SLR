#ifndef Kmer_H_INCLUDED 
#define Kmer_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "contig.h"

using namespace std;

typedef struct KmerHashNode{
    int position;
    int index;
	KmerHashNode * next;
}KmerHashNode;

typedef struct KmerHashNodeHead{
	KmerHashNode * kmerHashNode;
}KmerHashNodeHead;

typedef struct KmerHashTableHead{
    KmerHashNodeHead * kmerHashNodeHead;
	long int kmerHashNodeCount;
    long int kmerHashTableCount;
}KmerHashTableHead;


long int Hash(char * str, unsigned int len, unsigned long int max);

int SearchKmerInHashTable(ContigSetHead * contigSetHead, KmerHashNodeHead * kmerHashNodeHead, long int kmerHashTableCount, char * kmer, int kmerLength, char * tempKmer);

void InsertKmerToHashTable(ContigSetHead * contigSetHead, KmerHashTableHead * kmerHashTableHead, char * kmer, int kmerLength, long int index, long int position, char * tempKmer);

KmerHashTableHead * GetKmerHashTableHead(ContigSetHead * contigSetHead, int kmerLength);

void ReverseComplementKmer(char * kmer, int kmerLength);







#endif