#ifndef Aligning_H_INCLUDED 
#define Aligning_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "contig.h"
#include "read.h"

using namespace std;

typedef struct KmerAligningSet{
	int kmerPosition;
	int contigIndex;
	int contigPosition;
	bool orientation;
}KmerAligningSet;

typedef struct KmerAligningSetHead{
	KmerAligningSet * kmerAligningSet;
	int readIndex;
	int readLength;
	int kmerAligningSetCount;
	int allocateKmerAligningSetCount;
}KmerAligningSetHead;

typedef struct AligningResult{
	int readStartPosition;
	int readEndPosition;
	int contigStartPosition;
	int contigEndPosition;
	int contigIndex;
	int overlapLength;
	bool orientation;
	AligningResult * next;
}AligningResult;

typedef struct AligningResultHead{
	AligningResult * aligningResult;
	int readIndex;
	int readLength;
	int aligningResultCount;
	int allocateAligningResultCount;
}AligningResultHead;

typedef struct AligningResultSet{
	AligningResult * aligningResult;
	AligningResult * rcAligningResult;
	int readIndex;
	int readLength;
	int * sortNode;
	int sortNodeCount;
	bool sortNodeOrientation;
	AligningResultSet * next;
}AligningResultSet;

typedef struct ArrayList{
	int index;
	ArrayList * next;
}ArrayList;

typedef struct AligningGraphEdge{
	int nodeIndex;
	double weight;
	AligningGraphEdge * next;
}AligningGraphEdge;

typedef struct AligningGraphNode{
	int contigPosition;
	int kmerPosition;
	bool orientation;
}AligningGraphNode;

typedef struct AligningGraph{
	AligningGraphNode * aligningGraphNode;
	int aligningGraphNodeCount;
	int allocateAligningGraphNodeCount;
	int * subSet;
	bool * visited;
	AligningResult * subOverlap;
	int subOverlapCount;
	int allocateSubOverlapCount;
}AligningGraph;

typedef struct ContigGraph{
	int * out;
	bool inCount;
	int outCount;
	int allocateOutCount;
}ContigGraph;

typedef struct ContigGraphHead{
	ContigGraph * contigGraph;
	int contigGraphNodeCount;
	int * sortNode;
	int sortNodeCount;
	int * distance;
	int * previousNode;
	int * longestPath;
}ContigGraphHead;

bool InsertNodeSet(AligningGraph * aligningGraph, int leftNodeIndex, int rightNodeIndex);

void GetConnectedSubNodeSet(AligningGraph * aligningGraph);

void * SortAligningGraphAligningResult(AligningGraph * aligningGraph);

int GetOverlapRegion(AligningGraph * aligningGraph, int readLength, int contigLength);

AligningResult * ScoreOverlapRegion(AligningGraph * aligningGraph, int readLength, int readIndex, int contigIndex, int contigLength, FILE * fp, bool & token, bool orientation);

void * SortAligningResult(AligningResult * aligningResult);

void IncreaseKmerAligningSet(KmerAligningSetHead * kmerAligningSetHead, int count);

int GetAligningGraph(KmerAligningSetHead * kmerAligningSetHead, AligningGraph * aligningGraph, int contigIndex, bool rc);

AligningResult * AligningReadContig(KmerAligningSetHead * kmerAligningSetHead, ContigSetHead * contigSetHead, AligningGraph * aligningGraph, bool & token, bool rc, FILE * fp);

AligningResultSet * GetContigReadAligning(ContigSetHead * contigSetHead, ReadSetHead * readSetHead, AligningGraph * aligningGraph, int readStart, int readLengthCutOff, KmerAligningSetHead * kmerAligningSetHead, char * file, char * file1, char * line, int maxSize);

int GetAligningResultCount(AligningResult * aligningResult);

AligningResultSet * GetAligningResultSetFromFile(char * file);

void OptimizeALigningResultSet(AligningResultHead * aligningResultHead, ContigSetHead * contigSetHead, ContigGraphHead * contigGraphHead);

bool InsertSingleEdgeToContigGraph(ContigGraphHead * contigGraphHead, AligningResult * left, bool leftOrientation, AligningResult * right, bool rightOrientation, int readLength, ContigSetHead * contigSetHead);

void GetLongestPathInRead(ContigSetHead * contigSetHead, ContigGraphHead * contigGraphHead);

int * IncreaseAllocateNode(int * out, int outCount, int allocateOutCount, int count);

void InsertSortNode(int * sortNode, int * distance, int contigIndex, int readStartPosition, int nodeCount);

void OutputAligningResultSet(AligningResultSet * aligningResultSet, char * file);

#endif
