#ifndef aligningFromBam_H_INCLUDED 
#define aligningFromBam_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "contig.h"

using namespace std;
using namespace BamTools;

typedef struct AligningResult{
	int readStartPosition;
	int readEndPosition;
	int contigStartPosition;
	int contigEndPosition;
	int contigIndex;
	int overlapLength;
	bool orientation;
	long int * leftSoftClip;
	long int * rightSoftClip;
}AligningResult;


typedef struct AligningResultHead{
	AligningResult * aligningResult;
	int allocateAligningResultCount;
	int aligningResultCount;
	int aligningShortContigResultCount;
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



int GetAligningResut(ContigSetHead * contigSetHead, char * bamFileName, char * result, char * result1, long int readLengthCutOff);

bool GetAligningResultOneLine(AligningResultHead * aligningResultHead, BamAlignment alignment, ContigSetHead * contigSetHead, long int index);

void OutputAligningResultOneLine(AligningResultHead * aligningResultHead, ContigSetHead * contigSetHead, FILE * fp, FILE * fp1, long int readIndex, long int readLength);







#endif