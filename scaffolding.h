#ifndef SCAFFOLDING_H_INCLUDED 
#define SCAFFOLDING_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "contig.h"
#include "read.h"
#include "aligning.h"
#include "scaffoldgraph.h"

typedef struct ContigSequence{
    int index;
    int orientation;
    int gapDistance;
    ContigSequence * next;
}ContigSequence;

typedef struct ScaffoldSet{
    int length;
    ContigSequence * contigSequence;
	int * contigIndex;
	int contigNum;
    ScaffoldSet * next;
}ScaffoldSet;

typedef struct ScaffoldSetHead{
    ScaffoldSet * scaffoldSet;
}ScaffoldSetHead;




ScaffoldSetHead * GetScaffoldSet(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead);

ScaffoldGraphNode * GetOptimizeNodeIndex(ScaffoldGraphHead * scaffoldGraphHead, int nodeIndex, bool orientation, ContigSequence * tempContigSequence, bool right, bool & uniq);

int GetOverlapEdgeNum(ScaffoldGraphHead * scaffoldGraphHead, ContigSequence * contigSequence, int readIndex, bool right);

void OptimizeScaffoldSetCongtigSequence(ScaffoldSet * scaffoldSet, int contigNum);

bool DeterminOverlapInScaffolds(int * leftIndex, int leftNum, int * rightIndex, int rightNum, int * overlapResult);

bool GetOverlapBetweenArray(int * leftIndex, int leftNum, int * rightIndex, int rightNum, int * overlapResult);

void MergeContigSequence(ScaffoldSet * leftScaffoldSet, ScaffoldSet * rightScaffoldSet, int * overlapResult);

ContigSequence * GetContigSequenceIndex(int index, ContigSequence * contigSequence);

void InverseArray(int * array, int num);

int SearchScaffoldEdge(int index, ContigSequence * contigSequence);

void OutPutScaffoldTag(ScaffoldSet * scaffoldSet, char * result);

void OutPutScaffoldSet(ScaffoldSet * scaffoldSet, ContigSetHead * contigSetHead, char * result);




#endif