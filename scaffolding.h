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
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "lp/lp_lib.h"

using namespace std;

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
	int sequenceCount;
    ScaffoldSet * next;
}ScaffoldSet;

typedef struct ScaffoldSetHead{
    ScaffoldSet * scaffoldSet;
}ScaffoldSetHead;

long int DetermineOrientationOfContigs(ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation);

long int * IterativeDetermineOrderOfContigs(ContigSetHead * contigSetHead, ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation, long int * contigPosition, long int allContigLength);

ScaffoldSetHead * GetScaffoldSet(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, char * file, char * line, int maxSize, char * dir, char * longContigfile);

ScaffoldGraphNode * GetOptimizeNodeIndex(ScaffoldGraphHead * scaffoldGraphHead, int nodeIndex, bool orientation, ContigSequence * tempContigSequence, bool right, bool * printIndex);

void OptimizeScaffoldSetCongtigSequence(ScaffoldSet * scaffoldSet, int contigNum);

bool DeterminOverlapInScaffolds(int * leftIndex, int leftNum, int * rightIndex, int rightNum, int * overlapResult);

bool GetOverlapBetweenArray(int * leftIndex, int leftNum, int * rightIndex, int rightNum, int * overlapResult);

void MergeContigSequence(ScaffoldSet * leftScaffoldSet, ScaffoldSet * rightScaffoldSet, int * overlapResult);

ContigSequence * GetContigSequenceIndex(int index, ContigSequence * contigSequence);

void InverseArray(int * array, int num);

int SearchScaffoldEdge(int index, ContigSequence * contigSequence);

void OutPutScaffoldTag(ScaffoldSet * scaffoldSet, char * result);

void OutPutScaffoldSet(ScaffoldSet * scaffoldSet, ContigSetHead * contigSetHead, char * result);

void InsertRepeatContigToSequence(ContigSetHead * contigSetHead, ScaffoldSetHead * scaffoldSetHead, bool * printIndex, char * file, char * line, int maxSize);

bool SearchInsertContigSequence(ScaffoldSetHead * tempInsertSequenceHead, int startIndex, int endIndex, int * contigIndex, int * distance, bool * orientation, int * overlapLength, int count);

void InsertShortContigToSequence(ContigSequence * tempContigSequence, ContigSequence * preContigSequence, char * file, char * line, int maxSize, bool * printIndex, bool * lineIndex, long int lineCount);

void InsertRepeatContigToSequence1(ContigSetHead * contigSetHead, ScaffoldSetHead * scaffoldSetHead, bool * printIndex, char * file, char * line, int maxSize);

int GetContigSequenceNum(ContigSequence * tempContigSequence);

void OverlapHeadAndTail(ContigSetHead * contigSetHead, ScaffoldSetHead * scaffoldSetHead, bool * printIndex, LocalScaffoldSetHead * localScaffoldSetHead);

void OptimizeScaffoldSetCongtigSequence1(ScaffoldSet * scaffoldSet, int contigNum);

ScaffoldSetHead *  LocalScaffoldToScaffoldSet(LocalScaffoldSetHead * localScaffoldSetHead, ContigSetHead * contigSetHead, char * dir);

void * OutputLocalScaffoldNonUnique(ContigSetHead * contigSetHead, LocalScaffoldSetHead * localScaffoldSetHead, char * outputFile);

void OverlapHeadAndTailGraph(ContigSetHead * contigSetHead, ScaffoldSetHead * scaffoldSetHead, LocalScaffoldSetHead * localScaffoldSetHead, ScaffoldGraphHead * scaffoldGraphHead, char * shortGraph, bool * printIndex);

void OutputUnUsedLocalScaffold(ContigSetHead * contigSetHead, ScaffoldSetHead * scaffoldSetHead, LocalScaffoldSetHead * localScaffoldSetHead, ScaffoldGraphHead * scaffoldGraphHead, char * file, char * line, long int maxSize);

void sort(int * sortNode, int * length, long int left, long int right);

#endif
