#ifndef Scaffoldgraph_H_INCLUDED 
#define Scaffoldgraph_H_INCLUDED 
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
#include "aligningFromBam.h"

using namespace std;

typedef struct ScaffoldGraphNode{
	int nodeIndex;
	int readIndex;
	int gapDistance;
	int aligningReadCount;
	bool orientation;
	int overlapLength;
	int * readIndexArray;
	ScaffoldGraphNode * next;
}ScaffoldGraphNode;


typedef struct ScaffoldGraph{
	ScaffoldGraphNode * outEdge;
	ScaffoldGraphNode * inEdge;
}ScaffoldGraph;


typedef struct ScaffoldGraphHead{
	ScaffoldGraph * scaffoldGraph;
	int scaffoldGraphNodeCount;
}ScaffoldGraphHead;

typedef struct SimpleGraph{
	int leftIndex;
	int rightIndex;
	int count1;
	int count2;
	int count3;
	int count4;
}SimpleGraph;



void InsertOutOrInEdge(ScaffoldGraphHead * scaffoldGraphHead, int readIndex, int leftNodeIndex, bool leftOrientation, int rightNodeIndex, bool rightOrientation, int gapDistance, int overlapLength);

bool InsertSingleEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, int readIndex, AligningResult * left, bool leftOrientation, AligningResult * right, bool rightOrientation, int readLength, ContigSetHead * contigSetHead);

void InsertEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, AligningResult * aligningResult, AligningResult * rcAligningResult, int readIndex, int readLength, ContigSetHead * contigSetHead);

void InsertEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, AligningResultHead * aligningResultHead, ContigGraphHead * contigGraphHead, ContigSetHead * contigSetHead);

void IncreaseAligningResultHead(AligningResultHead * aligningResultHead, int count);

void GetPathInLongRead(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, ContigGraphHead * contigGraphHead, char * file, char * line, int maxSize, FILE * fp1, FILE * originalPathFileFP);

void GetScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, ContigGraphHead * contigGraphHead, char * file, char * line, int maxSize);

void GetScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, char * file, char * line, int maxSize);

void DeleteEdgeWithMinReadCount(ScaffoldGraphHead * scaffoldGraphHead, int minReadCount);

void OutputScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead);

void DeleteScaffoldGraphNode(ScaffoldGraphNode * edge);

ScaffoldGraphNode * DeleteScaffoldGraphSpecialNode(ScaffoldGraphNode * edge, int nodeIndex, bool orientation);

int DeleteSpecialScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, long int index1);

ScaffoldGraphNode * MergeMultipleEdges(ScaffoldGraphNode * edge, int contigLength, int contigLength0);

ScaffoldGraphNode * OptimizeEdgeInScaffoldGraph(ScaffoldGraphNode * edge, ContigSetHead * contigSetHead, int nodeIndex);

void OptimizeScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead);

void MergeBubbleInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead);

bool VisitOutNode(ScaffoldGraphHead * scaffoldGraphHead, int leftIndex, bool out, int rightIndex, bool * visited);

bool VisitInNode(ScaffoldGraphHead * scaffoldGraphHead, int rightIndex, bool in, int leftIndex, bool * visited);

void RemoveCycleInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead);

int GetEdgeNumber(ScaffoldGraphHead * scaffoldGraphHead, int index, bool out);

bool GetOverlapEdgeIndex(ScaffoldGraphNode * left, ScaffoldGraphNode * right);

void SortNodeOptimize(ContigSetHead * contigSetHead, char * file, char * line, int maxSize, FILE * longReadFileFP);

bool * GetLineIndex(ContigSetHead * contigSetHead, char * file, char * line, int maxSize, long int & lineCount);

void AddEdgeInSimpleGraph(SimpleGraph * simpleGraph, long int simpleGraphNodeCount, int leftIndex, bool leftOrientation, int rightIndex, bool rightOrientation);

#endif