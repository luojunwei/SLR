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
#include "aligning.h"

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



void InsertOutOrInEdge(ScaffoldGraphHead * scaffoldGraphHead, int readIndex, int leftNodeIndex, bool leftOrientation, int rightNodeIndex, bool rightOrientation, int gapDistance, int overlapLength);

bool InsertSingleEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, int readIndex, AligningResult * left, bool leftOrientation, AligningResult * right, bool rightOrientation, int readLength, ContigSetHead * contigSetHead);

void InsertEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, AligningResult * aligningResult, AligningResult * rcAligningResult, int readIndex, int readLength, ContigSetHead * contigSetHead);

void InsertEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, AligningResultHead * aligningResultHead, ContigGraphHead * contigGraphHead, ContigSetHead * contigSetHead);

void IncreaseAligningResultHead(AligningResultHead * aligningResultHead, int count);

void GetScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, ContigGraphHead * contigGraphHead, char * file, char * line, int maxSize);

void DeleteEdgeWithMinReadCount(ScaffoldGraphHead * scaffoldGraphHead, int minReadCount);

void OutputScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead);

void DeleteScaffoldGraphNode(ScaffoldGraphNode * edge);

ScaffoldGraphNode * DeleteScaffoldGraphSpecialNode(ScaffoldGraphNode * edge, int nodeIndex, bool orientation);

ScaffoldGraphNode * MergeMultipleEdges(ScaffoldGraphNode * edge, int contigLength, int contigLength0);

ScaffoldGraphNode * OptimizeEdgeInScaffoldGraph(ScaffoldGraphNode * edge, ContigSetHead * contigSetHead, int nodeIndex);

void OptimizeScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead);

void MergeBubbleInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead);

bool VisitOutNode(ScaffoldGraphHead * scaffoldGraphHead, int leftIndex, bool out, int rightIndex, bool * visited);

bool VisitInNode(ScaffoldGraphHead * scaffoldGraphHead, int rightIndex, bool in, int leftIndex, bool * visited);

void RemoveCycleInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead);

#endif