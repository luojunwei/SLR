#ifndef aligningFromBamOfScaffold_H_INCLUDED 
#define aligningFromBamOfScaffold_H_INCLUDED 
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

#include "scaffolding.h"

using namespace std;
using namespace BamTools;

typedef struct UniqueContigToScaffold{
	int scaffoldPosition;
	bool orientation;
	int contigIndex;
	int contigCount;
	UniqueContigToScaffold * next;
}UniqueContigToScaffold;

typedef struct UniqueScaffoldSetHead{
	UniqueContigToScaffold * uniqueScaffoldSet;
	long int uniqueScaffoldNum;
}UniqueScaffoldSetHead;


int ContigNameToIndex(ContigSetHead * contigSetHead, char * name);

void OutPutUniqueScaffoldSetHead(ContigSetHead * tempContigHead, ContigSetHead * uniqueContigSetHead, char * bamFile, char * outPutFile, char * line, int maxSize);

ScaffoldSetHead * GetUniqueScaffoldSetHead(char * file, char * line, int maxSize);

#endif