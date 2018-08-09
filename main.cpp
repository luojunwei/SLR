#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <malloc.h>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "scaffolding.h"

using namespace std;


int main(int argc,char** argv) {
	int maxSize = 2000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	
	char * resultOutPutDirectory = (char *)malloc(sizeof(char)*50);
	strcpy(resultOutPutDirectory, "./SLR-OutPut-Directory/");
	long int contigLengthThreshold = 2000;
	long int readCount = 2000;
	long int readLengthCutOff = 3000;
	int readStart = 0;
	int contigLengthIgnore = 0;
	int minAlignmentScore = 20;
	int minOverlapLength = 0;
	int overlapContigCount = 2;
	int minAlignmentRevised = 150;
	int weightType = 1;
	
	char * contigSetFile = NULL;
	char * longReadSetFile = NULL;
	char * alignSefBamFile = NULL;
	int ch = 0;
	while ((ch = getopt(argc, argv, "c:r:p:m:n:d:t:v:w:x:y:z:")) != -1) {
		switch (ch) {
			case 'c': contigSetFile = (char *)(optarg); break;
			case 'r': longReadSetFile = (char *)(optarg); break;
			case 'p': resultOutPutDirectory = (char *)optarg; break;
			case 'm': contigLengthThreshold = atoi(optarg); break;
			case 'n': readLengthCutOff = atoi(optarg); break;
			case 'd': alignSefBamFile = (char *)(optarg); break;
			case 't': weightType = atoi(optarg); break;
			case 'v': contigLengthIgnore = atoi(optarg); break;
			case 'w': minAlignmentScore = atoi(optarg); break;
			case 'x': minOverlapLength = atoi(optarg); break;
			case 'y': overlapContigCount = atoi(optarg); break;
			case 'z': minAlignmentRevised = atoi(optarg); break;
			default: break; 
		}
	}
	
	if(opendir(resultOutPutDirectory) == NULL){  
		mkdir(resultOutPutDirectory, 0777);       
    }
	
	int fileNameLen = strlen(resultOutPutDirectory);

	ContigSetHead * contigSetHead = GetContigSet(contigSetFile, contigLengthThreshold);
	contigSetHead->contigLengthIgnore = contigLengthIgnore;
	contigSetHead->minAlignmentScore = minAlignmentScore;
	contigSetHead->minOverlapLength = minOverlapLength;
	contigSetHead->overlapContigCount = overlapContigCount;
	contigSetHead->minAlignmentRevised = minAlignmentRevised;
	contigSetHead->weightType = weightType;
	
	if(alignSefBamFile != NULL){
		DetectRepeativeContigInSet(contigSetHead, alignSefBamFile, 1);
	}
	
	ScaffoldGraphHead * scaffoldGraphHead = (ScaffoldGraphHead *)malloc(sizeof(ScaffoldGraphHead));
	scaffoldGraphHead->scaffoldGraph = (ScaffoldGraph *)malloc(sizeof(ScaffoldGraph)*contigSetHead->contigCount);
	scaffoldGraphHead->scaffoldGraphNodeCount = contigSetHead->contigCount;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		scaffoldGraphHead->scaffoldGraph[i].outEdge = NULL;
		scaffoldGraphHead->scaffoldGraph[i].inEdge = NULL;
	}
	
	
	
	char * file2 = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(file2, resultOutPutDirectory);
	strcat(file2, "/shortContigPathInLongRead.fa");
	
	char * file3 = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(file3, resultOutPutDirectory);
	strcat(file3, "/originalPathInLongRead.fa");
	
	char * file4 = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(file4, resultOutPutDirectory);
	strcat(file4, "/optimizePathInLongRead.fa");
	
	FILE * fp4; 
    if((fp4 = fopen(file4, "w")) == NULL){
        printf("%s, does not exist!", file4);
        exit(0);
    }
	
	char * file5 = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(file5, resultOutPutDirectory);
	strcat(file5, "/shortContigOptimizePathInLongRead.fa");
	
	FILE * fp5; 
    if((fp5 = fopen(file5, "w")) == NULL){
        printf("%s, does not exist!", file5);
        exit(0);
    }
	
	char * file6 = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(file6, resultOutPutDirectory);
	strcat(file6, "/unique-contig-set.fa");
	
	FILE * fpUnique; 
    if((fpUnique = fopen(file6, "w")) == NULL){
        printf("%s, does not exist!", file6);
        exit(0);
    }
	
	char * file7 = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(file7, resultOutPutDirectory);
	strcat(file7, "/ambiguous-contig-set.fa");
	
	FILE * fpAmbiguous; 
    if((fpAmbiguous = fopen(file7, "w")) == NULL){
        printf("%s, does not exist!", file7);
        exit(0);
    }
	
	
	GetAligningResut(contigSetHead, longReadSetFile, file3, file2, readLengthCutOff);
	
	SortNodeOptimize(contigSetHead, file3, line, maxSize, fp4);

	SortNodeOptimize(contigSetHead, file2, line, maxSize, fp5);
	
	fflush(fp4);
	
	fflush(fp5);

	GetScaffoldGraph(scaffoldGraphHead, contigSetHead, file4, line, maxSize, fpUnique, fpAmbiguous);
	
	OptimizeScaffoldGraph(scaffoldGraphHead, contigSetHead);
	
	OutputScaffoldGraph(scaffoldGraphHead);

	ScaffoldSetHead * scaffoldSetHead = GetScaffoldSet(scaffoldGraphHead, contigSetHead, file5, line, maxSize, resultOutPutDirectory, file4);
	
	char * scaffoldTagFileName = (char *)malloc(sizeof(char)*fileNameLen + 50);
	
	strcpy(scaffoldTagFileName, resultOutPutDirectory);
	strcat(scaffoldTagFileName, "/scaffold");
	OutPutScaffoldSet(scaffoldSetHead->scaffoldSet, contigSetHead, scaffoldTagFileName);
		
	
}
