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
#include "read.h"
#include "kmer.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "scaffolding.h"

using namespace std;


int main(int argc,char** argv) {
	
	char * resultOutPutDirectory = (char *)malloc(sizeof(char)*50);
	strcpy(resultOutPutDirectory, "./SLR-OutPut-Directory/");
	long int contigLengthThreshold = 2000;
	long int readCount = 2000;
	long int readLengthCutOff = 3000;
	int readStart = 0;
	
	char * contigSetFile = NULL;
	char * longReadSetFile = NULL;
	int ch = 0;
	while ((ch = getopt(argc, argv, "c:r:p:m:n:")) != -1) {
		switch (ch) {
			case 'c': contigSetFile = (char *)(optarg); break;
			case 'r': longReadSetFile = (char *)(optarg); break;
			case 'p': resultOutPutDirectory = (char *)optarg; break;
			case 'm': contigLengthThreshold = atoi(optarg); break;
			case 'n': readLengthCutOff = atoi(optarg); break;
			default: break; 
		}
	}
	
	cout<<contigSetFile<<endl;
	cout<<longReadSetFile<<endl;
	cout<<resultOutPutDirectory<<endl;
	cout<<contigLengthThreshold<<endl;
	cout<<readLengthCutOff<<endl;
	
	if(opendir(resultOutPutDirectory) == NULL){  
		mkdir(resultOutPutDirectory, 0777);       
    }
	

	int fileNameLen = strlen(resultOutPutDirectory);

	ContigSetHead * contigSetHead = GetContigSet(contigSetFile, contigLengthThreshold);
	
	ScaffoldGraphHead * scaffoldGraphHead = (ScaffoldGraphHead *)malloc(sizeof(ScaffoldGraphHead));
	scaffoldGraphHead->scaffoldGraph = (ScaffoldGraph *)malloc(sizeof(ScaffoldGraph)*contigSetHead->contigCount);
	scaffoldGraphHead->scaffoldGraphNodeCount = contigSetHead->contigCount;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		scaffoldGraphHead->scaffoldGraph[i].outEdge = NULL;
		scaffoldGraphHead->scaffoldGraph[i].inEdge = NULL;
	}
	
	int maxSize = 2000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	
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
	
	
	GetAligningResut(contigSetHead, longReadSetFile, file3, file2, readLengthCutOff);
	
	SortNodeOptimize(contigSetHead, file3, line, maxSize, fp4);
	
	SortNodeOptimize(contigSetHead, file2, line, maxSize, fp5);
	
	fflush(fp4);
	
	fflush(fp5);

	GetScaffoldGraph(scaffoldGraphHead, contigSetHead, file4, line, maxSize);
	
	OptimizeScaffoldGraph(scaffoldGraphHead, contigSetHead);
	
	OutputScaffoldGraph(scaffoldGraphHead);
	
	ScaffoldSetHead * scaffoldSetHead = GetScaffoldSet(scaffoldGraphHead, contigSetHead, file5, line, maxSize, resultOutPutDirectory);
	
	char * scaffoldTagFileName = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(scaffoldTagFileName, resultOutPutDirectory);
	strcat(scaffoldTagFileName, "/scaffoldSet_tag.fa");
	OutPutScaffoldTag(scaffoldSetHead->scaffoldSet, scaffoldTagFileName);
	strcpy(scaffoldTagFileName, resultOutPutDirectory);
	strcat(scaffoldTagFileName, "/scaffoldSetFinal");
	OutPutScaffoldSet(scaffoldSetHead->scaffoldSet, contigSetHead, scaffoldTagFileName);
		
	
}
