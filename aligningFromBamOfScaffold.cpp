#ifndef aligningFromBamOfScaffold_CPP_INCLUDED 
#define aligningFromBamOfScaffold_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include<vector>


#include "aligningFromBamOfScaffold.h"

using namespace std;

int ContigNameToIndex(ContigSetHead * contigSetHead, char * name){
	for(int i = 0; i < contigSetHead->contigCount; i++){
		if(strcmp(contigSetHead->contigSet[i].name, name) == 0){
			return i;
			break;
		}
	}
}


void OutPutUniqueScaffoldSetHead(ContigSetHead * tempScaffoldHead, ContigSetHead * uniqueContigSetHead, char * bamFile, char * outPutFile, char * line, int maxSize){
	FILE * fp; 
	
	UniqueScaffoldSetHead * uniqueScaffoldSetHead = (UniqueScaffoldSetHead *)malloc(sizeof(UniqueScaffoldSetHead));
	uniqueScaffoldSetHead->uniqueScaffoldNum = tempScaffoldHead->contigCount;

	uniqueScaffoldSetHead->uniqueScaffoldSet = (UniqueContigToScaffold *)malloc(sizeof(UniqueContigToScaffold)*uniqueScaffoldSetHead->uniqueScaffoldNum);	

	for(int i = 0; i < uniqueScaffoldSetHead->uniqueScaffoldNum; i++){
		uniqueScaffoldSetHead->uniqueScaffoldSet[i].contigIndex = -1;
		uniqueScaffoldSetHead->uniqueScaffoldSet[i].orientation = false;
		uniqueScaffoldSetHead->uniqueScaffoldSet[i].scaffoldPosition = -1;
		uniqueScaffoldSetHead->uniqueScaffoldSet[i].contigCount = 0;
		uniqueScaffoldSetHead->uniqueScaffoldSet[i].next = NULL;
	}
	
	string bamFileName = bamFile;
	BamReader bamReader;
    bamReader.Open(bamFileName);
    BamAlignment alignment;
	
	
	std::vector< int > clipSizes;
	std::vector< int > readPositions;
	std::vector< int > genomePositions;
	
	int contigStartPosition = -1;
	int contigEndPosition = -1;
	
	
    long int contigIndex = -1;
	long int previousContigIndex = -1;
	int readNameLength = 2000;
	string readName;
	string previousReadName = "a";
   
    while(bamReader.GetNextAlignment(alignment)){
        
		
		readName = alignment.Name;
		if(!alignment.IsMapped() || alignment.MapQuality < 60){                                    
            previousReadName = readName;
			continue;
        }
		
		if(alignment.GetSoftClips(clipSizes, readPositions, genomePositions) == true){
			continue;
		}
		
		contigStartPosition = 0;
		contigEndPosition = alignment.Length - 1;
		
		int contigIndex = ContigNameToIndex(uniqueContigSetHead, const_cast<char*>(alignment.Name.c_str()));
		if(contigEndPosition != uniqueContigSetHead->contigSet[contigIndex].contigLength - 1){
			continue;
		}
		
		if(uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].contigIndex != -1){
			UniqueContigToScaffold * tempUniqueContigToScaffold = (UniqueContigToScaffold *)malloc(sizeof(UniqueContigToScaffold));
			tempUniqueContigToScaffold->contigIndex = atoi(const_cast<char*>(alignment.Name.c_str()));
			tempUniqueContigToScaffold->orientation = !alignment.IsReverseStrand();
			tempUniqueContigToScaffold->scaffoldPosition = alignment.Position;
			
			tempUniqueContigToScaffold->next = uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].next;
			uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].next = tempUniqueContigToScaffold;
			
		}else{
			uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].contigIndex = atoi(const_cast<char*>(alignment.Name.c_str()));
			uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].orientation = !alignment.IsReverseStrand();
			uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].scaffoldPosition = alignment.Position;
			uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].next = NULL;
		}
		uniqueScaffoldSetHead->uniqueScaffoldSet[alignment.RefID].contigCount++;
		
    }
	
	
	for(int i = 0; i < uniqueScaffoldSetHead->uniqueScaffoldNum; i++){
		UniqueContigToScaffold * tempUniqueContigToScaffold = & uniqueScaffoldSetHead->uniqueScaffoldSet[i];
		while(tempUniqueContigToScaffold->next != NULL){
			UniqueContigToScaffold * tempUniqueContigToScaffold1 = tempUniqueContigToScaffold->next;
			while(tempUniqueContigToScaffold1 != NULL){
				if(tempUniqueContigToScaffold1->scaffoldPosition < tempUniqueContigToScaffold->scaffoldPosition){
					int temp = tempUniqueContigToScaffold->contigIndex;
					tempUniqueContigToScaffold->contigIndex = tempUniqueContigToScaffold1->contigIndex;
					tempUniqueContigToScaffold1->contigIndex = temp;
					
					temp = tempUniqueContigToScaffold->scaffoldPosition;
					tempUniqueContigToScaffold->scaffoldPosition = tempUniqueContigToScaffold1->scaffoldPosition;
					tempUniqueContigToScaffold1->scaffoldPosition = temp;
					
					temp = tempUniqueContigToScaffold->orientation;
					tempUniqueContigToScaffold->orientation = tempUniqueContigToScaffold1->orientation;
					tempUniqueContigToScaffold1->orientation = temp;
					
				}
				tempUniqueContigToScaffold1 = tempUniqueContigToScaffold1->next;
			}
			tempUniqueContigToScaffold = tempUniqueContigToScaffold->next;
		}
	}
	
	
	if((fp = fopen(outPutFile, "w")) == NULL){
        printf("%s, does not exist!", outPutFile);
        exit(0);
    }
	int gapDistance = 0;
	
	for(int i = 0; i < uniqueScaffoldSetHead->uniqueScaffoldNum; i++){
		
		UniqueContigToScaffold * tempUniqueContigToScaffold = & uniqueScaffoldSetHead->uniqueScaffoldSet[i];
		if(tempUniqueContigToScaffold != NULL){
			fprintf(fp, "%d,", tempUniqueContigToScaffold->contigCount);
		}
		
		while(tempUniqueContigToScaffold != NULL){
			char string[30];
			sprintf(string, "%d", tempUniqueContigToScaffold->contigIndex);
			int tempIndex = ContigNameToIndex(uniqueContigSetHead, string);
			if(tempUniqueContigToScaffold->next != NULL){
				gapDistance = tempUniqueContigToScaffold->next->scaffoldPosition - tempUniqueContigToScaffold->scaffoldPosition 
					- uniqueContigSetHead->contigSet[tempIndex].contigLength;
			}else{
				gapDistance = 0;
			}
			fprintf(fp, "%d,%d,%d,", tempUniqueContigToScaffold->contigIndex, gapDistance, tempUniqueContigToScaffold->orientation);
			
			tempUniqueContigToScaffold = tempUniqueContigToScaffold->next;
		}
		fprintf(fp, "\n");
	}
	
	fflush(fp);
	fclose(fp);
	
}

ScaffoldSetHead * GetUniqueScaffoldSetHead(char * file, char * line, int maxSize){
	ScaffoldSetHead * scaffoldSetHead = (ScaffoldSetHead *)malloc(sizeof(ScaffoldSetHead));
	scaffoldSetHead->scaffoldSet = NULL;
	
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	char * p;
	const char * split = ","; 
	int contigCount = 0;
	int i = 0;
	ContigSequence * previous = NULL;
	while((fgets(line, maxSize, fp)) != NULL){ 
		ScaffoldSet * scaffoldSet = (ScaffoldSet *)malloc(sizeof(ScaffoldSet));
		scaffoldSet->contigSequence = NULL;
		p = strtok(line,split);
		contigCount = atoi(p);
		previous = NULL;
		i = 0;
		while(i < contigCount){
			ContigSequence * tempContigSequence = (ContigSequence *)malloc(sizeof(ContigSequence));
			p = strtok(NULL,split);
			tempContigSequence->index = atoi(p);
			
			p = strtok(NULL,split);
			tempContigSequence->gapDistance = atoi(p);
			p = strtok(NULL,split);
			tempContigSequence->orientation = atoi(p);
			tempContigSequence->next = NULL;
			if(i == 0){
				scaffoldSet->contigSequence = tempContigSequence;
			}
			if(previous != NULL){
				previous->next = tempContigSequence;
			}
			previous = tempContigSequence;
			i++;
		}
		
		scaffoldSet->next = scaffoldSetHead->scaffoldSet;
		scaffoldSetHead->scaffoldSet = scaffoldSet;
	}
	return scaffoldSetHead;
}









#endif