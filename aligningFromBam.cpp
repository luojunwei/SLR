#ifndef aligningFromBam_CPP_INCLUDED 
#define aligningFromBam_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include<vector>


#include "aligningFromBam.h"

using namespace std;

int GetAligningResut(ContigSetHead * contigSetHead, char * fileName, char * result, char * result1, long int readLengthCutOff){
    
	long int i = 0;
    long int j = 0;
	long int readLength = 0;
	long int previousReadLength = 0;
	long int readIndex = 0;
	
	FILE * fp; 
    if((fp = fopen(result, "w")) == NULL){
        printf("%s, does not exist!", result);
        exit(0);
    }
	
	FILE * fp1; 
    if((fp1 = fopen(result1, "w")) == NULL){
        printf("%s, does not exist!", result1);
        exit(0);
    }
	
	
	string bamFileName = fileName;
	
    long int contigIndex = -1;
	long int previousContigIndex = -1;
	int readNameLength = 2000;
	string readName;
	string previousReadName = "a";
    
    BamReader bamReader;
    bamReader.Open(bamFileName);
    BamAlignment alignment;
    
    long int contigCount = bamReader.GetReferenceCount();
	
    
    
	AligningResultHead * aligningResultHead = (AligningResultHead *)malloc(sizeof(AligningResultHead));
	aligningResultHead->allocateAligningResultCount = 100;
	aligningResultHead->aligningResultCount = 0;
	aligningResultHead->aligningShortContigResultCount = 0;
	aligningResultHead->aligningResult = (AligningResult *)malloc(sizeof(AligningResult)*aligningResultHead->allocateAligningResultCount);
	for(i = 0; i < aligningResultHead->allocateAligningResultCount; i++){
		aligningResultHead->aligningResult[i].readStartPosition = -1;
		aligningResultHead->aligningResult[i].readEndPosition = -1;
		aligningResultHead->aligningResult[i].contigStartPosition = -1;
		aligningResultHead->aligningResult[i].contigEndPosition = -1;
		aligningResultHead->aligningResult[i].contigIndex = -1;
		aligningResultHead->aligningResult[i].overlapLength = -1;
		aligningResultHead->aligningResult[i].orientation = false;
		aligningResultHead->aligningResult[i].leftSoftClip = (long int *)malloc(sizeof(long int)*2);
		aligningResultHead->aligningResult[i].rightSoftClip = (long int *)malloc(sizeof(long int)*2);
		for(j = 0; j < 2; j ++){
			aligningResultHead->aligningResult[i].leftSoftClip[j] = -1;
			aligningResultHead->aligningResult[i].rightSoftClip[j] = -1;
		}
	}
	i = 0;
    while(bamReader.GetNextAlignment(alignment)){
        if(!alignment.IsMapped() || alignment.MapQuality < 30){                                    
            continue;
        }
		
		
		contigIndex = alignment.RefID;
		readLength = alignment.Length;
		if(contigSetHead->contigSet[contigIndex].contigLength < 500 || readLength < readLengthCutOff){
			continue;
		}
		
		
		readName = alignment.Name;
		
		if(previousReadName != "a" && readName != previousReadName){
			
			if(aligningResultHead->aligningResultCount > 1){

				OutputAligningResultOneLine(aligningResultHead, contigSetHead, fp, fp1, readIndex, previousReadLength);
				readIndex++;
			}
			i = 0;
			aligningResultHead->aligningResultCount = 0;
			aligningResultHead->aligningShortContigResultCount = 0;
		}
		
		if(GetAligningResultOneLine(aligningResultHead, alignment, contigSetHead, i) != false){
			i++;
		}
		previousReadLength = readLength;
		previousReadName = readName;
    }
	
	fflush(fp);
	fflush(fp1);
	fclose(fp);
	fclose(fp1);
    
    
}

bool GetAligningResultOneLine(AligningResultHead * aligningResultHead, BamAlignment alignment, ContigSetHead * contigSetHead, long int index){
	
	int maxAlignmentLength = 500;
	
	std::vector< int > clipSizes;
	std::vector< int > readPositions;
	std::vector< int > genomePositions;
	
	
	if(alignment.GetSoftClips(clipSizes, readPositions, genomePositions) != true){
		return false;
	}
	
	int readStartPosition = -1;
	int readEndPosition = -1;
	int contigStartPosition = -1;
	int contigEndPosition = -1;

	long int contigIndex = alignment.RefID;
	long int contigLength = contigSetHead->contigSet[contigIndex].contigLength;
	long int readLength = alignment.Length;
	
	if(clipSizes.size() == 1){
		if(clipSizes[0] == readPositions[0]){
			readStartPosition = readPositions[0];
			readEndPosition = alignment.Length - 1;
			contigStartPosition = genomePositions[0];
			contigEndPosition = alignment.GetEndPosition();
		}else{
			readStartPosition = 0;
			readEndPosition = alignment.Length - clipSizes[0] - 1;
			contigStartPosition = alignment.Position;
			contigEndPosition = alignment.GetEndPosition() - 1;
		}
	}else{
		readStartPosition = readPositions[0];
		readEndPosition = alignment.Length - clipSizes[1] - 1;
		contigStartPosition = genomePositions[0];
		contigEndPosition = alignment.GetEndPosition() -1;
	}
	
	
	
	
	if(readStartPosition < contigStartPosition){
		if(readStartPosition > maxAlignmentLength){
			return false;
		}
		aligningResultHead->aligningResult[index].readStartPosition = 0;
		aligningResultHead->aligningResult[index].contigStartPosition = contigStartPosition - readStartPosition;
	}else{
		if(contigStartPosition > maxAlignmentLength){
			return false;
		}
		aligningResultHead->aligningResult[index].readStartPosition = readStartPosition - contigStartPosition;
		aligningResultHead->aligningResult[index].contigStartPosition = 0;
	}
					
	if(readLength - readEndPosition < contigLength - contigEndPosition){
		if(readLength - readEndPosition > maxAlignmentLength){
			return false;
		}
		aligningResultHead->aligningResult[index].readEndPosition = readLength - 1;
		aligningResultHead->aligningResult[index].contigEndPosition = contigEndPosition + (readLength - readEndPosition - 1);
	}else{
		if(contigLength - contigEndPosition > maxAlignmentLength){
			return false;
		}
		aligningResultHead->aligningResult[index].readEndPosition = readEndPosition + contigLength - contigEndPosition - 1;
		aligningResultHead->aligningResult[index].contigEndPosition = contigLength - 1;
	}
	
	if(alignment.IsReverseStrand() != false){
		int temp = aligningResultHead->aligningResult[index].readStartPosition;
		aligningResultHead->aligningResult[index].readStartPosition = readLength - aligningResultHead->aligningResult[index].readEndPosition - 1;
		aligningResultHead->aligningResult[index].readEndPosition = readLength - temp - 1;
	}
	
	aligningResultHead->aligningResult[index].contigIndex = contigIndex;
	aligningResultHead->aligningResult[index].overlapLength = aligningResultHead->aligningResult[index].contigEndPosition - aligningResultHead->aligningResult[index].contigStartPosition + 1;
	aligningResultHead->aligningResult[index].orientation = !alignment.IsReverseStrand();
	
	if(aligningResultHead->aligningResult[index].overlapLength < 500){
		return false;
	}
	
	clipSizes.clear();
	readPositions.clear();
	genomePositions.clear();
	
	if(contigSetHead->contigSet[contigIndex].shortContig == true){
		if(aligningResultHead->aligningResult[index].contigStartPosition == 0 && aligningResultHead->aligningResult[index].contigEndPosition == contigLength - 1){
			aligningResultHead->aligningShortContigResultCount++;
		}else{
			return false;
		}
	}else{
		aligningResultHead->aligningResultCount++;
	}
	
	
	return true;
}


void OutputAligningResultOneLine(AligningResultHead * aligningResultHead, ContigSetHead * contigSetHead, FILE * fp, FILE * fp1, long int readIndex, long int readLength){
	
	for(long int i = 0; i < aligningResultHead->aligningResultCount + aligningResultHead->aligningShortContigResultCount - 1; i++){
		for(long int j = i + 1; j < aligningResultHead->aligningResultCount + aligningResultHead->aligningShortContigResultCount; j++){
			if(aligningResultHead->aligningResult[i].contigIndex == aligningResultHead->aligningResult[j].contigIndex){
				return;
			}
		}
	}
	
	for(long int i = 0; i < aligningResultHead->aligningResultCount + aligningResultHead->aligningShortContigResultCount - 1; i++){
		for(long int j = i + 1; j < aligningResultHead->aligningResultCount + aligningResultHead->aligningShortContigResultCount; j++){
			
			if(aligningResultHead->aligningResult[i].readStartPosition > aligningResultHead->aligningResult[j].readStartPosition){
				int temp = aligningResultHead->aligningResult[i].readStartPosition;
				aligningResultHead->aligningResult[i].readStartPosition = aligningResultHead->aligningResult[j].readStartPosition;
				aligningResultHead->aligningResult[j].readStartPosition = temp;
				
				temp = aligningResultHead->aligningResult[i].readEndPosition;
				aligningResultHead->aligningResult[i].readEndPosition = aligningResultHead->aligningResult[j].readEndPosition;
				aligningResultHead->aligningResult[j].readEndPosition = temp;
				
				temp = aligningResultHead->aligningResult[i].contigStartPosition;
				aligningResultHead->aligningResult[i].contigStartPosition = aligningResultHead->aligningResult[j].contigStartPosition;
				aligningResultHead->aligningResult[j].contigStartPosition = temp;
				
				temp = aligningResultHead->aligningResult[i].contigEndPosition;
				aligningResultHead->aligningResult[i].contigEndPosition = aligningResultHead->aligningResult[j].contigEndPosition;
				aligningResultHead->aligningResult[j].contigEndPosition = temp;
				
				temp = aligningResultHead->aligningResult[i].contigIndex;
				aligningResultHead->aligningResult[i].contigIndex = aligningResultHead->aligningResult[j].contigIndex;
				aligningResultHead->aligningResult[j].contigIndex = temp;
				
				temp = aligningResultHead->aligningResult[i].overlapLength;
				aligningResultHead->aligningResult[i].overlapLength = aligningResultHead->aligningResult[j].overlapLength;
				aligningResultHead->aligningResult[j].overlapLength = temp;
				
				temp = aligningResultHead->aligningResult[i].orientation;
				aligningResultHead->aligningResult[i].orientation = aligningResultHead->aligningResult[j].orientation;
				aligningResultHead->aligningResult[j].orientation = temp;
			}
		}

	}
	
	fprintf(fp, "%d,%ld", aligningResultHead->aligningResultCount, readLength);
	fprintf(fp1, "%d,%ld", aligningResultHead->aligningResultCount + aligningResultHead->aligningShortContigResultCount, readLength);
	for(long int i = 0; i < aligningResultHead->aligningResultCount + aligningResultHead->aligningShortContigResultCount; i++){
		if(contigSetHead->contigSet[aligningResultHead->aligningResult[i].contigIndex].shortContig != true){
			fprintf(fp, ",%d,%d,%d,%d,%d,%d,%d", aligningResultHead->aligningResult[i].contigIndex, aligningResultHead->aligningResult[i].readStartPosition, 
				aligningResultHead->aligningResult[i].readEndPosition, aligningResultHead->aligningResult[i].contigStartPosition, 
				aligningResultHead->aligningResult[i].contigEndPosition, aligningResultHead->aligningResult[i].overlapLength,
				aligningResultHead->aligningResult[i].orientation);
		}
		fprintf(fp1, ",%d,%d,%d,%d,%d,%d,%d", aligningResultHead->aligningResult[i].contigIndex, aligningResultHead->aligningResult[i].readStartPosition, 
			aligningResultHead->aligningResult[i].readEndPosition, aligningResultHead->aligningResult[i].contigStartPosition, 
			aligningResultHead->aligningResult[i].contigEndPosition, aligningResultHead->aligningResult[i].overlapLength,
			aligningResultHead->aligningResult[i].orientation);
		
	}
	fprintf(fp, ",\n");
	fprintf(fp1, ",\n");
	
}








#endif