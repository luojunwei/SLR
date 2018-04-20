#ifndef Aligning_CPP_INCLUDED 
#define Aligning_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "aligning.h"

using namespace std;

int minAligningLength = 200;
double minOverlapRatioThreshold = 0.85;

int resultCountMax = 0;

bool InsertNodeSet(AligningGraph * aligningGraph, int leftNodeIndex, int rightNodeIndex){
	double minThreshold = 0.3;
	if(aligningGraph->aligningGraphNode[leftNodeIndex].contigPosition > aligningGraph->aligningGraphNode[rightNodeIndex].contigPosition
		|| aligningGraph->aligningGraphNode[leftNodeIndex].kmerPosition == aligningGraph->aligningGraphNode[rightNodeIndex].kmerPosition 
		|| aligningGraph->aligningGraphNode[leftNodeIndex].contigPosition == aligningGraph->aligningGraphNode[rightNodeIndex].contigPosition){
		return 0;
	}
	int readDistance = 0;
	int contigDistance = 0;
	
	readDistance = aligningGraph->aligningGraphNode[rightNodeIndex].kmerPosition - aligningGraph->aligningGraphNode[leftNodeIndex].kmerPosition;
	contigDistance = aligningGraph->aligningGraphNode[rightNodeIndex].contigPosition - aligningGraph->aligningGraphNode[leftNodeIndex].contigPosition;
	if(readDistance > minAligningLength || contigDistance > minAligningLength){
		return 0;
	}
	
	int minDistance = 0;
	if(readDistance < contigDistance){
		minDistance = readDistance;
	}else{
		minDistance = contigDistance;
	}
	
	double var = (double)(labs(readDistance - contigDistance))/(double)minDistance;
	if(var < minThreshold){
		return 1;
	}
	return 0;
}

void GetConnectedSubNodeSet(AligningGraph * aligningGraph){

	int * subSet = aligningGraph->subSet;
	
	int setIndex = 0;
	
	for(int i = 0; i < aligningGraph->aligningGraphNodeCount - 1; i++){
		if(subSet[i] == -1){
			subSet[i] = setIndex;
			setIndex++;
		}
		
		for(int j = i + 1; j < aligningGraph->aligningGraphNodeCount; j++){
			bool temp = InsertNodeSet(aligningGraph, i, j);
			if(temp == true){
				
				if(subSet[i] < 0 && subSet[j] < 0){
					subSet[i] = setIndex;
					subSet[j] = setIndex;
					setIndex++;
				}else{
					if(subSet[i] >= 0){
						subSet[j] = subSet[i];
					}else if(subSet[j] >= 0){
						subSet[i] = subSet[j];
					}
				}
			}
		}
		
	}
	
}



int GetOverlapRegion(AligningGraph * aligningGraph, int readLength, int contigLength){
	
	aligningGraph->subOverlapCount = 0;
	
	int * subSet = aligningGraph->subSet;
	
	bool * visited = aligningGraph->visited;
	
	int setIndex = -1;
	long int subSetCount = 0;
	
	int kmerPositionMin = readLength;
	int kmerPositionMax = -1;
	int readPositionMin = contigLength;
	int readPositionMax = -1;
	
	AligningResult * subOverlap = NULL;
	
	for(int i = 0; i < aligningGraph->aligningGraphNodeCount; i++){
		if(visited[i] == true || subSet[i] == -1){
			continue;
		}
		visited[i] = true;
		setIndex = subSet[i];
		kmerPositionMin = readLength;
		kmerPositionMax = -1;
		readPositionMin = contigLength;
		readPositionMax = -1;
		subSetCount = 0;
		for(int j = 0; j < aligningGraph->aligningGraphNodeCount; j++){
			if(i == j){
				//continue;
			}
			if(subSet[j] == setIndex){
				if(aligningGraph->aligningGraphNode[j].kmerPosition < kmerPositionMin){
					kmerPositionMin = aligningGraph->aligningGraphNode[j].kmerPosition;
					readPositionMin = aligningGraph->aligningGraphNode[j].contigPosition;
				}
				if(aligningGraph->aligningGraphNode[j].kmerPosition > kmerPositionMax){
					kmerPositionMax = aligningGraph->aligningGraphNode[j].kmerPosition;
					readPositionMax = aligningGraph->aligningGraphNode[j].contigPosition;
				}
				if(aligningGraph->aligningGraphNode[j].contigPosition < readPositionMin){
					//readPositionMin = aligningGraph->aligningGraphNode[j].contigPosition;
				}
				if(aligningGraph->aligningGraphNode[j].contigPosition > readPositionMax){
					//readPositionMax = aligningGraph->aligningGraphNode[j].contigPosition;
				}
				subSetCount++;
				visited[j] = true;
			}
		}
		if(subSetCount < 10){
			continue;
		}
		if(aligningGraph->subOverlapCount >= aligningGraph->allocateSubOverlapCount){
			aligningGraph->allocateSubOverlapCount = aligningGraph->allocateSubOverlapCount + 100;
			AligningResult * tempSubOverlap = (AligningResult *)malloc(sizeof(AligningResult)*aligningGraph->allocateSubOverlapCount);
			for(int i = 0; i < aligningGraph->subOverlapCount; i++){
				tempSubOverlap[i].readStartPosition = aligningGraph->subOverlap[i].readStartPosition;
				tempSubOverlap[i].readEndPosition = aligningGraph->subOverlap[i].readEndPosition;
				tempSubOverlap[i].contigStartPosition = aligningGraph->subOverlap[i].contigStartPosition;
				tempSubOverlap[i].contigEndPosition = aligningGraph->subOverlap[i].contigEndPosition;
			}
			free(aligningGraph->subOverlap);
			aligningGraph->subOverlap = tempSubOverlap;
			
		}
		
		aligningGraph->subOverlap[aligningGraph->subOverlapCount].readStartPosition = kmerPositionMin;
		aligningGraph->subOverlap[aligningGraph->subOverlapCount].readEndPosition = kmerPositionMax;
		aligningGraph->subOverlap[aligningGraph->subOverlapCount].contigStartPosition = readPositionMin;
		aligningGraph->subOverlap[aligningGraph->subOverlapCount].contigEndPosition = readPositionMax;
		
		aligningGraph->subOverlapCount ++;
	}
	
	return aligningGraph->subOverlapCount;
	
}

void * SortAligningGraphAligningResult(AligningGraph * aligningGraph){
	
	for(int i = 0; i < aligningGraph->subOverlapCount - 1; i++){
		for(int j = i + 1; j < aligningGraph->subOverlapCount; j++){
			if(aligningGraph->subOverlap[i].readStartPosition > aligningGraph->subOverlap[j].readStartPosition){
				
				int a = aligningGraph->subOverlap[i].readStartPosition;
				aligningGraph->subOverlap[i].readStartPosition = aligningGraph->subOverlap[j].readStartPosition;
				aligningGraph->subOverlap[j].readStartPosition = a;
		
				a = aligningGraph->subOverlap[i].readEndPosition;
				aligningGraph->subOverlap[i].readEndPosition = aligningGraph->subOverlap[j].readEndPosition;
				aligningGraph->subOverlap[j].readEndPosition = a;
		
				a = aligningGraph->subOverlap[i].contigStartPosition;
				aligningGraph->subOverlap[i].contigStartPosition = aligningGraph->subOverlap[j].contigStartPosition;
				aligningGraph->subOverlap[j].contigStartPosition = a;
		
				a = aligningGraph->subOverlap[i].contigEndPosition;
				aligningGraph->subOverlap[i].contigEndPosition = aligningGraph->subOverlap[j].contigEndPosition;
				aligningGraph->subOverlap[j].contigEndPosition = a;
				
			}
			
		}
		
		
	}	
}

AligningResult * ScoreOverlapRegion(AligningGraph * aligningGraph, int readLength, int readIndex, int contigIndex, int contigLength, FILE * fp, bool & token, bool orientation){
	
	SortAligningGraphAligningResult(aligningGraph);
	AligningResult * temp = aligningGraph->subOverlap;
	int subOverlapCount = aligningGraph->subOverlapCount;
	
	int readOverlapLength = 0;
	int contigOverlapLength = 0;
	
	int maxReadOverlapLength = 0;
	int maxContigOverlapLength = 0;
	int maxReadIntervalDistance = 0;
	int maxContigIntervalDistance = 0;
	
	int readOverlapStart = -1;
	int readOverlapEnd = -1;
	int contigOverlapStart = -1;
	int contigOverlapEnd = -1;
	
	int readIntervalLength = readLength;
	int contigIntervalLength = contigLength;
	
	int index = 0;
	int tt = 0;
	while(index >=0 ){
		int a = 0;
		temp = aligningGraph->subOverlap;
		if(index >= aligningGraph->subOverlapCount){
			break;
		}
		
		readOverlapLength = 0;
		contigOverlapLength = 0;
		
		int tempReadOverlapStart = temp[index].readStartPosition;
		int tempReadOverlapEnd = temp[index].readEndPosition;
		int tempContigOverlapStart = temp[index].contigStartPosition;
		int tempContigOverlapEnd = temp[index].contigEndPosition;
		
		readIntervalLength = 0;
		contigIntervalLength = 0;
		
		int readEnd = -1;
		int contigEnd = -1;
		bool b = false;
		
		tt++;
		int p = index;
		while(p < aligningGraph->subOverlapCount){
			
			if(temp[p].readStartPosition > readEnd && b == false){
				index++;
				b = true;
			}
			
			if(temp[p].readStartPosition > readEnd){
				readOverlapLength = readOverlapLength + temp[p].readEndPosition - temp[p].readStartPosition + 1;
				contigOverlapLength = contigOverlapLength + temp[p].contigEndPosition - temp[p].contigStartPosition + 1;
				readEnd = temp[p].readEndPosition;
				contigEnd = temp[p].contigEndPosition;
				tempReadOverlapEnd = temp[p].readEndPosition;
				tempContigOverlapEnd = temp[p].contigEndPosition;
				if(readIntervalLength < temp[p].readStartPosition - readEnd && readEnd != -1){
					readIntervalLength = temp[p].readStartPosition - readEnd;
				}
				if(contigIntervalLength < temp[p].contigStartPosition - contigEnd && contigEnd != -1){
					contigIntervalLength = temp[p].contigStartPosition - contigEnd;
				}
			}
			
			p ++;
		}
		
		if(readOverlapLength > maxReadOverlapLength){
			maxReadOverlapLength = readOverlapLength;
			maxContigOverlapLength = contigOverlapLength;
			
			readOverlapStart = tempReadOverlapStart;
			readOverlapEnd = tempReadOverlapEnd;
			contigOverlapStart = tempContigOverlapStart;
			contigOverlapEnd = tempContigOverlapEnd;
	
			maxReadIntervalDistance = readIntervalLength;
			maxContigIntervalDistance = contigIntervalLength;
		}
		
		if(a == index){
			break;
		}
		
		
	}
	
	int realReadOverlapStart = -1;
	int realReadOverlapEnd = -1;
	int realContigOverlapStart = -1;
	int realContigOverlapEnd = -1;
	
	if(readOverlapStart < contigOverlapStart){
		realReadOverlapStart = 0;
		realContigOverlapStart = contigOverlapStart - readOverlapStart;
	}else{
		realReadOverlapStart = readOverlapStart - contigOverlapStart;
		realContigOverlapStart = 0;
	}
					
	if(readLength - readOverlapEnd < contigLength - contigOverlapEnd){
		realReadOverlapEnd = readLength - 1;
		realContigOverlapEnd = contigOverlapEnd + (readLength - readOverlapEnd);
	}else{
		realReadOverlapEnd = readOverlapEnd + contigLength - contigOverlapEnd;
		realContigOverlapEnd = contigLength - 1;
	}
	
	
	double readOvelapRatio = (double)maxReadOverlapLength/(double)(realReadOverlapEnd - realReadOverlapStart);
	double contigOverlapRatio = (double)maxContigOverlapLength/(double)(realContigOverlapEnd - realContigOverlapStart);
	
	long int minOverlapLengthThreshold = 300;
	if(readOvelapRatio > minOverlapRatioThreshold && contigOverlapRatio > minOverlapRatioThreshold && readOverlapEnd - readOverlapStart > minOverlapLengthThreshold){
		if(token == false){
			fprintf(fp, "readIndex--%d--%d\n", readIndex, readLength);
			token = true;
		}
		fprintf(fp, "%d--%d--%d--%d--%d--%d--%d\n", orientation, realReadOverlapStart, realReadOverlapEnd, realContigOverlapStart, realContigOverlapEnd, contigIndex, maxReadOverlapLength);
		
	}

	return NULL;
	
}






void * SortAligningResult(AligningResult * aligningResult){
	if(aligningResult == NULL){
		return NULL;
	}
	
	AligningResult * next = aligningResult;
	AligningResult * temp = NULL;
	AligningResult * minTemp = NULL;
	int min = 9999999;
	
	while(next != NULL){
		temp = next;
		min = 9999999;
		while(temp != NULL){
			if(temp->readStartPosition < min){
				min = temp->readStartPosition;
				minTemp = temp;
			}
			temp = temp->next;
		}
		
		int a = minTemp->readStartPosition;
		minTemp->readStartPosition = next->readStartPosition;
		next->readStartPosition = a;
		
		a = minTemp->readEndPosition;
		minTemp->readEndPosition = next->readEndPosition;
		next->readEndPosition = a;
		
		a = minTemp->contigStartPosition;
		minTemp->contigStartPosition = next->contigStartPosition;
		next->contigStartPosition = a;
		
		a = minTemp->contigEndPosition;
		minTemp->contigEndPosition = next->contigEndPosition;
		next->contigEndPosition = a;
		
		a = minTemp->contigIndex;
		minTemp->contigIndex = next->contigIndex;
		next->contigIndex = a;
		
		a = minTemp->overlapLength;
		minTemp->overlapLength = next->overlapLength;
		next->overlapLength = a;
		
		next = next->next;
	}
	
	
	
	
}

int GetAligningGraph(KmerAligningSetHead * kmerAligningSetHead, AligningGraph * aligningGraph, int contigIndex, bool rc){
	int aligningGraphNodeCount = 0;
	aligningGraph->aligningGraphNodeCount = 0;
	for(int j = 0; j < kmerAligningSetHead->kmerAligningSetCount; j++){
		if(kmerAligningSetHead->kmerAligningSet[j].contigIndex == contigIndex && kmerAligningSetHead->kmerAligningSet[j].orientation == rc){
			aligningGraphNodeCount++;
		}
	}
	
	if(aligningGraphNodeCount < 10){
		return 0;
	}
	
	if(aligningGraphNodeCount == 0){
		return 0;
	}else{
		aligningGraph->aligningGraphNodeCount = aligningGraphNodeCount;
		if(aligningGraph->allocateAligningGraphNodeCount < aligningGraphNodeCount){
			aligningGraph->allocateAligningGraphNodeCount = aligningGraphNodeCount;
			free(aligningGraph->aligningGraphNode);
			free(aligningGraph->subSet);
			free(aligningGraph->visited);
			aligningGraph->aligningGraphNode = (AligningGraphNode *)malloc(sizeof(AligningGraphNode)*aligningGraph->aligningGraphNodeCount);
			aligningGraph->subSet = (int *)malloc(sizeof(int)*aligningGraph->aligningGraphNodeCount);
			aligningGraph->visited = (bool *)malloc(sizeof(bool)*aligningGraph->aligningGraphNodeCount);
		}
		for(int i = 0; i < aligningGraph->aligningGraphNodeCount; i++){
			aligningGraph->aligningGraphNode[i].contigPosition = -1;
			aligningGraph->aligningGraphNode[i].kmerPosition = -1;
			aligningGraph->aligningGraphNode[i].orientation = 0;
			
			aligningGraph->subSet[i] = -1;
			aligningGraph->visited[i] = false;
		}	
	}
	aligningGraphNodeCount = 0;
	
	for(int j = 0; j < kmerAligningSetHead->kmerAligningSetCount; j++){
		if(kmerAligningSetHead->kmerAligningSet[j].contigIndex == contigIndex && kmerAligningSetHead->kmerAligningSet[j].orientation == rc){
			aligningGraph->aligningGraphNode[aligningGraphNodeCount].kmerPosition = kmerAligningSetHead->kmerAligningSet[j].kmerPosition;
			aligningGraph->aligningGraphNode[aligningGraphNodeCount].contigPosition = kmerAligningSetHead->kmerAligningSet[j].contigPosition;
			aligningGraph->aligningGraphNode[aligningGraphNodeCount].orientation = kmerAligningSetHead->kmerAligningSet[j].orientation;
			
			aligningGraphNodeCount++;
		}
	}
	
	return aligningGraphNodeCount;

 
}

AligningResult * AligningReadContig(KmerAligningSetHead * kmerAligningSetHead, ContigSetHead * contigSetHead, AligningGraph * aligningGraph, bool & token, bool rc, FILE * fp){
	
	AligningResultSet * tempAligningResultSet = NULL;
	for(int i = 0; i < contigSetHead->contigCount; i++){
		contigSetHead->visited[i] = false;
	}
	
	int readIndex = kmerAligningSetHead->readIndex;
	int readLength = kmerAligningSetHead->readLength;
	
	for(int p = 0; p < kmerAligningSetHead->kmerAligningSetCount; p++){
		int contigIndex = kmerAligningSetHead->kmerAligningSet[p].contigIndex;
		if(contigSetHead->contigSet[contigIndex].shortContig != false){
			continue;
		}	
		if(contigSetHead->visited[contigIndex] == false && kmerAligningSetHead->kmerAligningSet[p].orientation == rc){
			contigSetHead->visited[contigIndex] = true;
		}else{
			continue;
		}
		
		int count = GetAligningGraph(kmerAligningSetHead, aligningGraph, contigIndex, rc);
		
		
		int contigLength = contigSetHead->contigSet[contigIndex].contigLength;
		if(count > 0){
				GetConnectedSubNodeSet(aligningGraph);

				int subOverlapCount = GetOverlapRegion(aligningGraph, readLength, contigLength);

				if(subOverlapCount > 0){
					ScoreOverlapRegion(aligningGraph, readLength, readIndex, contigIndex, contigLength, fp, token, rc);
				}

		}

	}


}

void IncreaseKmerAligningSet(KmerAligningSetHead * kmerAligningSetHead, int count){
	
	KmerAligningSet * tempKmerAligningSet = (KmerAligningSet *)malloc(sizeof(KmerAligningSet)*(kmerAligningSetHead->allocateKmerAligningSetCount + count));
	for(int i = 0; i < kmerAligningSetHead->allocateKmerAligningSetCount; i++){
		tempKmerAligningSet[i].contigIndex = kmerAligningSetHead->kmerAligningSet[i].contigIndex;
		tempKmerAligningSet[i].contigPosition = kmerAligningSetHead->kmerAligningSet[i].contigPosition;
		tempKmerAligningSet[i].kmerPosition = kmerAligningSetHead->kmerAligningSet[i].kmerPosition;
		tempKmerAligningSet[i].orientation = kmerAligningSetHead->kmerAligningSet[i].orientation;
	}
	for(int i = kmerAligningSetHead->allocateKmerAligningSetCount; i < kmerAligningSetHead->allocateKmerAligningSetCount + count; i++){
		tempKmerAligningSet[i].contigIndex = -1;
		tempKmerAligningSet[i].contigPosition = -1;
		tempKmerAligningSet[i].kmerPosition = -1;
		tempKmerAligningSet[i].orientation = false;
	}
	kmerAligningSetHead->allocateKmerAligningSetCount = kmerAligningSetHead->allocateKmerAligningSetCount + count;
	free(kmerAligningSetHead->kmerAligningSet);
	kmerAligningSetHead->kmerAligningSet = tempKmerAligningSet;
}

AligningResultSet * GetContigReadAligning(ContigSetHead * contigSetHead, ReadSetHead * readSetHead, AligningGraph * aligningGraph, int readStart, int readLengthCutOff, KmerAligningSetHead * kmerAligningSetHead, char * file, char * file1, char * line, int maxSize){
	
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	FILE * fp1; 
    if((fp1 = fopen(file1, "w")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }

	const char * split = "--"; 
	char * p;
	bool token = false;
	kmerAligningSetHead->kmerAligningSetCount = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		if(line[0] == 'r'){
			if(kmerAligningSetHead->kmerAligningSetCount >= 10){
				
				token = false;
				
				AligningReadContig(kmerAligningSetHead, contigSetHead, aligningGraph, token, 1, fp1);
				
				AligningReadContig(kmerAligningSetHead, contigSetHead, aligningGraph, token, 0, fp1);
				
			}
			p = strtok(line,split);
			p = strtok(NULL,split);
			kmerAligningSetHead->readIndex = atoi(p);
			p = strtok(NULL,split);
			kmerAligningSetHead->readLength = atoi(p);

			kmerAligningSetHead->kmerAligningSetCount = 0;
			continue;
		}
		
		p = strtok(line,split);
		kmerAligningSetHead->kmerAligningSet[kmerAligningSetHead->kmerAligningSetCount].kmerPosition = atoi(p);
		p = strtok(NULL,split);
		kmerAligningSetHead->kmerAligningSet[kmerAligningSetHead->kmerAligningSetCount].contigIndex = atoi(p);
		p = strtok(NULL,split);
		kmerAligningSetHead->kmerAligningSet[kmerAligningSetHead->kmerAligningSetCount].contigPosition = atoi(p);
		p = strtok(NULL,split);
		kmerAligningSetHead->kmerAligningSet[kmerAligningSetHead->kmerAligningSetCount].orientation = atoi(p);
		kmerAligningSetHead->kmerAligningSetCount++;
		if(kmerAligningSetHead->kmerAligningSetCount >= kmerAligningSetHead->allocateKmerAligningSetCount - 1){
			IncreaseKmerAligningSet(kmerAligningSetHead, 1000);
		}
		
    }
	if(kmerAligningSetHead->kmerAligningSetCount >= 10){
		token = false;
		AligningReadContig(kmerAligningSetHead, contigSetHead, aligningGraph, token, 1, fp1);
		AligningReadContig(kmerAligningSetHead, contigSetHead, aligningGraph, token, 0, fp1);
	}
	fflush(fp1);
	fclose(fp);
	fclose(fp1);
	return NULL;
}


void OptimizeALigningResultSet(AligningResultHead * aligningResultHead, ContigSetHead * contigSetHead, ContigGraphHead * contigGraphHead){
	
	for(int i = 0; i < contigGraphHead->contigGraphNodeCount; i++){
		contigGraphHead->contigGraph[i].outCount = 0;
		contigGraphHead->contigGraph[i].inCount = false;
		for(int j = 0; j < contigGraphHead->contigGraph[i].allocateOutCount; j++){
			contigGraphHead->contigGraph[i].out[j] = -1;
		}
		contigGraphHead->sortNode[i] = -1;
		contigGraphHead->distance[i] = -1;
		contigGraphHead->previousNode[i] = -1;
		contigGraphHead->longestPath[i] = -1;
		contigGraphHead->sortNodeCount = 0;
	}
	
	int readLength = aligningResultHead->readLength;
	int index = 0;
	
	for(int i = 0; i < aligningResultHead->aligningResultCount; i++){
		InsertSortNode(contigGraphHead->sortNode, contigGraphHead->distance, aligningResultHead->aligningResult[i].contigIndex, aligningResultHead->aligningResult[i].readStartPosition, contigGraphHead->contigGraphNodeCount);
		contigGraphHead->sortNodeCount++;
	}
	
	for(int i = 0; i < aligningResultHead->aligningResultCount - 1; i++){
		for(int j = i + 1; j < aligningResultHead->aligningResultCount; j++){
			InsertSingleEdgeToContigGraph(contigGraphHead, &aligningResultHead->aligningResult[i], aligningResultHead->aligningResult[i].orientation, 
				&aligningResultHead->aligningResult[j], aligningResultHead->aligningResult[j].orientation, readLength, contigSetHead);
		}
		
	}
	
	GetLongestPathInRead(contigSetHead, contigGraphHead);
		
}

void GetLongestPathInRead(ContigSetHead * contigSetHead, ContigGraphHead * contigGraphHead){
	int i = 0;
	int max = -1;
	int repeatMax = -1;
	while(contigGraphHead->sortNode[i] != -1){
		if(contigGraphHead->contigGraph[contigGraphHead->sortNode[i]].inCount != false){
			i++;
			continue;
		}
		
		for(int j = 0; j < contigGraphHead->contigGraphNodeCount; j++){
			contigGraphHead->distance[j] = -1;
		}
		contigGraphHead->distance[contigGraphHead->sortNode[i]] = 0;
		
		for(int j = 0; j < contigGraphHead->sortNodeCount; j++){
			
			if(contigGraphHead->distance[contigGraphHead->sortNode[j]] != -1){
				for(int p = 0; p < contigGraphHead->contigGraph[contigGraphHead->sortNode[j]].outCount; p++){
					int nodeIndex = contigGraphHead->sortNode[j];
					if(contigGraphHead->distance[nodeIndex] + 1 > contigGraphHead->distance[contigGraphHead->contigGraph[nodeIndex].out[p]]){
						contigGraphHead->distance[contigGraphHead->contigGraph[nodeIndex].out[p]] = contigGraphHead->distance[nodeIndex] + 1;
						contigGraphHead->previousNode[contigGraphHead->contigGraph[nodeIndex].out[p]] = nodeIndex;
					}
				}
			}
		}

		int tempMax = -1;
		int maxIndex = -1;
		
		for(int j = 0; j < contigGraphHead->sortNodeCount; j++){
			if(contigGraphHead->distance[contigGraphHead->sortNode[j]] > tempMax){
				tempMax = contigGraphHead->distance[contigGraphHead->sortNode[j]];
				maxIndex = contigGraphHead->sortNode[j];
			}
		}
		
		int maxCount = 0;
		
		for(int j = 0; j < contigGraphHead->sortNodeCount; j++){
			if(contigGraphHead->distance[contigGraphHead->sortNode[j]] == tempMax){
				maxCount++;
			}
		}
		
		
		
		if(tempMax > max && maxCount == 1){
			for(int j = 0; j < contigGraphHead->contigGraphNodeCount; j++){
				contigGraphHead->longestPath[j] = -1;
			}
			
			max = tempMax;
			int p = 0;
			while(maxIndex != -1){
				contigGraphHead->longestPath[p] = maxIndex;
				p++;
				maxIndex = contigGraphHead->previousNode[maxIndex];
			}
		}else if(tempMax == max && tempMax != -1){
			repeatMax = tempMax;
		}
		
		i++;
	}
	
	if(repeatMax == max){
		for(int j = 0; j < contigGraphHead->sortNodeCount; j++){
			contigGraphHead->longestPath[j] = -1;
		}
	}

}

int * IncreaseAllocateNode(int * out, int outCount, int allocateOutCount, int count){

	int * tempOut = (int *)malloc(sizeof(int)*(allocateOutCount + count));
	for(int i = 0; i < outCount; i++){
		tempOut[i] = out[i];
	}
	for(int i = outCount; i < allocateOutCount + count; i++){
		tempOut[i] = -1;
	}
	free(out);
	return tempOut;
}

void InsertSortNode(int * sortNode, int * distance, int contigIndex, int readStartPosition, int nodeCount){
	int tempPosition = -1;
	int tempDistance = -1;
	for(int i = 0; i < nodeCount; i++){
		if(readStartPosition <= distance[i] || sortNode[i] == -1){
			tempPosition = sortNode[i];
			tempDistance = distance[i];
			sortNode[i] = contigIndex;
			distance[i] = readStartPosition;
			
			int j = i + 1;
			while(tempPosition != -1){
				int temp = sortNode[j];
				sortNode[j] = tempPosition;
				tempPosition = temp;
				
				temp = distance[j];
				distance[j] = tempDistance;
				tempDistance = temp;
				
				j++;
			}
			break;
		}
	}

}




bool InsertSingleEdgeToContigGraph(ContigGraphHead * contigGraphHead, AligningResult * left, bool leftOrientation, AligningResult * right, bool rightOrientation, int readLength, ContigSetHead * contigSetHead){
	
	if(left->contigIndex == right->contigIndex || contigSetHead->repeatContigIndex[left->contigIndex] == true || contigSetHead->repeatContigIndex[right->contigIndex] == true){
		return false;
	}

	long int minOverlapLength = 500;
	float minRatio = 0.8;
	
	if(left->overlapLength < minOverlapLength && (float)left->overlapLength/(float)contigSetHead->contigSet[left->contigIndex].contigLength < minRatio){
		return false;
	}
	if(right->overlapLength < minOverlapLength && (float)right->overlapLength/(float)contigSetHead->contigSet[right->contigIndex].contigLength < minRatio){
		return false;
	}

	int gapDistance = 0;
	int gapDistance1 = 0;
	AligningResult * temp = NULL;
	
	int min = left->overlapLength;
	if(min > right->overlapLength){
		min = right->overlapLength;
	}

	if(leftOrientation == rightOrientation){
		
		gapDistance = right->readStartPosition - left->readEndPosition;
		gapDistance1 = left->readStartPosition - right->readEndPosition;
		if(gapDistance >= 0){
			int index = left->contigIndex;
			contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = right->contigIndex;
			contigGraphHead->contigGraph[right->contigIndex].inCount = true;
			contigGraphHead->contigGraph[index].outCount++;
			if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
				contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
				contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
			}
			return true;
		}else if(gapDistance1 >= 0){
			
			int index = right->contigIndex;
			contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = left->contigIndex;
			contigGraphHead->contigGraph[left->contigIndex].inCount = true;
			contigGraphHead->contigGraph[index].outCount++;
			if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
				contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
				contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
			}
			return true;
		}else{
			
			if(left->readStartPosition < right->readStartPosition && left->readEndPosition < right->readEndPosition){
				
				int index = left->contigIndex;
				contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = right->contigIndex;
				contigGraphHead->contigGraph[right->contigIndex].inCount = true;
				contigGraphHead->contigGraph[index].outCount++;
				if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
					contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
					contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
				}
				
				return true;
			}else if(left->readStartPosition > right->readStartPosition && left->readEndPosition > right->readEndPosition){
				
				int index = right->contigIndex;
				contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = left->contigIndex;
				contigGraphHead->contigGraph[left->contigIndex].inCount = true;
				contigGraphHead->contigGraph[index].outCount++;
				if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
					contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
					contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
				}
				return true;
			}
		}
		
	}
	
	if(leftOrientation != rightOrientation){
		
		if(leftOrientation != 1){
			temp = left;
			left = right;
			right = temp;
			leftOrientation = !leftOrientation;
			rightOrientation = !rightOrientation;
		}
		
		gapDistance = right->readStartPosition - left->readEndPosition;
		gapDistance1 = left->readStartPosition - right->readEndPosition;
		if(gapDistance >= 0){

			int index = left->contigIndex;
			contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = right->contigIndex;
			contigGraphHead->contigGraph[right->contigIndex].inCount = true;
			contigGraphHead->contigGraph[index].outCount++;
			if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
				contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
				contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
			}
			return true;
		}else if(gapDistance1 >= 0){

			int index = right->contigIndex;
			contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = left->contigIndex;
			contigGraphHead->contigGraph[left->contigIndex].inCount = true;
			contigGraphHead->contigGraph[index].outCount++;
			if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
				contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
				contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
			}
			return true;
		}else{
			if(left->readStartPosition < right->readStartPosition && left->readEndPosition < right->readEndPosition){

				int index = left->contigIndex;
				contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = right->contigIndex;
				contigGraphHead->contigGraph[right->contigIndex].inCount = true;
				contigGraphHead->contigGraph[index].outCount++;
				if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
					contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
					contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
				}
				return true;
			}else if(left->readStartPosition > right->readStartPosition && left->readEndPosition > right->readEndPosition){

				int index = right->contigIndex;
				contigGraphHead->contigGraph[index].out[contigGraphHead->contigGraph[index].outCount] = left->contigIndex;
				contigGraphHead->contigGraph[left->contigIndex].inCount = true;
				contigGraphHead->contigGraph[index].outCount++;
				if(contigGraphHead->contigGraph[index].outCount >= contigGraphHead->contigGraph[index].allocateOutCount){
					contigGraphHead->contigGraph[index].out = IncreaseAllocateNode(contigGraphHead->contigGraph[index].out, contigGraphHead->contigGraph[index].outCount, contigGraphHead->contigGraph[index].allocateOutCount, 10);
					contigGraphHead->contigGraph[index].allocateOutCount = contigGraphHead->contigGraph[index].allocateOutCount + 10;
				}
				return true;
			}
		}
		
	}
	return false;
	
	
}


void OutputAligningResultSet(AligningResultSet * aligningResultSet, char * file){
	ofstream ocout;
    ocout.open(file,ios::app);
	ocout<<"resultCountMax:"<<resultCountMax<<endl;
	while(aligningResultSet != NULL){
		ocout<<"readIndex:"<<aligningResultSet->readIndex<<"--"<<aligningResultSet->readLength<<endl;
		AligningResult * aligningResult1 = aligningResultSet->aligningResult;
		while(aligningResult1 != NULL){
			ocout<<"1--"<<aligningResult1->readStartPosition<<"--"<<aligningResult1->readEndPosition<<"--"<<aligningResult1->contigStartPosition<<"--"<<aligningResult1->contigEndPosition<<"--"<<aligningResult1->contigIndex<<"--"<<aligningResult1->overlapLength<<endl;
			aligningResult1 = aligningResult1->next;
		}
		aligningResult1 = aligningResultSet->rcAligningResult;
		while(aligningResult1 != NULL){
			ocout<<"0--"<<aligningResult1->readStartPosition<<"--"<<aligningResult1->readEndPosition<<"--"<<aligningResult1->contigStartPosition<<"--"<<aligningResult1->contigEndPosition<<"--"<<aligningResult1->contigIndex<<"--"<<aligningResult1->overlapLength<<endl;
			aligningResult1 = aligningResult1->next;
		}
		aligningResultSet = aligningResultSet->next;
	}

}


AligningResultSet * GetAligningResultSetFromFile(char * file){
	
	AligningResultSet * aligningResultSet = NULL;
	AligningResultSet * lastSet = NULL;
	AligningResultSet * tempAligningResultSet = NULL;
	AligningResult * last = NULL;
	AligningResult * lastRc = NULL;
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	char * line = NULL;
    if(NULL == (line = (char*)malloc(sizeof(char)*1000))){
        perror("malloc error!");
        exit(1);
    }
	
	int lineNum = 0;
	const char * split1 = "--"; 
	const char * split2 = "--"; 
	while((fgets(line, 1000, fp)) != NULL){
		if(line[0] == 'r'){
			if(tempAligningResultSet != NULL){
				if(lastSet != NULL){
					lastSet->next = tempAligningResultSet;
				}else{
					aligningResultSet = tempAligningResultSet;
				}
				lastSet = tempAligningResultSet;
			}
			lineNum ++;
			tempAligningResultSet = (AligningResultSet *)malloc(sizeof(AligningResultSet));
			tempAligningResultSet->aligningResult = NULL;
			tempAligningResultSet->rcAligningResult = NULL;
			tempAligningResultSet->next = NULL;
			last = NULL;
			lastRc = NULL;
			
			char * p; 
			p = strtok(line,split1);
			p = strtok(NULL,split1);
			char * p1;
			p1 = strtok(p,split2);
			tempAligningResultSet->readIndex = atoi(p1);
			//cout<<tempAligningResultSet->readIndex<<endl;
			p1 = strtok(NULL,split2);
			tempAligningResultSet->readLength = atoi(p1);
		}else if(line[0] == '1'){
			AligningResult * tempAligningResult = (AligningResult *)malloc(sizeof(AligningResult));
			tempAligningResult->next = NULL;
				
			char * p; 
			p = strtok(line,split2);
			p = strtok(NULL,split2);
			tempAligningResult->readStartPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->readEndPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->contigStartPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->contigEndPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->contigIndex = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->overlapLength = atoi(p);
			if(last != NULL){
				last->next = tempAligningResult;
			}else{
				tempAligningResultSet->aligningResult = tempAligningResult;
			}
			last = tempAligningResult;
			
		}else{
			AligningResult * tempAligningResult = (AligningResult *)malloc(sizeof(AligningResult));
			tempAligningResult->next = NULL;
				
			char * p; 
			p = strtok(line,split2);
			p = strtok(NULL,split2);
			tempAligningResult->readStartPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->readEndPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->contigStartPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->contigEndPosition = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->contigIndex = atoi(p);
			p = strtok(NULL,split2);
			tempAligningResult->overlapLength = atoi(p);
			if(lastRc != NULL){
				lastRc->next = tempAligningResult;
			}else{
				tempAligningResultSet->rcAligningResult = tempAligningResult;
			}
			lastRc = tempAligningResult;
		
		}
		lineNum++;
	}
	
	if(tempAligningResultSet != NULL){
		if(lastSet != NULL){
			lastSet->next = tempAligningResultSet;
		}else{
			aligningResultSet = tempAligningResultSet;
		}
		
		lastSet = tempAligningResultSet;
	}
	
	return aligningResultSet;
	
}


int GetAligningResultCount(AligningResult * aligningResult){
	int count = 0;
	while(aligningResult != NULL){
		count++;
		aligningResult = aligningResult->next;
	}
	return count;
}


























#endif
