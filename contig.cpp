#ifndef CONTIG_CPP_INCLUDED 
#define CONTIG_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "contig.h"

using namespace std;

ContigSetHead * GetContigSet(char * contigSetFile, int contigLengthThreshold){
    
    ContigSetHead * contigSetHead = (ContigSetHead *)malloc(sizeof(ContigSetHead));
    contigSetHead->contigSet = NULL;
    contigSetHead->contigCount = 0;
	contigSetHead->allContigLength =  0;
	contigSetHead->minContigSet = NULL;
    contigSetHead->minContigCount = 0;
	contigSetHead->minAllContigLength = 0;
    
    int maxSize = 90000;
    char * contig = NULL;
    if(NULL == (contig = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    FILE * fp; 
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    while((fgets(contig, maxSize, fp)) != NULL){ 
       if(contig[0] == '>'){  
           contigSetHead->contigCount++; 
       }  
    }  
    fclose(fp);
    
    contigSetHead->contigSet = (Contig *)malloc(sizeof(Contig)*contigSetHead->contigCount);
	contigSetHead->repeatContigIndex = (bool *)malloc(sizeof(bool)*contigSetHead->contigCount);
    for(int i = 0; i < contigSetHead->contigCount; i++){
        contigSetHead->contigSet[i].contig = NULL;
        contigSetHead->contigSet[i].contigLength = 0;
		contigSetHead->contigSet[i].shortContig = false;
		contigSetHead->contigSet[i].realContigIndex = -1;
		contigSetHead->repeatContigIndex[i] = false;
    }
    
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    
    int allocateLength = 0;
    int contigIndex = -1;
    while((fgets(contig, maxSize, fp)) != NULL){ 
       
       if(contig[0] == '>'){  
           if(strlen(contig) == maxSize-1){              
               while((fgets(contig, maxSize, fp)) != NULL){
                   if(strlen(contig) != maxSize-1){
                       break;
                   }
               }        
           }
		   if(contigIndex != -1){
			   contigSetHead->allContigLength = contigSetHead->allContigLength + contigSetHead->contigSet[contigIndex].contigLength;
		   }
		   contigIndex++;
           continue;
       }
       
       
       int extendLength = strlen(contig);
       if(contig[extendLength-1] == '\n'){
           extendLength--;
       }
       int contigLength = 0;
       char * tempContig = NULL;
       if(contigSetHead->contigSet[contigIndex].contig != NULL){
           if(contigSetHead->contigSet[contigIndex].contigLength + extendLength >= allocateLength){
               contigLength = contigSetHead->contigSet[contigIndex].contigLength;    
               contigSetHead->contigSet[contigIndex].contig = (char *)realloc(contigSetHead->contigSet[contigIndex].contig, allocateLength + maxSize + 1);

               allocateLength = allocateLength + maxSize + 1;
               
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigLength + extendLength] = '\0';    
               contigSetHead->contigSet[contigIndex].contigLength = contigLength + extendLength;
                       
           }else{
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigSetHead->contigSet[contigIndex].contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigSetHead->contigSet[contigIndex].contigLength + extendLength] = '\0';
               contigSetHead->contigSet[contigIndex].contigLength = contigSetHead->contigSet[contigIndex].contigLength + extendLength;
           }   
           
       }else{
           contigSetHead->contigSet[contigIndex].contig = (char *)malloc(sizeof(char)*(maxSize+1));
           strncpy(contigSetHead->contigSet[contigIndex].contig, contig, extendLength);
           contigSetHead->contigSet[contigIndex].contig[extendLength] = '\0';
           contigSetHead->contigSet[contigIndex].contigLength = extendLength;
           allocateLength = maxSize + 1;
       }  
    }  
	contigSetHead->allContigLength = contigSetHead->allContigLength + contigSetHead->contigSet[contigIndex].contigLength;
    fflush(fp);
    fclose(fp);
	
	contigSetHead->visited = (bool *)malloc(sizeof(bool)*contigSetHead->contigCount);
	
	for(int i = 0; i < contigSetHead->contigCount; i++){
		contigSetHead->visited[i] = false;
		if(contigSetHead->contigSet[i].contigLength < contigLengthThreshold){
			contigSetHead->minContigCount++;
			contigSetHead->minAllContigLength =  contigSetHead->minAllContigLength + contigSetHead->contigSet[i].contigLength;
		}
	}
	
	if(contigSetHead->minContigCount > 0){
		contigSetHead->minContigSet = (Contig *)malloc(sizeof(Contig)*contigSetHead->minContigCount);
		for(int i = 0; i < contigSetHead->minContigCount; i++){
			contigSetHead->minContigSet[i].contig = NULL;
			contigSetHead->minContigSet[i].contigLength = 0;
			contigSetHead->minContigSet[i].shortContig = false;
			contigSetHead->minContigSet[i].realContigIndex = -1;
		}
		
		Contig * tempContigSet = (Contig *)malloc(sizeof(Contig)*(contigSetHead->contigCount - contigSetHead->minContigCount));
		for(int i = 0; i < contigSetHead->contigCount - contigSetHead->minContigCount; i++){
			tempContigSet[i].contig = NULL;
			tempContigSet[i].contigLength = 0;
			tempContigSet[i].shortContig = false;
			tempContigSet[i].realContigIndex = -1;
		}
		
		int j = 0;
		int t = 0;
		for(int i = 0; i < contigSetHead->contigCount; i++){
			if(contigSetHead->contigSet[i].contigLength < contigLengthThreshold){
				contigSetHead->minContigSet[j].contig = contigSetHead->contigSet[i].contig;
				contigSetHead->minContigSet[j].realContigIndex = i;
				contigSetHead->contigSet[i].contig = NULL;
				contigSetHead->minContigSet[j].contigLength =  contigSetHead->contigSet[i].contigLength;
				j++;
			}else{
				tempContigSet[t].contig = contigSetHead->contigSet[i].contig;
				tempContigSet[t].realContigIndex = i;
				contigSetHead->contigSet[i].contig = NULL;
				tempContigSet[t].contigLength =  contigSetHead->contigSet[i].contigLength;
				t++;
			}
		}
		
		free(contigSetHead->contigSet);
		contigSetHead->contigSet = tempContigSet;
		contigSetHead->contigCount = t;
		contigSetHead->allContigLength = contigSetHead->allContigLength - contigSetHead->minAllContigLength;
		
	}
	
	
    
    return contigSetHead;
}

void GetContigIDandPosition(ContigSetHead * contigSetHead, int tempPosition, int * P){
	int tempLength = 0;
	bool last = false;
	P[0] = -1;
	P[1] = -1;
	P[2] = -1;
	int i = 0;
	for(i = 0; i < contigSetHead->contigCount; i++){
		if(tempPosition < tempLength + contigSetHead->contigSet[i].contigLength + 1){
			P[0] = i;
			P[1] = tempPosition - tempLength;
			P[2] = 1;
			return;
		}else{
			tempLength = tempLength + contigSetHead->contigSet[i].contigLength + 1;
			last = true;
		}
	}
	
	if(last == true && i == contigSetHead->contigCount){
		for(int j = 0; j < contigSetHead->contigCount; j++){
			if(tempPosition < tempLength + contigSetHead->contigSet[j].contigLength + 1){
				P[0] = j;
				P[1] = tempPosition - tempLength;
				P[2] = 0;
				return;
			}else{
				tempLength = tempLength + contigSetHead->contigSet[j].contigLength + 1;
			}
		}
	}
	return;
}


void DeleteTailOfContigSet(ContigSetHead * contigSetHead, int tailLength){
	for(int i = 0; i < contigSetHead->contigCount; i++){
		if(contigSetHead->contigSet[i].contigLength > 10*tailLength){
			/*
			char * tempContig = (char *)malloc(sizeof(char)*(contigSetHead->contigSet[i].contigLength - 2*tailLength + 1));
			strncpy(tempContig, contigSetHead->contigSet[i].contig + tailLength, contigSetHead->contigSet[i].contigLength - 2*tailLength);
			tempContig[contigSetHead->contigSet[i].contigLength - 2*tailLength] = '\0';
			free(contigSetHead->contigSet[i].contig);
			contigSetHead->contigSet[i].contig = tempContig;
			tempContig = NULL;
			contigSetHead->contigSet[i].contigLength = contigSetHead->contigSet[i].contigLength - 2*tailLength;
			contigSetHead->allContigLength = contigSetHead->allContigLength - 2*tailLength;
			*/
			contigSetHead->contigSet[i].shortContig = false;
		}else{
			contigSetHead->contigSet[i].shortContig = true;
		}
	}
}


void DeleteArrayTwoElementList(ArrayTwoElementList * arrayTwoElementList){
	ArrayTwoElementList * last = NULL;
	ArrayTwoElementList * first = arrayTwoElementList;
	while(first != NULL){
		last = first->next;
		free(first);
		first = last;
	}
}
















































#endif