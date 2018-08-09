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
    
    long int maxSize = 90000;
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
    for(long int i = 0; i < contigSetHead->contigCount; i++){
        contigSetHead->contigSet[i].contig = NULL;
        contigSetHead->contigSet[i].contigLength = 0;
		contigSetHead->contigSet[i].shortContig = false;
		contigSetHead->contigSet[i].uniqueContig = false;
		contigSetHead->contigSet[i].realContigIndex = -1;
		contigSetHead->contigSet[i].repeativeContig = false;
		contigSetHead->contigSet[i].name = NULL;
		contigSetHead->repeatContigIndex[i] = false;
    }
    
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    
    long int allocateLength = 0;
    long int contigIndex = -1;
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
		   int len = strlen(contig);
		   contigSetHead->contigSet[contigIndex].name = (char *)malloc(sizeof(char)*len);
		   strncpy(contigSetHead->contigSet[contigIndex].name, contig + 1, len - 1);
		   contigSetHead->contigSet[contigIndex].name[len - 2] = '\0';
           continue;
       }
       
       
       long int extendLength = strlen(contig);
       if(contig[extendLength-1] == '\n'){
           extendLength--;
       }
       long int contigLength = 0;
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
	int num = 0;
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		contigSetHead->visited[i] = false;
		if(contigSetHead->contigSet[i].contigLength < contigLengthThreshold){
			contigSetHead->contigSet[i].shortContig = true;
			num++;
		}else{
			contigSetHead->contigSet[i].shortContig = false;
		}
	}
	
    return contigSetHead;
}

void SortContigSet(char * contigSetFile, char * sortContigSetFile){
    
	ContigSetHead * contigSetHead = GetContigSet(contigSetFile, 0);
	
	long int * sortIndex = (long int *)malloc(sizeof(long int)*contigSetHead->contigCount);
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		sortIndex[i] = -1;
	}
	int maxIndex = -1;
	int maxLength = -1;
	bool token = false;
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		maxIndex = -1;
		maxLength = -1;
		for(long int j = 0; j < contigSetHead->contigCount; j++){
			token = false;
			for(long int t = 0; t < contigSetHead->contigCount; t++){
				if(sortIndex[t] == j){
					token = true;
					break;
				}
				if(sortIndex[t] == -1){
					break;
				}
			}
			if(token == true){
				continue;
			}
			if(contigSetHead->contigSet[j].contigLength > maxLength){
				maxLength = contigSetHead->contigSet[j].contigLength;
				maxIndex = j;
			}
		}
		sortIndex[i] = maxIndex;
	}
	
	FILE * fp; 
    if((fp = fopen(sortContigSetFile, "w")) == NULL){
        printf("%s, does not exist!", sortContigSetFile);
        exit(0);
    }
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		fprintf(fp, ">%ld\n", i);
		fprintf(fp, "%s\n", contigSetHead->contigSet[sortIndex[i]].contig);
	}
	
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
		if(contigSetHead->contigSet[i].shortContig == false){
			
			char * tempContig = (char *)malloc(sizeof(char)*(contigSetHead->contigSet[i].contigLength - 2*tailLength + 1));
			strncpy(tempContig, contigSetHead->contigSet[i].contig + tailLength, contigSetHead->contigSet[i].contigLength - 2*tailLength);
			tempContig[contigSetHead->contigSet[i].contigLength - 2*tailLength] = '\0';
			free(contigSetHead->contigSet[i].contig);
			contigSetHead->contigSet[i].contig = tempContig;
			tempContig = NULL;
			contigSetHead->contigSet[i].contigLength = contigSetHead->contigSet[i].contigLength - 2*tailLength;
			contigSetHead->allContigLength = contigSetHead->allContigLength - 2*tailLength;
			//contigSetHead->contigSet[i].shortContig = false;
		}else{
			//contigSetHead->contigSet[i].shortContig = true;
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




char * ReverseComplement(char * temp){
    int len = strlen(temp);
	char * rcTemp = (char *)malloc(sizeof(char)*(len + 1));
    for(int i = 0; i < len; i++){
        if(temp[i] == 'A' || temp[i] == 'a'){
            rcTemp[len - 1 - i] = 'T';
        }else if(temp[i] == 'T' || temp[i] == 't'){
            rcTemp[len - 1 - i] = 'A';
        }else if(temp[i] == 'G' || temp[i] == 'g'){
            rcTemp[len - 1 - i] = 'C';
        }else if(temp[i] == 'C' || temp[i] == 'c'){
            rcTemp[len - 1 - i] = 'G';
        }else if(temp[i] == 'N' || temp[i] == 'n'){
            rcTemp[len - 1 - i] = 'N';
        } 
    }
    rcTemp[len]='\0';
    return rcTemp;
}


void DetectRepeativeContigInSet(ContigSetHead * contigSetHead, char * bamFileName, float ratio){

	
    long int contigIndex = 0;
	long int referenceIndex = -1;
	string readName;
	string previousReadName = "a";
    
    BamReader bamReader;
    bamReader.Open(bamFileName);
    BamAlignment alignment;
	
	long int contigCount = bamReader.GetReferenceCount();

    while(bamReader.GetNextAlignment(alignment)){
		readName = alignment.Name;
		if(previousReadName != "a" && readName != previousReadName){
			contigIndex++;
		}
		
		referenceIndex = alignment.RefID;
		if(referenceIndex == contigIndex){
			previousReadName = readName;
			continue;
		}
		const vector<CigarOp>& cigarData = alignment.CigarData;
		vector<CigarOp>::const_iterator cigarBegin = cigarData.begin();
		vector<CigarOp>::const_iterator cigarIter  = cigarBegin;
		vector<CigarOp>::const_iterator cigarEnd   = cigarData.end();

		double matchCount = 0;
		double unMatchCount = 0;
		for ( ; cigarIter != cigarEnd; ++cigarIter ) {
			const CigarOp& op = (*cigarIter);
			if(op.Type == 'M'){
				matchCount = matchCount + op.Length;
			}else{
				unMatchCount = unMatchCount + op.Length;
			}
		}
		
		uint16_t value;
		if(alignment.HasTag("NM")){
			if(alignment.GetTag("NM", value)){
				if(value != 0){
					previousReadName = readName;
					continue;
				}
			}
		}
		
		if(matchCount/(matchCount + unMatchCount) >= ratio && matchCount + unMatchCount > 0){
			contigSetHead->contigSet[contigIndex].repeativeContig = true;
		}
		previousReadName = readName;
		
    }


}











































#endif