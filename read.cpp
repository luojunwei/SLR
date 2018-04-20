#ifndef READ_CPP_INCLUDED 
#define READ_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "read.h"
#include "kmer.h"

using namespace std;


void GetReadToKmer(ContigSetHead * contigSetHead, KmerHashTableHead * kmerHashTableHead, Read * read, int kmerLength, int step, char * kmer, char * tempKmer, FILE * fp, int readStart){
	int readLength = strlen(read->read);
	char * readS = read->read;
	
	int token = false;

	for(int i = 0; i < readLength - kmerLength + 1; i = i + step){
		strncpy(kmer, readS + i, kmerLength);
		kmer[kmerLength] = '\0';
		int locs = SearchKmerInHashTable(contigSetHead, kmerHashTableHead->kmerHashNodeHead, kmerHashTableHead->kmerHashTableCount, kmer, kmerLength, tempKmer);

		if(locs != -1){
			KmerHashNode * tempKmerHashNode = kmerHashTableHead->kmerHashNodeHead[locs].kmerHashNode;
			if(token == false){
				fprintf(fp, "read--%d--%d\n", readStart, readLength);
				token = true;
			}
			while(tempKmerHashNode != NULL){
				fprintf(fp, "%d--%d--%d--1\n", i, tempKmerHashNode->index, tempKmerHashNode->position);
				tempKmerHashNode = tempKmerHashNode->next;
			}
			
			
		}
		
		ReverseComplementKmer(kmer, kmerLength);
		locs = SearchKmerInHashTable(contigSetHead, kmerHashTableHead->kmerHashNodeHead, kmerHashTableHead->kmerHashTableCount, kmer, kmerLength, tempKmer);	
		if(locs != -1){	
			KmerHashNode * tempKmerHashNode = kmerHashTableHead->kmerHashNodeHead[locs].kmerHashNode;
			if(token == false){
				fprintf(fp, "read--%d--%d\n", readStart, readLength);
				token = true;
			}
			while(tempKmerHashNode != NULL){
				fprintf(fp, "%d--%d--%d--0\n", i, tempKmerHashNode->index, contigSetHead->contigSet[tempKmerHashNode->index].contigLength - tempKmerHashNode->position - kmerLength);
				tempKmerHashNode = tempKmerHashNode->next;
			}
			
		}
	}

}

void GetReadSetToKmer(Read * readSet, int readCount, ContigSetHead * contigSetHead, KmerHashTableHead * kmerHashTableHead, int kmerLength, int step, int readLengthCutOff, char * kmer, char * tempKmer, int readStart, char * file, char * line){
	
	FILE * fp; 
    if((fp = fopen(file, "w")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	for(int i = 0; i < readCount; i++){
		if(readSet[i].read != NULL){
			if(strlen(readSet[i].read) > readLengthCutOff){
				GetReadToKmer(contigSetHead, kmerHashTableHead, &readSet[i], kmerLength, step, kmer, tempKmer, fp, readStart + i);
			}
		}
	}
	
	
	fflush(fp);
}


void GetReadSet(char * readLibraryFileName, int startIndex, ReadSetHead * readSetHead, char * read, int maxSize){
	
	
	startIndex = startIndex*(readSetHead->readCount);
    
    FILE * fp; 
    if(NULL == (fp = fopen(readLibraryFileName, "r"))){
        printf("%s, does not exist!", readLibraryFileName);
        exit(0);
    }
   
    int readIndex = -1;
	
	while((fgets(read, maxSize, fp)) != NULL){ 
       
		if(read[0] == '>'){  
			if(strlen(read) == maxSize-1){              
               while((fgets(read, maxSize, fp)) != NULL){
                   if(strlen(read) != maxSize-1){
                       break;
                   }
               }        
			}
		   
			readIndex++;
			continue;
		}
		
		
		if(readIndex < startIndex){
			continue;
		}else if(readIndex >= startIndex + readSetHead->readCount){
			break;
		}
       
		int extendLength = strlen(read);
		if(read[extendLength-1] == '\n'){
			extendLength--;
		}
       
		if(readSetHead->readSet[readIndex - startIndex].readLength + extendLength >= readSetHead->readSet[readIndex - startIndex].allocateReadLength){
			readSetHead->readSet[readIndex - startIndex].read = (char *)realloc(readSetHead->readSet[readIndex - startIndex].read, readSetHead->readSet[readIndex - startIndex].allocateReadLength + maxSize + 1);
			readSetHead->readSet[readIndex - startIndex].allocateReadLength = readSetHead->readSet[readIndex - startIndex].allocateReadLength + maxSize + 1;
               
			strncpy(readSetHead->readSet[readIndex - startIndex].read + readSetHead->readSet[readIndex - startIndex].readLength, read, extendLength);
			readSetHead->readSet[readIndex - startIndex].read[readSetHead->readSet[readIndex - startIndex].readLength + extendLength] = '\0';    
			readSetHead->readSet[readIndex - startIndex].readLength = readSetHead->readSet[readIndex - startIndex].readLength + extendLength;
                       
		}else{
			strncpy(readSetHead->readSet[readIndex - startIndex].read + readSetHead->readSet[readIndex - startIndex].readLength, read, extendLength);
			readSetHead->readSet[readIndex - startIndex].read[readSetHead->readSet[readIndex - startIndex].readLength + extendLength] = '\0';
			readSetHead->readSet[readIndex - startIndex].readLength = readSetHead->readSet[readIndex - startIndex].readLength + extendLength;
		}   
	}      
       
	fflush(fp);
	fclose(fp);

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



























#endif
