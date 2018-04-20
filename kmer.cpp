#ifndef Kmer_CPP_INCLUDED 
#define Kmer_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "contig.h"
#include "kmer.h"

long int Hash(char * str, unsigned int len, unsigned long int max)  
{  
   unsigned int hash = 0;  
   unsigned int i = 0;  
  
   for(i = 0; i < len; str++, i++) {  
      hash = (*str) + (hash << 6) + (hash << 16) - hash;  
   }  
  
   return hash % max;  
} 

int SearchKmerInHashTable(ContigSetHead * contigSetHead, KmerHashNodeHead * kmerHashNodeHead, long int kmerHashTableCount, char * kmer, int kmerLength, char * tempKmer){
	long int hashIndex = Hash(kmer, kmerLength, kmerHashTableCount);
	while(true){
		if(kmerHashNodeHead[hashIndex].kmerHashNode == NULL){
			return -1;
		}else{
			strncpy(tempKmer, contigSetHead->contigSet[kmerHashNodeHead[hashIndex].kmerHashNode->index].contig + kmerHashNodeHead[hashIndex].kmerHashNode->position, kmerLength);
			tempKmer[kmerLength] = '\0';
			if(strcmp(tempKmer, kmer) == 0){
				return hashIndex;
			}
			hashIndex = (hashIndex + 1) % kmerHashTableCount;
		}
	}
}


void InsertKmerToHashTable(ContigSetHead * contigSetHead, KmerHashTableHead * kmerHashTableHead, char * kmer, int kmerLength, long int index, long int position, char * tempKmer){
	int kmerHashTableCount = kmerHashTableHead->kmerHashTableCount;
	KmerHashNodeHead * kmerHashNodeHead = kmerHashTableHead->kmerHashNodeHead;
	long int hashIndex = Hash(kmer, kmerLength, kmerHashTableCount);
	KmerHashNode * tempKmerHashNode = (KmerHashNode *)malloc(sizeof(KmerHashNode));
	tempKmerHashNode->index = index;
	tempKmerHashNode->position = position;
	tempKmerHashNode->next = NULL;
	while(true){
		if(kmerHashNodeHead[hashIndex].kmerHashNode == NULL){
			kmerHashNodeHead[hashIndex].kmerHashNode = tempKmerHashNode;
			break;
		}else{
			strncpy(tempKmer, contigSetHead->contigSet[kmerHashNodeHead[hashIndex].kmerHashNode->index].contig + kmerHashNodeHead[hashIndex].kmerHashNode->position, kmerLength);
			tempKmer[kmerLength] = '\0';
			if(strcmp(tempKmer, kmer) == 0){
				tempKmerHashNode->next = kmerHashNodeHead[hashIndex].kmerHashNode;
				kmerHashNodeHead[hashIndex].kmerHashNode = tempKmerHashNode;
				break;
			}
			hashIndex = (hashIndex + 1) % kmerHashTableCount;
		}
	}

}


KmerHashTableHead * GetKmerHashTableHead(ContigSetHead * contigSetHead, int kmerLength){
	
	KmerHashTableHead * kmerHashTableHead = (KmerHashTableHead * )malloc(sizeof(KmerHashTableHead));
	
	kmerHashTableHead->kmerHashTableCount = (contigSetHead->allContigLength - (kmerLength - 1)*contigSetHead->contigCount)*1.2;
	kmerHashTableHead->kmerHashNodeHead = (KmerHashNodeHead *)malloc(sizeof(KmerHashNodeHead)*kmerHashTableHead->kmerHashTableCount);
	kmerHashTableHead->kmerHashNodeCount = 0;
	for(int i = 0; i < kmerHashTableHead->kmerHashTableCount; i++){
		kmerHashTableHead->kmerHashNodeHead[i].kmerHashNode = NULL;
	}
	
	char * kmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * tempKmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
	
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		for(long int j = 0; j < contigSetHead->contigSet[i].contigLength - kmerLength + 1; j++){
			strncpy(kmer, contigSetHead->contigSet[i].contig + j, kmerLength);
			kmer[kmerLength] = '\0';
			InsertKmerToHashTable(contigSetHead, kmerHashTableHead, kmer, kmerLength, i, j, tempKmer);
		}
	}
	
	return kmerHashTableHead;
}


void ReverseComplementKmer(char * kmer, int kmerLength){
	for(int i = 0; i < kmerLength/2; i++){
		char temp = kmer[i];
		kmer[i] = kmer[kmerLength -1 - i];
		kmer[kmerLength -1 - i] = temp;
	}
	
	for(int i = 0; i < kmerLength; i++){
		if(kmer[i] == 'A' || kmer[i] == 'a'){
			kmer[i] = 'T';
		}else if(kmer[i] == 'T' || kmer[i] == 't'){
			kmer[i] = 'A';
		}else if(kmer[i] == 'G' || kmer[i] == 'g'){
			kmer[i] = 'C';
		}else if(kmer[i] == 'C' || kmer[i] == 'c'){
			kmer[i] = 'G';
		}else if(kmer[i] == 'N' || kmer[i] == 'n'){
			kmer[i] = 'N';
		}
	}
	
	
	
}

































#endif