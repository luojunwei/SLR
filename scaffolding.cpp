#ifndef SCAFFOLDING_CPP_INCLUDED 
#define SCAFFOLDING_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <typeinfo>

#include "scaffolding.h"
#include "scaffoldgraph.h"


ScaffoldSetHead * GetScaffoldSet(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead){
    
    int i = 0;
    int j = 0;
    
    bool * printIndex = (bool *)malloc(sizeof(bool)*scaffoldGraphHead->scaffoldGraphNodeCount);
	int * sortNode = (int *)malloc(sizeof(int)*scaffoldGraphHead->scaffoldGraphNodeCount);
    
    for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
        sortNode[i] = i;
		printIndex[i] = false;
    }
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
        for(int j = i; j < scaffoldGraphHead->scaffoldGraphNodeCount; j++){
			if(contigSetHead->contigSet[sortNode[i]].contigLength < contigSetHead->contigSet[sortNode[j]].contigLength){
				int temp = sortNode[i];
				sortNode[i] = sortNode[j];
				sortNode[j] = temp;
			}
    	}
    }
    
	ScaffoldGraph * scaffoldGraph = scaffoldGraphHead->scaffoldGraph;
    
    ScaffoldSetHead * scaffoldSetHead = (ScaffoldSetHead *)malloc(sizeof(ScaffoldSetHead));
	scaffoldSetHead->scaffoldSet = NULL;
    
    bool orientation = true;
    
    for(int p = 0; p < scaffoldGraphHead->scaffoldGraphNodeCount; p++){
        
		
		
		i = sortNode[p];
		
        if(printIndex[i] == true){
            continue;
        }
        
        ScaffoldGraphNode * temp = scaffoldGraph[i].outEdge;    
        
        ScaffoldSet * tempScaffoldSet = (ScaffoldSet *)malloc(sizeof(ScaffoldSet));
		tempScaffoldSet->next = NULL;
        if(scaffoldSetHead->scaffoldSet != NULL){
			tempScaffoldSet->next = scaffoldSetHead->scaffoldSet;
		}
		scaffoldSetHead->scaffoldSet = tempScaffoldSet;
        
        ContigSequence * tempContigSequence = (ContigSequence *)malloc(sizeof(ContigSequence));
        tempContigSequence->index = i;
        tempContigSequence->orientation = 1;
		tempContigSequence->gapDistance = 0;
		tempContigSequence->next = NULL;
        scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequence;
        j = i;

		
		if(temp == NULL){
			printIndex[i] = true;
		}
		orientation = 1;
        while(temp != NULL){
            printIndex[j] = true;
			
			bool uniq = false;
			
			temp = GetOptimizeNodeIndex(scaffoldGraphHead, j, orientation, scaffoldSetHead->scaffoldSet->contigSequence, 1, uniq);
			if(temp == NULL){
				break;
			}
			j = temp->nodeIndex; 
			if(printIndex[j] ==true){
				//break;
			}
			
			if((orientation == 1 && temp->orientation == 1) || (orientation == 0 && temp->orientation == 0)){
				orientation = 1;
			}else{
				orientation = 0;
			}
			
			
			int cc = SearchScaffoldEdge(temp->nodeIndex, scaffoldSetHead->scaffoldSet->contigSequence);
			if(cc > 0){
                printIndex[j] = true;
                break;
            }

            printIndex[j] = true;  
            
                        
            ContigSequence * tempContigSequence1 = (ContigSequence *)malloc(sizeof(ContigSequence));
            tempContigSequence1->index = temp->nodeIndex;
			tempContigSequence1->gapDistance = 0;
            tempContigSequence->gapDistance = temp->gapDistance;
            tempContigSequence->next = tempContigSequence1;
			tempContigSequence1->next = NULL;
            tempContigSequence = tempContigSequence1;
            tempContigSequence1 = NULL;    
			
			tempContigSequence->orientation = orientation;
			
			//cout<<"node3:"<<endl;
                                                               
        }
        
        temp = scaffoldGraph[i].inEdge;
        
        if(temp == NULL){
            continue;
        }
        

		
        orientation = 1;
        j = i;
        printIndex[j] = false;
		
        while(temp != NULL){

			printIndex[j] = true;
			
			bool uniq = false;
			
			temp = GetOptimizeNodeIndex(scaffoldGraphHead, j, orientation, scaffoldSetHead->scaffoldSet->contigSequence, 0, uniq);
			
			if(temp == NULL){
				break;
			}
			j = temp->nodeIndex;
			if(printIndex[j] ==true){
				//break;
			}
			
			if((orientation == 1 && temp->orientation == 1) || (orientation == 0 && temp->orientation == 0)){
				orientation = 1;
			}else{
				orientation = 0;
			}
			
            int cc = SearchScaffoldEdge(temp->nodeIndex, scaffoldSetHead->scaffoldSet->contigSequence);
            if(cc > 0){
                printIndex[j] = true;
                break;
            }
			
            ContigSequence * tempContigSequence1 = (ContigSequence *)malloc(sizeof(ContigSequence));
            tempContigSequence1->index = temp->nodeIndex;
            tempContigSequence1->gapDistance = temp->gapDistance;
            tempContigSequence1->next = scaffoldSetHead->scaffoldSet->contigSequence;
            scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequence1;
			
            scaffoldSetHead->scaffoldSet->contigSequence->orientation = orientation;
                                                     
        }
		
                
    }
    
    for(i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
        if(printIndex[i] != 1){
            ScaffoldSet * tempScaffoldSet = (ScaffoldSet *)malloc(sizeof(ScaffoldSet));
            tempScaffoldSet->next = scaffoldSetHead->scaffoldSet;
            scaffoldSetHead->scaffoldSet = tempScaffoldSet;
            
            ContigSequence * tempContigSequence = (ContigSequence * )malloc(sizeof(ContigSequence));
            tempContigSequence->index = i;
			tempContigSequence->gapDistance = 0;
            tempContigSequence->orientation = 1;
			tempContigSequence->next = NULL;
            scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequence;  
            
        }
    }
	char * tag = new char[20];
	strcpy(tag, "tag.fa");
	OutPutScaffoldTag(scaffoldSetHead->scaffoldSet, tag);
	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);
	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);
	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);
	strcpy(tag, "tag1.fa");
	OutPutScaffoldTag(scaffoldSetHead->scaffoldSet, tag);
	
    return scaffoldSetHead;
   
}


ScaffoldGraphNode * GetOptimizeNodeIndex(ScaffoldGraphHead * scaffoldGraphHead, int nodeIndex, bool orientation, ContigSequence * tempContigSequence, bool right, bool & uniq){
	ScaffoldGraphNode * temp = NULL;
	int edgeNum = 0;
	if((orientation == 1 && right == 1) || (orientation == 0 && right == 0)){
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge;
		while(temp != NULL){
			edgeNum++;
			temp = temp->next;
		}
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge;
	}else{
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge;
		while(temp != NULL){
			edgeNum++;
			temp = temp->next;
		}
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge;
		
	}
	
	
	if(edgeNum == 1){
		uniq = true;
		return temp;
	}
	
	int overlapEdgeNum[edgeNum];
	int index[edgeNum];
	for(int i = 0; i < edgeNum; i++){
		overlapEdgeNum[i] = -1;
		index[i] = -1;
	}
	for(int i = 0; i < edgeNum; i++){
		index[i] = temp->nodeIndex;
		
		for(int j = 0; j < temp->aligningReadCount; j++){
			int num = GetOverlapEdgeNum(scaffoldGraphHead, tempContigSequence, temp->readIndexArray[j], right);
			if(num > overlapEdgeNum[i]){
				overlapEdgeNum[i] = num;
			}
		}
		temp = temp->next;
	}
	
	int max = -1;
	int maxCount = 0;
	for(int i = 0; i < edgeNum; i++){
		if(max < overlapEdgeNum[i]){
			max = overlapEdgeNum[i];
		}
	}
	int returnIndex = -1;
	for(int i = 0; i < edgeNum; i++){
		if(max == overlapEdgeNum[i]){
			maxCount++;
			returnIndex = index[i];
		}
	}
	
	
	
	if(maxCount == 1){
		if((orientation == 1 && right == 1) || (orientation == 0 && right == 0)){
			temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge;
		}else{
			temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge;
		}
		
		while(temp != NULL){
			if(temp->nodeIndex == returnIndex){
				break;
			}
			temp = temp->next;
		}
		return temp;
	}else{

		return NULL;
	}
}

int GetOverlapEdgeNum(ScaffoldGraphHead * scaffoldGraphHead, ContigSequence * contigSequence, int readIndex, bool right){
	int resultEdgeNum = 0;
	ContigSequence * tempContigSequence = contigSequence;
	int contigSequenceNum = 0;
	while(tempContigSequence != NULL){
		contigSequenceNum++;
		tempContigSequence = tempContigSequence->next;
	}
	if(contigSequenceNum <=0){
		return 0;
	}
	if(contigSequenceNum == 1){
		return 1;
	}
	

	
	int contigIndex[contigSequenceNum];
	bool contigOrientation[contigSequenceNum];
	for(int i = 0; i < contigSequenceNum; i++){
		contigIndex[i] = -1;
		contigOrientation[i] = false;
	}

	tempContigSequence = contigSequence;
	for(int i = 0; i < contigSequenceNum; i++){
		contigIndex[i] = tempContigSequence->index;
		contigOrientation[i] = tempContigSequence->orientation;
		tempContigSequence = tempContigSequence->next;
	}
	
	ScaffoldGraphNode * temp = NULL;
	bool orientation = 1;

	if(right == false){
		for(int i = 1; i < contigSequenceNum; i++){

			orientation = contigOrientation[i];
			if(orientation == 1){
				temp = scaffoldGraphHead->scaffoldGraph[contigIndex[i]].inEdge;

			}else{
				temp = scaffoldGraphHead->scaffoldGraph[contigIndex[i]].outEdge;
			}
			if(temp == NULL){
				break;
			}
			while(temp != NULL){
				if(temp->nodeIndex == contigIndex[i - 1]){

					break;
				}
				temp = temp->next;
			}
			if(temp == NULL){
				break;
			}

			int * readIndexArray = temp->readIndexArray;
			bool a = false;
			for(int j = 0; j < temp->aligningReadCount; j++){
				if(readIndexArray[j] == readIndex){
					resultEdgeNum++;
					a = true;
					break;
				}
			}

			if(a == false){
				return resultEdgeNum;
			}
			
		}
	}
	if(right == true){
		for(int i = contigSequenceNum - 2; i >= 0; i--){
			orientation = contigOrientation[i];
			if(orientation == 1){
				temp = scaffoldGraphHead->scaffoldGraph[contigIndex[i]].outEdge;
			}else{
				temp = scaffoldGraphHead->scaffoldGraph[contigIndex[i]].inEdge;
			}
			if(temp == NULL){
				break;
			}
			while(temp != NULL){
				if(temp->nodeIndex == contigIndex[i + 1]){
					break;
				}
				temp = temp->next;
			}
			if(temp == NULL){
				break;
			}
			int * readIndexArray = temp->readIndexArray;
			bool a = false;
			for(int j = 0; j < temp->aligningReadCount; j++){
				if(readIndexArray[j] == readIndex){
					resultEdgeNum++;
					a = true;
					break;
				}
			}
			if(a == false){
				return resultEdgeNum;
			}
			
		}
	
	}
	
	return resultEdgeNum;

}

int GetContigSequenceNum(ContigSequence * tempContigSequence){

	int contigSequenceNum = 0;
	while(tempContigSequence != NULL){
		contigSequenceNum++;
		tempContigSequence = tempContigSequence->next;
	}
	return contigSequenceNum;

}

void OptimizeScaffoldSetCongtigSequence(ScaffoldSet * scaffoldSet, int contigNum){
	
	
	ScaffoldSet * tempScaffoldSet = scaffoldSet;
	ContigSequence * tempContigSequence = NULL;
	
	
	tempScaffoldSet = scaffoldSet;
	while(tempScaffoldSet != NULL){
		int contigSequenceNum = GetContigSequenceNum(tempScaffoldSet->contigSequence);
		tempScaffoldSet->contigNum = contigSequenceNum;
		if(contigSequenceNum > 0){
			tempScaffoldSet->contigIndex = (int *)malloc(sizeof(int)*contigSequenceNum);
		}else{
			tempScaffoldSet->contigIndex = NULL;
		}
		
		tempScaffoldSet = tempScaffoldSet->next;
	}
	
	tempScaffoldSet = scaffoldSet;
	while(tempScaffoldSet != NULL){
		tempContigSequence = tempScaffoldSet->contigSequence;
		int p = 0;
		while(tempContigSequence != NULL){
			tempScaffoldSet->contigIndex[p] = tempContigSequence->index;
			p++;
			tempContigSequence = tempContigSequence->next;
		}
		tempScaffoldSet = tempScaffoldSet->next;
	}
	
	tempScaffoldSet = scaffoldSet;
	int i = 0;
	int * overlapResult = (int *)malloc(sizeof(int)*6);//0:left-start-position, 1:left-end-position, 2:right-start-position,3:right-end-position,4:left-orientation,5:previous two node is right?
	while(tempScaffoldSet != NULL){
		ScaffoldSet * tempScaffoldSet1 = tempScaffoldSet->next;
		int j = 0;
		while(tempScaffoldSet1 != NULL){
			if(i == j){
				tempScaffoldSet1 = tempScaffoldSet1->next;
				j++;
				continue;
			}
			bool overlap = DeterminOverlapInScaffolds(tempScaffoldSet->contigIndex, tempScaffoldSet->contigNum, tempScaffoldSet1->contigIndex, tempScaffoldSet1->contigNum, overlapResult);
			if(overlap == true){
				if(overlapResult[0] == -2){
					tempScaffoldSet->contigSequence = NULL;
					tempScaffoldSet->contigIndex = NULL;
					tempScaffoldSet->contigNum = 0;
				}else if(overlapResult[2] == -2){
					tempScaffoldSet1->contigSequence = NULL;
					tempScaffoldSet1->contigIndex = NULL;
					tempScaffoldSet1->contigNum = 0;
				}else{
					MergeContigSequence(tempScaffoldSet, tempScaffoldSet1, overlapResult);
				}
				break;
			}
			tempScaffoldSet1 = tempScaffoldSet1->next;
			j++;
		}
		tempScaffoldSet = tempScaffoldSet->next;
		i++;
	}

}


bool DeterminOverlapInScaffolds(int * leftIndex, int leftNum, int * rightIndex, int rightNum, int * overlapResult){
	if(leftIndex == NULL || rightIndex == NULL){
		return false;
	}
	
	int overlapCount = 0;	
	overlapResult[0] = -1;
	
	for(int i = 0; i < leftNum; i++){
		for(int j = 0; j < rightNum; j++){
			if(leftIndex[i] == rightIndex[j]){
				overlapCount++;
				break;
			}
		}
	}
	if(overlapCount <= 0){
		overlapResult[0] = -1;
		return false;
	}
	
	GetOverlapBetweenArray(leftIndex, leftNum, rightIndex, rightNum, overlapResult);
	if(overlapResult[0] != -1){
		overlapResult[4] = 1;
		overlapResult[5] = 0;
		return true;
	}
	
	
	GetOverlapBetweenArray(rightIndex, rightNum, leftIndex, leftNum, overlapResult);
	if(overlapResult[0] != -1){
		if(overlapResult[0] == -2){
			overlapResult[0] = 0;
			overlapResult[2] = -2;
		}
		overlapResult[4] = 1;
		overlapResult[5] = 1;
		return true;
	}
	
	InverseArray(rightIndex, rightNum);
	GetOverlapBetweenArray(leftIndex, leftNum, rightIndex, rightNum, overlapResult);
	if(overlapResult[0] != -1){
		overlapResult[4] = 0;
		overlapResult[5] = 0;
		InverseArray(rightIndex, rightNum);
		return true;
	}
	
	GetOverlapBetweenArray(rightIndex, rightNum, leftIndex, leftNum, overlapResult);
	if(overlapResult[0] != -1){
		if(overlapResult[0] == -2){
			overlapResult[0] = 0;
			overlapResult[2] = -2;
		}
		overlapResult[4] = 0;
		overlapResult[5] = 1;
		InverseArray(rightIndex, rightNum);
		return true;
	}
	InverseArray(rightIndex, rightNum);
	return false;


}

bool GetOverlapBetweenArray(int * leftIndex, int leftNum, int * rightIndex, int rightNum, int * overlapResult){
	for(int i = 0; i < 6; i++){
		overlapResult[i] = -1;
	}
	
	if(leftNum > rightNum){
		for(int i = 0; i < leftNum; i++){
			if(leftIndex[i] == rightIndex[0]){
				int p = i;
				int j = 0; 
				for(j = 0; j < rightNum && p < leftNum; j++){
					if(rightIndex[j] != leftIndex[p]){
						break;
					}
					p++;
				}
				if(j == rightNum){
					overlapResult[2] = -2;
					return true;		
				}
			}
		}
	}else{
		for(int i = 0; i < rightNum; i++){
			if(rightIndex[i] == leftIndex[0]){
				int p = i;
				int j = 0; 
				for(j = 0; j < leftNum && p < rightNum; j++){
					if(leftIndex[j] != rightIndex[p]){
						break;
					}
					p++;
				}
				if(j == leftNum){
					overlapResult[0] = -2;
					return true;		
				}
			}
		}
	}
	
	for(int i = 0; i < leftNum; i++){
		if(leftIndex[i] == rightIndex[0]){
			for(int j = 0; j < rightNum; j++){
				if(rightIndex[j] == leftIndex[leftNum - 1]){
					if(leftNum - i == j + 1 && j + 1 >= 2){
						overlapResult[0] = i;
						overlapResult[1] = leftNum - 1;
						overlapResult[2] = 0;
						overlapResult[3] = j;
						return true;
					}
				}
			}
		}
	}
	
	
	
	
	
	return false;
}

void MergeContigSequence(ScaffoldSet * leftScaffoldSet, ScaffoldSet * rightScaffoldSet, int * overlapResult){
	
	if(overlapResult[5] == 1){
		if(overlapResult[4] == 0){
			ContigSequence * contigSequence = NULL;
			int gapDistance = 0;
			int i = 0;
			for(i = 0; i < overlapResult[0]; i++){
				contigSequence = GetContigSequenceIndex(rightScaffoldSet->contigNum - overlapResult[0] + i, rightScaffoldSet->contigSequence);
				if(i != 0){
					leftScaffoldSet->contigSequence->gapDistance = contigSequence->gapDistance;
				}
				if(contigSequence == NULL){
					return;
				}
				ContigSequence * tempContigSequence = (ContigSequence * )malloc(sizeof(ContigSequence));
            	tempContigSequence->index = contigSequence->index;
				tempContigSequence->gapDistance = 0;
            	tempContigSequence->orientation = !contigSequence->orientation;
				tempContigSequence->next = leftScaffoldSet->contigSequence;
				leftScaffoldSet->contigSequence = tempContigSequence;
			}
			contigSequence = GetContigSequenceIndex(i, rightScaffoldSet->contigSequence);
			leftScaffoldSet->contigSequence->gapDistance = contigSequence->gapDistance;
			rightScaffoldSet->contigSequence = NULL;
			rightScaffoldSet->contigIndex = NULL;
			rightScaffoldSet->contigNum = 0;
			
		}else{
			ContigSequence * contigSequence = GetContigSequenceIndex(overlapResult[0] - 1, rightScaffoldSet->contigSequence);
			if(contigSequence == NULL){
				return;
			}
			contigSequence->next = leftScaffoldSet->contigSequence;
			leftScaffoldSet->contigSequence = rightScaffoldSet->contigSequence;
			rightScaffoldSet->contigSequence = NULL;
			rightScaffoldSet->contigIndex = NULL;
			rightScaffoldSet->contigNum = 0;
		
		}
	
	}else{
		if(overlapResult[4] == 0){
			ContigSequence * contigSequence = NULL;
			int gapDistance = 0;
			
			int i = 0;
			for(i = 0; i < overlapResult[0]; i++){
				contigSequence = GetContigSequenceIndex(overlapResult[0] - i - 1, leftScaffoldSet->contigSequence);
				if(i != 0){
					rightScaffoldSet->contigSequence->gapDistance = contigSequence->gapDistance;
				}
				if(contigSequence == NULL){
					return;
				}
				ContigSequence * tempContigSequence = (ContigSequence * )malloc(sizeof(ContigSequence));
            	tempContigSequence->index = contigSequence->index;
				tempContigSequence->gapDistance = 0;
            	tempContigSequence->orientation = !contigSequence->orientation;
				tempContigSequence->next = rightScaffoldSet->contigSequence;
				rightScaffoldSet->contigSequence = tempContigSequence;
			}
			contigSequence = GetContigSequenceIndex(i, leftScaffoldSet->contigSequence);
			rightScaffoldSet->contigSequence->gapDistance = contigSequence->gapDistance;
			leftScaffoldSet->contigSequence = NULL;
			leftScaffoldSet->contigIndex = NULL;
			leftScaffoldSet->contigNum = 0;
			
		}else{
			ContigSequence * contigSequence = GetContigSequenceIndex(overlapResult[0] - 1, leftScaffoldSet->contigSequence);
			if(contigSequence == NULL){
				return;
			}
			contigSequence->next = rightScaffoldSet->contigSequence;
			rightScaffoldSet->contigSequence = NULL;
			rightScaffoldSet->contigIndex = NULL;
			rightScaffoldSet->contigNum = 0;
		
		}
	
	}
	
	
}

void InverseArray(int * array, int num){
	for(int i = 0; i < num/2; i++){
		int temp = array[i];
		array[i] = array[num -1 - i];
		array[num -1 - i] = temp;
	}
}

ContigSequence * GetContigSequenceIndex(int index, ContigSequence * contigSequence){
    int i = 0;
    while(contigSequence != NULL){
        if(i == index){
            return contigSequence;
        }
		i++;
        contigSequence = contigSequence->next;
    }
    
    return NULL;
    
}


int SearchScaffoldEdge(int index, ContigSequence * contigSequence){
    int num = 0;
    while(contigSequence != NULL){
        if(contigSequence->index == index){
            num ++;
        }
        contigSequence = contigSequence->next;
    }
    
    return num;
    
}





void OutPutScaffoldTag(ScaffoldSet * scaffoldSet, char * result){
    
    ofstream ocout;
    ocout.open(result);
    
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    int i = 0;
    int j = 0;
    while(tempScaffoldSet != NULL){
        
        ContigSequence * temp = tempScaffoldSet->contigSequence;
        while(temp != NULL){
            
            ocout<<temp->index<<",";
            temp = temp->next;
        }
        if(tempScaffoldSet->contigSequence != NULL){
            ocout<<endl;
        }
        tempScaffoldSet = tempScaffoldSet->next;
    }
    
    ocout.close();
    
}

void OutPutScaffoldSet(ScaffoldSet * scaffoldSet, ContigSetHead * contigSetHead, char * result){    
    int i = 0;
    int j = 0;

	Contig * contigSet = contigSetHead->contigSet;
	int contigCount = contigSetHead->contigCount;
	
	fflush(stdout);  
    setvbuf(stdout,NULL,_IONBF,0);
	
    bool * printContigIndex = new bool[contigCount];
    for(i = 0; i < contigCount; i++){
        printContigIndex[i] = false;
    }
	
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    scaffoldSet = tempScaffoldSet;
    
    char * scaffoldTagFileName = new char[1000];
    strcpy(scaffoldTagFileName, result);
    strcat(scaffoldTagFileName, "_Scaffold_Tag.fa");
    ofstream ocoutTag;
    ocoutTag.open(scaffoldTagFileName);
    
    char * scaffoldSetFileName = new char[1000];   
    strcpy(scaffoldSetFileName, result);
    strcat(scaffoldSetFileName, "_ScaffoldSet.fa");
    ofstream ocout1;
    ocout1.open(scaffoldSetFileName);
    j = 0;
    while(tempScaffoldSet != NULL){
        
		ContigSequence * tempContigSequence = tempScaffoldSet->contigSequence;
        if(tempContigSequence == NULL){
            ocoutTag<<endl;
			tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        ocout1<<">"<<j<<endl;
        
        int allLength = 0;
        
        while(tempContigSequence != NULL){
            
            if(printContigIndex[tempContigSequence->index] == true){
                ocoutTag<<"--";
            }
            
            printContigIndex[tempContigSequence->index] = true;
            ocoutTag<<tempContigSequence->index<<"("<<tempContigSequence->gapDistance<<"--"<<tempContigSequence->orientation<<"),";
            
            if(tempContigSequence->orientation==0){
                char * tempChar1 = ReverseComplement(contigSet[tempContigSequence->index].contig);
				if(tempContigSequence->gapDistance<=0 && tempContigSequence->next!=NULL){
					
					int contigLength = strlen(tempChar1);
					if(contigLength + tempContigSequence->gapDistance < 0){
						tempContigSequence = tempContigSequence->next;
						continue;
					}
					
                	char * tempChar2 = (char *)malloc(sizeof(char)*(contigLength + tempContigSequence->gapDistance + 1));
					
					strncpy(tempChar2, tempChar1, contigLength + tempContigSequence->gapDistance);
					
					tempChar2[contigLength + tempContigSequence->gapDistance] = '\0';
					
                	free(tempChar1);
					
					tempChar1 = tempChar2;
            	}
                ocout1<<tempChar1;
				free(tempChar1);
            }else{
				if(tempContigSequence->gapDistance<=0 && tempContigSequence->next!=NULL){
                	
					int contigLength = strlen(contigSet[tempContigSequence->index].contig);
					
					if(contigLength + tempContigSequence->gapDistance < 0){
						tempContigSequence = tempContigSequence->next;
						continue;	
					}
					
					char * tempChar2 = (char *)malloc(sizeof(char)*(contigLength + tempContigSequence->gapDistance + 1));
					strncpy(tempChar2, contigSet[tempContigSequence->index].contig, contigLength + tempContigSequence->gapDistance);
					
					tempChar2[contigLength + tempContigSequence->gapDistance] = '\0';
					
					ocout1<<tempChar2;
            	}else{
					ocout1<<contigSet[tempContigSequence->index].contig;
				}
                
            }
            
            if(tempContigSequence->gapDistance>0 && tempContigSequence->next!=NULL){
                for(int tt = 0; tt<tempContigSequence->gapDistance; tt++){
                    ocout1<<"N";
                }
            }
            tempContigSequence = tempContigSequence->next;
        }
        
        ocoutTag<<endl;
        ocout1<<endl;
        j++;
        tempScaffoldSet = tempScaffoldSet->next;
    }
    
    ocoutTag<<"----------------------------------------------------------------"<<endl;
    for(i = 0; i<contigCount; i++){
        if(printContigIndex[i] == false && contigSet[i].contig != NULL){
            ocout1<<">"<<j<<endl;
            ocout1<<contigSet[i].contig<<endl;
            ocoutTag<<i<<","<<endl;
            j++;
        }
    }
	
	Contig * minContigSet = contigSetHead->minContigSet;
	int minContigCount = contigSetHead->minContigCount;
	
	for(i = 0; i < minContigCount; i++){
        ocout1<<">min-"<<j<<endl;
        ocout1<<minContigSet[i].contig<<endl;
        ocoutTag<<i<<","<<endl;
        j++;
    }

}





#endif
