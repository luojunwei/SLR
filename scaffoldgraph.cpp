#ifndef Scaffoldgraph_CPP_INCLUDED 
#define Scaffoldgraph_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>



#include "scaffoldgraph.h"


using namespace std;


void InsertOutOrInEdge(ScaffoldGraphHead * scaffoldGraphHead, int readIndex, int leftNodeIndex, bool leftOrientation, int rightNodeIndex, bool rightOrientation, int gapDistance, int overlapLength){
	
	ScaffoldGraphNode * leftNode = (ScaffoldGraphNode *)malloc(sizeof(ScaffoldGraphNode));
	leftNode->nodeIndex = rightNodeIndex;
	leftNode->gapDistance = gapDistance;
	leftNode->overlapLength = overlapLength;
	leftNode->readIndex = readIndex;
	leftNode->readIndexArray = NULL;
	leftNode->aligningReadCount = 0;
	leftNode->next = NULL;
	ScaffoldGraphNode * rightNode = (ScaffoldGraphNode *)malloc(sizeof(ScaffoldGraphNode));
	rightNode->nodeIndex = leftNodeIndex;
	rightNode->gapDistance = gapDistance;
	rightNode->overlapLength = overlapLength;
	rightNode->readIndex = readIndex;
	rightNode->readIndexArray = NULL;
	rightNode->aligningReadCount = 0;
	rightNode->next = NULL;
	
	if(leftOrientation != rightOrientation){
		leftNode->orientation = 0;
		rightNode->orientation = 0;
	}else{
		leftNode->orientation = 1;
		rightNode->orientation = 1;
	}
	
	if(leftOrientation == true){
		if(scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge != NULL){
			leftNode->next = scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge;
		}
		scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge = leftNode;
	}else{
		if(scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge != NULL){
			leftNode->next = scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge;
		}
		scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge = leftNode;
	}
	
	if(rightOrientation == true){
		if(scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge != NULL){
			rightNode->next = scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge;
		}
		scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge = rightNode;
	}else{
		if(scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge != NULL){
			rightNode->next = scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge;
		}
		scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge = rightNode;
	}
	
}


bool InsertSingleEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, int readIndex, AligningResult * left, bool leftOrientation, AligningResult * right, bool rightOrientation, int readLength, ContigSetHead * contigSetHead){
	
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
			InsertOutOrInEdge(scaffoldGraphHead, readIndex, left->contigIndex, leftOrientation, right->contigIndex, rightOrientation, gapDistance, min);
			return true;
		}else if(gapDistance1 >= 0){
			InsertOutOrInEdge(scaffoldGraphHead, readIndex, right->contigIndex, rightOrientation, left->contigIndex, leftOrientation, gapDistance1, min);
			return true;
		}else{
			
			if(left->readStartPosition < right->readStartPosition && left->readEndPosition < right->readEndPosition){
				InsertOutOrInEdge(scaffoldGraphHead, readIndex, left->contigIndex, leftOrientation, right->contigIndex, rightOrientation, gapDistance, min);
				return true;
			}else if(left->readStartPosition > right->readStartPosition && left->readEndPosition > right->readEndPosition){
				InsertOutOrInEdge(scaffoldGraphHead, readIndex, right->contigIndex, rightOrientation, left->contigIndex, leftOrientation, gapDistance1, min);
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
			InsertOutOrInEdge(scaffoldGraphHead, readIndex, left->contigIndex, leftOrientation, right->contigIndex, rightOrientation, gapDistance, min);
			return true;
		}else if(gapDistance1 >= 0){
			InsertOutOrInEdge(scaffoldGraphHead, readIndex, right->contigIndex, rightOrientation, left->contigIndex, leftOrientation, gapDistance1, min);
			return true;
		}else{
			if(left->readStartPosition < right->readStartPosition && left->readEndPosition < right->readEndPosition){
				InsertOutOrInEdge(scaffoldGraphHead, readIndex, left->contigIndex, leftOrientation, right->contigIndex, rightOrientation, gapDistance, min);
				return true;
			}else if(left->readStartPosition > right->readStartPosition && left->readEndPosition > right->readEndPosition){
				InsertOutOrInEdge(scaffoldGraphHead, readIndex, right->contigIndex, rightOrientation, left->contigIndex, leftOrientation, gapDistance1, min);
				return true;
			}
		}
		
	}
	return false;
	
	
}

void InsertEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, AligningResult * aligningResult, AligningResult * rcAligningResult, int readIndex, int readLength, ContigSetHead * contigSetHead){
	
	
	AligningResult * tempLeft = aligningResult;
	AligningResult * tempRight = NULL;
	
	int startIndex = 0;
	
	while(tempLeft != NULL){
		tempRight = tempLeft->next;
		if(tempRight != NULL){		

			InsertSingleEdgeToScaffoldGraph(scaffoldGraphHead, readIndex, tempLeft, tempLeft->orientation, tempRight, tempRight->orientation, readLength, contigSetHead);

		}
		tempLeft = tempLeft->next;
	}
}


void OutputScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead){
	ScaffoldGraphNode * temp = NULL;
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		if(scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL){
			temp = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			while(temp != NULL){
				cout<<"outnodeIndex:"<<i<<",index:"<<temp->nodeIndex<<",gapDistance:"<<temp->gapDistance<<",ori:"<<temp->orientation<<",count:"<<temp->aligningReadCount<<",overlapLength:"<<temp->overlapLength<<endl;
				cout<<"readIndex:";
				for(int j = 0; j < temp->aligningReadCount; j++){
					cout<<"--"<<temp->readIndexArray[j];
				}
				cout<<endl;
				
				temp = temp->next;
			}
		}
		
		if(scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL){
			temp = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			while(temp != NULL){
				cout<<"innodeIndex:"<<i<<",index:"<<temp->nodeIndex<<",gapDistance:"<<temp->gapDistance<<",ori:"<<temp->orientation<<",count:"<<temp->aligningReadCount<<",overlapLength:"<<temp->overlapLength<<endl;
				cout<<"readIndex:";
				for(int j = 0; j < temp->aligningReadCount; j++){
					cout<<"--"<<temp->readIndexArray[j];
				}
				cout<<endl;
				
				temp = temp->next;
			}
		}
	}
}


void DeleteScaffoldGraphNode(ScaffoldGraphNode * edge){
	ScaffoldGraphNode * tempEdge = edge;
	while(edge != NULL){
		tempEdge = edge->next;
		edge->next = NULL;
		free(edge);
		edge = tempEdge;
	}
}

ScaffoldGraphNode * DeleteScaffoldGraphSpecialNode(ScaffoldGraphNode * edge, int nodeIndex, bool orientation){
	ScaffoldGraphNode * firstEdge = edge;
	ScaffoldGraphNode * previousEdge = NULL;

	while(edge != NULL){
		if(edge->nodeIndex == nodeIndex && edge->orientation == orientation){
			if(previousEdge == NULL){
				firstEdge = edge->next;
			}else{
				previousEdge->next = edge->next;
			}
			edge->next = NULL;
			free(edge);
			break;
		}

		previousEdge = edge;
		edge = edge->next;
	}
	return firstEdge;
}

ScaffoldGraphNode * MergeMultipleEdges(ScaffoldGraphNode * edge, int contigLength, int contigLength0){
	
	int forwardCount0 = 0;
	double forwardGapDistance0 = 0;
	int reverseCount0 = 0;
	int reverseGapDistance0 = 0;
	
	int forwardCount1 = 0;
	double forwardGapDistance1 = 0;
	int reverseCount1 = 0;
	int reverseGapDistance1 = 0;
	
	int forwardOverlapLength = 0;
	int reverseOverlapLength = 0;
	
	ScaffoldGraphNode * result = NULL;
	
	ScaffoldGraphNode * tempEdge = edge;
	while(tempEdge != NULL){
		if(tempEdge->orientation == 1){
			forwardCount1++;
			forwardGapDistance1 = forwardGapDistance1 + tempEdge->gapDistance;
			if(tempEdge->overlapLength > forwardOverlapLength){
				forwardOverlapLength = tempEdge->overlapLength;
			}
			
		}else{
			reverseCount1++;
			reverseGapDistance1 = reverseGapDistance1 + tempEdge->gapDistance;
			if(tempEdge->overlapLength > reverseOverlapLength){
				reverseOverlapLength = tempEdge->overlapLength;
			}

		}
		tempEdge = tempEdge->next;
	}
	
	if(forwardCount0 > 0){
		forwardGapDistance0 = forwardGapDistance0/forwardCount0;
	}
	if(reverseCount0 > 0){
		reverseGapDistance0 = reverseGapDistance0/reverseCount0;
	}
	if(forwardCount1 > 0){
		forwardGapDistance1 = forwardGapDistance1/forwardCount1;
	}
	if(reverseCount1 > 0){
		reverseGapDistance1 = reverseGapDistance1/reverseCount1;
	}
	
	int max = -1;
	int maxCount = -1;
	
	if(forwardCount1 > 0){
		ScaffoldGraphNode * result1 = (ScaffoldGraphNode *)malloc(sizeof(ScaffoldGraphNode));
		result1->nodeIndex = edge->nodeIndex;
		result1->gapDistance = forwardGapDistance1;
		result1->aligningReadCount = forwardCount1;
		result1->orientation = true;
		result1->overlapLength = forwardOverlapLength;
		result1->next = NULL;
		result1->readIndexArray = (int *)malloc(sizeof(int)*forwardCount1);
		tempEdge = edge;
		forwardCount1 = 0;
		while(tempEdge != NULL){
			if(tempEdge->orientation == 1){
				result1->readIndexArray[forwardCount1] = tempEdge->readIndex;
				forwardCount1++;
			}
			tempEdge = tempEdge->next;
		}
		result = result1;
	}
	if(reverseCount1 > 0){
		ScaffoldGraphNode * result1 = (ScaffoldGraphNode *)malloc(sizeof(ScaffoldGraphNode));
		result1->nodeIndex = edge->nodeIndex;
		result1->gapDistance = reverseGapDistance1;
		result1->aligningReadCount = reverseCount1;
		result1->orientation = false;
		result1->overlapLength = reverseOverlapLength;
		result1->next = NULL;
		result1->readIndexArray = (int *)malloc(sizeof(int)*reverseCount1);
		tempEdge = edge;
		reverseCount1 = 0;
		while(tempEdge != NULL){
			if(tempEdge->orientation == 0){
				result1->readIndexArray[reverseCount1] = tempEdge->readIndex;
				reverseCount1++;
			}
			tempEdge = tempEdge->next;
		}
		if(result != NULL){
			result1->next = result;
		}
		result = result1;
	}
	DeleteScaffoldGraphNode(edge);
	return result;
	
}




ScaffoldGraphNode *  OptimizeEdgeInScaffoldGraph(ScaffoldGraphNode * edge, ContigSetHead * cotigSetHead, int nodeIndex){
	
	
	ScaffoldGraphNode * previousEdge = NULL;
	ScaffoldGraphNode * previousEdge1 = NULL;
	
	ScaffoldGraphNode * result = NULL;

	
	ScaffoldGraphNode * tempStart = edge;
	
	while(tempStart != NULL){
		
		ScaffoldGraphNode * tempEdge = tempStart;
		tempStart = NULL;
		previousEdge1 = NULL;
		ScaffoldGraphNode * tempEdge1 = tempEdge->next;
		previousEdge = tempEdge;
		bool a = false;
		while(tempEdge1 != NULL){
			if(tempEdge1->nodeIndex == tempEdge->nodeIndex){
				previousEdge->next = tempEdge1;
				previousEdge = tempEdge1;
				
			}else if(a != true){
				tempStart = tempEdge1;
				previousEdge1 = tempStart;
				a = true;
			}else{
				previousEdge1->next = tempEdge1;
				previousEdge1 = tempEdge1;
			}
			tempEdge1 = tempEdge1->next;
		}
		if(previousEdge1 != NULL){
			previousEdge1->next = NULL;
		}
		if(previousEdge != NULL){
			previousEdge->next = NULL;
		}

		tempEdge = MergeMultipleEdges(tempEdge, cotigSetHead->contigSet[nodeIndex].contigLength, cotigSetHead->contigSet[tempEdge->nodeIndex].contigLength);

		if(tempEdge != NULL){
			if(result != NULL){
				ScaffoldGraphNode * tempEdge2 = tempEdge->next;
				while(tempEdge2 != NULL){
				    if(tempEdge2->next == NULL){
						break;
					}
				}
				if(tempEdge2 != NULL){
					tempEdge2->next = result;
				}else{
					tempEdge->next = result;
				}
				
			}
			result = tempEdge;
		}
		
	}
	return result;
	
}

bool VisitOutNode(ScaffoldGraphHead * scaffoldGraphHead, int leftIndex, bool out, int rightIndex, bool * visited){
    if(visited[leftIndex] == true){
		return false;
	}
	visited[leftIndex] = true;
	ScaffoldGraphNode * tempEdge = NULL;
	if(out == true){
		tempEdge = scaffoldGraphHead->scaffoldGraph[leftIndex].outEdge;
	}else{
		tempEdge = scaffoldGraphHead->scaffoldGraph[leftIndex].inEdge;
	}
	while(tempEdge != NULL){
		if(tempEdge->nodeIndex == rightIndex){
			return true;
		}
		bool tempOut = false;
		if(out == true){
			tempOut = tempEdge->orientation;
		}else{
			tempOut = !tempEdge->orientation;
		}
		bool result = VisitOutNode(scaffoldGraphHead, tempEdge->nodeIndex, tempOut, rightIndex, visited);
		if(result == true){
			return true;
		}
		tempEdge = tempEdge->next;
	}
	return false;   
}


bool VisitInNode(ScaffoldGraphHead * scaffoldGraphHead, int rightIndex, bool in, int leftIndex, bool * visited){
    if(visited[rightIndex] == true){
		return false;
	}
	visited[rightIndex] = true;
	ScaffoldGraphNode * tempEdge = NULL;
	if(in == true){
		tempEdge = scaffoldGraphHead->scaffoldGraph[rightIndex].inEdge;
	}else{
		tempEdge = scaffoldGraphHead->scaffoldGraph[rightIndex].outEdge;
	}
	while(tempEdge != NULL){
		if(tempEdge->nodeIndex == leftIndex){
			return true;
		}
		bool tempIn = false;
		if(in == true){
			tempIn = tempEdge->orientation;
		}else{
			tempIn = !tempEdge->orientation;
		}
		bool result = VisitOutNode(scaffoldGraphHead, tempEdge->nodeIndex, tempIn, leftIndex, visited);
		if(result == true){
			return true;
		}
		tempEdge = tempEdge->next;
	}
	return false;   
}


void MergeBubbleInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead){
	ScaffoldGraphNode * tempEdge = NULL;
	ScaffoldGraphNode * tempEdge1 = NULL;
	ScaffoldGraphNode * tempEdgeNext = NULL;
	ScaffoldGraphNode * tempEdgeNext1 = NULL;
	
	bool * visited = (bool *)malloc(sizeof(bool)*scaffoldGraphHead->scaffoldGraphNodeCount);
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
		while(tempEdge != NULL){
			tempEdge1 = tempEdge->next;
			tempEdgeNext = tempEdge->next;
			while(tempEdge1 != NULL){
				for(int j = 0; j < scaffoldGraphHead->scaffoldGraphNodeCount; j++){
					visited[j] = false;
				}
				bool visit = VisitOutNode(scaffoldGraphHead, tempEdge->nodeIndex, tempEdge->orientation, tempEdge1->nodeIndex, visited);
				
				tempEdgeNext1 = tempEdge1->next;
				if(visit == true){
					int nodeIndex = tempEdge->nodeIndex;
					bool orientation = tempEdge->orientation;
					scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge->nodeIndex, tempEdge->orientation);
					if(orientation == true){
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}else{
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}
					tempEdge = tempEdgeNext;
					break;
				}else{
					for(int j = 0; j < scaffoldGraphHead->scaffoldGraphNodeCount; j++){
						visited[j] = false;
					}
					visit = VisitOutNode(scaffoldGraphHead, tempEdge1->nodeIndex, tempEdge1->orientation, tempEdge->nodeIndex, visited);
					if(visit == true){
						int nodeIndex = tempEdge1->nodeIndex;
						bool orientation = tempEdge1->orientation;
						scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge1->nodeIndex, tempEdge1->orientation);
						if(orientation == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
						}
						tempEdge1 = tempEdgeNext1;
						continue;
					}
				}
				tempEdge1 = tempEdge1->next;
			}
			if(tempEdge1 == NULL){
				tempEdge = tempEdge->next;
			}
			
		}
		
		
		
		tempEdge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
		while(tempEdge != NULL){
			tempEdge1 = tempEdge->next;
			tempEdgeNext = tempEdge->next;
			while(tempEdge1 != NULL){
				for(int j = 0; j < scaffoldGraphHead->scaffoldGraphNodeCount; j++){
					visited[j] = false;
				}
				bool visit = VisitInNode(scaffoldGraphHead, tempEdge->nodeIndex, tempEdge->orientation, tempEdge1->nodeIndex, visited);
				tempEdgeNext1 = tempEdge1->next;
				if(visit == true){
					int nodeIndex = tempEdge->nodeIndex;
					bool orientation = tempEdge->orientation;
					scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge->nodeIndex, tempEdge->orientation);
					if(orientation == true){
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}else{
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}
					tempEdge = tempEdgeNext;
					break;
				}else{
					for(int j = 0; j < scaffoldGraphHead->scaffoldGraphNodeCount; j++){
						visited[j] = false;
					}
					visit = VisitInNode(scaffoldGraphHead, tempEdge1->nodeIndex, tempEdge1->orientation, tempEdge->nodeIndex, visited);
					if(visit == true){
						int nodeIndex = tempEdge1->nodeIndex;
						bool orientation = tempEdge1->orientation;
						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge1->nodeIndex, tempEdge1->orientation);
						if(orientation == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
						}
						tempEdge1 = tempEdgeNext1;
						continue;
					}
				}
				tempEdge1 = tempEdge1->next;
			}
			if(tempEdge1 == NULL){
				tempEdge = tempEdge->next;
			}
			
		}
		
		
		
		
	}

}

void RemoveCycleInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead){
	ScaffoldGraphNode * tempEdge = NULL;
	ScaffoldGraphNode * tempEdge1 = NULL;
	ScaffoldGraphNode * tempEdgeNext = NULL;
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		
		tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
		while(tempEdge != NULL){
			tempEdgeNext = tempEdge->next;
			tempEdge1 = tempEdge->next;
			bool del = false;
			while(tempEdge1 != NULL){
				del = false;
				if(tempEdge->nodeIndex == tempEdge1->nodeIndex){
					int nodeIndex = tempEdge1->nodeIndex;
					bool orientation = tempEdge->orientation;
					bool orientation1 = tempEdge1->orientation;
					if(tempEdge->aligningReadCount > tempEdge1->aligningReadCount){
						scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge1->nodeIndex, tempEdge1->orientation);
						if(orientation1 == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation1);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation1);
						}
					}else if(tempEdge->aligningReadCount < tempEdge1->aligningReadCount){
						scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge->nodeIndex, tempEdge->orientation);
						if(orientation == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
						}
					}else{
						scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge1->nodeIndex, tempEdge1->orientation);
						if(orientation1 == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation1);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation1);
						}
						scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge->nodeIndex, tempEdge->orientation);
						if(orientation == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
						}
					}
					del = true;
					break;
				}
				tempEdge1 = tempEdge1->next;
			}
			if(del != false){
				tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			}else{
				tempEdge = tempEdge->next;
			}
		}
		
		
		
		
		tempEdge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
		while(tempEdge != NULL){
			tempEdge1 = tempEdge->next;
			bool del = false;
			while(tempEdge1 != NULL){
				del = false;
				if(tempEdge->nodeIndex == tempEdge1->nodeIndex){
					int nodeIndex = tempEdge1->nodeIndex;
					bool orientation = tempEdge->orientation;
					bool orientation1 = tempEdge1->orientation;
					if(tempEdge->aligningReadCount > tempEdge1->aligningReadCount){
						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge1->nodeIndex, tempEdge1->orientation);
						if(orientation1 == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation1);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation1);
						}
					}else if(tempEdge->aligningReadCount < tempEdge1->aligningReadCount){
						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge->nodeIndex, tempEdge->orientation);
						if(orientation == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
						}
					}else{
						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge1->nodeIndex, tempEdge1->orientation);
						
						if(orientation1 == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation1);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation1);
						}
						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge->nodeIndex, tempEdge->orientation);
						if(orientation == true){
							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
						}
					}
					del = true;
					break;
				}
				tempEdge1 = tempEdge1->next;
			}
			if(del != false){
				tempEdge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			}else{
				tempEdge = tempEdge->next;
			}
		}
		
		
		
	}
	
}

AligningResult * FindIndexOfAligningResultHead(AligningResultHead * aligningResultHead, int index){
	
	for(int i = 0; i < aligningResultHead->aligningResultCount; i++){
		if(aligningResultHead->aligningResult[i].contigIndex == index){
			return & aligningResultHead->aligningResult[i];
		}
	}
	return NULL;
}

void InsertEdgeToScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, AligningResultHead * aligningResultHead, ContigGraphHead * contigGraphHead, ContigSetHead * contigSetHead){
	
	int readIndex = aligningResultHead->readIndex;
	int readLength= aligningResultHead->readLength;
	
	AligningResult * tempLeft = NULL;
	AligningResult * tempRight = NULL;
	
	for(int i = 0; i < contigGraphHead->contigGraphNodeCount - 1; i++){
		if(contigGraphHead->longestPath[i + 1] == -1){
			break;
		}
		tempLeft = FindIndexOfAligningResultHead(aligningResultHead, contigGraphHead->longestPath[i]);
		tempRight = FindIndexOfAligningResultHead(aligningResultHead, contigGraphHead->longestPath[i + 1]);
		InsertSingleEdgeToScaffoldGraph(scaffoldGraphHead, readIndex, tempLeft, tempLeft->orientation, tempRight, tempRight->orientation, readLength, contigSetHead);
	}	
	
}


void IncreaseAligningResultHead(AligningResultHead * aligningResultHead, int count){
	
	
	
}

void GetScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, ContigGraphHead * contigGraphHead, char * file, char * line, int maxSize){
	
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	const char * split = "--"; 
	char * p; 
	aligningResultHead->aligningResultCount = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		if(line[0] == 'r'){
			if(aligningResultHead->aligningResultCount > 1){
				OptimizeALigningResultSet(aligningResultHead, contigSetHead, contigGraphHead);
				InsertEdgeToScaffoldGraph(scaffoldGraphHead, aligningResultHead, contigGraphHead, contigSetHead);
			}
			p = strtok(line,split);
			p = strtok(NULL,split);
			aligningResultHead->readIndex = atoi(p);
			p = strtok(NULL,split);
			aligningResultHead->readLength = atoi(p);

			aligningResultHead->aligningResultCount = 0;
			continue;
		}
		
		p = strtok(line,split);
		aligningResultHead->aligningResult[aligningResultHead->aligningResultCount].orientation = atoi(p);
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[aligningResultHead->aligningResultCount].readStartPosition = atoi(p);
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[aligningResultHead->aligningResultCount].readEndPosition = atoi(p);
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[aligningResultHead->aligningResultCount].contigStartPosition = atoi(p);
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[aligningResultHead->aligningResultCount].contigEndPosition = atoi(p);
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[aligningResultHead->aligningResultCount].contigIndex = atoi(p);
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[aligningResultHead->aligningResultCount].overlapLength = atoi(p);
		
		
		
		aligningResultHead->aligningResultCount++;
		if(aligningResultHead->aligningResultCount >= aligningResultHead->allocateAligningResultCount - 1){
			IncreaseAligningResultHead(aligningResultHead, 100);
		}
		
    } 

	if(aligningResultHead->aligningResultCount > 1){
		OptimizeALigningResultSet(aligningResultHead, contigSetHead, contigGraphHead);
		InsertEdgeToScaffoldGraph(scaffoldGraphHead, aligningResultHead, contigGraphHead, contigSetHead);
	}
	fflush(fp);
	fclose(fp);
	
}


void DeleteEdgeWithMinReadCount(ScaffoldGraphHead * scaffoldGraphHead, int minReadCount){
	ScaffoldGraphNode * tempEdge = NULL;
	ScaffoldGraphNode * tempEdge1 = NULL;
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		
		if(scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL){
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			while(tempEdge != NULL){
				int nodeIndex = tempEdge->nodeIndex;
				bool orientation = tempEdge->orientation;
				
				if(tempEdge->aligningReadCount < minReadCount){
					tempEdge1 = tempEdge->next;
					scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge->nodeIndex, orientation);
					if(orientation == true){
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}else{
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}
					tempEdge = tempEdge1;
					continue;
				}
				
				tempEdge = tempEdge->next;
			}
		}
		
		if(scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL){
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			while(tempEdge != NULL){
				int nodeIndex = tempEdge->nodeIndex;
				bool orientation = tempEdge->orientation;
				
				if(tempEdge->aligningReadCount < minReadCount){
					tempEdge1 = tempEdge->next;
					scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge->nodeIndex, orientation);
					if(orientation == true){
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}else{
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}
					tempEdge = tempEdge1;
					continue;
				}
				tempEdge = tempEdge->next;
			}
		}
		
	}


}

int CountAverageReadCount(ScaffoldGraphHead * scaffoldGraphHead){
	ScaffoldGraphNode * tempEdge = NULL;
	ScaffoldGraphNode * tempEdge1 = NULL;
	int averageCount = 0;
	int edgeCount = 0;
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		
		if(scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL){
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			while(tempEdge != NULL){
				averageCount = tempEdge->aligningReadCount + averageCount;
				edgeCount++;
				tempEdge = tempEdge->next;
			}
		}
		
		
	}
	cout<<"averageCount:"<<averageCount<<"--edgeCount:"<<edgeCount<<"--ration:"<<(double)averageCount/edgeCount<<endl;
	return averageCount/(3*edgeCount);


}

void OptimizeScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * cotigSetHead){
	
	ScaffoldGraphNode * tempEdge = NULL;
	ScaffoldGraphNode * tempEdge1 = NULL;
	int maxCount = -1;
	int maxIndex = -1;
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		
		if(scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL){
			scaffoldGraphHead->scaffoldGraph[i].outEdge = OptimizeEdgeInScaffoldGraph(scaffoldGraphHead->scaffoldGraph[i].outEdge, cotigSetHead, i);
		}
		
		if(scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL){
			scaffoldGraphHead->scaffoldGraph[i].inEdge = OptimizeEdgeInScaffoldGraph(scaffoldGraphHead->scaffoldGraph[i].inEdge, cotigSetHead, i);
		}
		
	}
	
	RemoveCycleInScaffoldGraph(scaffoldGraphHead);
	int readAverageCount = CountAverageReadCount(scaffoldGraphHead);
	DeleteEdgeWithMinReadCount(scaffoldGraphHead, readAverageCount);
	
}























#endif