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


void RemoveCycleInScaffoldGraphTwoNode(ScaffoldGraphHead * scaffoldGraphHead){
	ScaffoldGraphNode * tempEdge = NULL;
	ScaffoldGraphNode * tempEdge1 = NULL;
	ScaffoldGraphNode * tempEdgeNext = NULL;
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		
		tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
		
		while(tempEdge != NULL){
			tempEdge1 = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			bool del = false;
			while(tempEdge1 != NULL){
				if(tempEdge1->nodeIndex == tempEdge->nodeIndex){
					int index = tempEdge1->nodeIndex;
					bool orientation = tempEdge->orientation;
					bool orientation1 = tempEdge1->orientation;
					scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, index, orientation1);
					scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, index, orientation);
					del = true;
					break;
				}
				tempEdge1 = tempEdge1->next;
			}
			if(del == true){
				tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			}else{
				tempEdge = tempEdge->next;
			}
			
			
		}
	}
	
}


void AddEdgeInSimpleGraph(SimpleGraph * simpleGraph, long int simpleGraphNodeCount, int leftIndex, bool leftOrientation, int rightIndex, bool rightOrientation){
	for(long int i = 0; i < simpleGraphNodeCount; i++){
		if(simpleGraph[i].leftIndex == -1){
			simpleGraph[i].leftIndex = leftIndex;
			simpleGraph[i].rightIndex = rightIndex;
			if(leftOrientation == 1 && rightOrientation == 1){
				simpleGraph[i].count1 ++;
			}
			if(leftOrientation == 1 && rightOrientation == 0){
				simpleGraph[i].count2 ++;
			}
			if(leftOrientation == 0 && rightOrientation == 0){
				simpleGraph[i].count3 ++;
			}
			if(leftOrientation == 0 && rightOrientation == 1){
				simpleGraph[i].count4 ++;
			}
			break;
		}
		if(simpleGraph[i].leftIndex == leftIndex && simpleGraph[i].rightIndex == rightIndex){
			if(leftOrientation == 1 && rightOrientation == 1){
				simpleGraph[i].count1 ++;
			}
			if(leftOrientation == 1 && rightOrientation == 0){
				simpleGraph[i].count2 ++;
			}
			if(leftOrientation == 0 && rightOrientation == 0){
				simpleGraph[i].count3 ++;
			}
			if(leftOrientation == 0 && rightOrientation == 1){
				simpleGraph[i].count4 ++;
			}
			break;
		}
		if(simpleGraph[i].leftIndex == rightIndex && simpleGraph[i].rightIndex == leftIndex){
			if(leftOrientation == 1 && rightOrientation == 1){
				simpleGraph[i].count3 ++;
			}
			if(leftOrientation == 0 && rightOrientation == 1){
				simpleGraph[i].count4 ++;
			}
			if(leftOrientation == 0 && rightOrientation == 0){
				simpleGraph[i].count1 ++;
			}
			if(leftOrientation == 1 && rightOrientation == 0){
				simpleGraph[i].count2 ++;
			}
			break;
		}
	}

}

SimpleGraph * GetLineIndex(ContigSetHead * contigSetHead, bool * lineIndex, char * file, char * line, int maxSize, long int & lineCount, long int & simpleGraphNodeCount){//the long reads support the two neigbor is large than a threshold, else line[i] = true; 
	
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	char * p;
	const char * split = ","; 
	simpleGraphNodeCount = 0;
	lineCount = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		lineCount++;
		p = strtok(line,split);
		
		int count = atoi(p);
		if(count <= 1){
			continue;
		}else{
			simpleGraphNodeCount = simpleGraphNodeCount + count - 1;
		}
	}
	
	
	SimpleGraph * simpleGraph = (SimpleGraph *)malloc(sizeof(SimpleGraph)*simpleGraphNodeCount);
	for(long int i = 0; i < simpleGraphNodeCount; i++){
		simpleGraph[i].leftIndex = -1;
		simpleGraph[i].rightIndex = -1;
		simpleGraph[i].count1 = 0;
		simpleGraph[i].count2 = 0;
		simpleGraph[i].count3 = 0;
		simpleGraph[i].count4 = 0;
	}
	
	fclose(fp);
	
	
	if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	while((fgets(line, maxSize, fp)) != NULL){ 
		p = strtok(line,split);
		
		int count = atoi(p);
		if(count <= 1){
			continue;
		}
		p = strtok(NULL,split);
		int startContigIndex = atoi(p);
		p = strtok(NULL,split);
		int distance = atoi(p);
		p = strtok(NULL,split);
		bool startOrientation = atoi(p);
		p = strtok(NULL,split);
	 	int startOverlapLength = atoi(p);
		
		int a = 2;
		while(a <= count){
			p = strtok(NULL,split);
			int endContigIndex = atoi(p);
			p = strtok(NULL,split);
			int distance1 = atoi(p);
			p = strtok(NULL,split);
			bool endOrientation = atoi(p);
			p = strtok(NULL,split);
	 		int endOverlapLength = atoi(p);
			
			AddEdgeInSimpleGraph(simpleGraph, simpleGraphNodeCount, startContigIndex, startOrientation, endContigIndex, endOrientation);
			
			startContigIndex = endContigIndex;
			startOrientation = endOrientation;
			startOverlapLength = endOverlapLength;
			a++;
		}
	}
	
	fclose(fp);
	
	long int edgeCount = 0;
	long int alignReadCount = 0;
	for(long int i = 0; i < simpleGraphNodeCount; i++){
		if(simpleGraph[i].leftIndex == -1){
			break;
		}
		bool a = false;
		if(simpleGraph[i].count1 >= 1){
			edgeCount =  edgeCount + simpleGraph[i].count1;
			a = true;
		}
		if(simpleGraph[i].count2 >= 1){
			edgeCount =  edgeCount + simpleGraph[i].count2;
			a = true;
		}
		if(simpleGraph[i].count3 >= 1){
			edgeCount =  edgeCount + simpleGraph[i].count3;
			a = true;
		}
		if(simpleGraph[i].count4 >= 1){
			edgeCount =  edgeCount + simpleGraph[i].count4;
			a = true;
		}
		if(a == true){
			alignReadCount++;
		}
		
	}

	int minEdgeCout = edgeCount/(10*alignReadCount);
	if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	lineCount = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		lineCount++;
		p = strtok(line,split);
		
		int count = atoi(p);
		if(count <= 1){
			lineIndex[lineCount - 1] = true;
			continue;
		}
		p = strtok(NULL,split);
		int startContigIndex = atoi(p);
		p = strtok(NULL,split);
		int distance = atoi(p);
		p = strtok(NULL,split);
		bool startOrientation = atoi(p);
		p = strtok(NULL,split);
	 	int startOverlapLength = atoi(p);
		
		int a = 2;
		while(a <= count){
			p = strtok(NULL,split);
			int endContigIndex = atoi(p);
			p = strtok(NULL,split);
			int distance1 = atoi(p);
			p = strtok(NULL,split);
			bool endOrientation = atoi(p);
			p = strtok(NULL,split);
	 		int endOverlapLength = atoi(p);
			
			for(long int i = 0; i < simpleGraphNodeCount; i++){
				if(simpleGraph[i].leftIndex == -1){
					break;
				}
				if(simpleGraph[i].leftIndex == startContigIndex && simpleGraph[i].rightIndex == endContigIndex){
					if(startOrientation == 1 && endOrientation == 1){
						if(simpleGraph[i].count1 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					if(startOrientation == 1 && endOrientation == 0){
						if(simpleGraph[i].count2 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					if(startOrientation == 0 && endOrientation == 0){
						if(simpleGraph[i].count3 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					if(startOrientation == 0 && endOrientation == 1){
						if(simpleGraph[i].count4 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					break;
				}
				if(simpleGraph[i].leftIndex == endContigIndex && simpleGraph[i].rightIndex == startContigIndex){
					if(startOrientation == 1 && endOrientation == 1){
						if(simpleGraph[i].count3 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					if(startOrientation == 0 && endOrientation == 1){
						if(simpleGraph[i].count4 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					if(startOrientation == 0 && endOrientation == 0){
						if(simpleGraph[i].count1 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					if(startOrientation == 1 && endOrientation == 0){
						if(simpleGraph[i].count2 <= minEdgeCout){
							lineIndex[lineCount - 1] = true;
						}
					}
					break;
				}
			}
			
			startContigIndex = endContigIndex;
			startOrientation = endOrientation;
			startOverlapLength = endOverlapLength;
			a++;
		}
	}
	
	fclose(fp);
	
	return simpleGraph;
}

void GetGraphNeighbor(SimpleGraph * simpleGraph, ContigSetHead * contigSetHead, long int simpleGraphNodeCount){
	for(int j = 0; j < contigSetHead->contigCount; j++){
		int leftNeighborIndex = -1;
		int rightNeighborIndex = -1;
		int leftMax = 0;
		int rightMax = 0;
		int leftSecondMax = 0;
		int rightSecondMax = 0;
		int leftCount = 0;;
		int rightCount = 0;
		for(int i = 0; i < simpleGraphNodeCount; i++){
			if(simpleGraph[i].leftIndex == -1){
					break;
			}
			if(simpleGraph[i].leftIndex == j){
				if(simpleGraph[i].count1 > 0){
					if(rightMax < simpleGraph[i].count1){
						rightMax = simpleGraph[i].count1;
						rightNeighborIndex = simpleGraph[i].rightIndex;
					}
					rightCount++;
				}
				if(simpleGraph[i].count2 > 0){
					if(rightMax < simpleGraph[i].count2){
						rightMax = simpleGraph[i].count2;
						rightNeighborIndex = simpleGraph[i].rightIndex;
					}
					rightCount++;
				}
				if(simpleGraph[i].count3 > 0){
					if(leftMax < simpleGraph[i].count3){
						leftMax = simpleGraph[i].count3;
						leftNeighborIndex = simpleGraph[i].rightIndex;
					}
					leftCount++;
				}
				if(simpleGraph[i].count4 > 0){
					if(leftMax < simpleGraph[i].count4){
						leftMax = simpleGraph[i].count4;
						leftNeighborIndex = simpleGraph[i].rightIndex;
					}
					leftCount++;
				}
			}
			if(simpleGraph[i].rightIndex == j){
				if(simpleGraph[i].count3 > 0){
					if(rightMax < simpleGraph[i].count3){
						rightMax = simpleGraph[i].count3;
						rightNeighborIndex = simpleGraph[i].leftIndex;
					}
					rightCount++;
				}
				if(simpleGraph[i].count4 > 0){
					if(leftMax < simpleGraph[i].count4){
						leftMax = simpleGraph[i].count4;
						leftNeighborIndex = simpleGraph[i].leftIndex;
					}
					leftCount++;
				}
				if(simpleGraph[i].count1 > 0){
					if(leftMax < simpleGraph[i].count1){
						leftMax = simpleGraph[i].count1;
						leftNeighborIndex = simpleGraph[i].leftIndex;
					}
					leftCount++;
				}
				if(simpleGraph[i].count2 > 0){
					if(rightMax < simpleGraph[i].count2){
						rightMax = simpleGraph[i].count2;
						rightNeighborIndex = simpleGraph[i].leftIndex;
					}
					rightCount++;
				}
			}
			
		}
		
		if(leftCount <=1 && rightCount <=1){
			continue;	
		}
		
		
		for(int i = 0; i < simpleGraphNodeCount; i++){
			if(simpleGraph[i].leftIndex == -1){
					break;
			}
			if(simpleGraph[i].leftIndex == j){
				if(simpleGraph[i].count1 > 0){
					if(rightSecondMax < simpleGraph[i].count1 && rightNeighborIndex != simpleGraph[i].rightIndex){
						rightSecondMax = simpleGraph[i].count1;
					}
				}
				if(simpleGraph[i].count2 > 0){
					if(rightSecondMax < simpleGraph[i].count2 && rightNeighborIndex != simpleGraph[i].rightIndex){
						rightSecondMax = simpleGraph[i].count2;
					}
				}
				if(simpleGraph[i].count3 > 0){
					if(leftSecondMax < simpleGraph[i].count3 && leftNeighborIndex != simpleGraph[i].rightIndex){
						leftSecondMax = simpleGraph[i].count3;
					}
				}
				if(simpleGraph[i].count4 > 0){
					if(leftSecondMax < simpleGraph[i].count4 && leftNeighborIndex != simpleGraph[i].rightIndex){
						leftSecondMax = simpleGraph[i].count4;
					}
				}
			}
			if(simpleGraph[i].rightIndex == j){
				if(simpleGraph[i].count3 > 0){
					if(rightSecondMax < simpleGraph[i].count3 && rightNeighborIndex != simpleGraph[i].leftIndex){
						rightSecondMax = simpleGraph[i].count3;
					}
				}
				if(simpleGraph[i].count4 > 0){
					if(leftSecondMax < simpleGraph[i].count4 && leftNeighborIndex != simpleGraph[i].leftIndex){
						leftSecondMax = simpleGraph[i].count4;
					}
				}
				if(simpleGraph[i].count1 > 0){
					if(leftSecondMax < simpleGraph[i].count1 && leftNeighborIndex != simpleGraph[i].leftIndex){
						leftSecondMax = simpleGraph[i].count1;
					}
				}
				if(simpleGraph[i].count2 > 0){
					if(rightSecondMax < simpleGraph[i].count2 && rightNeighborIndex != simpleGraph[i].leftIndex){
						rightSecondMax = simpleGraph[i].count2;
					}
				}
			}
			
		}

		if(leftSecondMax > 0){
			if(double(leftMax)/leftSecondMax < 6){
				contigSetHead->visited[j] = true;
			}
		}
		if(rightSecondMax > 0){
			if(double(rightMax)/rightSecondMax < 6){
				contigSetHead->visited[j] = true;
			}
		}
		
	}
	
}

void GetScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, char * file, char * line, int maxSize, FILE * fpUnique, FILE * fpAmbiguous){
	long int lineCount = 0;
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	char * p;
	const char * split = ","; 

	while((fgets(line, maxSize, fp)) != NULL){ 
		lineCount++;
		
	}
	fclose(fp);
	bool * lineIndex = (bool *)malloc(sizeof(bool)*lineCount);
	for(long int i = 0; i < lineCount; i++){
		lineIndex[i] = false;
		
	}
	for(int i = 0; i < contigSetHead->contigCount; i++){
		contigSetHead->visited[i] = false;
		contigSetHead->repeatContigIndex[i] = false;
	}
	long int simpleGraphNodeCount = 0;
	SimpleGraph * simpleGraph = GetLineIndex(contigSetHead, lineIndex, file, line, maxSize, lineCount, simpleGraphNodeCount);
	
	GetGraphNeighbor(simpleGraph, contigSetHead, simpleGraphNodeCount);
	 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }

	
	int leftNeighbor[contigSetHead->contigCount];
	int rightNeighbor[contigSetHead->contigCount];
	for(int i = 0; i < contigSetHead->contigCount; i++){
		//contigSetHead->visited[i] = false;
		//contigSetHead->repeatContigIndex[i] = false;
		leftNeighbor[i] = -1;
		rightNeighbor[i] = -1;
	}
	
	int contigIndex = -1;
	int previousContigIndex = -1;
	int orientation = -1;
	int previousOrientation = -1;
	lineCount = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		lineCount ++;
		if(lineIndex[lineCount - 1] == true){
			continue;
		}
		
		p = strtok(line,split);
		
		int count = atoi(p);
		if(count <= 1){
			continue;
		}
		int a = 1;
		contigIndex = -1;
		while(a <= count * 4){
			p = strtok(NULL,split);
			if(a%4 == 1){
				previousContigIndex = contigIndex;
				contigIndex = atoi(p);
				
				if(a/4 < count-1 && a/4 > 0){
					contigSetHead->repeatContigIndex[atoi(p)] = true;
				}
			}
			
			a++;
		}
	}
	fclose(fp);
	
	for(int i = 0; i < contigSetHead->contigCount; i++){
		if(contigSetHead->contigSet[i].shortContig == true){
			fprintf(fpAmbiguous, ">%d\n", i);
			fprintf(fpAmbiguous, "%s\n", contigSetHead->contigSet[i].contig);
			continue;
		}
		
		if(contigSetHead->visited[i] == true && contigSetHead->repeatContigIndex[i] == true){
			fprintf(fpAmbiguous, ">%d\n", i);
			fprintf(fpAmbiguous, "%s\n", contigSetHead->contigSet[i].contig);
		}else{
			fprintf(fpUnique, ">%d\n", i);
			fprintf(fpUnique, "%s\n", contigSetHead->contigSet[i].contig);
		}
	}
	fflush(fpAmbiguous);
	fflush(fpUnique);
	
	
	if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	int readIndex = 0;
	lineCount = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		
		
		lineCount ++;
		if(lineIndex[lineCount - 1] == true){
			continue;
		}
		p = strtok(line,split);
		
		int count = atoi(p);
		
		if(count <= 1){
			readIndex++;
			continue;
		}

		p = strtok(NULL,split);
		int startContigIndex = atoi(p);
		p = strtok(NULL,split);
		int distance = atoi(p);
		p = strtok(NULL,split);
		bool startOrientation = atoi(p);
		p = strtok(NULL,split);
	 	int startOverlapLength = atoi(p);

		if(contigSetHead->visited[startContigIndex] == true && contigSetHead->repeatContigIndex[startContigIndex] == true){
			startContigIndex = -1;
		}
		
		int a = 2;
		while(a <= count){
			p = strtok(NULL,split);
			int endContigIndex = atoi(p);
			p = strtok(NULL,split);
			int distance1 = atoi(p);
			p = strtok(NULL,split);
			bool endOrientation = atoi(p);
			p = strtok(NULL,split);
	 		int endOverlapLength = atoi(p);
	
			int realEndContigIndex = endContigIndex;
			if(contigSetHead->visited[endContigIndex] == true && contigSetHead->repeatContigIndex[endContigIndex] == true){
				endContigIndex = -1;
			}

			if(startContigIndex != -1){
				if(endContigIndex == -1){
					distance = distance + distance1 + contigSetHead->contigSet[realEndContigIndex].contigLength;
				}else{
					int min = 0;
					if(startOverlapLength > endOverlapLength){
						min = endOverlapLength;
					}else{
						min = startOverlapLength;
					}
					contigSetHead->contigSet[startContigIndex].uniqueContig = true;
					contigSetHead->contigSet[endContigIndex].uniqueContig = true;
					InsertOutOrInEdge(scaffoldGraphHead, readIndex, startContigIndex, startOrientation, endContigIndex, endOrientation, distance, min);
					distance = distance1;
					startContigIndex = endContigIndex;
					startOrientation = endOrientation;
					startOverlapLength = endOverlapLength;
					
				}
			}else{
				distance = distance1;
				startContigIndex = endContigIndex;
				startOrientation = endOrientation;
				startOverlapLength = endOverlapLength;
			}
			
			
			a++;
		}
		
		
		readIndex++;
		
	}
	fclose(fp);
}


void GetScaffoldGraphNonUniqueLocalScaffold(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, char * file, char * line, int maxSize){
	
	
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	const char * split = ","; 
	char * p; 
	
	if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	int readIndex = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		p = strtok(line,split);
		
		int count = atoi(p);
		
		if(count <= 1){
			readIndex++;
			continue;
		}

		p = strtok(NULL,split);
		int startContigIndex = atoi(p);
		p = strtok(NULL,split);
		int distance = atoi(p);
		p = strtok(NULL,split);
		bool startOrientation = atoi(p);
		p = strtok(NULL,split);
	 	int startOverlapLength = atoi(p);
		
		int a = 2;
		while(a <= count){
			p = strtok(NULL,split);
			int endContigIndex = atoi(p);
			p = strtok(NULL,split);
			int distance1 = atoi(p);
			p = strtok(NULL,split);
			bool endOrientation = atoi(p);
			p = strtok(NULL,split);
	 		int endOverlapLength = atoi(p);
			
			
			int min = 0;
			if(startOverlapLength > endOverlapLength){
				min = endOverlapLength;
			}else{
				min = startOverlapLength;
			}
			
			
			InsertOutOrInEdge(scaffoldGraphHead, readIndex, startContigIndex, startOrientation, endContigIndex, endOrientation, distance, min);
	
			int realEndContigIndex = endContigIndex;
			
			distance = distance1;
			startContigIndex = endContigIndex;
			startOrientation = endOrientation;
			startOverlapLength = endOverlapLength;
			
			a++;
		}
		
		
		readIndex++;
		
	}
	fclose(fp);
}


void OutputUniqueContigSet(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, char * file){
	

}


void SortNodeOptimize(ContigSetHead * contigSetHead, char * file, char * line, int maxSize, FILE * longReadFileFP){
	
	
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	const char * split = ","; 
	char * p;
	
	
	int maxCount = 100;
	int aligningSingle[7*maxCount];
	
	int lineNum = -1;
	while((fgets(line, maxSize, fp)) != NULL){ 
		
		lineNum++;
		p = strtok(line,split);
		int count = atoi(p);
		if(count <= 1){
			continue;
		}
		p = strtok(NULL,split);
		int readLength = atoi(p);
		
		int a = 1;
		
		while(a <= count*7){
			p = strtok(NULL,split);
			aligningSingle[a - 1] = atoi(p);
			a++;
		}
		
		for(int i = 0; i < count - 1; i++){
			for(int j = i + 1; j < count; j++){
				if(aligningSingle[i*7 + 1] == aligningSingle[j*7 + 1]
			  		&& aligningSingle[i*7 + 1] == 0){
					if(aligningSingle[i*7 + 5]  > aligningSingle[j*7 + 5]){
						aligningSingle[j*7] = -1;
					}
					if(aligningSingle[i*7 + 5]  < aligningSingle[j*7 + 5]){
						aligningSingle[i*7] = -1;
					}
					if(aligningSingle[i*7 + 5] == aligningSingle[j*7 + 5]){
						aligningSingle[j*7] = -1;
						aligningSingle[i*7] = -1;
					}
				}
				if(aligningSingle[i*7 + 2] == aligningSingle[j*7 + 2]
			 	 	&& aligningSingle[i*7 + 2] == readLength - 1){
					if(aligningSingle[i*7 + 5]  > aligningSingle[j*7 + 5]){
						aligningSingle[j*7] = -1;
					}
					if(aligningSingle[i*7 + 5]  < aligningSingle[j*7 + 5]){
						aligningSingle[i*7] = -1;
					}
					if(aligningSingle[i*7 + 5] == aligningSingle[j*7 + 5]){
						aligningSingle[j*7] = -1;
						aligningSingle[i*7] = -1;
					}
				}
			}
		}
		
		
		int realCount = count;
		int interval = 0;
		for(int i = 0; i < count; i++){
			if(aligningSingle[i*7] == -1){
				interval++;
			}else{
				for(int j = 0; j < 7; j++){
					aligningSingle[(i - interval)*7 + j] = aligningSingle[i*7 + j];
				}
			}
		
		}
		count = count - interval;
		realCount = count;
		
		int previousIndex = -1;
		fprintf(longReadFileFP, "%d,", realCount); 
		int gapDistance = 0;
		for(int i = 0; i < count; i++){
			
			if(aligningSingle[i*7] == -1){
				break;
			}
			if(i == count - 1){
				gapDistance = 0;
			}else{
				gapDistance =  aligningSingle[(i+1)*7 + 1] - aligningSingle[i*7 + 2];
			}

			int j = i + 1;
			bool token = false;
			for(j = i + 1; j < count - 1; j++){
				if(aligningSingle[j*7] == -1){
					gapDistance = gapDistance + contigSetHead->contigSet[aligningSingle[j*7]].contigLength + aligningSingle[(j+1)*7 + 1] - aligningSingle[j*7 + 2]; 
				}else{
					token = true;
					break;
				}
			}
			
			if(token == false){
				if(aligningSingle[(count - 1)*7 + 3] == -1){
					gapDistance = 0;
					fprintf(longReadFileFP, "%d,%d,%d,%d,", aligningSingle[i*7], gapDistance, aligningSingle[i*7 + 6], aligningSingle[i*7 + 5]);
					break;
				}
			}
			
			
			fprintf(longReadFileFP, "%d,%d,%d,%d,", aligningSingle[i*7], gapDistance, aligningSingle[i*7 + 6], aligningSingle[i*7 + 5]);
			gapDistance = 0;
		}
		
		fprintf(longReadFileFP, "\n"); 
		
	}
	
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
				
				if(tempEdge->aligningReadCount <= minReadCount){
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
				
				if(tempEdge->aligningReadCount <= minReadCount){
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
		if(scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL){
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			while(tempEdge != NULL){
				averageCount = tempEdge->aligningReadCount + averageCount;
				edgeCount++;
				tempEdge = tempEdge->next;
			}
		}
		
		
	}
	return averageCount/edgeCount;


}

void OptimizeScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * cotigSetHead){
	
	
	
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
	int minReadCount = readAverageCount/10;
	DeleteEdgeWithMinReadCount(scaffoldGraphHead, minReadCount);

}

int GetEdgeNumber(ScaffoldGraphHead * scaffoldGraphHead, int index, bool out){
	ScaffoldGraphNode * tempEdge = NULL;
	if(out == true){
		tempEdge = scaffoldGraphHead->scaffoldGraph[index].outEdge;
	}else{
		tempEdge = scaffoldGraphHead->scaffoldGraph[index].inEdge;
	}
	int edgeNum = 0;
	while(tempEdge != NULL){
		edgeNum++;
		tempEdge = tempEdge->next;
	}
	return edgeNum;
}


bool GetOverlapEdgeIndex(ScaffoldGraphNode * left, ScaffoldGraphNode * right){
	int * leftArray = left->readIndexArray;
	int * rightArray = right->readIndexArray;
	
	for(int i = 0; i < left->aligningReadCount; i++){
		for(int j = 0; j < right->aligningReadCount; j++){
			if(leftArray[i] == rightArray[j]){
				return true;
			}
		}
	}
	return false;
}

int DeleteSpecialScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, long int index1){
    
	scaffoldGraph[index].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraph[index].outEdge, index1, 0);
    scaffoldGraph[index].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraph[index].outEdge, index1, 1);
	scaffoldGraph[index].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraph[index].inEdge, index1, 0);
	scaffoldGraph[index].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraph[index].inEdge, index1, 1);
    
}

int ScaffoldGraphNullEdge(ScaffoldGraphHead * scaffoldGraphHead){
    
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		scaffoldGraphHead->scaffoldGraph[i].outEdge = NULL;
		scaffoldGraphHead->scaffoldGraph[i].inEdge = NULL;
	}
    
}




#endif