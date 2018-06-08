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
#include "lp/lp_lib.h"



long int DetermineOrientationOfContigs(ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation){
    
    long int i = 0;
    long int j = 0;
    long int p = 1;
    long int c = 1000;
    
    bool * contigIndex = new bool[contigCount];
    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
    }
    long int edgeNumber = 0;
    long int constraintNumber = 0;
    ScaffoldGraphNode * tempEdge = NULL;
    
    for(i=0;i<contigCount;i++){
        tempEdge = scaffoldGraph[i].outEdge;
        while(tempEdge!=NULL){
			edgeNumber++;
            tempEdge = tempEdge->next;
        }
        tempEdge = scaffoldGraph[i].inEdge;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
        contigIndex[i] = false;
    }

    edgeNumber = edgeNumber/2;
    constraintNumber = 0;
	
    lprec *lp;
    int Ncol, *colno = NULL, ret = 0;
    REAL *row = NULL;
    
    Ncol = edgeNumber + contigCount;
    lp = make_lp(0, Ncol);
    if(lp == NULL){
        printf("couldn't construct a new linear programming model");
        exit(0);
    }
    double * weight = new double[edgeNumber]; 
    long int * edgeLeftNode = new long int[edgeNumber];
    long int * edgeRightNode = new long int[edgeNumber];
    
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL)){
        printf("couldn't new colno and row");
        exit(0);
    }
    
    set_add_rowmode(lp, TRUE);
    
    int ttt = 0;
    for(i=0;i<contigCount;i++){
        for(int q=0;q<2;q++){
            if(q==0){
                tempEdge = scaffoldGraph[i].outEdge;
            }else{
                tempEdge = scaffoldGraph[i].inEdge;
            }
            while(tempEdge!=NULL){
                if(index[i][tempEdge->nodeIndex] == false && index[tempEdge->nodeIndex][i] == false){
                    if(tempEdge->orientation == false){
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = 1;
                        colno[j] = tempEdge->nodeIndex+1; 
                        row[j++] = 1;
                        
                        colno[j] = contigCount + p; 
                        row[j++] = c;
                        if(!add_constraintex(lp, j, row, colno, LE, c+1)){
                            printf("couldn't add_constraintex");
                            exit(0);
                        }
                        

                        j = 0;
                        colno[j] = i+1;
                        row[j++] = -1;
                        colno[j] = tempEdge->nodeIndex+1; 
                        row[j++] = -1;
                        
                        colno[j] = contigCount + p; 
                        row[j++] = c;
                        if(!add_constraintex(lp, j, row, colno, LE, c-1)){
                            printf("couldn't add_constraintex");
                            exit(0);
                        }
                        
                        ttt++;
                        constraintNumber = constraintNumber+2;
                        
                    }else{
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = 1;
                        colno[j] = tempEdge->nodeIndex+1; 
                        row[j++] = -1;
                        
                        colno[j] = contigCount + p; 
                        row[j++] = c;
                        if(!add_constraintex(lp, j, row, colno, LE, c)){
                            printf("couldn't add_constraintex");
                        	exit(0);
                        }
					
                       j = 0;
                       colno[j] = tempEdge->nodeIndex+1;
                       row[j++] = 1;
                       colno[j] = i+1; 
                       row[j++] = -1; 
                        
                       colno[j] = contigCount + p; 
                       row[j++] = c;
                       if(!add_constraintex(lp, j, row, colno, LE, c)){
                            printf("couldn't add_constraintex");
                            exit(0);
                        }
                        
                        constraintNumber = constraintNumber+2;
                        ttt++;
                    }
                    
                    
                    j = 0;
                    colno[j] = contigCount + p;
                    row[j++] = 1;
                    add_constraintex(lp, j, row, colno, LE, 1);
                    j = 0;
                    colno[j] = contigCount + p;
                    row[j++] = 1;
                    add_constraintex(lp, j, row, colno, GE, 0);
                    constraintNumber = constraintNumber+2;
                      
                    edgeLeftNode[p-1] = i;
                    edgeRightNode[p-1] = tempEdge->nodeIndex;
                    p++;
                    
                    weight[p-2] = tempEdge->overlapLength;
					
                    index[i][tempEdge->nodeIndex] = true;
                    index[tempEdge->nodeIndex][i] = true;

                }
                tempEdge = tempEdge->next;
            }    
        }
    }

    for(i=0;i< contigCount;i++){
        set_binary(lp, i+1, TRUE);
    }
	cout<<"cc"<<endl;
    p--;
    
    
    j=0;
    for(i=0;i<p;i++){
        colno[j] = contigCount + i + 1; 
        row[j] = weight[j];
        j++;
    }
    if(!set_obj_fnex(lp, j, row, colno)){
        printf("couldn't set_obj_fnex");
        exit(0);
    }
	set_add_rowmode(lp, FALSE);
	set_timeout(lp, 1200);

    set_maxim(lp);
    set_scaling(lp, 128); 

	

    ret = solve(lp);

    if(!(ret==0 || ret ==1)){
		cout<<"ee:"<<ret<<endl;
    }

    REAL * pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
    get_primal_solution(lp,pv);

    double temp = 1;
    long int result = 0;
    for(i=contigCount+constraintNumber+1;i<constraintNumber+contigCount+p+1;i++){
        if(pv[i] == 0){
           
            DeleteSpecialScaffoldEdge(scaffoldGraph, edgeLeftNode[i-contigCount-constraintNumber-1], edgeRightNode[i-contigCount-constraintNumber-1]);
            DeleteSpecialScaffoldEdge(scaffoldGraph, edgeRightNode[i-contigCount-constraintNumber-1], edgeLeftNode[i-contigCount-constraintNumber-1]);
            result++;
            
        }
		
    }

	for(i=constraintNumber+1;i<contigCount+constraintNumber+1;i++){
        contigOrientation[i-constraintNumber-1] = pv[i];
    }
    
    delete [] contigIndex;
    for(i=0;i<contigCount;i++){
        delete [] index[i];
    }
    delete [] index;
    delete [] weight;
    delete [] edgeLeftNode;
    delete [] edgeRightNode;
    delete [] pv;
    
    delete_lp(lp);
    
    return result;
    
}


long int * IterativeDetermineOrderOfContigs(ContigSetHead * contigSetHead, ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation, long int * contigPosition, long int allContigLength){
    
    
    long int i = 0;
    long int j = 0;
    long int p = 1;
    long int c = 2*allContigLength;
    
    bool * contigVisited = new bool[contigCount];

    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
        contigVisited[i] = false;
    }

    long int edgeNumber = 0;
    ScaffoldGraphNode * tempEdge = NULL;
    
    for(i=0;i<contigCount;i++){
		tempEdge = scaffoldGraph[i].outEdge;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
        tempEdge = scaffoldGraph[i].inEdge;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
    }

    edgeNumber = edgeNumber/2;

    lprec *lp;
    int Ncol, *colno = NULL, ret = 0;
    REAL *row = NULL;
    
    Ncol = contigCount + edgeNumber;
    lp = make_lp(0, Ncol);
    if(lp == NULL){
        printf("couldn't construct a new linear programming model");
        exit(0);
    }
    double * weight = new double[edgeNumber]; 
    long int * edgeLeftNode = new long int[edgeNumber];
    long int * edgeRightNode = new long int[edgeNumber];
    long int * gapDistance = new long int[edgeNumber];
    
    
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL)){
        printf("couldn't new colno and row");
        exit(0);
    }
    
    set_add_rowmode(lp, TRUE);

    long int constraintNumber = 0;
    
    for(i=0;i<contigCount;i++){
        for(int q=0;q<2;q++){
            if(q==0){
                tempEdge = scaffoldGraph[i].outEdge;
            }else{
                tempEdge = scaffoldGraph[i].inEdge;
            }
            long int start = 0;
            while(tempEdge!=NULL){  
                
                if(index[i][tempEdge->nodeIndex] == false && index[tempEdge->nodeIndex][i] == false){
                    
                    if((contigOrientation[i] == 1 && q ==0) || (contigOrientation[i]==0 && q==1)){
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = -1;
                        colno[j] = tempEdge->nodeIndex+1; 
                        row[j++] = 1;

                           
                        colno[j] = contigCount + p; 
                        row[j++] = c;
                        if(!add_constraintex(lp, j, row, colno, LE, c+contigSetHead->contigSet[i].contigLength+tempEdge->gapDistance)){
                            printf("couldn't add_constraintex");
                            exit(0);
                        }  
                        
                        
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = 1;
                        colno[j] = tempEdge->nodeIndex+1; 
                        row[j++] = -1;
                        colno[j] = contigCount + p; 
                        row[j++] = c;
                        if(!add_constraintex(lp, j, row, colno, LE, c-contigSetHead->contigSet[i].contigLength-tempEdge->gapDistance)){
                            printf("couldn't add_constraintex");
                            exit(0);
                        }
						
                        edgeLeftNode[p-1] = i;
                        edgeRightNode[p-1] = tempEdge->nodeIndex;  
                        
                        contigVisited[i] = true;
                        contigVisited[tempEdge->nodeIndex] = true;
                    }
                    if((contigOrientation[i] == 0 && q ==0) || (contigOrientation[i]==1 && q==1)){
                        
						j = 0;
                        colno[j] = i+1;
                        row[j++] = 1;
                        colno[j] = tempEdge->nodeIndex+1; 
                        row[j++] = -1;
                        
                        
                        colno[j] = contigCount + p; 
                        row[j++] = c;
                        if(!add_constraintex(lp, j, row, colno, LE, c+contigSetHead->contigSet[tempEdge->nodeIndex].contigLength+tempEdge->gapDistance)){
                            printf("couldn't add_constraintex");
                            exit(0);
                        }
                        
                        
                        
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = -1;
                        colno[j] = tempEdge->nodeIndex+1; 
                        row[j++] = 1;
                        
                        colno[j] = contigCount + p; 
                        row[j++] = c;
                        if(!add_constraintex(lp, j, row, colno, LE, c-contigSetHead->contigSet[tempEdge->nodeIndex].contigLength-tempEdge->gapDistance)){
                            printf("couldn't add_constraintex");
                            exit(0);
                        }
                        edgeLeftNode[p-1] = tempEdge->nodeIndex;
                        edgeRightNode[p-1] = i;
                        
                        contigVisited[i] = true;
                        contigVisited[tempEdge->nodeIndex] = true;
                    }
                   
                    j = 0;
                    colno[j] = contigCount + p;
                    row[j++] = 1;
                    add_constraintex(lp, j, row, colno, LE, 1);
                    j = 0;
                    colno[j] = contigCount + p;
                    row[j++] = 1;
                    add_constraintex(lp, j, row, colno, GE, 0);
                    weight[p-1] = tempEdge->overlapLength;
                    gapDistance[p-1] = tempEdge->gapDistance;
                    p++;
                    constraintNumber = constraintNumber + 4; 
                    index[i][tempEdge->nodeIndex] = true;
                    index[tempEdge->nodeIndex][i] = true;
                    

                }
                tempEdge = tempEdge->next;
            }
                                     
        }
    }

    for(i=0;i<contigCount;i++){
        set_int(lp, i+1, TRUE);
    }

    p--;
    
    
    
    j=0;
    for(i=0;i<p;i++){
        colno[j] = contigCount + i + 1; 
        row[j] = weight[j];
        j++;
    }

    if(!set_obj_fnex(lp, j, row, colno)){
        printf("couldn't set_obj_fnex");
        exit(0);
    }
	set_add_rowmode(lp, FALSE);

	
    set_maxim(lp);
    set_timeout(lp, 1200);

    ret = solve(lp);

    if(!(ret==0 || ret ==1)){

        set_break_at_first(lp, true);
        ret = solve(lp);
        if(!(ret==0 || ret ==1)){
            return NULL;
        }
    }
    
    REAL * pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
    get_primal_solution(lp,pv);
	
	
	
	for(i=contigCount+constraintNumber+1;i<contigCount+constraintNumber+p+1;i++){
        if(pv[i] < 1){
            long int d = pv[1+constraintNumber+edgeRightNode[i-contigCount-constraintNumber-1]] 
				- pv[1+constraintNumber+edgeLeftNode[i-contigCount-constraintNumber-1]] 
				- contigSetHead->contigSet[edgeLeftNode[i-contigCount-constraintNumber-1]].contigLength;
            double varD = double(labs(d - gapDistance[i-contigCount-constraintNumber-1]))/500;
			
            if(varD>3){
                DeleteSpecialScaffoldEdge(scaffoldGraph, edgeLeftNode[i-contigCount-constraintNumber-1], edgeRightNode[i-contigCount-constraintNumber-1]);
            	DeleteSpecialScaffoldEdge(scaffoldGraph, edgeRightNode[i-contigCount-constraintNumber-1], edgeLeftNode[i-contigCount-constraintNumber-1]);
				continue;
            }
        }

    }
	
    
    for(i=0;i<contigCount;i++){

        contigPosition[i] = pv[1+constraintNumber+i];
    }
    
    long int trueNumber = 0;
    long int realTrueNumber = 0;
    
    for(i=0;i<contigCount;i++){
        delete [] index[i];
    }
    
    
    delete [] contigVisited;
    delete [] index;
    delete [] weight;
    delete [] edgeLeftNode;
    delete [] edgeRightNode;
    delete [] gapDistance;
    delete [] pv;
    
    delete_lp(lp);
    
}


ScaffoldSetHead * GetScaffoldSet(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, char * file, char * line, int maxSize, char * dir){
    
	bool * contigOrientation = (bool *)malloc(sizeof(bool)*scaffoldGraphHead->scaffoldGraphNodeCount);
	long int * contigPosition = (long int *)malloc(sizeof(long int)*scaffoldGraphHead->scaffoldGraphNodeCount);
	long int allContigLength = contigSetHead->allContigLength;

	DetermineOrientationOfContigs(scaffoldGraphHead->scaffoldGraph, scaffoldGraphHead->scaffoldGraphNodeCount, contigOrientation);

	IterativeDetermineOrderOfContigs(contigSetHead, scaffoldGraphHead->scaffoldGraph, scaffoldGraphHead->scaffoldGraphNodeCount, contigOrientation, contigPosition, allContigLength);

	ScaffoldGraphNode * tempEdge = NULL;
	ScaffoldGraphNode * tempEdge1 = NULL;
	
	int maxCount = -1;
	int maxIndex = -1;
	
	for(int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		
		if(scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL){
			maxCount = -1;
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			maxIndex = tempEdge->nodeIndex;
			bool maxOrientation = false;
			while(tempEdge != NULL){
				int nodeIndex = tempEdge->nodeIndex;
				bool orientation = tempEdge->orientation;
				
				if(tempEdge->aligningReadCount < 0 || tempEdge->overlapLength <= maxCount){
					tempEdge1 = tempEdge->next;
					scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempEdge->nodeIndex, orientation);
					if(orientation == true){
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}else{
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}
					tempEdge = tempEdge1;
					continue;
				}else{
					if(maxCount != -1){
						scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, maxIndex, maxOrientation);
						if(maxOrientation == true){
							scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge, i, maxOrientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge, i, maxOrientation);
						}
					}
					maxCount = tempEdge->overlapLength;
					maxIndex = tempEdge->nodeIndex;
					maxOrientation = tempEdge->orientation;
				}
				tempEdge = tempEdge->next;
			}
		}
		
		if(scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL){
			maxCount = -1;
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			maxIndex = tempEdge->nodeIndex;
			bool maxOrientation = false;
			while(tempEdge != NULL){
				int nodeIndex = tempEdge->nodeIndex;
				bool orientation = tempEdge->orientation;
				
				if(tempEdge->aligningReadCount < 0 || tempEdge->overlapLength <= maxCount){
					tempEdge1 = tempEdge->next;
					scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempEdge->nodeIndex, orientation);
					if(orientation == true){
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}else{
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}
					tempEdge = tempEdge1;
					continue;
				}else{
					if(maxCount != -1){
						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, maxIndex, maxOrientation);
						if(maxOrientation == true){
							scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge, i, maxOrientation);
						}else{
							scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge = DeleteScaffoldGraphSpecialNode(scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge, i, maxOrientation);
						}
					}
					maxCount = tempEdge->overlapLength;
					maxIndex = tempEdge->nodeIndex;
					maxOrientation = tempEdge->orientation;
				}
				tempEdge = tempEdge->next;
			}
		}
		
	}
	
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
	
	for(int i = 0; i < contigSetHead->contigCount; i++){
		contigSetHead->visited[i] = false;
	}
    
    for(int p = 0; p < scaffoldGraphHead->scaffoldGraphNodeCount; p++){
        
		
		
		i = sortNode[p];
		
        if(printIndex[i] == true || contigSetHead->contigSet[i].contigLength < 3000){
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
			
			temp = GetOptimizeNodeIndex(scaffoldGraphHead, j, orientation, scaffoldSetHead->scaffoldSet->contigSequence, 1, printIndex);
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
			
			temp = GetOptimizeNodeIndex(scaffoldGraphHead, j, orientation, scaffoldSetHead->scaffoldSet->contigSequence, 0, printIndex);
			
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
	char * tag0 = new char[400];
	strcpy(tag0, dir);
	strcat(tag0, "/tag1.fa");
	OutPutScaffoldTag(scaffoldSetHead->scaffoldSet, tag0);
	
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
	
	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);
	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);
	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);
	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);
	
	char * tag1 = new char[400];
	strcpy(tag1, dir);
	strcat(tag1, "/tag2.fa");
	OutPutScaffoldTag(scaffoldSetHead->scaffoldSet, tag1);
	
	InsertRepeatContigToSequence1(contigSetHead, scaffoldSetHead, printIndex, file, line, maxSize);
	
	char * tag = new char[400];
	strcpy(tag, dir);
	strcat(tag, "/tag3.fa");

	OutPutScaffoldTag(scaffoldSetHead->scaffoldSet, tag);
	
    return scaffoldSetHead;
   
}

bool SearchInsertContigSequence(ScaffoldSetHead * tempInsertSequenceHead, int startIndex, int endIndex, int * contigIndex, int * distance, bool * orientation, int * overlapLength, int count){
	ScaffoldSet * tempScaffoldSet = tempInsertSequenceHead->scaffoldSet;
	while(tempScaffoldSet != NULL){
		if(tempScaffoldSet->contigSequence == NULL || tempScaffoldSet->sequenceCount != labs(endIndex - startIndex) + 1){
			tempScaffoldSet = tempScaffoldSet->next;
			continue;
		}
		ContigSequence * tempContigSequence = tempScaffoldSet->contigSequence;
		if(startIndex < endIndex){
			bool token = false;
			for(int i = startIndex; i <= endIndex; i++ ){
				if(tempContigSequence->index != contigIndex[i] || tempContigSequence->orientation != orientation[i]){
					token = true;
					break;
				}
				tempContigSequence = tempContigSequence->next;
				if(tempContigSequence == NULL){
					token = true;
					break;
				}
			}
			if(token == true){
				return false;
			}
			tempScaffoldSet->sequenceCount++;
			return true;
		}
		if(endIndex < startIndex){
			bool token = false;
			for(int i = startIndex; i <= endIndex; i-- ){
				if(tempContigSequence->index != contigIndex[i] || tempContigSequence->orientation == orientation[i]){
					token = true;
					break;
				}
				tempContigSequence = tempContigSequence->next;
				if(tempContigSequence == NULL){
					token = true;
					break;
				}
			}
			if(token == true){
				return false;
			}
			tempScaffoldSet->sequenceCount++;
			return true;
		}
		tempScaffoldSet = tempScaffoldSet->next;
	}
	return false;

}

void InsertShortContigToSequence(ContigSequence * tempContigSequence, ContigSequence * preContigSequence, char * file, char * line, int maxSize){
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	const char * split = ","; 
	char * p; 
	
	int maxCount = 100;
	int contigIndex[maxCount];
	int distance[maxCount];
	bool orientation[maxCount];
	int overlapLength[maxCount];
	
	ScaffoldSetHead * tempInsertSequenceHead = (ScaffoldSetHead *)malloc(sizeof(ScaffoldSetHead));
	tempInsertSequenceHead->scaffoldSet = NULL;
	cout<<"fff"<<endl;
	while((fgets(line, maxSize, fp)) != NULL){ 
		p = strtok(line,split);

		int count = atoi(p);
		
		if(count <= 2){
			continue;
		}
		
		if(count > maxCount){
			continue;
		}
		for(int i = 0; i < maxCount; i++){
			contigIndex[i] = -1;
			distance[i] = -1;
			orientation[i] = false;
			overlapLength[i] = -1;
		}
		
		int a = 1;
		while(a <= count){
			p = strtok(NULL,split);
			contigIndex[a - 1] = atoi(p);
			p = strtok(NULL,split);
			distance[a - 1] = atoi(p);
			p = strtok(NULL,split);
			orientation[a - 1] = atoi(p);
			p = strtok(NULL,split);
	 		overlapLength[a - 1] = atoi(p);
			a++;
		}
		int endIndex = -1;
		int startIndex = -1;
		for(int i = 0; i < count; i++){
			if(contigIndex[i] == tempContigSequence->index){
				endIndex = i;
			}
			if(contigIndex[i] == preContigSequence->index){
				startIndex = i;
			}
			
		}
		if(startIndex == -1 || endIndex == -1 || labs(startIndex - endIndex) == 1){
			continue;
		}

		bool token = SearchInsertContigSequence(tempInsertSequenceHead, startIndex, endIndex, contigIndex, distance, orientation, overlapLength, count);

		if(token == true){
			continue;
		}

		if(orientation[startIndex] == preContigSequence->orientation && orientation[endIndex] == tempContigSequence->orientation && endIndex > startIndex){
			ScaffoldSet * tempScaffoldSet = (ScaffoldSet *)malloc(sizeof(ScaffoldSet));
			tempScaffoldSet->contigSequence = NULL;
			tempScaffoldSet->sequenceCount = 1;
			tempScaffoldSet->contigNum = endIndex -startIndex + 1;
			tempScaffoldSet->next = NULL;
			ContigSequence * tempContigSequence0 = tempScaffoldSet->contigSequence;
			for(int i = startIndex; i <= endIndex; i++){
				ContigSequence * tempContigSequence1 = (ContigSequence *)malloc(sizeof(ContigSequence));
				tempContigSequence1->index = contigIndex[i];
				tempContigSequence1->orientation = orientation[i];
				tempContigSequence1->gapDistance = distance[i];
				tempContigSequence1->next = NULL;
				if(tempContigSequence0 != NULL){
					tempContigSequence0->next = tempContigSequence1;
				}else{
					tempScaffoldSet->contigSequence = tempContigSequence1;
				}
				tempContigSequence0 = tempContigSequence1;
			}
			tempScaffoldSet->next = tempInsertSequenceHead->scaffoldSet;
			tempInsertSequenceHead->scaffoldSet = tempScaffoldSet;
		}

		if(orientation[startIndex] != preContigSequence->orientation && orientation[endIndex] != tempContigSequence->orientation && startIndex > endIndex){
			ScaffoldSet * tempScaffoldSet = (ScaffoldSet *)malloc(sizeof(ScaffoldSet));
			tempScaffoldSet->contigSequence = NULL;
			tempScaffoldSet->sequenceCount = 1;
			tempScaffoldSet->contigNum = startIndex - endIndex + 1;
			tempScaffoldSet->next = NULL;
			ContigSequence * tempContigSequence0 = tempScaffoldSet->contigSequence;
			for(int i = startIndex; i <= endIndex; i--){
				ContigSequence * tempContigSequence1 = (ContigSequence *)malloc(sizeof(ContigSequence));
				tempContigSequence1->index = contigIndex[i];
				tempContigSequence1->orientation = !orientation[i];
				if(i == endIndex){
					tempContigSequence1->gapDistance = 0;
				}else{
					tempContigSequence1->gapDistance = distance[i - 1];
				}
				tempContigSequence1->next = NULL;
				if(tempContigSequence0 != NULL){
					tempContigSequence0->next = tempContigSequence1;
				}else{
					tempScaffoldSet->contigSequence = tempContigSequence1;
				}
				tempContigSequence0 = tempContigSequence1;
			}
			tempScaffoldSet->next = tempInsertSequenceHead->scaffoldSet;
			tempInsertSequenceHead->scaffoldSet = tempScaffoldSet;
		}
	}
	
	

	ScaffoldSet * tempScaffoldSet = tempInsertSequenceHead->scaffoldSet;
	maxCount = -1;
	while(tempScaffoldSet != NULL){
		if(tempScaffoldSet->sequenceCount > maxCount){
			maxCount = tempScaffoldSet->sequenceCount;
		}
		tempScaffoldSet = tempScaffoldSet->next;
	}
	
	tempScaffoldSet = tempInsertSequenceHead->scaffoldSet;
	ContigSequence * tempContigSequence1 = NULL;
	int largeNum = 0;
	while(tempScaffoldSet != NULL){
		if(tempScaffoldSet->sequenceCount == maxCount){
			tempContigSequence1 = tempScaffoldSet->contigSequence;
			largeNum++;
		}
		tempScaffoldSet = tempScaffoldSet->next;
	}

	if(largeNum == 1 && tempContigSequence1 != NULL){
		tempContigSequence1 = tempContigSequence1->next;
		ContigSequence * first = preContigSequence;
		while(tempContigSequence1 != NULL && tempContigSequence1->index != tempContigSequence->index){
			ContigSequence * tempContigSequence2 = (ContigSequence *)malloc(sizeof(ContigSequence));
			tempContigSequence2->index = tempContigSequence1->index;
			tempContigSequence2->orientation = tempContigSequence1->orientation;
			tempContigSequence2->gapDistance = tempContigSequence1->gapDistance;
			tempContigSequence2->next = tempContigSequence;
			
			first->next = tempContigSequence2;
			first = tempContigSequence2;
			tempContigSequence1 = tempContigSequence1->next;
		}
		first->next = tempContigSequence;
	}
	
}

void InsertRepeatContigToSequence1(ContigSetHead * contigSetHead, ScaffoldSetHead * scaffoldSetHead, bool * printIndex, char * file, char * line, int maxSize){
	
	ScaffoldSet * tempScaffoldSet = scaffoldSetHead->scaffoldSet;
	int i = 0;
	while(tempScaffoldSet != NULL){
		i++;
		if(tempScaffoldSet->contigNum < 2){
			tempScaffoldSet = tempScaffoldSet->next;
			continue;
		}
		ContigSequence * tempContigSequence = tempScaffoldSet->contigSequence;
		ContigSequence * preContigSequence = NULL;
        if(tempContigSequence == NULL){
			tempScaffoldSet = tempScaffoldSet->next;
           	 continue;
       	}
		
		while(tempContigSequence != NULL){
			if(preContigSequence != NULL){
				InsertShortContigToSequence(tempContigSequence, preContigSequence, file, line, maxSize);
			}
			preContigSequence = tempContigSequence;
			tempContigSequence = tempContigSequence->next;
		}
		
		tempScaffoldSet = tempScaffoldSet->next;
	}
	
	
}


void InsertRepeatContigToSequence(ContigSetHead * contigSetHead, ScaffoldSetHead * scaffoldSetHead, bool * printIndex, char * file, char * line, int maxSize){
	
	FILE * fp; 
    if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	const char * split = ","; 
	char * p; 
	
	int maxCount = 1000;
	int contigIndex[maxCount];
	int distance[maxCount];
	bool orientation[maxCount];
	int overlapLength[maxCount];
	
	for(int i = 0; i < contigSetHead->contigCount; i++){
		contigSetHead->visited[i] = false;
		contigSetHead->repeatContigIndex[i] = false;
	}
	int lineCount = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		p = strtok(line,split);
		lineCount++;
		int count = atoi(p);
		
		if(count <= 2){
			continue;
		}
		
		for(int i = 0; i < maxCount; i++){
			contigIndex[i] = -1;
			distance[i] = -1;
			orientation[i] = false;
			overlapLength[i] = -1;
		}
		
		if(count > maxCount){
			continue;
		}
		
		int a = 1;
		while(a <= count){
			p = strtok(NULL,split);
			contigIndex[a - 1] = atoi(p);
			p = strtok(NULL,split);
			distance[a - 1] = atoi(p);
			p = strtok(NULL,split);
			orientation[a - 1] = atoi(p);
			p = strtok(NULL,split);
	 		overlapLength[a - 1] = atoi(p);
			a++;
		}
		
		ScaffoldSet * tempScaffoldSet = scaffoldSetHead->scaffoldSet;
		while(tempScaffoldSet != NULL){
        
			ContigSequence * tempContigSequence = tempScaffoldSet->contigSequence;
			ContigSequence * preContigSequence = NULL;
        	if(tempContigSequence == NULL){
				tempScaffoldSet = tempScaffoldSet->next;
           	 	continue;
       	 	}
			
        	bool token = false;
        	while(tempContigSequence != NULL){
				token = false;
				if(tempContigSequence->index == contigIndex[0]){
					if(tempContigSequence->next != NULL && tempContigSequence->next->index == contigIndex[count - 1]){
						if((orientation[0] == orientation[count - 1] && tempContigSequence->orientation == tempContigSequence->next->orientation
						   && ((orientation[0] == 1 && (contigSetHead->visited[contigIndex[0]] == false || contigSetHead->repeatContigIndex[contigIndex[count - 1]] == false))
			 				|| (orientation[0] == 0 && (contigSetHead->repeatContigIndex[contigIndex[0]] == false || contigSetHead->visited[contigIndex[count - 1]] == false))))
						  	|| (orientation[0] != orientation[count - 1] && tempContigSequence->orientation != tempContigSequence->next->orientation 
							   && ((orientation[0] == 1 && (contigSetHead->visited[contigIndex[0]] == false || contigSetHead->visited[contigIndex[count - 1]] == false))
							|| (orientation[0] == 0 && (contigSetHead->repeatContigIndex[contigIndex[0]] == false || contigSetHead->repeatContigIndex[contigIndex[count - 1]] == false))))){
							ContigSequence * tempContigSequence2 = tempContigSequence;
							for(int i = 1; i < count - 1; i++){
								printIndex[contigIndex[i]] = true;
								ContigSequence * tempContigSequence1 = (ContigSequence * )malloc(sizeof(ContigSequence));
            					tempContigSequence1->index = contigIndex[i];
								tempContigSequence1->gapDistance = distance[i];
            					tempContigSequence1->orientation = orientation[i];
								tempContigSequence1->next = tempContigSequence2->next;
								tempContigSequence2->next = tempContigSequence1;
								tempContigSequence2 = tempContigSequence1;
							}
							tempContigSequence->gapDistance = distance[0];
							token = true;
							if(orientation[0] == 1){
								contigSetHead->visited[contigIndex[0]] = true;
							}else{
								contigSetHead->repeatContigIndex[contigIndex[0]] = true;
							}
							if(orientation[count - 1] == 1){
								contigSetHead->visited[contigIndex[count - 1]] = true;
							}else{
								contigSetHead->repeatContigIndex[contigIndex[count - 1]] = true;
							}
							
						}
						
					
					}
					if(preContigSequence != NULL && preContigSequence->index == contigIndex[count - 1]){
						if((orientation[0] == orientation[count - 1] && tempContigSequence->orientation == preContigSequence->orientation
						   && ((orientation[0] == 1 && (contigSetHead->visited[contigIndex[0]] == false || contigSetHead->repeatContigIndex[contigIndex[count - 1]] == false))
			 				|| (orientation[0] == 0 && (contigSetHead->repeatContigIndex[contigIndex[0]] == false || contigSetHead->visited[contigIndex[count - 1]] == false))))
						  	|| (orientation[0] != orientation[count - 1] && tempContigSequence->orientation != preContigSequence->orientation
							   && ((orientation[0] == 1 && (contigSetHead->visited[contigIndex[0]] == false || contigSetHead->visited[contigIndex[count - 1]] == false))
							|| (orientation[0] == 0 && (contigSetHead->repeatContigIndex[contigIndex[0]] == false || contigSetHead->repeatContigIndex[contigIndex[count - 1]] == false))))){
							ContigSequence * tempContigSequence2 = tempContigSequence;
							for(int i = 1; i < count - 1; i++){
								printIndex[contigIndex[i]] = true;
								ContigSequence * tempContigSequence1 = (ContigSequence * )malloc(sizeof(ContigSequence));
            					tempContigSequence1->index = contigIndex[i];
								tempContigSequence1->gapDistance = distance[i - 1];
            					tempContigSequence1->orientation = !orientation[i];
								tempContigSequence1->next = tempContigSequence2;
								preContigSequence->next = tempContigSequence1;
								tempContigSequence2 = tempContigSequence1;
							}
							preContigSequence->gapDistance = distance[count - 2];
							token = true;
							if(orientation[0] == 1){
								contigSetHead->visited[contigIndex[0]] = true;
							}else{
								contigSetHead->repeatContigIndex[contigIndex[0]] = true;
							}
							if(orientation[count - 1] == 1){
								contigSetHead->visited[contigIndex[count - 1]] = true;
							}else{
								contigSetHead->repeatContigIndex[contigIndex[count - 1]] = true;
							}
						}
					
					}
					break;
				}
				if(token == true){
					break;
				}
				preContigSequence = tempContigSequence;
				tempContigSequence = tempContigSequence->next;
			}
			
			tempScaffoldSet = tempScaffoldSet->next;
		}
	
	}
	
	
}


ScaffoldGraphNode * GetOptimizeNodeIndex(ScaffoldGraphHead * scaffoldGraphHead, int nodeIndex, bool orientation, ContigSequence * tempContigSequence, bool right, bool * printIndex){
	ScaffoldGraphNode * temp = NULL;
	ScaffoldGraphNode * tempOriginal = NULL;
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
	
	return temp;
	
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
    strcat(scaffoldTagFileName, "_tag.fa");
    ofstream ocoutTag;
    ocoutTag.open(scaffoldTagFileName);
    
    char * scaffoldSetFileName = new char[1000];   
    strcpy(scaffoldSetFileName, result);
    strcat(scaffoldSetFileName, "_set.fa");
    ofstream ocout1;
    ocout1.open(scaffoldSetFileName);
    j = 0;
	
	int tempLength = 0;
	int tempGapDis = 0;
	
    while(tempScaffoldSet != NULL){
        
		ContigSequence * tempContigSequence = tempScaffoldSet->contigSequence;
        if(tempContigSequence == NULL){
            ocoutTag<<endl;
			tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        ocout1<<">"<<j<<endl;
        
        int allLength = 0;
        tempLength = 0;
		tempGapDis = 0;
        while(tempContigSequence != NULL){
            
            if(printContigIndex[tempContigSequence->index] == true){
                ocoutTag<<"--";
            }
           	tempLength = tempGapDis + tempLength + strlen(contigSet[tempContigSequence->index].contig);
            printContigIndex[tempContigSequence->index] = true;
            ocoutTag<<tempContigSequence->index<<"("<<tempContigSequence->gapDistance<<"--"<<tempContigSequence->orientation<<"--"<<tempLength<<"),";
            tempGapDis = tempContigSequence->gapDistance;
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
