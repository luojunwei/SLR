#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <malloc.h>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include "contig.h"

using namespace std;


int main(int argc,char** argv) {
	char * contigSetFile = NULL;
	char * sortContigSetFile = NULL;
	int ch = 0;
	while ((ch = getopt(argc, argv, "c:s:")) != -1) {
		switch (ch) {
			case 'c': contigSetFile = (char *)(optarg); break;
			case 's': sortContigSetFile = (char *)(optarg); break;
			default: break; 
		}
	}
	SortContigSet(contigSetFile, sortContigSetFile);
}
