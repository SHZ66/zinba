#include <iostream>
#include <fstream>
#include <stdio.h>
#include "process_scores.h"
#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;

process_scores::process_scores(){
}

process_scores::~process_scores(){
}

int process_scores::adjustCoords(const char * alignFile,const char * outputFile,int aThresh,int adjustSize){
	
	unsigned int map_count;
	unsigned long int countBases = 0;
	unsigned short int * align_count = NULL;
	int cSize = 250000000;
	align_count = new unsigned short int[cSize];
	align_count[cSize] = 0;
	
	string line;
	unsigned short int * profile = NULL;

	cout << "Loading align count data......";
	ifstream seqfile(alignFile);
	while (getline(seqfile, line)){
		if (line[0] == 'f'){
		}else{
			countBases++;
			if(atoi(line.c_str()) <= aThresh && atoi(line.c_str()) > 0){
				align_count[countBases] = 1;
			}
		}
	}

	seqfile.close();
	profile = new unsigned short int[countBases];
	profile[countBases] = 0;
	
	cout << "\nCalculating adjustments......";
	for(int i = (adjustSize+1);i <= (countBases - adjustSize); i++){
		if(align_count[i] == 1){
			profile[(i-adjustSize)]++;
			profile[(i+adjustSize)]++;
		}
	}
	delete [] align_count;
	
	cout << "\nPrinting output......";
	FILE * fh;
	fh = fopen(outputFile,"w");
	for(int pos = 1; pos <= countBases;pos++){
		fprintf(fh,"%i\n",profile[pos]);
	}
	fclose (fh);
	delete [] profile;
	return 0;
}