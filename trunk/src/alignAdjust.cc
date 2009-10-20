#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "process_scores.h"
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

using namespace std;

/*void print_ops(){
		cout << "\n\nUSAGE: alignAdjust alignFile_LIST.txt align_threshold adjust_size\n\n";
		cout << "\tFile containing list of alignability files to process\n";
		cout << "\talign_threshold = max number of matches for a read\n";
		cout << "\tadjust_size = half the extension length for aligned reads\n";
		cout << "\n\tOUTPUT FORMAT is _ADJUST.wig file (values will be 0,1,2)\n\n\n";
}*/

//int main(int argc=0, char **argv=NULL){
extern "C" {
void alignAdjust(char **Rinputfile, int *RaThresh, int *RadjustSize){
 
	char* inputFile =Rinputfile[0];
	int* coordFile = RaThresh[0];
	int* outputFile = RadjustSize[0];

/*	if(argc < 4){
		print_ops();
		exit(1);
	}*/
	
	const char * inputFile = argv[1];
	int aThresh = atoi(argv[2]);
	int adjustSize = atoi(argv[3]);
	
//	cout << "\n\nInput file is " << inputFile << "\naThresh is " << aThresh << "\nadjustSize is " << adjustSize << "\n";
	
	string aFile;
	string outFile;
	ifstream listfile(inputFile);	
	while(getline(listfile,aFile)){
//		cout << "\n\nGetting data from " << aFile.c_str() << "\n";
		outFile = aFile;
		outFile.erase((outFile.length()-4),4);
		outFile.insert(outFile.length(), "_ADJUST.wig");
//		cout << "Printing output to " << outFile.c_str() << "\n";
		process_scores newProcess;
		newProcess.adjustCoords(aFile.c_str(),outFile.c_str(),aThresh,adjustSize);
	}
//	cout << "\nFINISHED all files, exiting\n\n\n";
	listfile.close();
//	return 0;
}
}