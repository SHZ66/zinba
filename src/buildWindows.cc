#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "calcCovs.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;
/*
void print_ops(){
		cout << "USAGE: buildWindows:\n";
		cout << "\tseq file\n\tinput file\n\talign directory\n\ttwoBit\n\tz win\n\tz offset\n\tcnv win size\n\tcnv offset\n\tparamfile\n";
}
*/
//int main(int argc=0, char **argv=NULL){
extern "C" {
void buildWindows(char **RexpSeqFile,char **RinSeqFile,char **RalignDir,char **RtwoBitFile,int *RzWinSize,int *RzOffsetSize,int *RcWinSize,int *RcOffsetSize,char **RparamFile){

	string expSeqFile = RexpSeqFile[0];
	string inSeqFile = RinSeqFile[0];
	string alignDir = RalignDir[0];
	string twoBitFile = RtwoBitFile[0];
	int zWinSize = RzWinSize[0];
	int zOffsetSize = RzOffsetSize[0];
	int cWinSize = RcWinSize[0];
	int cOffsetSize = RcOffsetSize[0];
	string paramFile = RparamFile[0];
	
/*	if(argc < 9){
		print_ops();
		exit(1);
	}else{
		string expSeqFile = argv[1];
		string inSeqFile = argv[2];
		string alignDir = argv[3];
		string twoBitFile = argv[4];
		int zWinSize = atoi(argv[5]);
		int zOffsetSize = atoi(argv[6]);
		int cWinSize = atoi(argv[7]);
		int cOffsetSize = atoi(argv[8]);
		string paramFile = argv[9];
*/
		int ret;
		calcCovs newAnalysis;// = new analysis;
		Rprintf("\nImporting reads from file %s \n", expSeqFile.c_str());
		ret=newAnalysis.importRawSignal(expSeqFile.c_str(),0);

		if(ret == 1){
			Rprintf("ERROR opening file %s \n",expSeqFile.c_str());
			exit(1);
		}

		Rprintf("\nBuilding window data\n");
		ret = newAnalysis.processSignals(zWinSize,zOffsetSize,cWinSize,cOffsetSize,alignDir,twoBitFile,paramFile.c_str(),inSeqFile.c_str());
		Rprintf("\n\n--------BUILD WINDOWS COMPLETE-------\n\n");
//	}
//	return 0;
}
}