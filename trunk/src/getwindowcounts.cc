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

extern "C" {
void getWindowCounts(char **RexpSeqFile,char **RtwoBitFile,int *RzWinSize,int *RzOffsetSize, char **Rfiletype,int *Rextension){

	string expSeqFile = RexpSeqFile[0];
	const char * twoBitFile = RtwoBitFile[0];
	int zWinSize = RzWinSize[0];
	int zOffsetSize = RzOffsetSize[0];
	const char * filetype = Rfiletype[0];
	int extension = Rextension[0];
	
	int ret;
	calcCovs newAnalysis;// = new analysis;
	Rprintf("\nImporting reads from file %s \n", expSeqFile.c_str());
	ret=newAnalysis.importRawSignal(expSeqFile.c_str(),extension,filetype,0,twoBitFile);

	if(ret == 1){
		Rprintf("ERROR opening file %s \n",expSeqFile.c_str());
	}else{
		Rprintf("\nBuilding window data\n");
		size_t found = expSeqFile.find_last_of(".");
		string outfile_prefix = expSeqFile.substr(0,found);
		ret = newAnalysis.processWinSignal(zWinSize,zOffsetSize,twoBitFile,outfile_prefix,extension,filetype);
		if(ret == 1){
			Rprintf("ERROR: building windows was unsuccssful\n");
		}
	}
	Rprintf("\n\n--------GET WINDOW COUNTS COMPLETE-------\n\n");
}
}
