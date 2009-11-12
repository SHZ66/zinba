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
void buildWindows(char **RexpSeqFile,char **RinSeqFile,char **RalignDir,char **RtwoBitFile,int *RzWinSize,int *RzOffsetSize,int *RcWinSize,int *RcOffsetSize,char **Rfilelist){

	string expSeqFile = RexpSeqFile[0];
	string inSeqFile = RinSeqFile[0];
	string alignDir = RalignDir[0];
	string twoBitFile = RtwoBitFile[0];
	int zWinSize = RzWinSize[0];
	int zOffsetSize = RzOffsetSize[0];
	int cWinSize = RcWinSize[0];
	int cOffsetSize = RcOffsetSize[0];
	string filelist = Rfilelist[0];
	
		int ret;
		calcCovs newAnalysis;// = new analysis;
		Rprintf("\nImporting reads from file %s \n", expSeqFile.c_str());
		ret=newAnalysis.importRawSignal(expSeqFile.c_str(),0);

		if(ret == 1){
			Rprintf("ERROR opening file %s \n",expSeqFile.c_str());
			exit(1);
		}

		Rprintf("\nBuilding window data\n");
		size_t found = expSeqFile.find_last_of(".");
		string outfile_prefix = expSeqFile.substr(0,found);
		ret = newAnalysis.processSignals(zWinSize,zOffsetSize,cWinSize,cOffsetSize,alignDir,twoBitFile,inSeqFile.c_str(),outfile_prefix.c_str(),filelist);
		Rprintf("\n\n--------BUILD WINDOWS COMPLETE-------\n\n");
}
}