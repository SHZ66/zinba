#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include "bcanalysis.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C"{
void baseAlignCounts(char **RinputFile,char **RoutputFile, char **Rtwobitfile,int *RextendLength,char **Rfiletype){
	const char * inputFile = RinputFile[0];
	const char * outputFile = RoutputFile[0];
	const char * twobitfile = Rtwobitfile[0];
	int extendLength = RextendLength[0];
	const char * filetype = Rfiletype[0];
	
	bcanalysis newAnalysis;// = new analysis;
	Rprintf("\nImporting reads from file %s ....\n",inputFile);
	Rprintf("Reads are formatted as %s ....\n",filetype);
	int ret=newAnalysis.importRawSignal(inputFile,filetype);
	if(ret == 0){
		Rprintf("Calculating counts at each base\nPrinting output to %s\n",outputFile);
		ret = newAnalysis.processSignals(outputFile,twobitfile,extendLength);
	}
	Rprintf("-------- BASE ALIGN COUNTS COMPLETE --------\n");

}
}
