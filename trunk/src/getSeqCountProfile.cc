#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include "analysis.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
void getSeqCountProfile(char **Rinputfile,char **Rwinlist,double *Rthreshold,char **Rmethod,int *Rwformat,char **Routputfile, char **Rtwobitfile,char **Rchromosome){
	const char*  inputFile = Rinputfile[0];
	const char* winlist = Rwinlist[0];
	const char* outputFile = Routputfile[0];
	double threshold = Rthreshold[0];
	const char* method = Rmethod[0];
	int wformat = Rwformat[0];
	const char* twobitfile = Rtwobitfile[0];
	const char* chromosome = Rchromosome[0];
	
	analysis newAnalysis;// = new analysis;
	int ret=newAnalysis.importCoords(winlist,threshold,method,wformat);
		
	Rprintf("Getting basecount data for %s\n",chromosome);
	ret = newAnalysis.processCoords(inputFile,outputFile,twobitfile,chromosome);
	if(ret == 1)
		Rprintf("\nERROR occurred in processing\n");
	Rprintf("\ngetSeqCountProfile COMPLETE\n");

}
}
