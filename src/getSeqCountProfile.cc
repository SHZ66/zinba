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
void getSeqCountProfile(char **Rinputfile, char **Rcoordfile, char **Routputfile, char **Rtwobitfile,char **Rchromosome){
		const char*  inputFile = Rinputfile[0];
		const char* coordFile = Rcoordfile[0];
		const char* outputFile = Routputfile[0];
		const char* twobitfile = Rtwobitfile[0];
		const char* chromosome = Rchromosome[0];
		
		analysis newAnalysis;// = new analysis;
		int ret=newAnalysis.importCoords(coordFile);

		if(ret == 1){
			Rprintf("ERROR opening file: %s\n",coordFile);
		}else{
			Rprintf("Getting basecount data for %s\n",chromosome);
			int retP = newAnalysis.processCoords(inputFile,outputFile,twobitfile,chromosome);
			if(retP == 1){
				Rprintf("\nError occurred in processing\n");
			}
			Rprintf("\ngetSeqCountProfile COMPLETE\n");
		}
}
}
