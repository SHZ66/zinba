#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "analysis.h"
#include <ext/slist>
#include <ext/slist>
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
	void getSeqCountProfile(char **Rinputfile, char **Rcoordfile, char **Routputfile, char **Rchromosome){
		char* inputFile =Rinputfile[0];
		char* coordFile = Rcoordfile[0];
		char* outputFile = Routputfile[0];
		char* chromosome = Rchromosome[0];	
		int ret;
		analysis newAnalysis;// = new analysis;
		Rprintf("\nImporting coordinates from file %s\n",coordFile);
		ret=newAnalysis.importCoords(coordFile);

		if(ret == 1){
			Rprintf("ERROR opening file: %s\n",inputFile);
			return;
		}else if(ret == 2){
			Rprintf("FILE FORMATTING ERROR- wrong number of columns, exiting");
			return;
		}

		Rprintf("FINISHED importing coordinates");
		Rprintf("\nGetting profiles for coordinates %s\n",inputFile);
		Rprintf("Printing output to %s \n",outputFile);
	 
		ret = newAnalysis.processCoords(inputFile,outputFile,chromosome);	
		if(ret == 1){
			Rprintf("\nError occurred in processing\n");
		}
		Rprintf("\nFinished all coordinates, COMPLETE\n");
	}
}
