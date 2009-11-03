#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "analysis.h"
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

using namespace std;

extern "C" {
	void getSeqCountProfile(char **Rinputfile, char **Rcoordfile, char **Routputfile, char **Rchromosome){
		string inputFile = Rinputfile[0];
		string coordFile = Rcoordfile[0];
		string outputFile = Routputfile[0];
		string chromosome = Rchromosome[0];
		
		analysis newAnalysis;// = new analysis;
		Rprintf("\nImporting coordinates from file %s\n",coordFile.c_str());
		int ret=newAnalysis.importCoords(coordFile.c_str());

		if(ret == 1){
			Rprintf("ERROR opening file: %s\n",inputFile.c_str());
		}else if(ret == 2){
			Rprintf("FILE FORMATTING ERROR- wrong number of columns, exiting");
		}else{
			Rprintf("FINISHED importing coordinates");
			Rprintf("\nGetting profiles for coordinates %s\n",inputFile.c_str());
			Rprintf("Printing output to %s \n",outputFile.c_str());
			
			int retP = newAnalysis.processCoords(inputFile.c_str(),outputFile.c_str(),chromosome.c_str());	
			if(retP == 1){
				Rprintf("\nError occurred in processing\n");
			}
			Rprintf("\ngetSeqCountProfile COMPLETE\n");
		}
	}
}
