#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include "ccanalysis.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C"{
void crossCor(char **RinputFile,char **Rtwobitfile,char **Rfiletype){
	string inputFile = RinputFile[0];
	string twobitfile = Rtwobitfile[0];
	string filetype = Rfiletype[0];
	
		ccanalysis newAnalysis;// = new analysis;
		Rprintf("\nImporting reads from file %s ....\n",inputFile.c_str());
		int ret=newAnalysis.importRawSignal(inputFile.c_str(),filetype.c_str());

		if(ret == 1){
			Rprintf("ERROR opening file %s \n",inputFile.c_str());
			exit(1);
		}

		Rprintf("Calculating cross correlation profile\n");
		size_t found = inputFile.find_last_of(".");
		string oPrefix = inputFile.substr(0,found);
		string outfile = oPrefix + ".crosscor";
		ret = newAnalysis.processSignals(outfile.c_str(),twobitfile.c_str());
		Rprintf("-------- CROSS CORRELATION COMPLETE --------\n");

}
}
