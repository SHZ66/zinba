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

/*void print_ops(){
		cout << "USAGE: baseAlignCounts inputFile.bed outputFile genomeBuild\n";
		cout << "\tInput file must be formatted -> CHROM START STOP\n";
		cout << "\t\tLines starting with track are ignored\n";
}

int main(int argc=0, char **argv=NULL){

	if(argc < 3){
		print_ops();
		exit(1);
	}else{
		string inputFile = argv[1];
		string outputFile = argv[2];
		int extendLength = argv[3];
*/
extern "C"{
void baseAlignCounts(char **RinputFile,char **RoutputFile, char **Rtwobitfile,int *RextendLength){
	string inputFile = RinputFile[0];
	string outputFile = RoutputFile[0];
	string twobitfile = Rtwobitfile[0];
	int extendLength = RextendLength[0];
	
		bcanalysis newAnalysis;// = new analysis;
		Rprintf("\nImporting reads from file %s ....\n",inputFile.c_str());
		int ret=newAnalysis.importRawSignal(inputFile.c_str());

		if(ret == 1){
			Rprintf("ERROR opening file %s \n",inputFile.c_str());
			exit(1);
		}

		Rprintf("Calculating counts at each base\nPrinting output to %s\n",outputFile.c_str());
		ret = newAnalysis.processSignals(outputFile.c_str(),twobitfile.c_str(),extendLength);
		Rprintf("-------- BASE ALIGN COUNTS COMPLETE --------\n");
//	}
//	return 0;
}
}
