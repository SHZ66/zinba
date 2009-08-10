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

/*void print_ops(){s
		cout << "\n\nUSAGE: getSeqCountProfile inputFile_COUNT.wig coordFile.tsv outputFile genomeBuild [chrm]\n";
		cout << "\n\tInput file must be formatted as a fixedStep wig file\n";
		cout << "\n\tCoordinate file must be TAB separated with the following columns->\n\t\tID CHRM START STOP STRAND\n";
		cout << "\tStrandedness is considered, use all '+' if you want it ignored\n";
		cout << "\n\tgenomeBuild is the UCSC build (and prefix of chromosome length file)\n";
		cout << "\n\tchrm is optional if left blank will do all chromosomes\n";
		cout << "\n\tHeader lines in coordFile.tsv should start with # and are ignored\n";
		cout << "\n\tLooks for envirnoment variable HG_BFADIR for path to chromosome length file\n";
		cout << "\tchromosome length file should be named hg18_c_chromLength.txt\n";
		cout << "\tFORMAT is tab-delimited CHRM START STOP\n\n\n";
}



int main(int argc=0, char **argv=NULL){

	string chrmSearch;
	if(argc < 5){
		print_ops();
		exit(1);
	}else if (argc == 5){
		cout << "\nUsing all chromosomes\n";
		chrmSearch = "all";
	}else if (argc == 6){
		chrmSearch = argv[5];
	}
	string inputFile = "testSeqDataset1_COUNT.wig";
	string coordFile = "testCoords1.tsv";
	string outputFile = "testCoords1_COUNT.tsv";
	string genomeBuild = "test";*/
extern "C" {
void getSeqCountProfile(char **Rinputfile, char **Rcoordfile, char **Routputfile, char **Rchromosome){
	char* inputFile =Rinputfile[0];
	char* coordFile = Rcoordfile[0];
	char* outputFile = Routputfile[0];
	char* chromosome = Rchromosome[0];	
	int ret;
	analysis newAnalysis;// = new analysis;
	/*cout << "\nImporting coordinates from file " << coordFile << "..." << endl;*/
	Rprintf("\nImporting coordinates from file %s\n",coordFile);
	ret=newAnalysis.importCoords(coordFile);

	if(ret == 1){
		/*cout << "ERROR opening file " << inputFile;*/
		Rprintf("ERROR opening file: %s\n",inputFile);
		exit(1);
	}else if(ret == 2){
		/* cout << "FILE FORMATTING ERROR- wrong number of columns, exiting" << endl;*/
		Rprintf("FILE FORMATTING ERROR- wrong number of columns, exiting");
		exit(1);
	}

	/*cout << "FINISHED importing coordinates" << endl;
	cout << "\nGetting profiles for coordinates " << inputFile << endl;
	cout << "Printing output to " << outputFile << endl;*/

	Rprintf("FINISHED importing coordinates");
	Rprintf("\nGetting profiles for coordinates %s\n",inputFile);
	Rprintf("Printing output to %s \n",outputFile);
 
	ret = newAnalysis.processCoords(inputFile,outputFile,chromosome);		
	if(ret == 1){
		/*cout <<"ERROR: Environment Variable HG_BFADIR not declared, exiting" << endl;*/
		Rprintf("ERROR: Environment Variable HG_BFADIR not declared, exiting");
		exit(1);	
	}else if (ret == 2){
		/*cout << "Unable to open chrom info file, exiting" << endl;*/
		Rprintf("Unable to open chrom info file, exiting");
		exit(1);
	}
	/*cout << "\nFinished all coordinates, COMPLETE\n";
	return 0;*/
	Rprintf("\nFinished all coordinates, COMPLETE\n");
}
}
