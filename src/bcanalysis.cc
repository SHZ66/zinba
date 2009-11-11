#include <iostream>
#include <fstream>
#include <stdio.h>
#include "bcanalysis.h"
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;

bcanalysis::bcanalysis(){
	chromCounter = 0;
}

bcanalysis::~bcanalysis(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int bcanalysis::processSignals(const char* outputFile,const char *twoBitFile,int extend){
	
	FILE * tempTB;
	const char * tInfo = "tempInfo.txt"; 
	const char * tChrSize = "tempChromSize.txt";
	tempTB = fopen(tInfo,"w");
	fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
	fclose (tempTB);
	
	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
	int s = system("R CMD BATCH tempInfo.txt /dev/null");
	if(s == 0){
		remove(tInfo);
	}else{
		cout << "twoBitInfo failed\n";
		exit(1);
	}
	
	tempTB = fopen(tChrSize,"r");
	char tbChrom[128];
	unsigned long int tbStart;
	while(!feof(tempTB)){
		int ret = fscanf(tempTB,"%s%lu",tbChrom,&tbStart);
		unsigned short int chromIntTB = getHashValue(tbChrom);
		chr_size[chromIntTB] = tbStart;
	}
	fclose(tempTB);
	remove(tChrSize);
	
	slist<aRead>::iterator i = signal_slist.begin();
	unsigned short int currchr = i->chrom;
	const char * chromReport;
	chromReport = getKey(currchr);
	cout << "\tProcessing " << chromReport << ".........." << endl;
	unsigned short int * basepair = NULL;
	basepair = new unsigned short int[chr_size[currchr]+1];
	for(int in = 0; in <= chr_size[currchr]; in++)
		basepair[in] = 0;
//	basepair[chr_size[currchr]] = 0;
	unsigned short int printFLAG = 0;
	
	while(!signal_slist.empty()){
		if(i == signal_slist.end()){
			if(outputData(outputFile,currchr,printFLAG,basepair) != 0){
				cout << "Error printing output to file, exiting" << endl;
				exit(1);
			}
			printFLAG = 1;
			delete [] basepair;
			i = signal_slist.begin();
			currchr = i->chrom;
			chromReport = getKey(currchr);
			basepair = new unsigned short int[chr_size[currchr]+1];
			for(int in = 0; in <= chr_size[currchr]; in++)
				basepair[in] = 0;
//			basepair[chr_size[currchr]+1] = 0;
			cout << "\tProcessing " << chromReport << ".........." << endl;
		}

		if(i->chrom==currchr){
			unsigned long int aStart = 1;
			unsigned long int aStop = chr_size[currchr];
			if(i->start > extend)
				aStart = i->start - extend;
			if((i->start+extend) < chr_size[currchr])
				aStop = i->start + extend;
			for(unsigned long int pos = aStart; pos <= aStop;pos++){
				basepair[pos]++;
			}
			signal_slist.erase(i++);
		}else{
			i++;
		}
	}

	if(outputData(outputFile,currchr,printFLAG,basepair) != 0){
		cout << "Error printing output to file, exiting" << endl;
		exit(1);
	}
	delete [] basepair;
	return 0;
}

int bcanalysis::outputData(const char * outputFile, unsigned short int currChr,unsigned short int pFLAG,unsigned short int basepair[]){
	FILE * fh;
	if(pFLAG == 0){
		fh = fopen(outputFile,"w");
		fprintf(fh,"track type=wiggle_0 name=\"%s\" desc=\"%s\" visibility=full\n",outputFile,outputFile);
	}else{
		fh = fopen(outputFile,"a");
	}
	const char * chrom = getKey(currChr);
	fprintf(fh,"fixedStep chrom=%s start=1 step=1\n",chrom);

	for(int posP = 1; posP <= chr_size[currChr];posP++){
		fprintf(fh,"%hu\n",basepair[posP]);
	}
	fclose (fh);
	return 0;
}

unsigned short int bcanalysis::getHashValue(char *currChrom){
	map<const char*, int>::iterator i;
	i = chroms.find(currChrom);
	if(i == chroms.end()){
		char * chromosome = new char[128];
		strcpy(chromosome,currChrom);
		chroms[chromosome] = chromCounter;
		intsToChrom[chromCounter] = chromosome;
		return(chromCounter++);
	}else{
		return i->second;	
	}
}

const char * bcanalysis::getKey(unsigned short int chrom){
	map<int, const char*>::iterator i;
	i = intsToChrom.find(chrom);
	if(i == intsToChrom.end()){
		cout << chrom << endl;
		cout << "REALLY REALLY BAD ERROR!" << endl;
		exit(1);
	}else{
		return i->second;	
	}
}

int bcanalysis::importRawSignal(const char * signalFile){
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){return 1;}
	
	unsigned long int lineCount = 0;
	char cChrom[128];
	unsigned long int iStart;
	slist<aRead>::iterator back =  signal_slist.previous(signal_slist.end());

	while(!feof(fh)){
		fscanf(fh,"%s%lu",cChrom,&iStart);
		unsigned short int chromInt = getHashValue(cChrom);
		aRead sig(chromInt,iStart);
		lineCount++;	
		back = signal_slist.insert_after(back,sig);
	}
	fclose(fh);
	cout << "\tLoaded " << lineCount << " reads" << endl;
	signal_slist.sort();
	return 0;
}
