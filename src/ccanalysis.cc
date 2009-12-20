#include <iostream>
#include <fstream>
#include <stdio.h>
#include "ccanalysis.h"
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

ccanalysis::ccanalysis(){
	chromCounter = 0;
}

ccanalysis::~ccanalysis(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int ccanalysis::processSignals(const char* outputFile,const char *twoBitFile){
	
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
	
	slist<bwRead>::iterator i = reads_slist.begin();
	unsigned short int currchr = i->chrom;
	int firstBP = 0;
	const char * chromReport;
	chromReport = getKey(currchr);
	cout << "\tProcessing " << chromReport << ".........." << endl;

	double corVals [501];
	for(int v = 0; v <= 500; v++)
		corVals[v] = 0;
	unsigned short int * plusReads = NULL;
	unsigned short int * minusReads = NULL;
	plusReads = new unsigned short int[chr_size[currchr]+1];
	minusReads = new unsigned short int[chr_size[currchr]+1];
cout << "Initializing arrays" << endl;
	for(int in = 0; in <= chr_size[currchr]; in++){
		plusReads[in] = 0;
		minusReads[in] = 0;
	}
cout << "Arrays initialized" << endl;
	double readcountPlus = 0;
	double readcountMinus = 0;
	double totalReads = reads_slist.size();
	
	while(!reads_slist.empty()){
		if(i->chrom==currchr){
			if(i->strand == 1){
				plusReads[i->pos]++;
				readcountPlus++;
				if(firstBP == 0)
					firstBP = i->pos;
			}else if(i->strand == 0){
				minusReads[i->pos]++;
				readcountMinus++;
			}
			reads_slist.erase(i++);
		}else{
			i++;
		}
		
		if(i == reads_slist.end()){
cout << "Analyzing data:" << endl;
			for(int c = 0; c <= 500; c++){
				double numerator = 0;
				double denX = 0;
				double denY = 0;
				for(int bp = firstBP; bp <= (chr_size[currchr]-c); bp++){
					numerator += ((double)(plusReads[bp]-(readcountPlus/chr_size[currchr]))*(double)(minusReads[bp+c]-(readcountMinus/chr_size[currchr])));
					denX += pow((double)(plusReads[bp]-(readcountPlus/chr_size[currchr])),2);
					denY += pow((double)(minusReads[bp+c]-(readcountMinus/chr_size[currchr])),2);
				}
				corVals[c] += ((readcountPlus+readcountMinus)/totalReads)*(numerator/(sqrt(denX)*sqrt(denY)));
cout << "num is " << numerator << " denX is " << denX << " denY is " << denY << endl;
cout << c << " " << corVals[c] << endl;

			}
			readcountPlus = 0;
			readcountMinus = 0;
			delete [] plusReads;plusReads = NULL;
			delete [] minusReads;minusReads = NULL;
			if(!reads_slist.empty()){
				firstBP = 0;
				i = reads_slist.begin();
				currchr = i->chrom;
				chromReport = getKey(currchr);
				plusReads = new unsigned short int[chr_size[currchr]+1];
				minusReads = new unsigned short int[chr_size[currchr]+1];
				for(int in = 0; in <= chr_size[currchr]; in++){
					plusReads[in] = 0;
					minusReads[in] = 0;
				}
				cout << "\tProcessing " << chromReport << ".........." << endl;
			}
		}
	}

	FILE * fh;
	fh = fopen(outputFile,"w");
	for(int c = 0; c <= 500; c++)
		fprintf(fh,"%i%f\n",c,corVals[c]);
	fclose (fh);	
	return 0;
}

unsigned short int ccanalysis::getHashValue(char *currChrom){
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

const char * ccanalysis::getKey(unsigned short int chrom){
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

int ccanalysis::importRawSignal(const char * signalFile,const char * filetype){
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){return 1;}

	char cChrom[128];
	unsigned long int pos;
	unsigned short int sval;
	char strand[1];
	string bowtie = "bowtie";
	string bed = "bed";
	char minus[] = "-";char plus[] = "+";
	slist<bwRead>::iterator back =  reads_slist.previous(reads_slist.end());
	unsigned long int start;unsigned long int stop;
	char line[512];char seq[128];
	
	while(!feof(fh)){
		if(strcmp(filetype,bed.c_str()) == 0){
			fscanf(fh,"%s%lu%lu%s",cChrom,&start,&stop,strand);
			if(strcmp(strand,minus) == 0){
				pos = stop;
				sval = 1;
			}else if(strcmp(strand,plus) == 0){
				pos = start;
				sval = 0;
			}
		}else if (strcmp(filetype,bowtie.c_str()) == 0){
			fscanf(fh,"%*s%s%s%lu%s%*s%*d",strand,cChrom,&pos,seq);
			fgets(line,512,fh);
			sval = 1;
			if(strcmp(strand,minus) == 0){
				pos = pos + strlen(seq);
				sval = 0;
			}
		}
		unsigned short int chromInt = getHashValue(cChrom);
		bwRead sig(chromInt,pos,sval);
		back = reads_slist.insert_after(back,sig);
	}
	fclose(fh);
	cout << "\tLoaded " << reads_slist.size() << " reads" << endl;
	reads_slist.sort();
	return 0;
}
