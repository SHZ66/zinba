#include <iostream>
#include <fstream>
#include <stdio.h>
#include "coord.h"
#include <sstream>
#include <string>
#include <string.h>
#include <cstring>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <list>
#include <ext/hash_map>

extern "C" {
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;
using namespace __gnu_cxx; 

extern "C" {
void collapse_windows(const char * winlist, const char * method, int wformat, double thresholds[]){
	hash_map<string, int, hash<string>,equal_to<string> > cind_map;
	vector<string> cnames;

	string winstring = string(winlist);
	size_t found = winstring.find_last_of(".");
	string outfile = winstring.substr(0,found);
	
	FILE * wlist;
	wlist = fopen(winlist,"r");
	if(wlist == NULL){error("Unable to open list file: %s\n", winlist);}
	char cChrom[128];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	double sigVal;
	const char *pscl = "pscl";
	const char *mixture = "mixture";
	char sigFile [256];
	int readResult = 0;
	char firstline [256];
	int rwline;
	///////////////////////////
	unsigned long int winSizeThresh;
	int winSize;
	//////////////////////////
	list<coord> coordIN_slist;
	list<coord> coordOUT_slist;
	
	int lengthThresholds = sizeof(thresholds)/sizeof(double);
	double highThresh;
	if(strcmp(method,pscl) == 0){
		highThresh = 0;
		for(int t = 0; t <= lengthThresholds; t++){
			if(thresholds[t] > highThresh)
				highThresh = thresholds[t];
		}
	}else if(strcmp(method,mixture) == 0){
		highThresh = 1;
		for(int t = 0; t <= lengthThresholds; t++){
			if(thresholds[t] < highThresh)
				highThresh = thresholds[t];
		}
	}
	
	while(!feof(wlist)){
		int rline = fscanf(wlist,"%s",sigFile);
		if(rline == 1){
			string signalFile(sigFile);
			FILE * fh;
			fh = fopen(signalFile.c_str(),"r");
			if(fh == NULL){error("Unable to open input file: %s\n", signalFile.c_str());}
			cout << "\tImporting windows from " << signalFile.c_str() << "..." << endl;
			fgets(firstline,256,fh);
			while(!feof(fh)){
				readResult = 0;
				if(strcmp(method,pscl) == 0){
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%lf%*lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal <= highThresh && rwline > 0)
							readResult = 1;
					}else if(wformat == 1){
						rwline = fscanf(fh,"%s%lu%lu%*d%*lf%*lf%*lf%*lf%hu%lf%*lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal <= highThresh && rwline > 0)
							readResult = 1;
					}
				}else if(strcmp(method,mixture) == 0){
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal >= highThresh && rwline > 0)
							readResult = 1;
					}else if(wformat == 1){
						int exp;double inp,gc,ap,ecl;
						rwline = fscanf(fh,"%s%lu%lu%*d%*lf%*lf%*lf%*lf%hu%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal >= highThresh && rwline > 0)
							readResult = 1;
					}
				}
				if(readResult == 1){
					if(qFlag == 1){
						// determine the chromosome index
						chr = string(cChrom);
						hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
						int chromInt=-1;
						if(li==cind_map.end()) {
							// register new chromosome
							chromInt=cnames.size();
							cnames.push_back(chr);
							cind_map[chr]=chromInt;
						}else{
							chromInt=li->second;
						}
						coord c(chromInt,iStart,iEnd,qFlag,sigVal);
						coordIN_slist.push_back(c);
					}
				}
			}
			fclose(fh);
		}
	}
	fclose(wlist);
	winSize = iEnd - iStart;
	winSizeThresh = winSize * 10;
	cout << "\nImported " << coordIN_slist.size() << " coordinates" << endl;
	coordIN_slist.sort();
	list<coord> tempCoords;	
	
	for(int t = 0; t <= lengthThresholds; t++){
		if(strcmp(method,pscl) == 0){
			while(back != coordIN_slist.end()){
				if(back->sigVal <= thresholds[t])
					tempCoords.push_back(*back);
				back++;
			}
		}else if(strcmp(method,mixture) == 0){
			while(back != coordIN_slist.end()){
				if(back->sigVal >= thresholds[t])
					tempCoords.push_back(*back);
				back++;
			}
		}
		tempCoords.sort();
		
		list<coord>::iterator back =  tempCoords.begin();
		coord tempCoord = *back;
		back++;
		while(back != tempCoords.end()){
			if(back->chrom == tempCoord.chrom && back->start <= tempCoord.end+1){
				tempCoord.end = back->end;
				if(strcmp(method,pscl) == 0 && tempCoord.sigVal > back->sigVal)
					tempCoord.sigVal = back->sigVal;
				else if(strcmp(method,mixture) == 0 && tempCoord.sigVal < back->sigVal)
					tempCoord.sigVal = back->sigVal;
			}else{
				if((tempCoord.end-tempCoord.start) <= winSizeThresh){
					coordOUT_slist.push_back(tempCoord);
				}else{
					//REMOVE ONCE C PEAKBOUND IS RUNNING
					string cName = cnames[tempCoord.chrom];
					cout << "Excluding " << cName.c_str() << ":" << tempCoord.start << "-" << tempCoord.end << " SIZE=" << (tempCoord.end-tempCoord.start) << endl;
				}
				tempCoord = *back;
				back++;
			}
		}
		coordOUT_slist.sort();
		cout << "\nCollapsed to " << coordOUT_slist.size() << " non-overlapping regions" << endl;
	
		///output data
		FILE * fh;
		string outputFile = output + "_" + thresholds[t] + ".coords";
		fh = fopen(outputFile.c_str(),"w");
		if(fh==NULL){error("Unable to open output file: %s\n", outputFile);}
		fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND\tSIGVAL\n");
		string chromName;
		char strand[] = "+";
		char winID[255];
		char start[10];
		char stop[10];
		back = coordOUT_slist.begin();
		while(back != coordOUT_slist.end()){
			chromName = cnames[back->chrom];
			sprintf( start,"%d", back->start);
			sprintf( stop,"%d", back->end);
			strcpy(winID,chromName.c_str());strcat(winID,":");strcat(winID,start);strcat(winID,"-");strcat(winID,stop);
			fprintf(fh,"%s\t%s\t%lu\t%lu\t%s\t%.14f\n",winID,chromName.c_str(),back->start,back->end,strand,back->sigVal);
			back++;
		}
		fclose (fh);
		tempCoords.clear();
		coordOUT_slist.clear();
	}
}
}