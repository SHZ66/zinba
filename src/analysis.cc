#include <iostream>
#include <fstream>
#include <stdio.h>
#include "analysis.h"
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

analysis::analysis(){
	chromCounter = 0;
}

analysis::~analysis(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int analysis::processCoords(const char* inputFile,const char* outputFile,const char* twoBitFile,const char* chrmSearch){

	FILE * tempTB;
	const char * tInfo = "tempInfo.txt"; 
	const char * tChrSize = "tempChromSize.txt";
	tempTB = fopen(tInfo,"w");
	fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
	fclose (tempTB);
	
//	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
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
		if(ret == 2){
			string sChr(tbChrom);
			unsigned short int chromIntTB = getHashValue(sChr.c_str());
			chr_size[chromIntTB] = tbStart;
		}
	}
	fclose(tempTB);
	remove(tChrSize);
	
	unsigned short int chromInt;
	string line;
	string field;
//////////////////////////////
	int profile_extend = 500;
//////////////////////////////
	int printFlag = 0;
	int getChrmData = 0;
	unsigned short int collectData = 0;
	list<coord>::iterator i = coord_slist.begin();
	unsigned short int * basepair = NULL;
	unsigned short int * profile = NULL;
	unsigned long int countBases = 250000000;
	unsigned long int startOffset = 0;
	cout << "Getting basecount data from " << inputFile << endl;
	ifstream seqfile(inputFile);
	
	while (getline(seqfile, line)){
		if (line[0] == 'f'){
			if(collectData == 1 && getChrmData == 1){
				i = coord_slist.begin();
				while(i!=coord_slist.end()){
					if(i->chrom == chromInt){
						int pIndex = 0;
						unsigned long int startPos = 1;
						if(i->start > profile_extend){
							startPos = i->start-profile_extend;
						}
						unsigned long int stopPos = i->end+profile_extend;
						profile = new unsigned short int[(stopPos-startPos)+1];
						profile[(stopPos-startPos)] = 0;
						for(int s = startPos; s <= stopPos; s++){
							profile[pIndex] = basepair[s];
							pIndex++;
						}
						if(outputData( outputFile,printFlag,i->chrom,startPos,stopPos,pIndex,profile) != 0){
							cout << "Error printing output to file, exiting" << endl;
							return 1;
						}
						delete [] profile;
						profile = NULL;
						printFlag = 1;
						coord_slist.erase(i++);
					}else{
						i++;
					}
				}
				delete [] basepair;
				basepair = NULL;
				
				if(coord_slist.empty()){
					seqfile.close();
					return 0;
				}
			}
			
			istringstream iss(line);
			while(iss >> field){
				if (field[0] == 'c'){
					string chrom;
					for ( int it= 6; it < field.length(); it++ ){
						chrom += field.at(it);
					}
					chromInt = getHashValue(chrom.c_str());					
					string allChrm = "all";
					if(strcmp(chrmSearch,allChrm.c_str()) == 0){
						getChrmData = 1;
					}else if (strcmp(chrmSearch,chrom.c_str()) == 0){
						getChrmData = 1;
					}else{
						getChrmData = 0;
						cout << "Skipping " << chrom.c_str() << endl;
					}
					
					if(getChrmData == 1){
						cout << "Loading data for " << chrom.c_str() << endl;
						basepair = new unsigned short int[chr_size[chromInt]+1];
						basepair[chr_size[chromInt]] = 0;
						countBases = 0;
					}
				}else if (field[0] == 's' && field[2] == 'a' && getChrmData == 1){
					string startVal;
					for ( int it= 6; it < field.length(); it++ ){
						startVal += field.at(it);
					}
					startOffset = atoi(startVal.c_str());
					countBases = startOffset - 1;
				}
			}
		}else if(getChrmData == 1 && line[0] != 't'){
				collectData = 1;
				countBases++;
				basepair[countBases] = atoi(line.c_str());
				if(countBases > chr_size[chromInt]){
					cout << "Print some error, adding more data than basepairs in chrom\n";
				}
		}
	}

	if(getChrmData == 1){
		i = coord_slist.begin();
		while(i!=coord_slist.end()){
			if(i->chrom == chromInt){
				int pIndex = 0;
				unsigned long int startPos = 1;
				if(i->start > profile_extend)
					startPos = i->start-profile_extend;
				unsigned long int stopPos = i->end+profile_extend;
				profile = new unsigned short int[(stopPos-startPos)+1];
				profile[(stopPos-startPos)] = 0;
				for(int s = startPos; s <= stopPos; s++){
					profile[pIndex] = basepair[s];
					pIndex++;
				}
				if(outputData( outputFile,printFlag,i->chrom,startPos,stopPos,pIndex,profile) != 0){
					cout << "Error printing output to file, exiting" << endl;
					return 1;
				}
				delete [] profile;
				profile = NULL;
				printFlag = 1;
				coord_slist.erase(i++);
			}else{
				i++;
			}
		}
	}
	seqfile.close();
	delete [] basepair;
	basepair = NULL;
	return 0;
}

int analysis::outputData(const char * outputFile,int pFlag,unsigned short int pChrom,unsigned long int pStart,unsigned long int pStop,int printStop,unsigned short int pProfile[]){
	FILE * fh;
	if(pFlag == 0){
		fh = fopen(outputFile,"w"); 
		fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND");
		for(int p = 1;p <= printStop;p++){
			fprintf(fh,"\tPosition%i",p);
		}
		fprintf(fh,"\n");
	}else if (pFlag == 1){
		fh = fopen(outputFile,"a");
	}
	const char * chromName = getKey(pChrom);
	char strand[] = "+";
	char winID[255];
	char start[10];
	char stop[10];
	sprintf( start,"%d", pStart);
	sprintf( stop,"%d", pStop);
	strcpy(winID,chromName);strcat(winID,":");strcat(winID,start);strcat(winID,"-");strcat(winID,stop);
	fprintf(fh,"%s\t%s\t%lu\t%lu\t%s",winID,chromName,pStart,pStop,strand);
	for(int posP = 0; posP < printStop;posP++){
		fprintf(fh,"\t%i",pProfile[posP]);
	}
	fprintf(fh,"\n");
	fclose (fh);
	return 0;
}

unsigned short int analysis::getHashValue(const char *currChrom){
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

const char * analysis::getKey(unsigned short int chrom){
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

int analysis::importCoords(const char * signalFile){

	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){return 1;}
	
	unsigned long int lineCount = 0;
	char cChrom[128];
	char id[128];
	char strand[1];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	list<coord>::iterator back =  coord_slist.begin();
	while(!feof(fh)){
		int readResult = fscanf(fh,"%s%s%lu%lu%hu%s",id,cChrom,&iStart,&iEnd,&qFlag,strand);
		if(readResult == 6){
			string chromIn(cChrom);
			unsigned short int chromInt = getHashValue(chromIn.c_str());
			coord c(chromInt,iStart,iEnd,qFlag);
			lineCount++;
			coord_slist.push_back(c);
		}
	}
	fclose(fh);
//	cout << lineCount << " coordinates imported" << endl;
	coord_slist.sort();
	back = coord_slist.begin();
	coord lastCoord = *back;
	back++;
	while(back != coord_slist.end()){
		if(back->chrom == lastCoord.chrom && back->start <= (lastCoord.end+1)){
			long unsigned int lcStop = back->end;
			coord_slist.erase(back--);
			back->end = lcStop;
			back->qFlag = 1;
		}
		lastCoord = *back;		
		back++;
	}
	
	back = coord_slist.begin();	
	while(back != coord_slist.end()){
		if(back->qFlag == 0){
			coord_slist.erase(back++);
		}else{
			back++;
		}
	}
	cout << "\nCollapsed windows down to " << coord_slist.size() << " non-overlapping regions" << endl;
	return 0;
}
