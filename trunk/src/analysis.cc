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
	list<coord>::iterator i = coordOUT_slist.begin();
	unsigned short int * basepair = NULL;
	unsigned short int * profile = NULL;
	unsigned long int countBases = 250000000;
	unsigned long int startOffset = 0;
	cout << "Getting basecount data from " << inputFile << endl;
	ifstream seqfile(inputFile);
	
	while (getline(seqfile, line)){
		if (line[0] == 'f'){
			if(collectData == 1 && getChrmData == 1){
				i = coordOUT_slist.begin();
				while(i!=coordOUT_slist.end()){
					if(i->chrom == chromInt){
						int pIndex = 0;
						unsigned long int startPos = 1;
						if(i->start > profile_extend)
							startPos = i->start-profile_extend;
						unsigned long int stopPos = chr_size[chromInt];
						if((i->end+profile_extend) < chr_size[chromInt])
							stopPos = i->end+profile_extend;
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
						coordOUT_slist.erase(i++);
					}else{
						i++;
					}
				}
				delete [] basepair;
				basepair = NULL;
				
				if(coordOUT_slist.empty()){
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
						for(int c = 0; c <= chr_size[chromInt];c++)
							basepair[c] = 0;
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
		i = coordOUT_slist.begin();
		while(i!=coordOUT_slist.end()){
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
				coordOUT_slist.erase(i++);
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
		cout << "printing output to " << outputFile << endl;
	}
	
	if(pFlag == 0){
		fh = fopen(outputFile,"w");
		if(fh==NULL){
			cout << "Unable to open output file: " << outputFile << endl;
			return 1;
		}
		fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND");
		for(int p = 1;p <= printStop;p++){
			fprintf(fh,"\tPosition%i",p);
		}
		fprintf(fh,"\n");
	}else if (pFlag == 1){
		fh = fopen(outputFile,"a");
		if(fh==NULL){
			cout << "Unable to open output file: " << outputFile << endl;
			return 1;
		}
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

int analysis::importCoords(const char *winlist,double threshold,const char *method,int wformat){
	FILE * wlist;
	wlist = fopen(winlist,"r");
	if(wlist == NULL){return 1;}
	char cChrom[128];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	long double sigVal;
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
	while(!feof(wlist)){
		int rline = fscanf(wlist,"%s",sigFile);
		if(rline == 1){
			string signalFile(sigFile);
			FILE * fh;
			fh = fopen(signalFile.c_str(),"r");
			if(fh == NULL){return 1;}
			cout << "\tImporting windows from " << signalFile.c_str() << "..." << endl;
			fgets(firstline,256,fh);
			while(!feof(fh)){
				readResult = 0;
				if(strcmp(method,pscl) == 0){
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%Lf%*f",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal <= threshold && rwline > 0){
							readResult = 1;
							winSize = iEnd - iStart;
						}
					}else if(wformat == 1){
						rwline = fscanf(fh,"%s%lu%lu%*d%*f%*f%*f%*f%hu%Lf%*f",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal <= threshold && rwline > 0){
							readResult = 1;
							winSize = iEnd - iStart;
						}
					}
				}else if(strcmp(method,mixture) == 0){
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%Lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal >= threshold && rwline > 0){
							readResult = 1;
							winSize = iEnd - iStart;
						}
					}else if(wformat == 1){
						int exp;double inp,gc,ap,ecl;
						rwline = fscanf(fh,"%s%lu%lu%*d%*f%*f%*f%*f%hu%Lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal >= threshold && rwline > 0){
							readResult = 1;
							winSize = iEnd - iStart;
						}
					}
				}
				if(readResult == 1){
					string chromIn(cChrom);
					unsigned short int chromInt = getHashValue(chromIn.c_str());
					coord c(chromInt,iStart,iEnd,qFlag);
					coordIN_slist.push_back(c);
				}
			}
			fclose(fh);
		}
	}
	fclose(wlist);
	winSizeThresh = winSize * 10;
	cout << "\nImported " << coordIN_slist.size() << " coordinates" << endl;
	coordIN_slist.sort();
	
	list<coord>::iterator back =  coordIN_slist.begin();	
	list<coord> tempcoord_list;
	while(back->qFlag == 0)
		coordIN_slist.erase(back++);
	tempcoord_list.push_front(*back);
	coordIN_slist.erase(back++);
	list<coord>::iterator tempIt = tempcoord_list.begin();
	int flagFinish = 0;
	
	while(flagFinish == 0){
		if(coordIN_slist.empty())
			flagFinish = 1;
		if(flagFinish == 0 && back->qFlag == 1 && back->chrom == tempIt->chrom && back->start <= tempIt->end+1){
			tempcoord_list.push_back(*back);
			tempIt = tempcoord_list.end();
			tempIt--;
			coordIN_slist.erase(back++);
		}else{
			tempIt = tempcoord_list.begin();
			unsigned short int tChrom = tempIt->chrom;
			unsigned long int start = tempIt->start;
			unsigned long int stop = tempIt->end;
			tempIt++;
			while(tempIt != tempcoord_list.end()){
				stop = tempIt->end;
				tempIt++;
			}
			if((stop-start) <= winSizeThresh){
				coord c(tChrom,start,stop,1);
				coordOUT_slist.push_back(c);
			}else{
				//REMOVE ONCE C PEAKBOUND IS RUNNING
				const char * cName = getKey(tChrom);
				cout << "Excluding " << cName << ":" << start << "-" << stop << " SIZE=" << (stop-start) << endl;
			}

			if(flagFinish == 0){
				tempcoord_list.clear();
				if(back->qFlag == 1){
					tempcoord_list.push_front(*back);
					coordIN_slist.erase(back++);
				}else{
					while(back->qFlag == 0)
						coordIN_slist.erase(back++);
					tempcoord_list.push_front(*back);
				}
				tempIt = tempcoord_list.begin();
			}
		}
	}
	coordOUT_slist.sort();
	cout << "\nCollapsed to " << coordOUT_slist.size() << " non-overlapping regions" << endl;
}