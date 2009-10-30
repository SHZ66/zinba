#include <iostream>
#include <fstream>
#include <stdio.h>
#include "analysis.h"
#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <list>
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

int analysis::processCoords(const char* inputFile,const char* outputFile,string chrmSearch){

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
	int cSize = 250000000;
	basepair = new unsigned short int[cSize];
	basepair[cSize] = 0;

	char *chChr = (char *) malloc(1024);
	unsigned long int sizeProfile = 0;
	unsigned long int countBases = cSize;
	unsigned long int startOffset = 0;
	ifstream seqfile(inputFile);

	while (getline(seqfile, line)){
		if(line[0] == 't'){
			cout << "\nIgnoring track definition line\n";
		}else if (line[0] == 'f'){
			if(collectData == 1 && getChrmData == 1){
				cout << "Finished, loaded " << countBases << "\nGetting data for coordinates and printing output...";
				for(i = coord_slist.begin();i!=coord_slist.end();++i){
					if(i->chrom == chromInt){
						profile = new unsigned short int[(profile_extend * 2)+1];
						int pIndex = 0;
						long int startPos = i->start-profile_extend+1;
						long int stopPos = i->end+profile_extend;
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
					}
				}
				cout << "Finished\n";
				
				if(coord_slist.empty()){
					cout << "\nFinished all coordinates, COMPLETE\n";
					delete [] basepair;
					basepair = NULL;
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
					strcpy(chChr,chrom.c_str());
					chromInt = getHashValue(chChr);
					
					string allChrm = "all";
					if(strcmp(chrmSearch.c_str(),allChrm.c_str()) == 0){
						getChrmData = 1;
					}else if (strcmp(chrmSearch.c_str(),chrom.c_str()) == 0){
						getChrmData = 1;
					}else{
						getChrmData = 0;
						cout << "\tSkipping " << chrom.c_str() << "\n";
					}
					
					if(getChrmData == 1){
						cout << "\nLoading data for " << chChr << endl;
						basepair[cSize] = 0;
						countBases = 0;
					}
				}else if (field[0] == 's' && field[2] == 'a' && getChrmData == 1){
					string startVal;
					for ( int it= 6; it < field.length(); it++ ){
						startVal += field.at(it);
					}
					startOffset = atoi(startVal.c_str());
					countBases = startOffset - 1;
					cout << "Setting start position to " << startOffset << endl;
					cout << "Getting base counts...";
				}
			}
		}else if(getChrmData == 1){
				collectData = 1;
				countBases++;
				basepair[countBases] = atoi(line.c_str());
				if(countBases > cSize){
					cout << "Print some error, adding more data than basepairs in chrom\n";
				}
		}
	}

	if(getChrmData == 1){
		cout << "Finished, loaded " << countBases << "\nGetting data for coordinates and printing output...";
		for(i = coord_slist.begin();i!=coord_slist.end();++i){
			if(i->chrom == chromInt){
				int pIndex = 0;
				long int startPos = i->start-profile_extend+1;
				long int stopPos = i->end+profile_extend;
				profile = new unsigned short int[(stopPos-startPos)+2];
				profile[(stopPos-startPos)+1] = 0;
				for(int s = startPos; s <= stopPos; s++){
					profile[pIndex] = basepair[s];
					pIndex++;
				}
				coord_slist.erase(i++);
				if(outputData( outputFile,printFlag,i->chrom,startPos,stopPos,pIndex,profile) != 0){
					cout << "Error printing output to file, exiting" << endl;
					return 1;
				}
				delete [] profile;
				printFlag = 1;
			}
		}
	}
	delete [] basepair;
	cout << "Finished\nFinished all coordinates, COMPLETE\n";
	seqfile.close();
	return 0;
}

int analysis::outputData(const char * outputFile,int pFlag,unsigned short int pChrom,long int pStart,long int pStop,int printStop,unsigned short int pProfile[]){
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
	strcpy(winID,chromName);strcat(winID,":");strcat(winID,start);strcat(winID,":");strcat(winID,stop);
	fprintf(fh,"%s\t%s\t%i\t%i\t%s",winID,chromName,pStart,pStop,strand);
	for(int posP = 0; posP < printStop;posP++){
		fprintf(fh,"\t%i",pProfile[posP]);
	}
	fprintf(fh,"\n");
	fclose (fh);
	return 0;
}

unsigned short int analysis::getHashValue(char *currChrom){
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
	char line[1024];

	if(!feof(fh)){
		fpos_t position;
		fgetpos (fh, &position);		
		fgets(line,1024,fh);
		while(strncmp(line, "#", 1)==0){
			fgetpos (fh, &position);
			fgets(line, 1024, fh);
		}
		fsetpos (fh, &position);

		unsigned short int countTabs = 0;
  		string tabCounter = line;
  		string::size_type i = tabCounter.find("\t", 0);
  		while(i!=string::npos){
  			countTabs++;
  			i = tabCounter.find("\t", i+1);
  		}
  		
  		if(countTabs!=5){
  			return 2;
  		}
	 }

	while(!feof(fh)){
		int readResult = fscanf(fh,"%s%s%lu%lu%hu%s",id,cChrom,&iStart,&iEnd,&qFlag,strand);
		if(readResult == 6){
			unsigned short int chromInt = getHashValue(cChrom);
			coord c(id,chromInt,iStart,iEnd,qFlag,strand);
			lineCount++;
			coord_slist.push_back(c);
		}
	}
	fclose(fh);
	cout << lineCount << " coordinates imported\n";
	coord_slist.sort();
	back = coord_slist.begin();
	coord lastCoord = *back;
	back++;
	while(back != coord_slist.end()){
		if(back->start <= (lastCoord.end+1)){
			long unsigned int lcStop = back->end;
			coord_slist.erase(back--);
			back->end = lcStop;
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
	return 0;
}
