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
	int printFlag = 0;
	int getChrmData = 0;
	unsigned short int collectData = 0;
	slist<coord>::iterator i = coord_slist.begin();
	int * basepair = NULL;
	int * profile = NULL;
	int cSize = 250000000;
	basepair = new int[cSize];
	basepair[cSize] = 0;

	char *chChr = (char *) malloc(1024);

	unsigned long int sizeProfile = 0;
////	unsigned long int countBases = 0;
	unsigned long int countBases = cSize;
	
	unsigned long int startOffset = 0;
	ifstream seqfile(inputFile);

	while (getline(seqfile, line)){
		if(line[0] == 't'){
			cout << "\nIgnoring track definition line\n";
		}else if (line[0] == 'f'){
			if(collectData == 1 && getChrmData == 1){
				cout << "Finished, loaded " << countBases << "\nGetting data for coordinates and printing output...";
				for(i = coord_slist.begin();i!=coord_slist.end();i++){
					if(i->chrom == chromInt){
						sizeProfile = i->end - i->start;
//						profile = new double[(sizeProfile + 1)];
						profile = new int[(sizeProfile + 1)];
						int pIndex = 0;
						for(int s = i->start; s <= i->end; s++){
							profile[pIndex] = basepair[s];
							pIndex++;
						}

						if(outputData( outputFile,printFlag,i->ident,i->chrom,i->start,i->end,i->strand,sizeProfile,profile) != 0){
							cout << "Error printing output to file, exiting" << endl;
							return 1;
						}
						delete [] profile;
						sizeProfile = NULL;
						printFlag = 1;
						coord_slist.erase(i);
					}
				}
////				delete [] basepair;
				cout << "Finished\n";
				
				if(coord_slist.empty()){
					cout << "\nFinished all coordinates, COMPLETE\n";
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
////					map<int, int>::iterator cIter;
////					cIter = chrom_stop.find(chromInt);
////					cEnd = cIter->second;
					
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
						cout << "\nLoading data for " << chChr << "\nInitializing array...";
//						basepair = new int[(cEnd+1)];
////						basepair = new double[(cEnd+1)];
////						for(int p = 0;p <= cEnd;p++){
						for(int p = startOffset;p <= countBases;p++){
////
							basepair[p] = 0;
						}
						cout << "Finished\n";
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
//				basepair[countBases] = atof(line.c_str());
////				if(countBases > cEnd){
				if(countBases > cSize){
					cout << "Print some error, adding more data than basepairs in chrom\n";
				}
////		}
		}
	}

	if(getChrmData == 1){
		cout << "Finished, loaded " << countBases << "\nGetting data for coordinates and printing output...";
		for(i = coord_slist.begin();i!=coord_slist.end();i++){
			if(i->chrom == chromInt){
				sizeProfile = i->end - i->start;
//				profile = new double[(sizeProfile + 1)];
				profile = new int[(sizeProfile + 1)];
				int pIndex = 0;
				for(int s = i->start; s <= i->end; s++){
					profile[pIndex] = basepair[s];
					pIndex++;
				}

				if(outputData( outputFile,printFlag,i->ident,i->chrom,i->start,i->end,i->strand,sizeProfile,profile) != 0){
					cout << "Error printing output to file, exiting" << endl;
//					exit(1);
					return 1;
				}

				delete [] profile;
				sizeProfile = NULL;
				printFlag = 1;
				coord_slist.erase(i);
			}
		}
		delete [] basepair;
	}
	cout << "Finished\nFinished all coordinates, COMPLETE\n";
	seqfile.close();
	return 0;
}

int analysis::outputData(const char * outputFile,int pFlag,const char * pId,unsigned short int pChrom,unsigned long int pStart, unsigned long int pStop,const char * pStrand,int printStop,int pProfile[]){
//int analysis::outputData(const char * outputFile,int pFlag,const char * pId,unsigned short int pChrom,unsigned long int pStart, unsigned long int pStop,const char * pStrand,int printStop,double pProfile[],const char * statoutputFile){
	FILE * fh;
//	FILE * fhstats;
	if(pFlag == 0){
		fh = fopen(outputFile,"w"); 
//		fhstats = fopen(statoutputFile,"w");
		fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND\t");
//		fprintf(fhstats,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND\tAVG\tSTDEV\n");
		for(int p = 1;p <= (printStop+1);p++){
			fprintf(fh,"Position%i\t",p);
		}
		fprintf(fh,"\n");
	}else if (pFlag == 1){
		fh = fopen(outputFile,"a");
//		fhstats = fopen(statoutputFile,"a");
	}

	const char * chromName = getKey(pChrom);
	fprintf(fh,"%s\t%s\t%i\t%i\t%s\t",pId,chromName,pStart,pStop,pStrand);
//	fprintf(fhstats,"%s\t%s\t%i\t%i\t%s\t",pId,chromName,pStart,pStop,pStrand);
	char plus[] = "+";
	char minus[] = "-";
//	double sumScores = 0;
//	double sumScoresX2 = 0;
	if(strcmp(pStrand,plus) == 0){
		for(int posP = 0; posP <= printStop;posP++){
			fprintf(fh,"%i\t",pProfile[posP]);
//			sumScores += pProfile[posP];
//			sumScoresX2 += pow((pProfile[posP]),2);
		}
	}else if (strcmp(pStrand,minus) == 0){
		for(int posP = printStop; posP >= 0;posP--){
			fprintf(fh,"%i\t",pProfile[posP]);
//			sumScores += pProfile[posP];
//			sumScoresX2 += pow((pProfile[posP]),2);
		}
	}
//	double avgScore = sumScores/printStop;
//	double stdevScore = sqrt((sumScoresX2 -((sumScores)*(sumScores)/printStop))/(printStop-1));
//	fprintf(fhstats,"%f\t%f\n", avgScore, stdevScore);
//	fprintf(fh,"%f\n", avgScore);
	fprintf(fh,"\n");
	fclose (fh);
//	fclose (fhstats);
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
	slist<coord>::iterator back =  coord_slist.previous(coord_slist.end());
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
  		
  		if(countTabs!=4){
  			return 2;
  		}
	 }

	while(!feof(fh)){
		fscanf(fh,"%s%s%lu%lu%s",id,cChrom,&iStart,&iEnd,strand);
		
		//cout << "Got coord " << id << " chrom " << cChrom << " start " << iStart << " stop " << iEnd << " strand " << strand << " line count is " << lineCount << endl;
		
		unsigned short int chromInt = getHashValue(cChrom);
		coord c(id,chromInt,iStart,iEnd,strand);
		lineCount++;
		back = coord_slist.insert_after(back,c);
	}
	fclose(fh);
	cout << lineCount << " or " << coord_slist.size() << " coordinates imported\n";
	return 0;
}
