#include <iostream>
#include <fstream>
#include <stdio.h>
#include "zCollapseWins.h"
#include <list>
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>

using namespace std;

zCollapseWins::zCollapseWins(){
	chromCounter = 0;
}

zCollapseWins::~zCollapseWins(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int zCollapseWins::outputData(const char * outputFile){
	FILE * fh;
	fh = fopen(outputFile,"w"); 
	fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tQFLAG\tSTRAND\n");

	char strand[] = "+";
	char winID[255];
	char start[10];
	char stop[10];
	
int countCoords = 0;
	
	list<coord>::iterator outIt =  coordOut_list.begin();
	while(outIt != coordOut_list.end()){

countCoords++;
//cout << countCoords << " " << outIt->chrom << endl;

		const char * chromName = getKey(outIt->chrom);
		const char * neg = "NONE";
		if(strcmp(chromName,neg) != 0){
			sprintf(start,"%d", outIt->start);
			sprintf(stop,"%d", outIt->end);
			strcpy(winID,chromName);strcat(winID,":");strcat(winID,start);strcat(winID,"-");strcat(winID,stop);
			fprintf(fh,"%s\t%s\t%lu\t%lu\t1\t+\n",winID,chromName,outIt->start,outIt->end);
			coordOut_list.erase(outIt++);
		}
	}
	fclose (fh);
	return 0;
}

int zCollapseWins::importCoords(const char * winlist,const char * coordfile,unsigned short int winSize,unsigned short int wingap,double threshold,string method,int wformat){

	FILE * wlist;
	wlist = fopen(winlist,"r");
	if(wlist == NULL){return 1;}
	char cChrom[128];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	long double sigVal;
	string pscl = "pscl";
	string mixture = "mixture";
	char sigFile [256];
	int readResult = 0;
	char firstline [256];
	int rwline;
	///////////////////////////
	unsigned long int winSizeThresh = winSize * 10;
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
				if(strcmp(method.c_str(),pscl.c_str()) == 0){
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%Lf%*f",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal <= threshold && rwline > 0){
							readResult = 1;
						}
					}else if(wformat == 1){
						rwline = fscanf(fh,"%s%lu%lu%*d%*f%*f%*f%*f%hu%Lf%*f",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal <= threshold && rwline > 0){
							readResult = 1;
						}
					}
				}else if(strcmp(method.c_str(),mixture.c_str()) == 0){
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%Lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal >= threshold && rwline > 0){
							readResult = 1;
						}
					}else if(wformat == 1){
						int exp;double inp,gc,ap,ecl;
						rwline = fscanf(fh,"%s%lu%lu%*d%*f%*f%*f%*f%hu%Lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal);
						if(sigVal >= threshold && rwline > 0){
							readResult = 1;
						}
					}
				}
				if(readResult == 1){
					string chromIn(cChrom);
					unsigned short int chromInt = getHashValue(chromIn.c_str());
					coord c(chromInt,iStart,iEnd,qFlag);
					coord_slist.push_back(c);
				}
			}
			fclose(fh);
		}
	}
	fclose(wlist);

	for(int ci = 0; ci < chromCounter;ci++){
		const char * chromName = getKey(ci);
		cout << ci << " is " << chromName << endl;
	}
	
	cout << "\nImported " << coord_slist.size() << " coordinates" << endl;
	coord_slist.sort();
	
	list<coord>::iterator back =  coord_slist.begin();
	list<coord> tempcoord_list;
	tempcoord_list.push_front(*back);
	list<coord>::iterator tempIt = tempcoord_list.begin();
	coord_slist.erase(back++);
	int flagFinish = 0;

	while(flagFinish == 0){
		if(coord_slist.empty())
			flagFinish = 1;
		if(flagFinish == 0 && back->chrom == tempIt->chrom && back->start <= tempIt->end+1){
			tempcoord_list.push_back(*back);
			tempIt = tempcoord_list.rbegin();
			coord_slist.erase(back++);
		}else if( (back->chrom == tempIt->chrom && back->start > tempIt->end+1) || back->chrom != tempIt->chrom || flagFinish == 1){
			if(tempcoord_list.size() == 1 && tempIt->qFlag == 1){

if(tempIt->chrom > 3 || tempIt->chrom < 0){
	cout << "Chrom1 is " << tempIt->chrom << endl;
}

				coordOut_list.push_back(*tempIt);
			}else if(tempcoord_list.size() > 1){
				// TRIM PRECEEDING WINDOWS BELOW QUARTILE
				tempIt = tempcoord_list.begin();
				while(tempIt != tempcoord_list.end()){
					if(tempIt->qFlag == 0)
						tempcoord_list.erase(tempIt++);
					else if(tempIt->qFlag == 1)
						tempIt = tempcoord_list.end();
					else
						tempIt++;
				}
				// TRIM TRAILING WINDOWS BELOW QUARTILE
				tempIt = tempcoord_list.rbegin();
				while(tempIt != tempcoord_list.begin()){
					if(tempIt->qFlag == 0)
						tempcoord_list.erase(tempIt--);
					else if(tempIt->qFlag == 1)
						tempIt = tempcoord_list.begin();
					else
						tempIt--;
				}
				if(tempcoord_list.size() <= wingap){
					tempIt = tempcoord_list.begin();
					unsigned long int start = tempIt->start;
					unsigned long int stop = tempIt->end;
					tempIt++;
					while(tempIt != tempcoord_list.end()){
						stop = tempIt->end;
						tempIt++;
					}
					
if(tempIt->chrom > 3 || tempIt->chrom < 0){
	cout << "Chrom2 is " << tempIt->chrom << endl;
}
					
					coord c(tempIt->chrom,start,stop,1);
					coordOut_list.push_back(c);
				}else{
					//BREAK UP WINDOWS WITH CONTIGUOUS WINDOWS BELOW QUARTILE IN MIDDLE
					tempIt = tempcoord_list.begin();
					int countLowQuart = 0;
					int countHighQuart = 1;
					coord lastCoord = *tempIt;			
					tempIt++;
					while(tempIt != tempcoord_list.end()){
						if(lastCoord.qFlag == 0 && tempIt->qFlag == 0){
							countLowQuart++;
						}else if(lastCoord.qFlag == 1 && tempIt->qFlag == 0){
							countLowQuart = 1;
						}else if(lastCoord.qFlag == 1 && tempIt->qFlag == 1){
							countHighQuart++;
						}else if(lastCoord.qFlag == 0 && tempIt->qFlag == 1){
							if(countLowQuart >= wingap){
								tempIt--;
								while(countLowQuart > 0){
									tempcoord_list.erase(tempIt--);
									countLowQuart--;
								}
								unsigned long int start = tempIt->start;
								unsigned long int stop = tempIt->end;
								countHighQuart--;
								tempcoord_list.erase(tempIt--);
								while(countHighQuart > 0){
									start = tempIt->start;
									tempcoord_list.erase(tempIt--);
									countHighQuart--;
								}
								if((stop-start) <= winSizeThresh){
									
if(tempIt->chrom > 3 || tempIt->chrom < 0){
	cout << "Chrom3 is " << tempIt->chrom << endl;
}
									
									coord c(tempIt->chrom,start,stop,1);
									coordOut_list.push_back(c);
								}else{
									//COULD FURTHER DROP OUT LOW QUART REGIONS TO POTENTIALLY SAVE THESE REGIONS
									cout << "Excluding " << back->chrom << ":" << start << "-" << stop << " SIZE=" << (stop-start) << endl;
								}
								tempIt = tempcoord_list.begin();
								countLowQuart = 0;
								countHighQuart = 1;
							}else{
								countLowQuart = 0;
								countHighQuart = countHighQuart + countLowQuart;
							}
						}
						lastCoord = *tempIt;
						tempIt++;
					}
					if(countLowQuart >= wingap){
						while(countLowQuart > 0){
							tempcoord_list.erase(tempIt--);
							countLowQuart--;
						}
					}
					tempIt = tempcoord_list.begin();
					unsigned long int start = tempIt->start;
					unsigned long int stop = tempIt->end;
					tempIt++;
					while(tempIt != tempcoord_list.end()){
						stop = tempIt->end;
						tempIt++;
					}
					if((stop-start) <= winSizeThresh){
if(tempIt->chrom > 3 || tempIt->chrom < 0){
	cout << "Chrom4 is " << tempIt->chrom << endl;
}
						coord c(tempIt->chrom,start,stop,1);
						coordOut_list.push_back(c);
					}else{
						//COULD FURTHER DROP OUT LOW QUART REGIONS TO POTENTIALLY SAVE THESE REGIONS
						cout << "Excluding " << back->chrom << ":" << start << "-" << stop << " SIZE=" << (stop-start) << endl;
					}
				}
			}
			tempcoord_list.clear();
			tempcoord_list.push_front(*back);
			tempIt = tempcoord_list.begin();
			coord_slist.erase(back++);
		}
	}
	coordOut_list.sort();
	cout << "\nCollapsed to " << coordOut_list.size() << " non-overlapping regions" << endl;
	outputData(coordfile);
}

unsigned short int zCollapseWins::getHashValue(const char *currChrom){
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

const char * zCollapseWins::getKey(unsigned short int chrom){
	map<int, const char*>::iterator i;
	i = intsToChrom.find(chrom);
	if(i == intsToChrom.end()){
		const char * neg = "NONE";
		return neg;
//		cout << chrom << endl;
//		cout << "REALLY REALLY BAD ERROR!" << endl;
//		exit(1);
	}else{
		return i->second;	
	}
}

