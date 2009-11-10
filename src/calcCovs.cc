#include <iostream>
#include <fstream>
#include <stdio.h>
#include "calcCovs.h"
#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "chisquaredistr.h"
#include "binomialdistr.h"


using namespace std;

calcCovs::calcCovs(){
	chromCounter = 0;
}

calcCovs::~calcCovs(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int calcCovs::processSignals(int zWinSize, int zOffsetSize, int cWinSize, int cOffsetSize, string alignDir, string twoBitFile, const char * inputFile,string outfile,string flist){
			
	FILE * tempTB;
	const char * tInfo = "tempInfo.txt"; 
	const char * tChrSize = "tempChromSize.txt";
	tempTB = fopen(tInfo,"w");
	fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile.c_str(),tChrSize);
	fclose (tempTB);
	
	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile.c_str() << endl;
	int s = system("R CMD BATCH tempInfo.txt /dev/null");
	if(s == 0){
		remove(tInfo);
	}else{
		cout << "twoBitInfo failed\n";
		exit(1);
	}
	
	tempTB = fopen(tChrSize,"r");
	char cChrom[128];
	unsigned long int cStart;
	while(!feof(tempTB)){
		int ret = fscanf(tempTB,"%s%lu",cChrom,&cStart);
		unsigned short int chromInt = getHashValue(cChrom);
		chr_size[chromInt] = cStart;
	}
	fclose(tempTB);
	remove(tChrSize);
	
	slist<aRead>::iterator i = signal_slist.begin();
	unsigned short int currchr = 999;
	unsigned short int * basepair = NULL;
	unsigned short int * ibasepair = NULL;
	unsigned short int * gcContent = NULL;	
	unsigned short int * alignability = NULL;
	int readInput = 0;
	
	while(!signal_slist.empty()){
		i = signal_slist.begin();
		currchr = i->chrom;
		const char * chromReport = getKey(currchr);
		cout << "\nProcessing " << chromReport << endl;
		basepair = new unsigned short int[chr_size[currchr]+1];
		basepair[chr_size[currchr]+1] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(i != signal_slist.end()){			
			if(i->chrom==currchr){
				if(i->start < chr_size[currchr]){
					basepair[i->start]++;
					signal_slist.erase(i++);
				}else{
					cout << "WARNING: Read maps off end of chromosome: " << i->start << endl;
					basepair[chr_size[currchr]]++;
					signal_slist.erase(i++);
				}
			}else{
				i++;
			}
		}

		alignability = new unsigned short int[chr_size[currchr] + 1];
		alignability[chr_size[currchr]+1] = 0;
		string alignFileS = alignDir + chromReport + ".wig";
		char * alignFile = new char[alignFileS.size() + 1];
		strcpy(alignFile, alignFileS.c_str());
		cout << "\tGetting alignability info from:\n\t\t" << alignFile << endl;
		tempTB = fopen(alignFile,"r");
		unsigned short int aScore;
		unsigned long int pos = 1;
		while(!feof(tempTB)){
			int ret = fscanf(tempTB,"%hu",&aScore);
			alignability[pos] = aScore;
			pos++;
		}
		fclose(tempTB);
		
		cout << "\tGetting sequence from .2bit file:\n\t\t" << twoBitFile.c_str() << endl;
		gcContent = new unsigned short int[chr_size[currchr] + 1];
		gcContent[chr_size[currchr]+1] = 0;
		const char * tInfo = "tempInfo.txt"; 
		const char * tSeq = "tempSeq.txt";
		tempTB = fopen(tInfo,"w");
		fprintf(tempTB,"library(zinba);\ntwobittofa(chrm=\"%s\",start=1,end=%lu,twoBitFile=\"%s\",gcSeq=\"%s\");\n",chromReport,chr_size[currchr],twoBitFile.c_str(),tSeq);
		fclose (tempTB);

		int s = system("R CMD BATCH tempInfo.txt /dev/null");
		if(s == 0){
			remove(tInfo);
		}else{
			cout << "twoBitToFa failed\n";
			exit(1);
		}
		ifstream seqfile(tSeq);
		string line;
		pos = 1;
		unsigned long int nStart = 0;
		unsigned long int nStop = 0;
		unsigned short int prevEnt = 0;
		while(getline(seqfile,line)){
			if(line[0] != '>'){
				for(int s = 0; s < line.length();s++){
					unsigned short int ent = 0;
					if(line[s] == 'G' || line[s] == 'C'){
						ent = 1;
					}else if (line[s] == 'N'){
						ent = 2;
					}
					gcContent[pos] = ent;
					pos++;
				}
			}
		}
		seqfile.close();
		remove(tSeq);

		cout << "\tGetting counts for " << cWinSize << "bp windows.........." << endl;
		int numOffsets = 1;
		if(cOffsetSize > 0){
			numOffsets = int(cWinSize/cOffsetSize);
		}
		for(int o = 0; o < numOffsets; o++){
			unsigned long int cWinStart = (cOffsetSize * o) + 1;
			unsigned long int cWinStop = cWinStart + cWinSize - 1;			
			while(cWinStop <= chr_size[currchr]){
				double cnvCount = 0;
				int alignCount = 0;
				double nCount = 0;
				unsigned long int gapStart = 0;
				for(int b = cWinStart; b <= cWinStop; b++){
					cnvCount += basepair[b];
					alignCount += alignability[b];
					if(gcContent[b] == 2){
						if(nCount == 0)
							gapStart = b;
						nCount++;
					}
				}
				double cScore = 0;
				if (nCount > (double)(0.1*zWinSize)){
					if((gapStart - cWinStart) > cOffsetSize){
						cnvCount = 0;
						alignCount = 0;
						for(int b = cWinStart; b < gapStart; b++){
							cnvCount += basepair[b];
							alignCount += alignability[b];
						}
						if(alignCount > 0)
							cScore = cnvCount/alignCount;
						cnvWins gcnv(currchr,cWinStart,(gapStart-1),cScore,0,0,0);	
						cnv_wins.push_back(gcnv);
					}
					while(gcContent[gapStart] == 2 && gapStart < cWinStop)
						gapStart++;
					if((cWinStop - gapStart) > cOffsetSize){
						cnvCount = 0;
						alignCount = 0;
						for(int b = gapStart; b <= cWinStop; b++){
							cnvCount += basepair[b];
							alignCount += alignability[b];
						}
						if(alignCount > 0)
							cScore = cnvCount/alignCount;
						cnvWins gcnv2(currchr,gapStart,cWinStop,cScore,0,0,0);	
						cnv_wins.push_back(gcnv2);
					}					
				}else{
					if(alignCount > 0)
						cScore = cnvCount/alignCount;
					cnvWins cnv(currchr,cWinStart,cWinStop,cScore,0,0,0);
					cnv_wins.push_back(cnv);
				}
				cWinStart += cWinSize;
				cWinStop += cWinSize;
			}
		}
			
		cout << "\t\tRefining boundaries...." << endl;
		cnv_wins.sort();
		int numCnvWins = cnv_wins.size();
		double slideWinSize = numOffsets * 2;
		double globalSum = 0;
		double globalSumX2 = 0;
		list<double> localSum;
		list<double> localSumX2;

		list<cnvWins>::iterator c = cnv_wins.begin();
		while(c != cnv_wins.end()){
			globalSum += c->cnvScore;
			globalSumX2 += pow((c->cnvScore),2);
			if(localSum.size() < (slideWinSize - 1)){
				localSum.push_back(c->cnvScore);
				localSumX2.push_back(pow((c->cnvScore),2));
			}else{
				localSum.push_back(c->cnvScore);
				localSumX2.push_back(pow((c->cnvScore),2));
				if(localSum.size() == (slideWinSize + 1)){
					localSum.pop_front();
					localSumX2.pop_front();
				}

				list<double>::iterator x = localSum.begin();
				list<double>::iterator x2 = localSumX2.begin();
				double sum = 0;
				double sumX2 = 0;
				while(x != localSum.end()){
					sum += *x;
					sumX2 += *x2;
					x++;x2++;
				}
				c->varWin = (sumX2 -((sum)*(sum)/slideWinSize))/(slideWinSize-1);
			}
			c++;
		}
		localSum.clear();
		localSumX2.clear();		
		double globalVar = (globalSumX2 -((globalSum)*(globalSum)/numCnvWins))/(numCnvWins-1);
		cout << "\t\t\tGlobal variance is " << globalVar << endl;
		list<cnvWins> sigBoundary; 
		c = cnv_wins.begin();
		while(c != cnv_wins.end()){
			double degfree = (slideWinSize - 1);
			c->chiSq = (degfree*pow(c->varWin,2))/pow(globalVar,2);
			c->pval = chisquarecdistribution(degfree,c->chiSq);
			if(c->pval <= 0.001){
				sigBoundary.push_back(*c);
			}
			c++;
		}
		sigBoundary.sort();
		list<cnvWins>::iterator sb = sigBoundary.begin();
		cnvWins lastCNV = *sb;
		sb++;
		while(sb != sigBoundary.end()){
			if(sb->start >= lastCNV.start && sb->start <= lastCNV.stop){
				long unsigned int lcStop = sb->stop;
				sigBoundary.erase(sb--);
				sb->stop = lcStop;
			}
			lastCNV = *sb;
			sb++;
		}
		cout << "\t\t\tRefining " << sigBoundary.size() << " boundaries" << endl;		
		sb = sigBoundary.begin();
		list<long unsigned int> transPts;
		long unsigned int refOffset = numOffsets * cOffsetSize;
		while(sb != sigBoundary.end()){
			long unsigned int leftWinStart = 1;
			if(sb->start > refOffset)
				leftWinStart = sb->start - refOffset;
			long unsigned int leftWinStop = leftWinStart + cWinSize;
			double lowPval = 1;
			double maxRatio = 0;
			long unsigned int transPoint = 0;
			int hitGap = 0;
			while((leftWinStop+1) <= sb->stop && (leftWinStop+cWinSize+1) <= chr_size[currchr]){
				int leftWinCount = 0;
				int rightWinCount = 0;
				unsigned long int gapStart = 0;
				for(int s = leftWinStart; s <= leftWinStop; s++){
					leftWinCount += basepair[s];
					rightWinCount += basepair[(s+cWinSize+1)];
					if(gcContent[s] == 2){
						if(hitGap == 0)
							gapStart = s;
						hitGap = 1;
					}
					if(gcContent[s+cWinSize+1] == 2){
						if(hitGap == 0)
							gapStart = s+cWinSize+1;
						hitGap = 1;
					}
				}
				if(hitGap == 1){
					transPoint = gapStart;
					leftWinStop = sb->stop;
				}else if (hitGap == 0){
					unsigned int maxCount = leftWinCount;
					if(rightWinCount > leftWinCount)
						maxCount = rightWinCount;
					unsigned int sumCount = leftWinCount + rightWinCount;
					double ratioDiff = (double) maxCount/sumCount;
					double bPval = binomialcdistribution(maxCount,sumCount,0.5);
					if(bPval <= lowPval){
						if(ratioDiff > maxRatio){
							lowPval = bPval;
							maxRatio = ratioDiff;
							transPoint = leftWinStop + 1;
						}
					}
					leftWinStart += zOffsetSize;
					leftWinStop += zOffsetSize;
				}
			}
			transPts.push_back(transPoint);
			if(hitGap == 1){
				while(gcContent[transPoint] == 2 && transPoint <= chr_size[currchr])
					transPoint++;
				transPts.push_back((transPoint-1));
			}
			sigBoundary.erase(sb++);
		}
		transPts.sort();
		list<long unsigned int>::iterator tp = transPts.begin();
		c = cnv_wins.begin();
		while(c != cnv_wins.end() && tp != transPts.end()){
			if(c->start >= *tp)
				tp++;
			if(c->start < *tp && c->stop > *tp && tp != transPts.end())
				cnv_wins.erase(c++);
			else
				c++;
		}
		tp = transPts.begin();
		list<long unsigned int>::iterator nexttp = transPts.begin();
		nexttp++;
		while(tp != transPts.end()){
			double cnvCount = 0;
			int alignCount = 0;
			unsigned long int tpStart = (*tp-cWinSize-1);
			unsigned long int tpStop = (*tp-1);
			for(int b = tpStart; b <= tpStop; b++){
				cnvCount += basepair[b];
				alignCount += alignability[b];
			}
			double cScore = cnvCount/alignCount;
			cnvWins cnv(currchr,tpStart,tpStop,cScore,0,0,0);
			cnv_wins.push_back(cnv);
			if(gcContent[*tp] == 2 && gcContent[*nexttp] == 2){
				tpStart = *nexttp+1;
				tpStop = (*nexttp+cWinSize);
			}else{
				tpStart = *tp;
				tpStop = (*tp+cWinSize);				
			}
			if(tpStop <= chr_size[currchr]){
				cnvCount = 0;
				alignCount = 0;
				for(int b = tpStart; b < tpStop; b++){
					cnvCount += basepair[b];
					alignCount += alignability[b];
				}
				cScore = cnvCount/alignCount;
				cnvWins cnv_two(currchr,tpStart,tpStop,cScore,0,0,0);
				cnv_wins.push_back(cnv_two);
				if(gcContent[*tp] == 2 && gcContent[*nexttp] == 2){
					tp++;
					nexttp++;
				}
			}else{
				tp++;
			}
			tp++;
			nexttp++;
		}
		cnv_wins.sort();

		const char* inVal = "none";
		if(strcmp(inputFile,inVal)!=0){
			if(readInput == 0){
				cout << "\tLoading reads from input file " << inputFile << "........." << endl;
				int rVal = importRawSignal(inputFile,1);
				readInput = 1;
			}
			ibasepair = new unsigned short int[chr_size[currchr]+1];
			ibasepair[chr_size[currchr]+1] = 0;
			slist<aRead>::iterator in = input_slist.begin();
			while(in != input_slist.end()){			
				if(in->chrom==currchr){					
					if(in->start < chr_size[currchr]){
						basepair[in->start]++;
						input_slist.erase(in++);
					}else{
						cout << "WARNING: Read maps off end of chromosome: " << in->start << endl;
						basepair[chr_size[currchr]]++;
						input_slist.erase(in++);
					}
				}else{
					in++;
				}
			}
		}
		
		cout << "\tGetting counts for zinba windows.........." << endl;
		tempTB = fopen(flist.c_str(),"a");
		numOffsets = 1;
		if(zOffsetSize > 0){
			numOffsets = int(zWinSize/zOffsetSize);
		}
		string outfileDATA;
		slist<dataWins>::iterator z;
		for(int o = 0; o < numOffsets; o++){
			z = peak_wins.previous(peak_wins.end());
			list<cnvWins>::iterator cnvBegin = cnv_wins.begin();
			list<cnvWins>::iterator cnvEnd = cnv_wins.begin();
			cout << "\t\tOffset " << (zOffsetSize * o) << "bp......" << endl;
			stringstream offset;
			offset << (zOffsetSize * o);
			outfileDATA = outfile + "_" + chromReport + "_offset" + offset.str() + "bp.txt";
			if(o == (numOffsets-1))
				fprintf(tempTB,"%s\n",outfileDATA.c_str());
			else
				fprintf(tempTB,"%s;",outfileDATA.c_str());
			unsigned long int zWinStart = (zOffsetSize * o) + 1;
			unsigned long int zWinStop = zWinStart + zWinSize - 1;			
			while(zWinStop <= chr_size[currchr]){
				int peakCount = 0;
				double alignCount = 0;
				double gcCount = 0;
				double nCount = 0;
				double inCount = 0;
				for(int b = zWinStart; b <= zWinStop; b++){
					peakCount += basepair[b];
					alignCount += alignability[b];
					if(gcContent[b] == 1){
						gcCount++;
					}else if(gcContent[b] == 2){
						nCount++;
					}
					if(strcmp(inputFile,inVal)!=0)
						inCount += ibasepair[b];
				}
				if((nCount/zWinSize) < 0.1){
					double cnvSum = 0;
					int cnvCount = 0;
					while(cnvBegin->stop < zWinStart && cnvBegin != cnv_wins.end())
						cnvBegin++;
					cnvEnd = cnvBegin;
					while( (cnvEnd->start <= zWinStart && cnvEnd->stop >= zWinStart) || (cnvEnd->start <= zWinStop && cnvEnd->stop >= zWinStop)){
						cnvSum += cnvEnd->cnvScore;
						cnvCount++;
						cnvEnd++;
					}
					
					double cnvLogScore = 0;
					if(cnvSum > 0)
						cnvLogScore = log((cnvSum/cnvCount)*(zWinSize*2));
					if(inCount > 0)
						inCount = log(inCount);
					double gcPerc = gcCount/zWinSize;
					double aPerc = alignCount/(zWinSize*2);
					dataWins zwin(currchr,zWinStart,zWinStop,peakCount,inCount,gcPerc,aPerc,cnvLogScore);
					z = peak_wins.insert_after(z,zwin);
				}
				zWinStart += zWinSize;
				zWinStop += zWinSize;
			}
			
			if(outputData(outfileDATA.c_str(),currchr) != 0){
				cout << "Error printing output to file, exiting" << endl;
				exit(1);
			}
			
		}
		fclose(tempTB);
		cnv_wins.clear();
		
cout << "cnv window size is " << cnv_wins.size() << endl;
		
		sigBoundary.clear();
		transPts.clear();
		delete [] basepair;
		delete [] ibasepair;
		delete [] gcContent;
		delete [] alignability;
	}
	return 0;
}

int calcCovs::outputData(const char * outputFile, unsigned short int currChr){
	FILE * fh;
	fh = fopen(outputFile,"w");
	fprintf(fh,"chromosome\tstart\tstop\texp_count\tinput_count\tgcPerc\talign_perc\texp_cnvwin_log\n");
	const char * chrom = getKey(currChr);
	slist<dataWins>::iterator c = peak_wins.begin();
	while(c != peak_wins.end()){
		fprintf(fh,"%s\t%lu\t%lu\t%i\t%f\t%f\t%f\t%f\n",chrom,c->start,c->stop,c->eCount,c->iCount,c->gcPerc,c->alignPerc,c->cnvScore);
		peak_wins.erase(c++);
	}
	fclose (fh);
	return 0;
}

unsigned short int calcCovs::getHashValue(char *currChrom){
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

const char * calcCovs::getKey(unsigned short int chrom){
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

int calcCovs::importRawSignal(const char * signalFile,int dataType){
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){return 1;}
	
	unsigned long int lineCount = 0;
	char cChrom[128];
	unsigned long int iStart;
	
	slist<aRead>::iterator back =  signal_slist.previous(signal_slist.end());
	slist<aRead>::iterator iback = input_slist.previous(input_slist.end());
	
	while(!feof(fh)){
		int ret = fscanf(fh,"%s%lu",cChrom,&iStart);
		unsigned short int chromInt = getHashValue(cChrom);
		aRead sig(chromInt,iStart);
		lineCount++;
		if(dataType == 0)
			back = signal_slist.insert_after(back,sig);
		else if(dataType == 1)
			iback = input_slist.insert_after(iback,sig);			
	}
	fclose(fh);
	if(dataType == 0){
		cout << "\tImported " << lineCount << " reads" << endl;
		cout << "\tSorting reads ...";
		signal_slist.sort();
	}else if(dataType == 1){
		cout << "\t\tImported " << lineCount << " reads" << endl;
		cout << "\t\tSorting reads ...";
		input_slist.sort();
	}
	cout << "COMPLETE" << endl;
	return 0;
}