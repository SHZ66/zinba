#include <iostream>
#include <fstream>
#include <stdio.h>
#include "calcCovs.h"
//#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <time.h>

using namespace std;

calcCovs::calcCovs(){
	chromCounter = 0;
	tbSizeFlag = 0;
}

calcCovs::~calcCovs(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int calcCovs::processSignals(int zWinSize, int zOffsetSize, int cWinSize, int cOffsetSize, string alignDir,const char * twoBitFile,const char * inputFile,string outfile,const char * flist,int extension,const char * filetype){

	FILE * tempTB;
	time_t rtime;
	struct tm *timeinfo;
	char tInfo[128];// = "tempInfo.txt";
	char sysCall[256];
	
	unsigned short int currchr = 999;
	unsigned short int * basepair = NULL;
	unsigned short int * ibasepair = NULL;
	unsigned char * gcContent = NULL;	
	unsigned char * alignability = NULL;
	int i;
	int printflag = 0;

	int readInput = 0;
	const char* noneVal = "none";
	
	while(!signal_slist.empty()){
		i = 0;
		currchr = signal_slist[0].chrom;
		const char * chromReport = getKey(currchr);
		cout << "\nProcessing " << chromReport << endl;
		basepair = new unsigned short int[chr_size[currchr]+1];
		for(int ch = chr_size[currchr]; ch--;)
			basepair[ch] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(signal_slist[i].chrom==currchr && i < (int) signal_slist.size()){		
			basepair[signal_slist[i].pos]++;
			i++;
		}
		signal_slist.erase(signal_slist.begin(),signal_slist.begin()+i);

		alignability = new unsigned char[chr_size[currchr] + 1];
		for(int ch = chr_size[currchr]; ch--;)
			alignability[ch] = (unsigned char) 0;

		string alignFileS = alignDir + chromReport + ".wig";
		char * alignFile = new char[alignFileS.size() + 1];
		strcpy(alignFile, alignFileS.c_str());
		cout << "\tGetting alignability info from:\n\t\t" << alignFile << endl;
		tempTB = fopen(alignFile,"r");
		if(tempTB == NULL){
			cout << "Unable to open alignability file " << alignFile << endl;
			return 1;
		}
		unsigned short int aScore;
		unsigned long int pos = 1;
		while(!feof(tempTB)){
			int ret = fscanf(tempTB,"%hu",&aScore);
			alignability[pos] = (unsigned char) aScore;
			pos++;
		}
		fclose(tempTB);
		
		cout << "\tGetting sequence from .2bit file:\n\t\t" << twoBitFile << endl;
		gcContent = new unsigned char[chr_size[currchr] + 1];
		for(int ch = chr_size[currchr]; ch--;)
			gcContent[ch] = (unsigned char) 0;

		char tSeq[128];
		time(&rtime);
		timeinfo=localtime(&rtime);
		strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
		strftime(tSeq,128,"tempSeq_%H_%M_%S.txt",timeinfo);

		tempTB = fopen(tInfo,"w");
		fprintf(tempTB,"library(zinba);\ntwobittofa(chrm=\"%s\",start=1,end=%lu,twoBitFile=\"%s\",gcSeq=\"%s\");\n",chromReport,chr_size[currchr],twoBitFile,tSeq);
		fclose (tempTB);
		sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		int s = 1;
		int twobitCount = 0;
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitToFa again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "twoBitToFa failed, exiting" << endl;
				return 1;
			}
		}
		remove(tInfo);
		
		ifstream seqfile(tSeq);
		string line;
		pos = 1;
		unsigned long int nStart = 0;
		unsigned long int nStop = 0;
		unsigned short int prevEnt = 0;
		while(getline(seqfile,line)){
			if(line[0] != '>'){
				for(int s = 0; s < line.length();s++){
					if(line[s] == 'G' || line[s] == 'C'){
						gcContent[pos] = (unsigned char) 1;
					}else if (line[s] == 'N'){
						gcContent[pos] = (unsigned char) 2;
					}else{
						gcContent[pos] = (unsigned char) 0;
					}
					pos++;
				}
			}
		}
		seqfile.close();
		remove(tSeq);

		cout << "\tGetting counts for " << cWinSize << "bp windows.........." << endl;
		int numOffsets = 1;
		if(cOffsetSize > 0){
			numOffsets = (int) cWinSize/cOffsetSize;
		}
		for(int o = 0; o < numOffsets; o++){
			unsigned long int cWinStart = (cOffsetSize * o) + 1;
			unsigned long int cWinStop = cWinStart + cWinSize - 1;
			if(cWinStop > chr_size[currchr])
				cWinStop = chr_size[currchr];
			while(cWinStop <= chr_size[currchr]){
				double cnvCount = 0.0;
				int alignCount = 0;
				double nCount = 0.0;
				for(int b = cWinStart; b <= cWinStop; b++){
					if((int) gcContent[b] == 2){
						nCount++;						
					}else{
						cnvCount += basepair[b];
						alignCount += (int) alignability[b];
					}
				}
				double cScore = 0.0;
				if(alignCount > 0)
					cScore = (double) cnvCount/alignCount;
				double percGap = (double) nCount/cWinSize;
				cnvWins cnv(cWinStart,cWinStop,cScore,0,0,percGap);
				cnv_wins.push_back(cnv);
				cWinStart += cWinSize;
				cWinStop += cWinSize;
			}
		}
			
		cout << "\t\tRefining boundaries...." << endl;
		cnv_wins.sort();
		
//const char * cWinFile = "cnv_wins.txt";
//tempTB = fopen(cWinFile,"w");
//list<cnvWins>::iterator cf = cnv_wins.begin();
//while(cf != cnv_wins.end()){
//	if(cf->pGap < 0.01){
//		long unsigned int cwpos = (long unsigned int) (cf->stop+cf->start)/2;
//		fprintf(tempTB,"%lu\t%f\n",cwpos,cf->cnvScore);
//	}
//	cf++;
//}
//fclose (tempTB);
				
		int numCnvWins = 0;
		double slideWinSize = numOffsets * 2.0;
		double globalSum = 0.0;
		double globalSumX2 = 0.0;
		double gapThresh = 0.0;
		int firstVar = cnv_wins.size();
		int countWin = 0;
		list<double> localSum;
		list<double> localSumX2;

//const char * varFile = "var_wins.txt";
//tempTB = fopen(varFile,"w");
		
		list<cnvWins>::iterator c = cnv_wins.begin();
		while(c != cnv_wins.end()){
			if(c->pGap <= gapThresh){
				globalSum += c->cnvScore;
				globalSumX2 += pow((c->cnvScore),2);
				numCnvWins++;
				if(localSum.size() < (slideWinSize - 1.0)){
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
					double sum = 0.0;
					double sumX2 = 0.0;
					while(x != localSum.end()){
						sum += *x;
						sumX2 += *x2;
						x++;x2++;
					}
					c->varWin = (sumX2 -((sum * sum)/slideWinSize))/(slideWinSize-1);
//long unsigned int vpos = (unsigned long int) ((c->start+c->stop)/2)-cWinSize;
//fprintf(tempTB,"%lu\t%e\n",vpos,c->varWin);
					if(countWin < firstVar)
						firstVar = countWin;
				}
				c++;
			}else{
				if(localSum.size() > 0){
					localSum.clear();
					localSumX2.clear();
				}
				while(c->pGap != 0 && c!= cnv_wins.end())
					c++;
			}
			countWin++;
		}
		
//fclose (tempTB);

//const char * pvalFile = "pval_wins.txt";
//tempTB = fopen(pvalFile,"w");
		
		localSum.clear();
		localSumX2.clear();
		double globalVar = (globalSumX2 -((globalSum * globalSum)/numCnvWins))/(numCnvWins-1);
		double degfree = (slideWinSize - 1);
		cout << "\t\t\tGlobal variance is " << globalVar << endl;
		list<cnvWins> sigBoundary;
		c = cnv_wins.begin();
		int cnvpos = 0;
		while(cnvpos < firstVar){
			c++;
			cnvpos++;
		}
		while(c != cnv_wins.end()){
			if(c->pGap <= gapThresh){
				c->chiSq = (degfree*pow(c->varWin,2))/pow(globalVar,2);
				double cPval = dchisq(c->chiSq, degfree, 0);
//				c->pval = chisquarecdistribution(degfree,c->chiSq);
				
//long unsigned int vpos = (unsigned long int) ((c->start+c->stop)/2)-cWinSize;
//fprintf(tempTB,"%lu\t%e\n",vpos,cPval);
				
				if(cPval >= 0.000001){
					sigBoundary.push_back(*c);
				}
				c++;
			}else{
				c++;
			}
		}
		
//fclose (tempTB);
		
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
		
//while(sb != sigBoundary.end()){
//	cout << "\t\t\tSIG BOUND is " << sb->start << " " << sb->stop << endl;
//	sb++;
//}
//sb = sigBoundary.begin();
		
		list<long unsigned int> transPts;
		long unsigned int refOffset = 2 * cWinSize;
		long unsigned int leftWinStart = sb->start;
		if(sb->start > refOffset)
			leftWinStart = sb->start - refOffset;

		while(!sigBoundary.empty()){
			
			int searchLen = cWinSize;
			leftWinStart = sb->start - refOffset;

			long unsigned int leftWinStop;
			if((chr_size[currchr] - leftWinStart) > refOffset){
				leftWinStop = leftWinStart + cWinSize;
			}else{
				searchLen = (int) (chr_size[currchr] - leftWinStart)/4;
				leftWinStop = leftWinStart + searchLen;
			}
			
			double lowPval = 1;
			double maxRatio = 0;
			long unsigned int transPoint = 0;
			
			while(leftWinStop <= (sb->stop-searchLen) && (leftWinStop+searchLen+1) <= chr_size[currchr]){
				int leftWinCount = 0;
				int rightWinCount = 0;

				for(int s = leftWinStart; s <= leftWinStop; s++){
					leftWinCount += basepair[s];
					rightWinCount += basepair[(s+searchLen+1)];
				}

				double maxCount = (double) leftWinCount;
				if(rightWinCount > leftWinCount)
					maxCount = (double) rightWinCount;
				double sumCount = (double) leftWinCount + rightWinCount;
				double ratioDiff = (double) maxCount/sumCount;
				double bPval = dbinom(maxCount,sumCount,0.5,0);
//				double bPval = binomialcdistribution(maxCount,sumCount,0.5);
				if(bPval <= lowPval){
					if(ratioDiff > maxRatio){
						lowPval = bPval;
						maxRatio = ratioDiff;
						transPoint = leftWinStop + 1;
					}
				}
				leftWinStart += zWinSize;
				leftWinStop += zWinSize;
//				leftWinStart += zOffsetSize;
//				leftWinStop += zOffsetSize;
			}
			transPts.push_back(transPoint);
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
		
////////////////////////////////////////////
//
// need to make sure that new cnv windows
// don't overlap other trans points
//
////////////////////////////////////////////
		
		while(!transPts.empty()){
//cout << "\t\t\t\tTransition point at " << *tp << endl;
			double leftCnvCount = 0;
			int leftAlignCount = 0;
			double rightCnvCount = 0;
			int rightAlignCount = 0;
			
			unsigned long int tpStart = 1;
			if(*tp > (cWinSize+1))
				tpStart = (*tp-cWinSize);
			unsigned long int tpStop = (unsigned long int) (chr_size[currchr]-tpStart)/2;
			if(chr_size[currchr] > (*tp+cWinSize)){
				tpStop = *tp;
			}
			for(int b = tpStart; b <= tpStop; b++){
				leftCnvCount += basepair[b];
				leftAlignCount += (int) alignability[b];
				rightCnvCount += basepair[b+cWinSize];
				rightAlignCount += (int) alignability[b+cWinSize];
			}
			double cScore = 0;
			if(leftAlignCount > 0)
				cScore = leftCnvCount/leftAlignCount;
			cnvWins cnv(tpStart,tpStop,cScore,0,0,0);
			cnv_wins.push_back(cnv);
			cScore = 0;
			if(rightAlignCount > 0)
				cScore = rightCnvCount/rightAlignCount;
			cnvWins cnvR((tpStart+cWinSize),(tpStop+cWinSize),cScore,0,0,0);
			cnv_wins.push_back(cnvR);
			transPts.erase(tp++);
		}
		cnv_wins.sort();

		if(strcmp(inputFile,noneVal)!=0){
			if(readInput == 0){
				cout << "\tLoading reads from input file " << inputFile << "........." << endl;
				int rVal = importRawSignal(inputFile,extension,filetype,1,twoBitFile);
				if(rVal == 1){
					cout << "Unable to open file with input reads" << endl;
					return 1;
				}
				readInput = 1;
			}
			cout << "\tMapping input tags to the genome........." << endl;
			i = 0;
			ibasepair = new unsigned short int[chr_size[currchr]+1];
			for(int ch = chr_size[currchr]; ch--;)
				ibasepair[ch] = 0;

			while(input_slist[i].chrom==currchr && i < (int) input_slist.size()){
				ibasepair[input_slist[i].pos]++;
				i++;
			}
			input_slist.erase(input_slist.begin(),input_slist.begin()+i);
		}
		
		cout << "\tGetting counts for zinba windows.........." << endl;
		if(printflag==0){
			tempTB = fopen(flist,"w");
			printflag = 1;
		}else{
			tempTB = fopen(flist,"a");
		}
		numOffsets = 1;
		if(zOffsetSize > 0){
			numOffsets = (int) zWinSize/zOffsetSize;
		}
		string outfileDATA;
		slist<dataWins>::iterator z;
		for(int o = 0; o < numOffsets; o++){		
			z = peak_wins.previous(peak_wins.end());	
			list<cnvWins>::iterator cnvBegin = cnv_wins.begin();
			list<cnvWins>::iterator cnvEnd = cnv_wins.begin();
			cout << "\t\tOffset " << (zOffsetSize * o) << "bp......" << endl;
			char offset[128];
			sprintf(offset,"%d",(zOffsetSize * o));
			char winsize[128];
			sprintf(winsize,"%d",zWinSize);			
			outfileDATA = outfile + "_" + chromReport + "_win" + string(winsize) + "bp_offset" + string(offset) + "bp.txt";
			
			if(o == (numOffsets-1))
				fprintf(tempTB,"%s\n",outfileDATA.c_str());
			else
				fprintf(tempTB,"%s;",outfileDATA.c_str());
			unsigned long int zWinStart = (zOffsetSize * o) + 1;
			unsigned long int zWinStop = zWinStart + zWinSize - 1;
			cout << zWinStart<< endl;
			while(zWinStop <= chr_size[currchr]){
				int peakCount = 0;
				double alignCount = 0;
				double gcCount = 0;
				double nCount = 0;
				double inCount = 0;
				for(int b = zWinStart; b <= zWinStop; b++){
					peakCount += basepair[b];
					alignCount += (int) alignability[b];
					
					if((int) gcContent[b] == 1)
						gcCount++;
					else if((int) gcContent[b] == 2)
						nCount++;
					
					if(readInput == 1)
						inCount += ibasepair[b];
				}
				if((nCount/zWinSize) < 0.1){
					double cnvSum = 0;
					int cnvCount = 0;
					while(cnvBegin->stop < zWinStart && cnvBegin != cnv_wins.end())
						cnvBegin++;
					cnvEnd = cnvBegin;
					while( ((cnvEnd->start <= zWinStart && cnvEnd->stop >= zWinStart) || (cnvEnd->start <= zWinStop && cnvEnd->stop >= zWinStop)) && cnvEnd != cnv_wins.end() ){
						cnvSum += cnvEnd->cnvScore;
						cnvCount++;
						cnvEnd++;
					}
					
					double cnvLogScore = 0.0;
					if(cnvCount > 0)
						cnvLogScore = log(((cnvSum/cnvCount)*(zWinSize*2.0))+1.0);
					inCount = log((inCount+1));
					double gcPerc = gcCount/zWinSize;
					double aPerc = alignCount/(zWinSize*2.0);
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
		delete [] basepair;
		basepair = NULL;
		delete [] gcContent;
		gcContent = NULL;
		delete [] alignability;
		alignability = NULL;
		if(readInput == 1){
			delete [] ibasepair;
			ibasepair = NULL;
		}
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

int calcCovs::importRawSignal(const char * signalFile,int extension,const char * filetype,int dataType,const char * twoBitFile){

	if(tbSizeFlag == 0){
		FILE * tempTB;
		time_t rtime;
		struct tm *timeinfo;

		char tInfo[128];// = "tempInfo.txt";
		char tChrSize[128];// = "tempChromSize.txt";
		char sysCall[256];
		time(&rtime);
		timeinfo=localtime(&rtime);
		strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
		strftime(tChrSize,128,"tempChromSize_%H_%M_%S.txt",timeinfo);

		tempTB = fopen(tInfo,"w");
		fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
		fclose (tempTB);

		cout << "\tGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
		int s = 1;
		int twobitCount = 0;
		sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitInfo again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "TwoBitInfo failed, exiting" << endl;
				return 1;
			}
		}
		remove(tInfo);

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
		tbSizeFlag = 1;
	}

	const char * bowtie = "bowtie";
	const char * bed = "bed";
	const char * tagAlign = "tagAlign";
	int rval = 0;

	if(strcmp(filetype,bed) == 0){
		rval = importBed(signalFile,extension,dataType);
	}else if (strcmp(filetype,bowtie) == 0){
		rval = importBowtie(signalFile,extension,dataType);
	}else if(strcmp(filetype,tagAlign) == 0){
		rval = importTagAlign(signalFile,extension,dataType);
	}else{
		cout << "Unrecognized type of file " << filetype << ", must be either bowtie, bed, or tagAlign" << endl;
		return 1;
	}

	if(rval == 0){
		if(dataType == 0){
			cout << "\tImported " << signal_slist.size() << " reads" << endl;
			cout << "\tSorting reads ...";
			sort (signal_slist.begin(), signal_slist.end());
		}else if(dataType == 1){
			cout << "\t\tImported " << input_slist.size() << " inpput reads" << endl;
			cout << "\t\tSorting reads ...";
			sort (input_slist.begin(), input_slist.end());
		}
		cout << "COMPLETE" << endl;
	}else{
		cout << "Stopping buildwindows" << endl;
		return 1;
	}
	return 0;
}

int calcCovs::importBowtie(const char * signalFile,int extension,int dataType){
	
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
		
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	char line[512];char seq[128];
	char name[128];char sscore[128];int ival;
	int extend = (int)(extension/2);
	int rval;

	while(!feof(fh)){
		rval = fscanf(fh,"%s%s%s%lu%s%s%i",name,strand,cChrom,&pos,seq,sscore,&ival);
		fgets(line,512,fh);
		if(rval == 7){
			if(strcmp(strand,minus) == 0){
				if((pos + strlen(seq)) >= extend)
					pos = (pos + strlen(seq)) - extend + 1;
				else
					pos = 1;
			}else{
				if((pos + extend-1) <= chr_size[getHashValue(cChrom)])
					pos += extend - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				signal_slist.push_back(sig);
			else if(dataType == 1)
				input_slist.push_back(sig);
		}
	}
	fclose(fh);
	return 0;
}

int calcCovs::importTagAlign(const char * signalFile,int extension,int dataType){
	
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char seq[128];int score;
	int extend = (int)(extension/2);
	int rval;
	
	while(!feof(fh)){
		rval = fscanf(fh,"%s%lu%lu%s%i%s",cChrom,&start,&stop,seq,&score,strand);
		if(rval == 6){
			if(strcmp(strand,minus) == 0){
				if(stop >= extend)
					pos = stop - extend + 1;
				else
					pos = 1;
			}else{
				if((start + extend-1) <= chr_size[getHashValue(cChrom)])
					pos = start + extend - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				signal_slist.push_back(sig);
			else if(dataType == 1)
				input_slist.push_back(sig);
		}
	}
	fclose(fh);
	return 0;
}

int calcCovs::importBed(const char * signalFile,int extension,int dataType){
	
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char name[128];int bscore;
	int extend = (int)(extension/2);
	int rval;
	
	while(!feof(fh)){
		rval = fscanf(fh,"%s%lu%lu%s%i%s",cChrom,&start,&stop,name,&bscore,strand);
		if(rval == 6){
			if(strcmp(strand,minus) == 0){
				if(stop >= extend)
					pos = stop - extend + 1;
				else
					pos = 1;
			}else{
				if((start + extend-1) <= chr_size[getHashValue(cChrom)])
					pos = start + extend - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				signal_slist.push_back(sig);
			else if(dataType == 1)
				input_slist.push_back(sig);
		}
	}
	fclose(fh);
	return 0;
}


int calcCovs::outputDataWinCount(const char * outputFile, unsigned short int currChr){
	FILE * fh;
	fh = fopen(outputFile,"w");
	fprintf(fh,"chromosome\tstart\tstop\texp_count\n");
	const char * chrom = getKey(currChr);
	slist<dataWinsCount>::iterator c = peak_wins2.begin();
	while(c != peak_wins2.end()){
		fprintf(fh,"%s\t%lu\t%lu\t%i\n",chrom,c->start,c->stop,c->eCount);
		peak_wins2.erase(c++);
	}
	fclose (fh);
	return 0;
}

int calcCovs::processWinSignal(int zWinSize, int zOffsetSize,const char * twoBitFile,string outfile,int extension,const char * filetype){

	time_t rtime;
	struct tm *timeinfo;
	FILE * tempTB;
	char tInfo[128];// = "tempInfo.txt";
	char sysCall[256];
	
	unsigned short int currchr = 999;
	unsigned short int * basepair = NULL;
	unsigned short int * ibasepair = NULL;
	unsigned char * gcContent = NULL;	
	unsigned char * alignability = NULL;
	int i;
	int printflag = 0;

	int readInput = 0;
	const char* noneVal = "none";
	
	while(!signal_slist.empty()){
		i = 0;
		currchr = signal_slist[0].chrom;
		const char * chromReport = getKey(currchr);
		cout << "\nProcessing " << chromReport << endl;
		basepair = new unsigned short int[chr_size[currchr]+1];
		for(int ch = chr_size[currchr]; ch--;)
			basepair[ch] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(signal_slist[i].chrom==currchr && i < (int) signal_slist.size()){
			basepair[signal_slist[i].pos]++;
			i++;
		}	
		signal_slist.erase(signal_slist.begin(),signal_slist.begin()+i);

		cout << "\tGetting sequence from .2bit file:\n\t\t" << twoBitFile << endl;
		gcContent = new unsigned char[chr_size[currchr] + 1];
		for(int ch = chr_size[currchr]; ch--;)
			gcContent[ch] = (unsigned char) 0;

		char tSeq[128];
		time(&rtime);
		timeinfo=localtime(&rtime);
		strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
		strftime(tSeq,128,"tempSeq_%H_%M_%S.txt",timeinfo);

		tempTB = fopen(tInfo,"w");
		fprintf(tempTB,"library(zinba);\ntwobittofa(chrm=\"%s\",start=1,end=%lu,twoBitFile=\"%s\",gcSeq=\"%s\");\n",chromReport,chr_size[currchr],twoBitFile,tSeq);
		fclose (tempTB);
		sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		int s = 1;
		int twobitCount = 0;
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitToFa again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "twoBitToFa failed, exiting" << endl;
				return 1;
			}
		}
		remove(tInfo);
		
		ifstream seqfile(tSeq);
		string line;
		int pos = 1;
		unsigned long int nStart = 0;
		unsigned long int nStop = 0;
		unsigned short int prevEnt = 0;
		while(getline(seqfile,line)){
			if(line[0] != '>'){
				for(int s = 0; s < line.length();s++){
					if(line[s] == 'G' || line[s] == 'C'){
						gcContent[pos] = (unsigned char) 1;
					}else if (line[s] == 'N'){
						gcContent[pos] = (unsigned char) 2;
					}else{
						gcContent[pos] = (unsigned char) 0;
					}
					pos++;
				}
			}
		}
		seqfile.close();
		remove(tSeq);


		cout << "\tGetting counts for zinba windows.........." << endl;
		int numOffsets = 1;
		if(zOffsetSize > 0){
			numOffsets = (int) zWinSize/zOffsetSize;
		}
		string outfileDATA;
		slist<dataWinsCount>::iterator z;
		for(int o = 0; o < numOffsets; o++){
			z = peak_wins2.previous(peak_wins2.end());
			cout << "\t\tOffset " << (zOffsetSize * o) << "bp......" << endl;
			char offset[128];
			sprintf(offset,"%d",(zOffsetSize * o));
			char winsize[128];
			sprintf(winsize,"%d",zWinSize);
			outfileDATA = outfile + "_" + chromReport + "_win" + string(winsize) + "bp_offset" + string(offset) + "bp.txt";
			unsigned long int zWinStart = (zOffsetSize * o) + 1;
			unsigned long int zWinStop = zWinStart + zWinSize - 1;
			while(zWinStop <= chr_size[currchr]){
				int peakCount = 0;
				double gcCount = 0;
				double nCount = 0;
				
				for(int b = zWinStart; b <= zWinStop; b++){
					peakCount += basepair[b];
					if((int) gcContent[b] == 1)
						gcCount++;
					else if((int) gcContent[b] == 2)
						nCount++;
				}
				if((nCount/zWinSize) < 0.1){
					dataWinsCount zwin(currchr,zWinStart,zWinStop,peakCount);
					z = peak_wins2.insert_after(z,zwin);
				}
				zWinStart += zWinSize;
				zWinStop += zWinSize;
			}
			
			if(outputDataWinCount(outfileDATA.c_str(),currchr) != 0){
				cout << "Error printing output to file, exiting" << endl;
				exit(1);
			}
			
		}
		delete [] basepair;
		basepair = NULL;
		}
	return 0;
}

