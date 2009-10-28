#include "dataWins.h"

dataWins::~dataWins(){}

dataWins::dataWins(){}

dataWins::dataWins(const dataWins& s){
	chrom = s.chrom;
	start = s.start;
	stop = s.stop;
	eCount = s.eCount;
	iCount = s.iCount;
	gcPerc = s.gcPerc;
	alignPerc = s.alignPerc;
	cnvScore = s.cnvScore;
}

dataWins::dataWins(unsigned short int _chrom, unsigned long int _start, unsigned long int _stop,int _eCount, int _iCount, double _gcPerc, double _alignPerc, double _cnvScore){
	chrom = _chrom;
	start = _start;
	stop = _stop;
	eCount = _eCount;
	iCount = _iCount;
	gcPerc = _gcPerc;
	alignPerc = _alignPerc;
	cnvScore = _cnvScore;
}

bool dataWins::operator <(dataWins other) const{
	if(chrom==other.chrom){
		if( start < other.start ){
			return true;
		}else{
			return false;	
		}
	}else{
		return (chrom<other.chrom);	
	}
}