#include "cnvWins.h"

cnvWins::~cnvWins(){}

cnvWins::cnvWins(){}

cnvWins::cnvWins(const cnvWins& s){
	chrom = s.chrom;
	start = s.start;
	stop = s.stop;
	cnvScore = s.cnvScore;
	varWin = s.varWin;
	chiSq = s.chiSq;
	pval = s.pval;
}

cnvWins::cnvWins(unsigned short int _chrom, unsigned long int _start, unsigned long int _stop,double _cnvScore,double _varWin, double _chiSq, double _pval){
	chrom = _chrom;
	start = _start;
	stop = _stop;
	cnvScore = _cnvScore;
	varWin = _varWin;
	chiSq = _chiSq;
	pval = _pval;
}

bool cnvWins::operator <(cnvWins other) const{
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