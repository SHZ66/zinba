#include "aRead.h"

aRead::aRead(){}

aRead::~aRead(){}

aRead::aRead(const aRead& s){
	chrom = s.chrom;
	start = s.start;
}

aRead::aRead(unsigned short int _chrom, unsigned long int _start){
	chrom = _chrom;
	start = _start;
}

bool aRead::operator <(aRead other) const{
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