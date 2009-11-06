#include "coord.h"

coord::coord(){}

coord::~coord(){}

coord::coord(const coord& c){
	chrom = c.chrom;
	start = c.start;
	end = c.end;
	qFlag = c.qFlag;
}

coord::coord(unsigned short int _chrom, unsigned long int _start, unsigned long int _end, unsigned short int _qFlag){
	chrom = _chrom;
	start = _start;
	end = _end;
	qFlag = _qFlag;
}

bool coord::operator <(coord other) const{
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