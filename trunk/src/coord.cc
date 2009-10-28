#include "coord.h"

coord::coord(){
	ident = new char[128];
	strand = new char[1];
}

coord::~coord(){
	delete [] ident;
	delete [] strand;
}

coord::coord(const coord& c){
	ident = new char[128];
	strcpy(ident,c.ident);
	chrom = c.chrom;
	start = c.start;
	end = c.end;
	qFlag = c.qFlag;
	strand = new char[1];
	strcpy(strand,c.strand);
}

coord::coord(char * _ident, unsigned short int _chrom, unsigned long int _start, unsigned long int _end, unsigned short int _qFlag,char * _strand){
	ident = new char[128];
	strcpy(ident,_ident);
	chrom = _chrom;
	start = _start;
	end = _end;
	qFlag = _qFlag;
	strand = new char[1];
	strcpy(strand,_strand);
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