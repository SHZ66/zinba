#include "bwRead.h"

bwRead::bwRead(){}

bwRead::~bwRead(){}

bwRead::bwRead(const bwRead& s){
	chrom = s.chrom;
	pos = s.pos;
//	strand = s.strand;
}

//bwRead::bwRead(unsigned short int _chrom, unsigned long int _pos,unsigned short int _strand){
bwRead::bwRead(unsigned short int _chrom, unsigned int _pos){
	chrom = _chrom;
	pos = _pos;
//	strand = _strand;
}

bool bwRead::operator <(bwRead other) const{
//	if(chrom==other.chrom){
//		if( pos < other.pos ){
//			return true;
//		}else{
//			return false;	
//		}
//	}else{
		return (chrom<other.chrom);	
//	}
}