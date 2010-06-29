#include "abp.h"

abp::abp(){}

abp::~abp(){}

abp::abp(const abp& s){
	score = s.score;
	ascore = s.ascore;
//	strand = s.strand;
}

//bwRead::bwRead(unsigned short int _chrom, unsigned long int _pos,unsigned short int _strand){
abp::abp(unsigned short int _score, unsigned long int _ascore){
	score = _score;
	ascore = _ascore;
//	strand = _strand;
}

bool abp::operator <(abp other) const{
		return (score<other.score);
}