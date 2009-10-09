#include "maxbp.h"

maxbp::maxbp(){
}

maxbp::~maxbp(){
}

maxbp::maxbp(const maxbp& c){
	chrom = c.chrom;
	position = c.position;
	score = c.score;
}

maxbp::maxbp(unsigned short int _chrom, unsigned long int _position, int _score){
	chrom = _chrom;
	position = _position;
	score = _score;
}

bool maxbp::operator<(const maxbp& other) const{
		return position < other.position;
}

bool maxbp::operator==(const maxbp& other) const{
	return position == other.position;
}

bool maxbp::operator>(const maxbp& other) const{
	return position > other.position;
}