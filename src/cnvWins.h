#ifndef CNVWINS_H_
#define CNVWINS_H_

#include <iostream>

class cnvWins {
	
	public:
			cnvWins();
			~cnvWins();
			cnvWins(unsigned short int, unsigned long int, unsigned long int,double,double,double,double);
			cnvWins(const cnvWins&);
			bool operator <(cnvWins) const;

			unsigned short int chrom;//2
			unsigned long int start;//8
			unsigned long int stop;
			double cnvScore;
			double varWin;
			double chiSq;
			double pval;
	
	private:
			
};

#endif /*CNVWINS_H_*/