#ifndef MAXBP_H_
#define MAXBP_H_

#include <iostream>
#include <cstring>

class maxbp{

	public:
			maxbp();
			~maxbp();
			maxbp(unsigned short int, long int, int);
			maxbp(const maxbp&);
			bool operator<(const maxbp& other) const;
			bool operator==(const maxbp& other) const;
			bool operator>(const maxbp& other) const;

			unsigned short int chrom;//2
			long int position;//8
			int score;
	
	private:
			
};

#endif /*MAXBP_H_*/