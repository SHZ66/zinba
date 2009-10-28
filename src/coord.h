#ifndef COORD_H_
#define COORD_H_

#include <iostream>
#include <cstring>

class coord {
	
	public:
			coord();
			~coord();
			coord(char[],unsigned short int, unsigned long int, unsigned long int,char[]);
			coord(const coord&);
			bool operator <(coord) const;

			char * ident;
			unsigned short int chrom;//2
			unsigned long int start;//8
			unsigned long int end;//8
			unsigned short int qFlag;
			char * strand;
	
	private:
			
};

#endif /*COORD_H_*/