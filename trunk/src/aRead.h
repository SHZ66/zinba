#ifndef AREAD_H_
#define AREAD_H_

#include <iostream>

class aRead {
	
	public:
			aRead();
			~aRead();
			aRead(unsigned short int, unsigned long int);
			aRead(const aRead&);
			bool operator <(aRead) const;

			unsigned short int chrom;//2
			unsigned long int start;//8
			
	private:
			
};

#endif /*AREAD_H_*/