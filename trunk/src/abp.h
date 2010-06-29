#ifndef ABP_H_
#define ABP_H_

#include <iostream>

class abp {
	
	public:
			abp();
			~abp();
//			bwRead(unsigned short int, unsigned long int,unsigned short int);
			abp(unsigned short int, unsigned long int);
			abp(const abp&);
			bool operator <(abp) const;

			unsigned short int score;//2
			unsigned long int ascore;//8
//			unsigned short int strand;//Plus = 1;Minus = 0
	
	private:
			
};

#endif /*ABP_H_*/