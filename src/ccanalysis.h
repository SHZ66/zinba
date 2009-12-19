#ifndef CCANALYSIS_H_
#define CCANALYSIS_H_

#include <string>
#include <cstring>
#include <vector>
#include "bwRead.h"
//#include "base.h"
#include <map>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class ccanalysis{
	
	public:
	
		ccanalysis(); //Implemented
		~ccanalysis(); //Implemented
		
		int importRawSignal(const char *,const char *);//Implemented
		int processSignals(const char *,const char *);//Implemented
		
		struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};
	
		struct ltint{
			bool operator()(const int i1, const int i2) const
			{
				return i1<i2;
		 	}
		};

	private:
		
		slist<bwRead> reads_slist;//Implemented	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;
};

#endif /*CCANALYSIS_H_*/
