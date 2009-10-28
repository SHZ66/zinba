#ifndef BCANALYSIS_H_
#define BCANALYSIS_H_

#include <string>
#include <vector>
#include "aRead.h"
//#include "base.h"
#include <map>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class bcanalysis{
	
	public:
	
		bcanalysis(); //Implemented
		~bcanalysis(); //Implemented
		
		int importRawSignal(const char *);//Implemented
		int processSignals(const char *,int);//Implemented
		int outputData(const char *, unsigned short int,unsigned short int,unsigned short int[]);//Implemented
		
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
		
		slist<aRead> signal_slist;//Implemented	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
};

#endif /*BCANALYSIS_H_*/
