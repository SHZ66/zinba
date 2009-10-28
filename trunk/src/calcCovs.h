#ifndef CALCCOVS_H_
#define CALCCOVS_H_

#include <string>
#include <vector>
#include "aRead.h"
#include "dataWins.h"
#include "cnvWins.h"
#include <cstring>
//#include "base.h"
#include <map>
#include <list>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class calcCovs{
	
	public:
	
		calcCovs(); //Implemented
		~calcCovs(); //Implemented
		
		int importRawSignal(const char *,int);//Implemented
		int processSignals(int,int,int,int,string,string,const char*,const char*);//Implemented
		int outputData(const char *, unsigned short int);//Implemented
		
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
		slist<aRead> input_slist;
		list<cnvWins> cnv_wins;
		slist<dataWins> peak_wins;
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;
};

#endif /*CALCCOVS_H_*/
