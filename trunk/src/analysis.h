#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <string>
#include <vector>
#include "coord.h"
#include <map>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class analysis{
	
	public:
	
		analysis(); //Implemented
		~analysis(); //Implemented
		int importCoords(const char *);//Implemented
		int processCoords(const char *,const char *,string);//Implemented
		int outputData(const char *,int,const char *,unsigned short int,unsigned long int,unsigned long int,const char *,int,int[]);//Implemented
//		int outputData(const char *,int,const char *,unsigned short int,unsigned long int,unsigned long int,const char *,int,double[],const char *);//Implemented
	
		struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};

	private:
		
		slist<coord> coord_slist;//Implemented
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<int, int> chrom_start;//Implemented
		map<int, int> chrom_stop;//Implemented
};

#endif /*ANALYSIS_H_*/
