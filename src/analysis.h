#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <string>
#include <vector>
#include "coord.h"
#include <map>
#include <list>

using namespace std;

class analysis{
	
	public:
	
		analysis(); //Implemented
		~analysis(); //Implemented
		int importCoords(const char *);//Implemented
		int processCoords(const char *,const char *,string);//Implemented
		int outputData(const char *,int,unsigned short int,long int,long int,int,unsigned short int[]);//Implemented

		struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};

	private:
		
		list<coord> coord_slist;//Implemented
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented

};

#endif /*ANALYSIS_H_*/
