#ifndef ZCOLLAPSEWINS_H_
#define ZCOLLAPSEWINS_H_

#include <string>
#include <cstring>
#include <vector>
#include "coord.h"
#include <map>
#include <list>

using namespace std;

class zCollapseWins{
	
	public:
	
		zCollapseWins(); //Implemented
		~zCollapseWins(); //Implemented
		int importCoords(const char *,const char *,unsigned short int,unsigned short int,double,string,int);//Implemented
		int outputData(const char *);//Implemented

		struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};

	private:
		
		list<coord> coord_slist;//Implemented
		list<coord> coordOut_list;
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(const char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
};

#endif /*ZCOLLAPSEWINS_H_*/
