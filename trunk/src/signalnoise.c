#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define min(a,b)  ((a) < (b) ? (a) : (b))
#define maxx(x,y) ((x) > (y) ? (x) : (y))
#define MAX_LEN 5000
#define MAX_LENB 250000000

int int_cmp(const void *a, const void *b)
{
    const int *ia = (const int *)a; // casting pointer types
    const int *ib = (const int *)b;
    return *ia  - *ib; 
	/* integer comparison: returns negative if b > a 
	and positive if a > b */
}


void signalnoise(char **RsbpcFile, char **Routputfile ,int *RwinSize){
  int i, j, lbasecount, winstart, winstop, winSize;
  winSize=0;
  FILE *FI, *FO;
  FILE * tempTB;
  char str_buf[MAX_LEN], line[MAX_LEN];
  char chr[127];
  int basecount[MAX_LENB];
  winSize=*RwinSize;
    //Rprintf("%d\n",winSize);
  *RwinSize=winSize+1;
  int window[5];
  size_t numbers_len = sizeof(window)/sizeof(int);
  const char *input = RsbpcFile[0];
  const char *output = Routputfile[0];

  tempTB = fopen(output,"w");
  FI = fopen(input, "r");	
  if(FI == NULL){		
    error("cannot open file %s\n", input);
  }
  if(tempTB == NULL){		
    error("cannot open file %s\n", input);
  }
  //skip the first line of the wig file, usually the "track" line
    if(fgets(str_buf, MAX_LEN, FI) == NULL){
      error("there are no lines in file\n");
    }


Rprintf("Begin Signal to Noise Analysis\n"); 
//read in a line, save the information in each, perform peakbounds, then print out, repeat for each line
i=0;
j=0;
winstart=0;
winstop=winstart+winSize;
int m=0;
 while(fgets(str_buf, MAX_LEN, FI) != NULL){

  if (str_buf[0] != 't' && str_buf[0] != 'f') {
  	basecount[m]=atol(str_buf);
	m++;
  }else if (i>0){
	//Done reading information into basecount
	//Now perform window analysis
	while(winstop<m){
		//copy basecount corresponding to window position
		for(j=winstart;j<winstop;j++){
			window[j-winstart]=basecount[j];		
		}
		qsort(window, numbers_len, sizeof(int), int_cmp); 
		fprintf(tempTB,"%i\t%i\n",window[(int) winSize/2],window[winSize-1]);
		winstart=winstop;
		winstop=winstart+winSize;
	}
	//Reset basecount vector
  	for(j=0; j<MAX_LENB;j++){
		basecount[j]=0;
  	}
	//reset basecount index to zero
	m=0;
	//reset winstart
	winstart=0;
  }else{
	//Reset basecount vector
	for(j=0; j<MAX_LENB;j++){
		basecount[j]=0;
  	}
	i++;
  }

  }
fclose(tempTB);
}

