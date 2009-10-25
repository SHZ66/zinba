#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

#define min(a,b)  ((a) < (b) ? (a) : (b))

void peakboundc(int basecount[],int *Rlbasecount, int maxvec[], int *Rlmaxvec, int *Rsearchlength, int maxvecbounds[])
{


int lbasecount = *Rlbasecount;
int lmaxvec = *Rlmaxvec;
int max;
int searchlength=*Rsearchlength;
int  peakstart=0;
int  peakend=0;
int  fit=0;
int  bound,h, i, j, n;
double val,sumxy,sumx, sumx2, sumy, sumy2;

for(j=0; j<lmaxvec; j++){
	max=maxvec[j];
	double hold[lbasecount];
	for(i=0; i<lbasecount; i++){
		hold[i]=0;
	}
	if(min(searchlength, max)<50){
		peakstart=1;
	}else{		
		for(bound=50; bound<=min(searchlength, max); bound++){
			sumxy=0;
			sumx=0;
			sumy=0;
			sumx2=0;
			sumy2=0;
			n=bound;
			for(h=max-bound; h<=max-1;h++){ 
				sumxy=sumxy+(h+1)*basecount[h];
				sumx=sumx+h+1;
				sumy=sumy+basecount[h];
				sumx2=sumx2+(h+1)*(h+1);
				sumy2=sumy2+basecount[h]*basecount[h];						
			}
			hold[max-bound]=pow(sumxy-sumy*sumx/n,2)/((sumy2-sumy*sumy/n)*(sumx2-sumx*sumx/n))	;
		}
		//find max R2 position
		val=0;
		fit=0;
		for(i=0; i<lbasecount; i++){
        		if(hold[i] > val){
            		val = hold[i];
			fit=i;
        		}
			hold[i]=0;
    		}
		peakstart=fit+1;
	}
	//right bound
	if(min(searchlength, lbasecount-max)<50){
		peakend=lbasecount;
	}else{		

		for(bound=50; bound<=min(searchlength, lbasecount-max); bound++){
			sumxy=0;
			sumx=0;			
			sumy=0;
			sumx2=0;
			sumy2=0;
			n=bound;
			for(h=max-1; h<max+bound-1;h++){ 
				sumxy=sumxy+(h+1)*basecount[h];
				sumx=sumx+h+1;
				sumy=sumy+basecount[h];
				sumx2=sumx2+(h+1)*(h+1);
				sumy2=sumy2+basecount[h]*basecount[h];
					
			}
			hold[max-1+bound]=pow(sumxy-sumy*sumx/n,2)/((sumy2-sumy*sumy/n)*(sumx2-sumx*sumx/n));

			
			
		}
		//find max R2 position
		val=0;
		fit=0;
		for(i=0; i<lbasecount; i++){
        		if(hold[i] > val){
            		val = hold[i];
			fit=i;
        		}
		//	hold[i]=0;
    		}
		peakend=fit;
	}
//save boundaries, adjust back to R index start at 1 
maxvecbounds[2*j]=peakstart;
maxvecbounds[2*j+1]=peakend;
}	


}
