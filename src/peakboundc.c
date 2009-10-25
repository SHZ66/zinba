#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

#define min(a,b)  ((a) < (b) ? (a) : (b))

void peakboundc(int **Rbasecount,int *Rlbasecount, int **Rmaxvec, int *Rlmaxvec, int *Rsearchlength, double maxvecbounds[])
{

int* basecount = Rbasecount[0];
int lbasecount = *Rlbasecount;
int* maxvec= Rmaxvec[0];
int lmaxvec = *Rlmaxvec;
int searchlength=*Rsearchlength;
int max, peakstart, peakend, fit, bound,val,h, i, n, sumxy,sumx, sumx2, sumy, sumy2, tmp;

/*
for(i=0; i<lmaxvec; i++){
	//adjust for C start index at 0
	max=maxvec[i]-1;
	double hold[min(searchlength, max+1)];
	if(min(searchlength, max+1)<50){
		peakstart=1;
	}else{		
		for(bound=min(searchlength, max-1); bound>=50; bound--){
			sumxy=0;
			sumx=0;
			sumy=0;
			sumx2=0;
			sumy2=0;
			n=bound+1;
			for(h=max-bound; h<=max;h++){ 
				sumxy=sumxy+h*basecount[h];
				sumx=sumx+h;
				sumy=sumy+basecount[h];
				sumx2=sumx2+pow(h,2);
				sumy2=sumy2+pow(basecount[h],2);					
			}
			hold[bound]=pow(sumxy-sumy*sumx/n,2)/((sumy2+sumy*sumy/n+sumy*sumy/(n*n))*(sumx2+sumx*sumx/n+sumx*sumx/(n*n)));	
		}
		//find max R2 position
		val=0;
		fit=0;
		for(h=0; h<min(searchlength, max+1); h++){
        		if(hold[h] > val){
            		val = hold[min(searchlength, max+1)];
			fit=h;
        		}
    		}
		peakstart=max-fit;
	}
	//right bound
	double hold2[min(searchlength, lbasecount-max+1)];

	if(min(searchlength, lbasecount-max+1)<50){
		peakend=lbasecount;
	}else{
		for(bound=50; bound<=min(500, searchlength-max); bound++){
			sumxy=0;
			sumx=0;
			sumy=0;
			sumx2=0;
			sumy2=0;
			n=bound+1;
			for(h=max; h<=max+bound;h++){ 
				sumxy=sumxy+h*basecount[h];
				sumx=sumx+h;
				sumy=sumy+basecount[h];
				sumx2=sumx2+pow(h,2);
				sumy2=sumy2+pow(basecount[h],2);					
			}
			hold2[bound]=pow(sumxy-sumy*sumx/n,2)/((sumy2+sumy*sumy/n+sumy*sumy/(n*n))*(sumx2+sumx*sumx/n+sumx*sumx/(n*n)));	
		}
		//find max R2 position
		val=0;
		fit=0;
		for(h=0; h<min(searchlength, lbasecount-max+1); h++){
        		if(hold2[h] > val){
            		val = hold2[h];
			fit=h;
        		}
    		}
		peakend=max+fit;
	}
//save boundaries, adjust back to R index start at 1
maxvecbounds[2*i]=peakstart+1;
maxvecbounds[2*i+1]=peakend+1;
	
} 
*/
maxvecbounds[0]=basecount[0]
}
