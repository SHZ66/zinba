#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

#define min(a,b)  ((a) < (b) ? (a) : (b))
#define maxx(x,y) ((x) > (y) ? (x) : (y))

void peakboundc(int basecount[],int *Rlbasecount, int maxvec[], int *Rlmaxvec, int *Rsearchlength, int results[])
{


int lbasecount = *Rlbasecount;
int lmaxvec = *Rlmaxvec;
int max;
int searchlength=*Rsearchlength;
int  peakstart=0;
int  peakend=0;
int  fit=0;
int  bound,h, i, j, n;
int maxvecbounds[lmaxvec*2];
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

int dist, l, mingap, indexmingap, sep;

//final outputted vector of bounds after merging, contains zeros
//intialize results vector
results[0]=maxvecbounds[0];
results[1]=maxvecbounds[1];

h=0;
i=1;
while(i<lmaxvec){
	//calculate metrics per round
	dist=maxvec[i]-maxvec[i-1];
	val=10000000;
	fit=0;
	for(l=maxvec[i-1]-1;l<maxvec[i]; l++){
      		if(basecount[l] < val){
       			val = basecount[l];
			fit=l;
		}
	}
	mingap=val;
	indexmingap=fit;
	sep=maxvecbounds[2*i]-results[2*h+1];

	if(dist<=100){
	//merge if too close
		results[2*h]=min(results[2*h],maxvecbounds[2*i]);
		results[2*h+1]=maxx(results[2*h+1],maxvecbounds[2*i+1]);		
	}else if(dist>100 & mingap>1.0*basecount[maxvec[i-1]]/2 & sep<40){
		//merge if not close enough but if bound overlap and no valley inbetween 
		//need to normalize with respect to min of whole basecount?
		results[2*h]=min(results[2*h],maxvecbounds[2*i]);
		results[2*h+1]=maxx(results[2*h+1],maxvecbounds[2*i+1])	;
	}else if(dist>100 & mingap<=1.0*basecount[maxvec[i-1]]/2 & sep<40){
		//separate if not close enough with some bound overlap with valley inbetween 
		//set new boundary in between at the site of the min in between
		results[2*h+1]=indexmingap;
		h=h+1;
		results[2*h]=indexmingap +1 ;
		results[2*h+1]=maxvecbounds[2*i+1];
		}else if(dist>100 & sep>=40){
		//separate peak here, no boundary overlap and too far 
		h=h+1;
		results[2*h]=maxvecbounds[2*i];
		results[2*h+1]=maxvecbounds[2*i+1]	;
	}
	i=i+1;
}


}
