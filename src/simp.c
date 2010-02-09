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
#define MAX_LEN 50000

SEXP mkans(int *, int);
//int feval(int *, SEXP ,SEXP );

SEXP peakboundc(SEXP f, SEXP bpprofile, SEXP outputfile, SEXP rho){

  int i, j, lbasecount, pstart, pstop;
  int basecount[MAX_LEN];
  FILE *FI, *FO;
  char str_buf[MAX_LEN], line[MAX_LEN];
  char chr[127], ID[127], strand[2];
	char *delim = "\t";
	char *ptr= NULL;
  const char *input =CHAR(STRING_ELT(bpprofile,0));
  const char *output=CHAR(STRING_ELT(outputfile,0));

/////////////////
double sig;
/////////////////
	
  for(j=0; j<MAX_LEN;j++){
	basecount[j]=0;
  }

  
  FI = fopen(input, "r");	
  if(FI == NULL){		
    error("cannot open file %s\n", input);
  }

  FO = fopen(output, "w");	
  if(FO == NULL){
    error("cannot open file %s\n", output);
  }
  
  

skip the first line of the coordinate file
    if(fgets(str_buf, MAX_LEN, FI) == NULL){
      error("there are only %d lines in file %s\n", i, input);
    }
//print headers to file
sprintf(line, "PEAKID\tChrom\tStart\tStop\tStrand\tSig\tMaxloc\tMax\tpStart\tpStop\n");
fputs(line, FO);
Rprintf("Begin Peak Refinement\n"); 
//read in a line, save the information in each, perform peakbounds, then print out, repeat for each line
int m=0;
 while(fgets(str_buf, MAX_LEN, FI) != NULL){
    i = 0;
     if ((ptr = strtok(str_buf, delim)) != NULL) {
    do {
      i++;
      if(i==1){
        strcpy(ID, ptr);
      }else if(i==2){
        strcpy(chr, ptr);
      }else if(i==3){
        pstart = atol(ptr);
      }else if(i==4){
        pstop = atol(ptr);
      }else if(i==5){
        strcpy(strand, ptr);
      }else if(i==6){
		  sig = atof(ptr);
	  }else if(i>6){
		  basecount[i-7]  = atoi(ptr);
//	  }else if(i>5){
//        basecount[i-6]  = atoi(ptr);
      }
//      lbasecount=i-5;
		lbasecount=i-6;
    } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", input);
  }
 
  defineVar(install("x"), mkans(basecount, lbasecount), rho);
    int l=round(INTEGER(eval(f, rho))[0]);
    int maxvec[l];
  for(i=1;i<=l;i++){
  	defineVar(install("x"), mkans(basecount, lbasecount), rho);
  	maxvec[i-1]=INTEGER(eval(f, rho))[i];
  }
////////////////////////////////////////////////////////////////////////

int lmaxvec = l;
int results[lmaxvec*2];
int max;
int searchlength=0;
int  peakstart=0;
int  peakend=0;
int  fit=0;
int  bound,h, n;
int maxvecbounds[lmaxvec*2];
double val,sumxy,sumx, sumx2, sumy, sumy2;

for(j=0; j<lmaxvec; j++){
	max=maxvec[j];
	double hold[lbasecount];
	for(i=0; i<lbasecount; i++){
		hold[i]=0;
	}
	searchlength=max-1;		
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
	//right bound
	searchlength=lbasecount-max-1;
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
//save boundaries, adjust back to R index start at 1 
maxvecbounds[2*j]=peakstart;
maxvecbounds[2*j+1]=peakend;
}	

int dist, sep;
int hmaxvec[lmaxvec];

//final outputted vector of bounds after merging, contains zeros
//intialize results vector
for(i=0; i<2*lmaxvec;i++){
	results[i]=0;
}
results[0]=maxvecbounds[0];
results[1]=maxvecbounds[1];
hmaxvec[0]=maxvec[0];
i=0;
Rprintf("i %d h %d \tmaxpos %d sep %d, dist %d\t maxvecleft %d maxvecright %d resultshleft %d resultshright %d hmaxvec[h] %d \n", i, h, maxvec[i],sep, dist, maxvecbounds[2*i], maxvecbounds[2*i+1], results[2*h], results[2*h+1], hmaxvec[h]);
h=0;
i=1;
while(i<lmaxvec){
	//calculate metrics per round
	dist=maxvec[i]-maxvec[i-1];
	sep=maxvecbounds[2*i]-results[2*h+1];
	
	if(dist<=250){
	//merge if too close
		results[2*h]=min(results[2*h],maxvecbounds[2*i]);
		results[2*h+1]=maxx(results[2*h+1],maxvecbounds[2*i+1]);			if(basecount[hmaxvec[h]]<basecount[maxvec[i]]){		
			hmaxvec[h]=maxvec[i];
		}
	}else if(dist>250 & sep<100){
		//merge if not close enough but if bound overlap and no valley inbetween 
		//need to normalize with respect to min of whole basecount?
		results[2*h]=min(results[2*h],maxvecbounds[2*i]);
		results[2*h+1]=maxx(results[2*h+1],maxvecbounds[2*i+1])	;
		if(basecount[hmaxvec[h]]<basecount[maxvec[i]]){		
			hmaxvec[h]=maxvec[i];
		}
	}else if(dist>250 & sep>=100){
		//separate peak here, no boundary overlap and too far 
		h=h+1;
		hmaxvec[h]=maxvec[i];
		results[2*h]=maxvecbounds[2*i];
		results[2*h+1]=maxvecbounds[2*i+1]	;
	}
	Rprintf("i %d h %d \tmaxpos %d sep %d, dist %d\t maxvecleft %d maxvecright %d resultshleft %d resultshright %d hmaxvec[h] %d \n", i, h, maxvec[i],sep, dist, maxvecbounds[2*i], maxvecbounds[2*i+1], results[2*h], results[2*h+1], hmaxvec[h]);
	i=i+1;
}
/////////////////////////////////////////////////////////////////////////
int pos=0;

for(i=0;i<=h;i++){
Rprintf("i %d hmaxvec[i] %d", i, hmaxvec[i]);
Rprintf("%s\t%s\t%d\t%d\t%s\t%.14f\t%d\t%d\t%d\t%d\n",ID, chr, pstart, pstop, strand,sig ,hmaxvec[i],basecount[hmaxvec[i]],results[2*i],results[2*i+1]);
	
sprintf(line, "%s\t%s\t%d\t%d\t%s\t%.14f\t%d\t%d\t%d\t%d\n",ID, chr, pstart, pstop, strand,sig ,pstart+maxvec[i],basecount[maxvec[i]],pstart+results[2*i],pstart+results[2*i+1]);
fputs(line, FO);
}



}
Rprintf("Peak Refinement Complete\n"); 
 fclose(FO);
 fclose(FI);
 return(bpprofile);
}



SEXP mkans(int *x, int len)
     {
         SEXP ans;
	 	
	 int l = len;
         PROTECT(ans = allocVector(INTSXP, l));
	 for(int i=0;i<l; i++){
         	INTEGER(ans)[i] = x[i];
	 }
         UNPROTECT(1);
         return ans;
     }
     
   /*  int feval(int *x, SEXP f,SEXP rho)
     {
         defineVar(install("x"), mkans(x), rho);
	 //PROTECT(ans = allocVector(INTSXP, l));
	 return(INTEGER(eval(f, rho))[1]);
     }
     
    SEXP zero(SEXP f, SEXP guesses, SEXP tol, SEXP rho)
     {
         int * x0 = INTEGER(guesses);
         int f0, f1, fc, xc;
     
        f0 = feval(x0, f, rho);
	defineVar(install("x"), mkans(x0), rho);
	int l=round(INTEGER(eval(f, rho))[0]);
	int maxvec[l];
	for(int i=1;i<=l;i++){
		defineVar(install("x"), mkans(x0), rho);
		maxvec[i-1]=INTEGER(eval(f, rho))[i];
	}

		
	return mkans(x0);
     }*/

