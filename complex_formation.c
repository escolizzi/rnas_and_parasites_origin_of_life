#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cash2.h"
#include "complex_formation.h"
#include "parameters.h"
#include "constants.h"
#include "mersenne.h"

void IntProbComplForm(int *s1, int *s2, double *k_asso, double *k_diss)
{  
  int i, maxmatch, lens2;
  int revs2[MAXLEN];
  
  lens2=MAXLEN;
  for(i=0;i<lens2;i++) revs2[i]= 1-s2[lens2 -1 - i];
  maxmatch=IntMaxMatch(s1,revs2);
  
  *k_asso = 1. - exp(par_factor_asso*(double)maxmatch);
  *k_diss = exp(par_factor_diss*(double)maxmatch);
}

/*s2 must have been reversed already!!! */
int IntMaxMatch(int *s1, int *s2)
{
  int i,j,len1,len2,len_diff;
  int nmatch; /* current */
  int pnmatch=0; /* previous */
  int *end_last,*end_middle; /* dont be confused by the name. */

  len1=MAXLEN;
  len2=MAXLEN;

  if(len1==0 || len2==0)
    return 0;

  if(len1>=len2) {
    end_last = s2 + len2;			//end_last points to 'one after the end' of end5
    len_diff = len1 - len2;
    end_middle = s1 + len_diff;	//end_middle points to the point of end3 which is at the length of the difference in length between end5 and end3
    /*
      Let the length of 3' be 'll' and 5' be 'dd'.
      Let l = ll-1 and d = dd-1.
      "end_last" is pointing to "5[dd]".
      (i=1)
      3'       012........l   -> 3[0] is the begining
      5' 012...d              -> 5[dd-i] is the begining

      (i=dd)
      3' 0123........l        -> 3[0]
      5' 0123.....d           -> 5[dd-dd]
    */
    
    /*
    starts from end3
    
    */
    for(i=1,j=len1; i<=j; i++){
      nmatch = nbitpair(s2,end_last-i,i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    
    /*
      (i=1)
      3' 0123.........l     -> 3[i]
      5'  0123....d         -> 5[0]

      (i=ll-dd)
      3' 0123.........l     -> 3[(ll-dd)]
      5'      0123....d     -> 5[0]
    */
    for(i=1,j=len_diff;i<=j;i++){
      nmatch = nbitpair(s1+i,s2,len2);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    /*
      Let "end_middle" point to "3[(ll-dd)]".
      (i=1)
      3' 0123.........l      -> 3[(ll-dd)+i]
      5'       0123....d     -> 5[0]

      (i=dd-1)
      3' 0123........l            -> 3[(ll-dd)+dd-1] = 3[l]
      5'             0123.....d   -> 5[0]
    */
    for(i=1,j=len2-1;i<=j;i++){
      nmatch = nbitpair(end_middle+i,s2,len2-i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    return pnmatch;
  }
  
  else {
    len_diff = len2 - len1;
    end_last = s1 + len1;
    end_middle = s2 + len_diff;
    for(i=1,j=len1;i<=j;i++){
      nmatch = nbitpair(s2,end_last-i,i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    for(i=1,j=len_diff;i<=j;i++){
      nmatch = nbitpair(s2+i,s1,len1);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    for(i=1,j=len1-1;i<=j;i++){
      nmatch = nbitpair(end_middle+i,s1,len1-i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
  }
  return pnmatch;
}

/* Sum over the bit pair between "this" and "that". Don't forget
   the length of the sequence we want to look. */
int nbitpair(int* this,int* that,int len)
{
  int i;
  int nmatch=0;

  for(i=0;i<len;i++){
    //adds 1 if the bits are different, 0 if they are the same
    nmatch+= abs( *(this+i) - *(that+i) );
  }
  return nmatch;
}


/*Probability of complex formation: 
  returns the probability associated to the least energetic 
  continuous match between two given dangling ends.
  The score is additive relative to single nucleotide matches:
  G-C = 3, A-U = 2, G-U = 1
*/

/*For now this is exactly Nobuto's function
  Corrected for 5'->3' direction, however, MaxMatch stil gets it wrong!!!
*/
void ProbComplForm(char *end5, char *end3,double *k_asso, double *k_diss)
{

  int maxmatch,len3,i;
  char revend3[MAXLEN];
  
  len3=strlen(end3);
  for(i=0;i<len3;i++) revend3[i]=end3[len3-i-1];
  revend3[len3]='\0';
  
  maxmatch = MaxMatch(end5, revend3 );
  

  *k_asso = 1. - exp(par_factor_asso*(double)maxmatch);
  *k_diss = exp(par_factor_diss*(double)maxmatch);
  
}


/* This will give back the maximam number of base pairs between
   "end3" and "end5". We count GC as 3, AU as 2, GU as 1. */
int MaxMatch(char *end3, char *end5)
{
  int i,j,len3,len5,len_diff;
  int nmatch; /* current */
  int pnmatch=0; /* previous */
  char *end_last,*end_middle; /* dont be confused by the name. */

  len3=strlen(end3);
  len5=strlen(end5);

  if(len3==0 || len5==0)
    return 0;

  if(len3>=len5) {
    end_last = end5 + len5;			//end_last points to 'one after the end' of end5
    len_diff = len3 - len5;
    end_middle = end3 + len_diff;	//end_middle points to the point of end3 which is at the length of the difference in length between end5 and end3
    /*
      Let the length of 3' be 'll' and 5' be 'dd'.
      Let l = ll-1 and d = dd-1.
      "end_last" is pointing to "5[dd]".
      (i=1)
      3'       012........l   -> 3[0] is the begining
      5' 012...d              -> 5[dd-i] is the begining

      (i=dd)
      3' 0123........l        -> 3[0]
      5' 0123.....d           -> 5[dd-dd]
    */
    
    /*
    starts from end3
    
    */
    for(i=1,j=len5;i<=j;i++){
      nmatch = nbasepair(end3,end_last-i,i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    /*
      (i=1)
      3' 0123.........l     -> 3[i]
      5'  0123....d         -> 5[0]

      (i=ll-dd)
      3' 0123.........l     -> 3[(ll-dd)]
      5'      0123....d     -> 5[0]
    */
    for(i=1,j=len_diff;i<=j;i++){
      nmatch = nbasepair(end3+i,end5,len5);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    /*
      Let "end_middle" point to "3[(ll-dd)]".
      (i=1)
      3' 0123.........l      -> 3[(ll-dd)+i]
      5'       0123....d     -> 5[0]

      (i=dd-1)
      3' 0123........l            -> 3[(ll-dd)+dd-1] = 3[l]
      5'             0123.....d   -> 5[0]
    */
    for(i=1,j=len5-1;i<=j;i++){
      nmatch = nbasepair(end_middle+i,end5,len5-i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    return pnmatch;
  }
  else {
    len_diff = len5 - len3;
    end_last = end3 + len3;
    end_middle = end5 + len_diff;
    for(i=1,j=len3;i<=j;i++){
      nmatch = nbasepair(end5,end_last-i,i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    for(i=1,j=len_diff;i<=j;i++){
      nmatch = nbasepair(end5+i,end3,len3);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
    for(i=1,j=len3-1;i<=j;i++){
      nmatch = nbasepair(end_middle+i,end3,len3-i);
      if(pnmatch<nmatch)
	pnmatch = nmatch;
    }
  }
  return pnmatch;
}


/* Sum over the base pair between "this" and "that". Don't forget
   the length of the sequence we want to look. */
int nbasepair(char* this,char* that,int len)
{
  int i;
  int nmatch=0;

  for(i=0;i<len;i++){
    nmatch+=ifbasepair(*(this+i),*(that+i));
  }

  return nmatch;
}
/* If bases can make a pair, return 1(GU) or 2(AU) or 3(GC),
   other wise 0 */
int ifbasepair(char this,char that)
{
  switch (this) {
  case 'A': {
    switch(that) {
    case 'A': return 0;
    case 'C': return 0;
    case 'G': return 0;
    case 'U': return 2;
    default : return 0;
    }
  }
  case 'C': {
    switch(that) {
    case 'A': return 0;
    case 'C': return 0;
    case 'G': return 3;
    case 'U': return 0;
    default : return 0;
    }
  }
  case 'G': {
    switch(that) {
    case 'A': return 0;
    case 'C': return 3;
    case 'G': return 0;
    case 'U': return 1;
    default : return 0;
    }
  }
  case 'U': {
    switch(that) {
    case 'A': return 2;
    case 'C': return 0;
    case 'G': return 1;
    case 'U': return 0;
    default : return 0;
    }
  }
  default:
    return 0;
  }

}


