#include <stdio.h>
#include <stdlib.h>		//remove me after testing
#include <string.h>
#include "cash2.h"
#include "replication.h"
#include "mersenne.h"


//Replication function, takes into account mutations
void Replication(TYPE2 *child,TYPE2 *parent,TYPE2 *cat){
  //first perfect replicate parent to child;
  PerfectReplication(child,parent);
  //then mutations happen
  Mutations(child,cat);
}

/*Replication with no errors, from source to target*/
void PerfectReplication(TYPE2 *child,TYPE2 *parent)
{
  *child=*parent;	//makes exact copy from child to parent, you should be careful with variable data!
  
  child->val=parent->val;
  
  /*
  child->kcat=parent->ckcat;
  child->kmut=parent->ckmut;
  child->kdelta=parent->ckdelta;
  
  child->ckcat=parent->kcat;
  child->ckmut=parent->kmut;
  child->ckdelta=parent->kdelta;
  
  //SeqRepl(child->end5,parent->end3);
  //SeqRepl(child->end3,parent->end5);
  
  IntSeqRepl(child->seq,parent->seq);	//makes complementary of array of integers
  
  if(parent->sign=='+') child->sign='-';
  else child->sign='+';  
  
  //If there, copy ancestor trace
  //for(j=0; j<50;j++) icel->anc[j]=parent->anc[j];
  */
}

/*Mutational operators apply here: so far 3 INDEPENDENT mutational events:
  1)-2) Point mutations to 5' and 3' dangling ends: binomially distributed
  3) Continuous mutations on kcat: continuously within [+delta/2,-delta/2]
*/
void Mutations(TYPE2 *child,TYPE2 *cat)
{
  //Mutate sequences
  //IntSeqSubstitutions(child->seq, cat->kmut);
  //SeqSubstitutions(child->end5, cat->kmut);
  //SeqSubstitutions(child->end3, cat->kmut);
  
  
  //Mutate krec
  if(genrand_real1() <  cat->kmut){
    child->krec += (cat->kdelta)*( genrand_real1()-0.5 );
    
    //krec for cat cannot be higher than 1, for par it should be
    //if(child->val==1){
    if(child->krec>2.) child->krec =  4. - child->krec;
    //}
    if(child->krec<0.) child->krec = -1. * child->krec;
  }
  
  
  
  /*
  //and ckcat
  if(genrand_real1() <  cat->kmut){
    child->ckcat += (cat->kdelta)*( genrand_real1()-0.5 );
  }
  if(child->ckcat>1.) child->ckcat =  2. - child->ckcat;
  if(child->ckcat<0.) child->ckcat = -1. * child->ckcat;
  */
  /*
  //Mutate kcat
  if(genrand_real1() <  cat->kmut){
    child->kcat += (cat->kdelta)*( genrand_real1()-0.5 );
  }
  if(child->kcat>1.) child->kcat =  2. - child->kcat;
  if(child->kcat<0.) child->kcat = -1. * child->kcat;
  //and ckcat
  if(genrand_real1() <  cat->kmut){
    child->ckcat += (cat->kdelta)*( genrand_real1()-0.5 );
  }
  if(child->ckcat>1.) child->ckcat =  2. - child->ckcat;
  if(child->ckcat<0.) child->ckcat = -1. * child->ckcat;
  */
}

//makes substitutions out of a seq of integers (bits, really)
void IntSeqSubstitutions(int *seq, double kmut)
{
  int i,nmut,childlen,index,end;
  int site[MAXLEN],positions[MAXLEN];		//site is for indexes, position contains the sites (on the sequence) where mutations occur
  int original, mutated;
  
  childlen = MAXLEN;
  nmut = (int)bnldev(kmut,childlen);	//returns an integer in the form of a double, telling how many mutations occur
  if(nmut!=0){
    /* Before choosing where mutations occur, initialize
     'rna_site'. We use rna_site to choose where mutations occur. */
    for(i=0;i<childlen;i++) site[i]=i;
    end=childlen;	//end is going to be decremented for each mutation
    for(i=0;i<nmut;i++){
      //find position: pos is int [0,end-1]
      index=(int)(  ((double)end) * genrand_real2()  );
      
      //we put the mutation in positions[]
      positions[i] = site[index];
      //we copy the last element of site into where the mutation has just happenend
      site[index]=site[end-1];
      end-- ;	//decrease end for next index generation
    }
    //now position contains the indexes of where mutations are
    for(i=0;i<nmut;i++){
      original=seq[positions[i]];
      mutated=1 - original;
      seq[positions[i]]=mutated;
    }
  }
}

void SeqSubstitutions(char *seq, double kmut)
{
  
  int i,nmut,childlen,index,end;
  int site[MAXLEN],positions[MAXLEN];		//site is for indexes, position contains the sites (on the sequence) where mutations occur
  char original, mutated;
  
  childlen = strlen(seq);
  nmut = (int)bnldev(kmut,childlen);	//returns an integer in the form of a double, telling how many mutations occur
  if(nmut!=0){
    /* Before choosing where mutations occur, initialize
     'rna_site'. We use this to chose where mutations occur. */
    for(i=0;i<childlen;i++) site[i]=i;
    end=childlen;	//end is going to be decremented for each mutation
    for(i=0;i<nmut;i++){
      //find position: pos is int [0,end-1]
      index=(int)(  ((double)end) * genrand_real2()  );
      
      //we put the mutation in positions[]
      positions[i] = site[index];
      //we copy the last element of site into where the mutation has just happenend
      site[index]=site[end-1];
      end-- ;	//decrease end for next index generation
    }
    //now position contains the indexes of where mutations are
    for(i=0;i<nmut;i++){
      original=seq[positions[i]];
      mutated=NucleotideMutation(original);
      seq[positions[i]]=mutated;
    }
  }
}

char NucleotideMutation(char nt)
{
  char A[3]="UGC";
  char U[3]="GCA";
  char G[3]="CAU";
  char C[3]="AUG";
  int rn;
  
  rn=(int)(3*genrand_real2());	//rn is now integer in [0,2]
  switch(nt){
    case 'A':
      return A[rn];
    case 'U':
      return U[rn];
    case 'G':
      return G[rn];
    case 'C':
      return C[rn];
    default:
      fprintf(stderr,"NucleotideMutation(): Error. I got a non RNA character\n");
      return '\0';
  }
  
}


//makes complementary of integer string
void IntSeqRepl(int *trgt, int *src)
{
  int seqlen,i;
  
  seqlen = MAXLEN;
  for(i=0;i<seqlen;i++) 
    trgt[i]= 1 - src[seqlen-1- i];
}

void SeqRepl(char *trgt, char *src)
{
  int seqlen,i;
  
  seqlen=strlen(src);
  for(i=0;i<seqlen;i++)
    trgt[i]=RNAc2cc(src[seqlen-1-i]);
}

char RNAc2cc(char rna)
{
  switch (rna) {
  case 'A':
    return 'U';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'U':
    return 'A';
  default:
    fprintf(stderr,"RNAcc2c: I got a non RNA character\n");
    return '\0';
  }
}
