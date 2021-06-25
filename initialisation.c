#include <stdio.h>
#include <string.h>

#include "cash.h"
#include "cash2.h"
#include "constants.h"
#include "parameters.h"
#include "replication.h"
#include "mersenne.h"
#include "initialisation.h"

#define BOUNDARY -1
#define BG 0	/*Background val*/
#define CAT 1	/*Catalyst val, notice that there can be cat that behave like parasites*/
#define PAR 2	/*Parasites val*/

extern int nrow;
extern int ncol;

int InitialPlane(TYPE2 **world, TYPE2 empty)
{
  int control;
  
  if(strlen(par_input_file_name)!=0){
    control=InitialiseFromBackup(world, empty);
  }else
    control=InitialiseFromScratch(world, empty);
  
  return control;
  
}

int InitialiseFromScratch(TYPE2 **world, TYPE2 empty)
{
  int nr,nc,i,j;
  TYPE2 init,c_init;
  
  nr=nrow; nc=ncol;
  
  //Initialise the first TYPE2;
  init=empty;	//this takes care especially of variable data (bonr, etc...)
  
  init.val = init_val;
  memcpy( init.seq, init_seq, MAXLEN*sizeof(int) );
  init.sign=init_sign;
  init.krec=init_krec;
  //fprintf(stderr,"Warning: Random initialisation for krec, about 0.5\n");
  init.kcat=init_kcat;
  init.kmut=init_kmut;
  init.kdelta=init_kdelta;
  init.ktime=init_ktime;
  /*
  init.ckrec=init_ckrec;
  init.ckcat=init_ckcat;
  init.ckmut=init_ckmut;
  init.ckdelta=init_ckdelta;
  */
  PerfectReplication(&c_init,&init);
  //printf("%d %s %s %f %f %f %d %d\n", c_init.val, c_init.end5, c_init.end3, c_init.kcat, c_init.kmut, c_init.kdelta, c_init.bonr, c_init.bonc);
   
  
  for(i=1;i<=nr;i++)for(j=1;j<=nc;j++){
    world[i][j]=empty;
    //initialise a square 50x50
    //if(i>nrow/2 -25 && i<nrow/2 +25 && j>ncol/2 -25 && j<nrow/2 +25 ){
      if(genrand_real1()<.1){
        if(genrand_real1()<0.5) world[i][j]=init;
        else world[i][j]=c_init;
	if(genrand_real1()<0.1){
	  world[i][j].val=PAR;
	  world[i][j].krec=par_beta;
	  //world[i][j].krec=init_beta;
	}
	//world[i][j].krec=0.5 + (genrand_real2()-0.5)/10.;
      }
    //}
  }
  
  /*
  //intialise 1,2,3 cell (or cell pairs) only to check whether diffusion is correct
  world[50][50].val=1;
  */
  /*
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    world[i][j]=empty;
  }
  
  world[50][50]=init;
  world[50][50].bonr=51;
  world[50][50].bonc=51;
  world[50][50].term=5;
  
  world[51][51]=c_init;
  world[51][51].bonr=50;
  world[51][51].bonc=50;
  world[51][51].term=3;
  */
  /*
  world[50][50].val=1;
  world[50][50].bonr=51;
  world[50][50].bonc=51;
  world[50][50].term=5;
  
  world[51][51].val=2;
  world[51][51].bonr=50;
  world[51][51].bonc=50;
  world[51][51].term=3;
  
  world[50][51].val=1;
  world[50][51].bonr=51;
  world[50][51].bonc=50;
  world[51][50].val=3;
  world[51][50].bonr=50;
  world[51][50].bonc=51;
  
  world[100][100].val=1;
  world[100][100].bonr=101;
  world[100][100].bonc=101;
  world[101][101].val=4;
  world[101][101].bonr=100;
  world[101][101].bonc=100;
  */
  
  return 0;
}


/* In backup reading:
 * it can happen that someone from the backup file we are using as input
 * is the edge of the plane and is in complex with someone 
 * on the "other side" (because plane is wrapped, no?). 
 * It can also happen that the new field is larger (or smaller) 
 * than the one we had, so we have to deal them being in complex carefully 
 */
int InitialiseFromBackup(TYPE2 **world, TYPE2 empty)
{
  int nr,nc,fr,fc,i,j,k,control=0;
  char instr[BUFSIZE],anc[PAR_ANCESTOR];
  char *sign;
  //stuff to assign to TYPE:
  
  FILE *fin;
  
  fin=fopen(par_input_file_name,"r");
  if(fin==NULL) {
    fprintf(stderr, "InitialiseFromBackup(): Error. Can't open %s\n",par_input_file_name);
    exit(1);
  }
  
  //read heading of backup file: "#Time nr X nc"
  fgets(instr,BUFSIZE,fin);
  strtok(instr," ");	//strip #Time
  fr = atoi(strtok(NULL," ")); //read # of row
  strtok(NULL," "); // strip X
  fc = atoi(strtok(NULL," \n")); //read # of col
  
  fprintf(stderr,"Initialise from backup file. #row: %d, #col: %d\n",fr,fc);
  
  if(nrow!=fr || ncol!=fc) control = 1;
  
  nr = (nrow>fr) ? nrow : fr;
  nc = (ncol>fc) ? ncol : fc;
  for(i=1;i<=nr;i++){
    for(j=1;j<=nc;j++){
      if(i<=fr && j<=fc){
	    world[i][j]=empty;
	    
	    //read a line  
	    if(fgets(instr,BUFSIZE,fin)==NULL){
	      fprintf(stderr, "InitialiseFromBackup(): Error at reading file\n");
	      exit(1);
	    }
	    
	    sign=strtok(instr," ");
	    if(sign[0]=='N') continue;
	    
	    world[i][j].val = atoi(sign);
	    sign=strtok(NULL," ");		//for now, all sign is +
	    world[i][j].sign = sign[0];
	    //world[i][j].bonr = atoi(strtok(NULL," "));	//This will come when I'll also save info on time span of replication
	    //world[i][j].bonc = atoi(strtok(NULL," "));
	    //world[i][j].term = atoi(strtok(NULL," "));
	    //world[i][j].kdiss= atof(strtok(NULL," "));
	    
	    (strtok(NULL," "));
	    (strtok(NULL," "));
	    (strtok(NULL," "));
	    (strtok(NULL," "));
	    //(strtok(NULL," "));
	    
	    world[i][j].krec = atof(strtok(NULL," "));
	    //world[i][j].krec = init_krec;			//We overwrite our own krec //we do not
	    world[i][j].kcat = atof(strtok(NULL," "));
	    world[i][j].kmut = atof(strtok(NULL," "));
	    world[i][j].kdelta = atof(strtok(NULL," "));
	    
	    /*
	    world[i][j].ckrec = atof(strtok(NULL," "));
	    world[i][j].ckcat = atof(strtok(NULL," "));
	    world[i][j].ckmut = atof(strtok(NULL," "));
	    world[i][j].ckdelta = atof(strtok(NULL," "));
	    */
	    
	    //sign=strtok(NULL," ");
	    //for(k=0;sign[k]!='\0';k++) world[i][j].seq[k]=sign[k]; //sequence (when needed)
	    //world[i][j].seq[k]='\0';
	    world[i][j].seq[0]=0;
	    //printf("string ---%s---\n", sign);
	    
	    
	    world[i][j].ktime=init_ktime;
	    
	    //length = strlen(sign);
	    //strncpy(world[i][j].seq,sign,length);
	    	    
	    //after sequence comes \n sign <- dealt with by default
	    //ancestor trace (when needed, remember that separator is -)
	    
	    // --- IMPORTANT --- //
	    //Here we check that bonr and bonc are less than frow and fcol
	    //and also that if you get that bonr and bonc are on the other
	    //side of the lattice, that nrow=frow and ncol=fcol
	    //otherwise just destroy bonr,bonc,kdiss values (to -1, -1, 0.)
	    if(world[i][j].bonr>fr || world[i][j].bonc>fc) {
	      world[i][j].bonr=-1;
	      world[i][j].bonc=-1;
	      world[i][j].kdiss=0.;
	    }
	    if(nrow!=fr || ncol!=fc)
	      if( abs(i-world[i][j].bonr)>1 || abs(j-world[i][j].bonc)>1 ){
	        world[i][j].bonr=-1;
	        world[i][j].bonc=-1;
	        world[i][j].kdiss=0.;
	      }
	    
	    //printf("i:%d j:%d val:%d krec:%f\n",i,j,world[i][j].val,world[i][j].krec);
	    
	    
      }
    }
  }
  fclose(fin);
  return control;
}




