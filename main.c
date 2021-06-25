/* Version 3.1
  In this version replication takes time (trepl). 
  
  We treat time to the end of replication as first order reaction for a complex
  In this way we can scale 
  
  so if replication 
  
  
  all this down here is bullshit
  Notice that we measure replication time in ATU's, however 1 ATU = (max_number_events/global_scale) 
  time steps so if trepl=4.3 ATU, this means we wait 4.3 * 4/0.5 = 34.4 time steps. 
  After 34 timesteps (decreasing the counter), we put the counter to zero
  with probability = 1 - 0.3, else we wait one more turn.
  
  Pseudocode:
  - If complex: 
      if time of replication is not finished increase time counter
      (during this time "falling off" can happen, as in complex dissociation)
    if empty space:
      if nei is complex and its replication time is passed:
        replication
  - Replication time is evolvable, but (eventually) the longer it takes, the lower mutation rate is      
    
    
    
*/

/* Version 1.3
  beta - the parasite advantage - can also evolve freely
*/

/*
  Introduced parameter for global probability scaling
  Technicals: 
 - Backup file reading fixed
 - Backup file writing adds who is at 3' and 5' of a complex, heading includes which parameters are saved
*/

/*
  Parasites are introduced, they have val=2
  - parasites can form complexes only in one way, i.e. they can only be at 5
  - All the times we are asking if the state of a grid point is non empty
    we ask if icel->val > 0: 
      0 is background, 
      -1 is boundary,
      1 is cat
      2 is par
  - technicals: in color.c, function added that prints colours to stderr
  - color map is fixed, now it yields colours as they should be.
  - backup maker and reader introduced: backup now also says with whom you're making complex, if anyone.
  -
  
*/
/*Version 1.0
 - Nobuto's RNA model
 
 Technical:
 - The probabilities of each event are normalised
 - Possibly some optimisation will be implemented as well... 
*/

/*Version 0.3: BACKBONE-MINIMAL MODEL:
 - one bit string (int seq) as both template and recognition seq 
 - replication makes complementary
 - kcat and ckcat
 
 Technical:
 - modify scaling parameters for complex formation probabilities to sensible values
 - save file dataMut.txt with sequences in random order
 - save ALL the information
 - make some video colouring of the tails 
   (maybe based on nucleotide homogeneity and preference, in this case it would 1 or 0)
*/

/*The program simulates the evolution of interacting particles on a grid.
  Essentially, this is a spatially extended, discrete, interaction-reaction-diffusion system
*/

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "cash.h"
#include "cash2.h"
#include "my.h"
#include "parameters.h"
#include "initialisation.h"
#include "replication.h"
#include "mersenne.h"
#include "graphical_output.h"
#include "complex_formation.h"

#define BOUNDARY -1
#define BG 0	/*Background val*/
#define CAT 1	/*Catalyst val, notice that there can be cat that behave like parasites*/
#define PAR 2	/*Parasites val*/

/*
TYPE2 global_empty={0,"N","N","N",0.,0.,0.,	//permanent data: val,seq,end5,end3,kcat,kmut,kdelta
                    0.,0.,0.,				//values of the complementary
                    '\0',					//sign/
		            -1,-1,0,0.};	//variable data: bonr,bonc,term,kdiss
*/
TYPE2 global_empty={
  0,{-1},0., 0.,0.,0.,	//permanent data: val,seq,krec,kcat,kmut,kdelta
  //0.,0.,0.,0.,		//complementary information
  '\0',-1.,			//sign, ktime
  -1,-1,0,0.,-1.,-1		//variable data: bonr,bonc,term,kdiss,trepl,to_be_zero
};

int main(int argc, char **argv)
{
  TYPE2 **world;
  struct updorder *updorder_p;
  struct point **neighbor[9];
  //  MARGOLUS ***margolus;
  int Time, maxtime;
  int stop_signal;
  int i;
  
  nrow = 2048; /* This does not include boundary */
  ncol = 2048;
  scale= 1;

  boundary=WRAP;	//FIXED
  boundaryvalue=0;
  boundaryvalue2.val=BOUNDARY;
  boundaryvalue2.seq[0]='B'; /* Boundary. Do not change to N! */
  boundaryvalue2.bonr = -1;
  
  /* Allocating memory to plane*/
  world = New2();
  /* Arguments */
  //Args(argc,argv);

  /* Initialize 1) planes 
  				2) boundary conditions 
  				3) random number generator 
  				4) output file(s) 
  				5) color map 
  				6) the window. */
  Initial(world);
  
  /*Maybe you'll have to initialise the "observer"
    as Nobuto likes to call it.
  */
  
  /* Set neighbor plane. "boundary" and boundary value variables has to be
     set before. */
  for(i=0;i<9;i++)
    neighbor[i]=NeighSet(i);
  
  fprintf(stderr,"Simulation start, have fun!\n");
  fprintf(stderr,"Warning: so far we mutate nothing!\n");
  /******************/
  /*** MAIN CYCLE ***/
  /******************/
  for (Time = par_inittime,maxtime=par_maxtime; Time <= maxtime; Time++) {
  
    //if(Time%1000==0)Meteorites(world);
    
    /*Save data if time*/
    if(Time%par_time_savemovie==0) SaveMovie(world,Time);
    if(Time%par_time_savedata==0){
      UpdOrdReset(&updorder_p);
      stop_signal=SaveDataRandomOrder(world,Time,updorder_p);	//saves data in random order (useful for random samples)
    }
    if(stop_signal==0) {
      fprintf(stderr,"Simulation finished due to extinction\n");
      break;
    }
    
    if(Time%par_time_backup==0 && Time!=0) OutputBackup(world,Time);
    
    /*All the dynamics passes through this function*/
    UpdOrdReset(&updorder_p);
    AllThatCanHappen(world, updorder_p, neighbor);
    
  }
  
  //Exit everything
  Finish(world,neighbor,updorder_p);
  fprintf(stderr, "Simulation run finished. Everything properly quit.\nYou won't see this message very often, you know?\n");
  return 0;
}


void Meteorites(TYPE2 **world)
{
  int i,j,k;
  double metr,metc;
  int n_meteorites;
  TYPE2 empty = global_empty;
  TYPE2 *icel;
  
  //n_meteorites=(int)(10*genrand_real2());
  n_meteorites=10;
  for(k=0;k<n_meteorites;k++){
    metr=1 + (int)(nrow*genrand_real2());
    metc=1 + (int)(ncol*genrand_real2());
    
    for(i=metr-49;i<metr+50;i++)for(j=metc-49;j<metc+50;j++){
      icel = &world[1 + ((nrow+i)%nrow)][1 + ((ncol+j)%ncol)];
      if(icel->bonr!=-1) world[icel->bonr][icel->bonc].bonr=-1;
      *icel=empty;
    }
  }
    
}


/******************************/
/*** --- INITIALISATION --- ***/
/******************************/
void Initial(TYPE2 **world)
{
  //int i,j;
  int control;
  TYPE2 empty = global_empty;
  
  Boundaries2(world);		//initialise boundaries
  init_genrand(par_ulseed);	//initialise random number generator with seed par_ulseed (see parameters.c)
  InitColormap();			//go make a colormap ()
  OpenPNG(par_movie_directory_name, nrow, par_nplane*ncol);
  
  //ColorPrint();		//print colormap
  
  control=InitialPlane(world, empty);
  switch(control){
    case -1: fprintf(stderr,"InitialPlane() failed to initialise the plane properly.\n");
             fprintf(stderr,"The program will exit badly, now.\n");
	     exit(1);
    case 0 : fprintf(stderr,"InitialPlane(): Initialisation successful.\n");
             break;
    case 1 : fprintf(stderr,"InitialPlane(): Warning. Nothing failed, but something is strange.\n");
             fprintf(stderr,"Did you change plane size? Is it correct? Check things out.\n");
	     break;
    default: fprintf(stderr,"Status from InitialPlane() not recognised\n");
             exit(1);
  }
  
  //for(i=1;i<=nrow;i++) for(j=1;j<=ncol;j++) if((world[i][j].ktime!=init_ktime && world[i][j].ktime!=-1.) || world[i][j].trepl !=-1.) printf("%f %f\n",world[i][j].ktime,world[i][j].trepl);
  //exit(1);
}

/*
 What can happen in the model? For now:
  icel is empty, inei is empty 		-> nothing
		 inei is free mol 	-> swapping (2nd ord)
		 inei is complex 	-> swapping (3rd ord), replication
  icel is free mol, inei is empty 	-> death icel, swapping (2nd ord)
		    inei is free mol	-> death icel, swapping (2nd ord), compl form i->n, compl form n->i.
		    inei is compelx	-> death icel, swapping (3rd)
  icel is complex, inei is empty	-> death icel, swapping (3rd), compl dissociation
		   inei is free mol	-> death icel, swapping (3rd), compl dissociation
		   inei is in complex 	-> death icel, swapping (4th ord), swapping (2nd) if with eachother
  The highest possible number of reactions is the case icel free, inei free.
  This means that what_happens={0,1,2,3} and 1 time step (aut) = 4 updates of the variable Time
*/

void AllThatCanHappen(TYPE2 **world, struct updorder* updorder_p,struct point** neighbor[])
{
  int i,length,irow,icol,idir,ineirow,ineicol;
  TYPE2 *icel,*inei;
  int what_happens;
  //events that never happen at the same time can be tagged by the same number
  int event_death=0, event_repl=0; 	
  int event_diff2=1, event_diff3=1, event_diff4=1;
  int event_in_compl=2, event_compl_diss=2; 
  int event_ni_compl=3, event_time_pass=3;
  double p_ithappens;
  double max_number_events=4.; //double to avoid conversions when used
  
  length = nrow*ncol;
  for(i=0;i<length;i++){
    irow=(updorder_p+i)->row;
    icol=(updorder_p+i)->col;
    
    icel=&world[irow][icol];
    
    idir = 1 + (int)(8.*genrand_real3());
    ineirow=neighbor[idir][irow][icol].row;
    ineicol=neighbor[idir][irow][icol].col;
    inei=&world[ineirow][ineicol];
    if(inei->val==BOUNDARY) {
      fprintf(stderr,"Test: event {inei->val = -1} i.e. neighbour, does happen.\n Delete this message from code\n");
      continue;
    }
    
    what_happens=(int)(max_number_events*genrand_real2());
    p_ithappens=genrand_real1()/global_scale;
    
    /////////////////////////////
    //  --- ICEL IS EMPTY ---  //
    /////////////////////////////
    if(icel->val==BG){
      //inei is empty -> nothing happens
      if(inei->val==BG) continue;
      //if inei is not empty
      else if(inei->val>BG){
	//if inei not in complex
	if(inei->bonr==-1){
	  //all that can happen is second order diffusion
	  if(what_happens==event_diff2 && p_ithappens<par_kdiff) SwapCells(2,icel,irow,icol, inei,ineirow,ineicol);
	}else if(inei->bonr!=-1){
	  //replication and 3rd order diffusion can happen
	  if(what_happens==event_repl) ToReplication(icel,inei,&world[inei->bonr][inei->bonc]);
	  else if(what_happens==event_diff3 && p_ithappens<par_kdiff) SwapCells(3,icel,irow,icol, inei,ineirow,ineicol, &world[inei->bonr][inei->bonc]);
	}
      }
    }else if(icel->val>BG){
      /////////////////////////////////////
      //  --- ICEL IS FREE MOLECULE ---  //
      /////////////////////////////////////      
      if(icel->bonr==-1){
	//if inei is empty
	if(inei->val==0){
	  //death, diffusion (2nd ord)
	  if(what_happens==event_death && p_ithappens<par_kdec) MolDecay(1,icel);
	  else if(what_happens==event_diff2 && p_ithappens<par_kdiff) SwapCells(2,icel,irow,icol, inei,ineirow,ineicol);
	}
	//if inei is not empty
	else if(inei->val>BG){
	  //if inei is a free molecule
	  if(inei->bonr==-1){
	    //icel and inei can be cat or par so
	    //icel==cat, inei==cat: death of icel, diffusion (2nd order), icel->inei compl, inei->icel compl
	    //		 inei==par: death of icel, diffusion (2nd order), icel->inei compl
	    //icel==par, inei==cat: death of icel, diffusion (2nd order), inei->icel compl
	    //		 inei==par: death of icel, diffusion (2nd order)
	    if(icel->val==CAT && inei->val==CAT){
	      if(what_happens==event_death && p_ithappens<par_kdec) MolDecay(1,icel);
	      else if(what_happens==event_diff2 && p_ithappens<par_kdiff) SwapCells(2,icel,irow,icol, inei,ineirow,ineicol);
	      else if(what_happens==event_in_compl && p_ithappens<icel->krec) Complex_Formation(par_kdiss,icel,irow,icol,inei,ineirow,ineicol);
	      else if(what_happens==event_ni_compl && p_ithappens<inei->krec) Complex_Formation(par_kdiss,inei,ineirow,ineicol,icel,irow,icol);
	    }else if(icel->val==CAT && inei->val==PAR){
	      if(what_happens==event_death && p_ithappens<par_kdec) MolDecay(1,icel);
	      else if(what_happens==event_diff2 && p_ithappens<par_kdiff) SwapCells(2,icel,irow,icol, inei,ineirow,ineicol);
	      else if(what_happens==event_in_compl && p_ithappens<inei->krec*icel->krec) Complex_Formation(par_kdiss,icel,irow,icol,inei,ineirow,ineicol);
	    }else if(icel->val==PAR && inei->val==CAT){
	      if(what_happens==event_death && p_ithappens<par_kdec) MolDecay(1,icel);
	      else if(what_happens==event_diff2 && p_ithappens<par_kdiff) SwapCells(2,icel,irow,icol, inei,ineirow,ineicol);
	      else if(what_happens==event_ni_compl && p_ithappens<icel->krec*inei->krec) Complex_Formation(par_kdiss,inei,ineirow,ineicol,icel,irow,icol);
	    }else if(icel->val==PAR && inei->val==PAR){
	      if(what_happens==event_death && p_ithappens<par_kdec) MolDecay(1,icel);
	      else if(what_happens==event_diff2 && p_ithappens<par_kdiff) SwapCells(2,icel,irow,icol, inei,ineirow,ineicol);
	    }
	  }else if(inei->bonr!=-1){
	    //diffusion (3rd order), death of icel
	    if(what_happens==event_diff3 && p_ithappens<par_kdiff) SwapCells(3,icel,irow,icol, inei,ineirow,ineicol, &world[inei->bonr][inei->bonc]);
	    else if(what_happens==event_death && p_ithappens<par_kdec) MolDecay(1,icel);
	  }
	}
      }
      //////////////////////////////////
      //  --- ICEL IS IN COMPLEX ---  //
      //////////////////////////////////
      else if(icel->bonr!=-1){
	if(what_happens==event_death && p_ithappens<par_kdec) MolDecay(2,icel, &world[icel->bonr][icel->bonc]); //degradation is always first order reaction
	
	//Here time can pass: first order reaction (depends only on cat)
	//so if icel is at the 3' of the complex, if event_time_pass
	//  if time left >1.: time --
	//  if time left >0.: 
	//    if it's first occasion (to_be_zero==1): chances of putting time left to zero = 1. - (time left)
	//    if second occasion: time left = 0.
	//  if time == 0.: nothing happens
	//  if time left<0.: we have a problem because time = -1. is the flag we use for molecules not in complexes...
	else if(what_happens==event_time_pass && icel->term==3){
	  
	  
	  //printf("\nTime left: %f\n",icel->trepl );
	  //if(icel->trepl==0.) printf("%d %d %d %f %f %d\n", irow, icol, icel->val, icel->krec, icel->kdiss, world[icel->bonr][icel->bonc].val);
	  
	  
	  if( icel->trepl >= 1. && p_ithappens<1.) icel->trepl -=1.;
	  else if(icel->trepl < 1. && icel->trepl > 0.){
	    if(icel->to_be_zero==1){
	      if ( p_ithappens<(1. - icel->trepl) ) icel->trepl=0.;
	      else icel->to_be_zero=2;
	    }else if( icel->to_be_zero==2 ) icel->trepl=0.;
	  }
	  //else if time left for replication is 0, nothing happens
	  else if( icel->trepl < 0. ) fprintf(stderr,"AllThatCanHappen(): Warning, trepl<0. for molecule in complex. Troubles ahead!\n");
	}else if(inei->val==BG || (inei->val>BG && inei->bonr==-1) ){
	  //complex can dissociate (2nd ord reaction), diffusion (3rd order, notice the order of the cel passed to the funct), death of icel (icel is in complex)
	  if(what_happens==event_diff3 && p_ithappens<par_kdiff) SwapCells(3,inei,ineirow,ineicol, icel,irow,icol, &world[icel->bonr][icel->bonc]);
	  else if(what_happens==event_compl_diss && p_ithappens<icel->kdiss) Complex_Dissociation(icel, &world[icel->bonr][icel->bonc]);
	}else if(inei->val>BG && inei->bonr!=-1){
	  //complex dissociation for icel, diffusion (4th ord), death of icel (icel is in complex!!!)
	  if(what_happens==event_compl_diss && p_ithappens<icel->kdiss) Complex_Dissociation(icel, &world[icel->bonr][icel->bonc]);
	  else if(what_happens==event_diff4 && p_ithappens<par_kdiff){
	    //we have 2 options, if icel is in complex with inei
	    //or if each is in complex with other molecules (Notice ORDER !!!)
	    if(icel->bonr==ineirow && icel->bonc==ineicol) SwapCells(2,inei,ineirow,ineicol,icel,irow,icol);
	    else SwapCells(4,icel,irow,icol,&world[icel->bonr][icel->bonc],inei,ineirow,ineicol,&world[inei->bonr][inei->bonc]);
	  }
	}
      }
    }
  }//END OF LOOP
}

//////////////////////////////////
// ---  MOLECULAR DYNAMICS  --- //
//////////////////////////////////

//the probability of complex formation has been already calculated
//Here we only set complex formation
void Complex_Formation(double kdiss, TYPE2 *cat,int irow, int icol,TYPE2 *tmpl,int ineirow,int ineicol)
{
  cat->bonr = ineirow; cat->bonc = ineicol;
  tmpl->bonr = irow; tmpl->bonc = icol;

  cat->term = 3;
  tmpl->term = 5;
  cat->kdiss = tmpl->kdiss = kdiss; 
  cat->trepl=cat->ktime;
  cat->to_be_zero=1;
}

void Complex_Dissociation(TYPE2 *x,TYPE2 *xn)
{
  x->bonr=-1;
  xn->bonr=-1;
}

// ToReplication();	you have to find who is replicase and who is template
//			also, you have to scale this as a 2nd order reaction
//			icel is the empty spot, inei and neinei are in complex.
void ToReplication(TYPE2 *icel,TYPE2 *inei,TYPE2 *neinei)
{
  double scaling_factor=0.5;
  TYPE2 *parent,*replicase;
  
  //printf("Inside ToReplication function\n");
  
  if(genrand_real1()<scaling_factor){
    if(inei==NULL || neinei==NULL){
      fprintf(stderr, "Replication(): Error. Cell points to NULL\n");
      exit(1);
    }
    if(inei->term==3){
      replicase=inei;
      parent=neinei;
    }else if(neinei->term==3){
      replicase=neinei;
      parent=inei;
    }else{
      fprintf(stderr,"ToReplication(): Error. Molecules are in complex but 3 and 5 terms are undefined.\n");
      exit(1);
    }
    //time has been dealt with elsewhere
    if( replicase->trepl==0.){
      
      if(genrand_real1() < replicase->kcat){
        Replication(icel,parent,replicase);
        
        icel->bonr = -1; /* Don't forget this */
        /* The complex molecule breaks apart */
        replicase->bonr = parent->bonr = -1;
	icel->trepl = replicase->trepl = parent->trepl = -1.;
	
      }
    }
  }
}

// MolDecay(TYPE2 *icel, ...)
// takes int and 1 or 2 TYPE2* depending on icel being free or in complex
void MolDecay(int howmany, ...)
{
  va_list arguments;
  TYPE2 *x, *xn;
  TYPE2 empty=global_empty;
  
  va_start( arguments, howmany);	//initialize variable length argument list
  
  if(howmany==1){
    x=va_arg( arguments, TYPE2 * );
    //free molecule, it just decays
    *x=empty;
  }else if(howmany==2){
    x=va_arg( arguments, TYPE2 * );
    xn=va_arg( arguments, TYPE2 * );
    
    xn->bonr=-1;
    *x=empty;
  }
}


/*******************************
 *          Diffusion          *
 *******************************/

// SwapCells(int, ...)	takes an integer and a variable number of TYPE2 *, 
//			and performs swapping, with probability scaled by the react order
//			Notice that it is exactly the same as MyDiffusion except it
//			does not deal with all molecules at once
void SwapCells(int order, ...)
{
  va_list arguments;
  TYPE2 *x, *xn, *y, *yn;
  TYPE2 tmp;
  double scaling_factor;
  int xrow,xcol,yrow,ycol,tmpbonr,tmpbonc;
  
  va_start ( arguments, order);	//initialize variable length argument list
  switch (order){
    case 2: 
      //fprintf(stderr,"2nd order diff\n");
      x=va_arg( arguments, TYPE2 * );
      xrow=va_arg( arguments, int);
      xcol=va_arg( arguments, int);
      y=va_arg( arguments, TYPE2 * );
      yrow=va_arg( arguments, int);
      ycol=va_arg( arguments, int);
      break;
    case 3: 
      //fprintf(stderr,"3rd order diff\n");
      y=va_arg( arguments, TYPE2 * ); 
      yrow=va_arg( arguments, int);
      ycol=va_arg( arguments, int);
      x=va_arg( arguments, TYPE2 * ); 
      xrow=va_arg( arguments, int);
      xcol=va_arg( arguments, int);
      xn=va_arg( arguments, TYPE2 * ); 
      break;
    case 4: 
      //fprintf(stderr,"4th order diff\n");
      x=va_arg( arguments, TYPE2 * );
      xrow=va_arg( arguments, int);
      xcol=va_arg( arguments, int);
      xn=va_arg( arguments, TYPE2 * ); 
      y=va_arg( arguments, TYPE2 * ); 
      yrow=va_arg( arguments, int);
      ycol=va_arg( arguments, int);
      yn=va_arg( arguments, TYPE2 * ); 
      break;
    default:
      fprintf(stderr, "SwapCells(): Error. Reaction order parameter is strange: %d.\n", order);
      exit(1);
  }
  va_end ( arguments ); 	// Cleans up the list
  
  //Now we got all the arguments
  //depending on the order of the reaction we perform the swapping
  //Of course dividing by order here as well is pure overhead
  //but to me it's clearer what's happening
  if(order==2){
    scaling_factor=0.5;
    if(genrand_real1() < scaling_factor) {
      //fprintf(stderr," diff 2");
      //if the molecules are in complex with each other
      //we have first to swap the location they point to for complex
      if(x->bonr==yrow){
        x->bonr=xrow; x->bonc=xcol;
	y->bonr=yrow; y->bonc=ycol;
      }
      
      tmp=*x;
      *x=*y;
      *y=tmp;
    }
    return;  
  }else if(order==3){
    scaling_factor=1/6.;
    if(genrand_real1() < scaling_factor) {
      //fprintf(stderr," diff 3");
      x->bonr=xrow;     x->bonc=xcol;	//x in complex with xn, which will be where x is now
      xn->bonr=yrow; xn->bonc=ycol;	//xn in complex with x, which will be where y is now
      
      tmp=*y;		/*put y in tmp*/
      *y=*x;		/*x -> y*/
      *x=*xn;		/*x' -> x*/
      *xn=tmp;		/*y -> x'*/
    }
    return;
  }else if(order==4){
    scaling_factor=1/24.;
    if(genrand_real1() < scaling_factor){
      //fprintf(stderr," diff 4");
      tmpbonr=x->bonr; tmpbonc=x->bonc;
      x->bonr= yrow; x->bonc=ycol;		//x in complex with xn, which will be where y is now
      xn->bonr=y->bonr; xn->bonc=y->bonc;	//xn in complex with x, which will be where yn is now
          
      y->bonr=xrow; y->bonc=xcol;		//y in complex with yn, which will be where x is now
      yn->bonr=tmpbonr; yn->bonc=tmpbonc;	//yn in complex with y, which will be where xn is now
      
      tmp=*x;
      *x=*yn;	//yn -> x
      *yn=tmp;	//x -> yn
        
      tmp=*y;
      *y=*xn;	//xn -> y
      *xn=tmp;	//y -> xn
    }  
    return;
  }
  /*ELSE IF WE GOT HERE THERE IS SOMETHING WRONG BECAUSE WE SHOULD HAVE RUN OUT POSSIBILITIES*/
  else fprintf(stderr, "DiffusionComplex(). Warning, I got to a possibility I did not consider\n xval: %d, yval: %d, irow: %d, icol: %d, xbonr: %d, xbonc: %d, ineirow: %d, ineicol %d, ybonr: %d, ybonc: %d\n",x->val,y->val,xrow,xcol,x->bonr,x->bonc,yrow,ycol,y->bonr,y->bonc);
}

/****************************/
/*** --- FINALISATION --- ***/
/****************************/
void Finish(TYPE2 **world,struct point **neighbor[],struct updorder *updorder_p)
{
  int i;
  
  //RemoveEndFile();
  PlaneFree2(world);
  for(i=0;i<9;i++)
    NeighFree(neighbor[i]);
  free(updorder_p);
}




