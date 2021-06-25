#include <stdio.h>
#include <string.h>
#include "cash.h"
#include "cash2.h"
#include "mersenne.h"

/********************************************** basic */
extern FILE *parfile;

extern int seed;
extern int seedset;
extern int graphics;
extern int nrow;
extern int ncol;
extern int first;
extern int last;
extern int boundary;
extern int boundaryvalue;
extern int scale;
extern int ranmargolus;
extern int psborder;
extern int psreverse;


TYPE2 boundaryvalue2;

TYPE2 **NewP2(void)
{
  TYPE2 **a;
  a = (TYPE2 **)calloc(nrow+2, sizeof(TYPE2 *));
  if (a == NULL) {
    fprintf(stderr,"NewP2: error in memory allocation\n");
    exit(1);
  }
  return a;
}

TYPE2 **New2(void)
{
  TYPE2 **a;
  int i,j;
  a = NewP2(); 
  a[0] = (TYPE2 *)calloc((nrow+2)*(ncol+2),sizeof(TYPE2));
  if (a[0] == NULL) {
    fprintf(stderr,"New2:error in memory allocation\n");
    exit(1);
  }
  for (i=1,j=nrow+2; i < j; i++)
    a[i] = a[i-1] + ncol + 2;

  if (first == 0) {
    first = ncol + 3;
    last = (nrow+2)*(ncol+2) - ncol - 4;
  }
  return a;
}

int PlaneFree2(TYPE2 **a)
{
  free(a[0]);
  free(a);
  return 0;
}

/* UpdOrdReset (Update Order Reset) will reset the update order of the
   cells in asynchronous updating. It will take "struct updorder" array and
   will shuffle it randomly. User updates the cell along the order this
   "struct updorder" array specifies so that User can update every cell
   only once in one time step in randomely chosen order. 
   User has to declear a pointer variable to "struct updorder". This
   pointer variable will be made to a array of the cell coodinate which
   specifies the order of updating by UpdOrdReset(). To use UpdOrdReset(),
   give the address of User's pointer variable to "struct updorder". The
   function will asign address to the pointer variable.*/

void UpdOrdReset(struct updorder **updorder_pp)
{
  static int initial_flag=0; /* If 0, "updorder_p" is not initialized */
  int nc1,nc2,length,location;
  int i;
  int target;
  struct updorder temp;

  nc1 = ncol+1;
  nc2 = nc1 + 1;
  length = ncol*nrow; /* the length of the "*updorder_pp" array */

  /* Initilization of "*updorder_pp". */
  if(initial_flag==0)
    {
      *updorder_pp = (struct updorder *) malloc(length*sizeof(struct updorder));

      location = nc1 + 2; /* a currunt location of the head along the length
			   of array */
      for(i=0;i<length;i++)
	{
	  if((location%nc2)==nc1)
	    location+=2; /* Skip boder picxels. */

	  /* Assign a initial coorinate */
	  (*updorder_pp+i)->row=location/nc2;
	  (*updorder_pp+i)->col=location%nc2;

	  location++;
	}
      initial_flag=1; /* *updorder_pp is initialized */
    }


  /* We shuffle the contents of "*updorder_pp" array. */
  for(i=0;i<length;i++)
    {

      /* A current position in "updorder" array is "i". We chose one of the
	 elements (target) in this array among [i,length-1] randomly. We
	 switch the contents of the current element and that of the
	 target */
      target = (double)i + ((double)(length-i))*genrand_real3();
      /* The avobe is the same as
	 target = floor((double)i + (double)(length-i)*genrand_real2()); 
	 because the fractional part is trancated when double values are
	 assigned to int variables */
      /* genrand_real2() is a random number on [0,1)-real-interval. Do not
	 use [0,1] or (0,1] here because if genrand_real==1, "target" will
	 be "length", but "*updorder_pp + length" is outside of the
	 limit!! */
      if(i!=target){
	temp=*(*updorder_pp+target);
	*(*updorder_pp+target)=*(*updorder_pp+i);
	*(*updorder_pp+i)= temp;
      }
      /* else do nothing. we don't have to swap. */
    }
}

/* ATTENTION: This includes the copying of boundary */
TYPE2 **Copy2(TYPE2 **a,TYPE2 **b)
{
  memcpy(*a,*b,sizeof(TYPE2)*(nrow+2)*(ncol+2));
  return a;
}

/* ATTENTION: It does not fill the boundary */
TYPE2 **Fill2(TYPE2 **a,TYPE2 c)
{
  int i,j;
  /*$dir force_vector*/
  for (i=first,j=last; i <= j; i++)
    *(*a+i)=c;
  return a;
}


/************************************************ shift */

/* This function set the boundaries of TYPE2. Actually, it does
   not make sense to call this function more than once, because
   NeighSet() puts intermediate layer in the neighbor retrieval,
   which makes user's rule uniform for different boundary
   conditions. For details, look at NeighSet() and get to know
   how to use it. */
TYPE2 **Boundaries2(TYPE2 **a)
{
  int i,j;
  int nc,nr,nr1,nc1;
  TYPE2 bound;
  nr = nrow;
  nc = ncol;
  nr1 = nrow + 1;
  nc1 = ncol + 1;
  if (boundary == WRAP) {
    /*$dir force_vector*/
    for (i=1; i <= nr; i++) {
      a[i][0] = a[i][nc];
      a[i][nc1] = a[i][1];
    }
    /*$dir force_vector*/
    for (j=0; j <= nc1; j++) {
      a[0][j] = a[nr][j];
      a[nr1][j] = a[1][j];
    }
  }
  else if (boundary == FIXED) {
    bound=boundaryvalue2;
    /*$dir force_vector*/
    for (i=1; i <= nr; i++){
      a[i][0] = bound;
      a[i][nc1] = bound;
    }
    /*$dir force_vector*/
    for (j=0; j <= nc1; j++){
      a[0][j] = bound;
      a[nr1][j] = bound;
    //      a[0][j] = a[nr1][j] = bound;
    }
  }
  else if (boundary == ECHO) {
    /*$dir force_vector*/
    for (i=1; i <= nr; i++) {
      a[i][0]=a[i][1];
      a[i][nc1]=a[i][nc];
      //      a[i][0] = a[i][1];
      //      a[i][nc1] = a[i][nc];
    }
    /*$dir force_vector*/
    for (j=0; j <= nc1; j++) {
      a[0][j] = a[1][j];
      a[nr1][j] = a[nr][j];
    }
  }
  return a;
}


/*****************************************neighbors*/

/* NeighSet(int direction)

   It will return a pointer to the plain of "struct
   point". "return[i][j]" will give a correct coordinate of the
   neighborhood of [i][j]. The direction of neighborhood is
   specified by "direction" argument. The function argument can
   take CENTRAL, NORTH, WEST, EAST, SOUTH, NORTHWEST, NORTHEAST,
   SOUTHWEST, SOUTHEAST (defined in cash.h). This is usefull when
   "boundary" is WRAP. For instance, you want to get the north
   neighbor of "plane[1][1]" (where "plane[][]" is the user's CA
   plane). Then, you can do
   
   North=NeighSet(NORTH);

   and the 'coordinate' (not the value) of the neighbor is
   [North[1][1].row][North[1][1].col]. If you want to get the
   value of the neighbor, you can do

   plane[North[1][1].row][North[1][1].col].

   If "bounday==WRAP", then 
   North[1][1].row = nrow
   North[1][1].col = 1. 
   If "boundary==FIXED", then
   North[1][1].row = 0
   North[1][1].col = 1.
   If "boundary==ECHO", then
   North[1][1].row = 1
   North[1][1].col = 1.

   Therefore, if the user makes rules by using above way of
   neighbor retrieving, the user does not have to change the rule
   according to the boundary condition (For "boundary==FIXED",
   the user must set boundary values of the plane by
   Boundaries2(); for the other boundaries, it is not
   necessary).

   Note that for "boundary==ECHO",
   NorthEast[1][1] = {1,2}
   NorthWest[1][1] = {1,1}
   West[1][1] = {1,1}
   SouthWest[1][1] = {2,1} etc.
   
   Additionally, it's easier to use this function as
   for(i=1;i<=8;i++)
     neighbors[i] = NeighSet(i)
   rathar than the above. Then make a random number between 1 and
   8, you put it in "neighbors[random]", you get a random
   neighbor from Moor neighborhood.

   NOTE: this function is relatively easy to extended to include
   higher order neighbors (like North-North, North-West-West
   etc.).
*/
struct point **NeighSet(int direction)
{
  int i,j;
  int nc,nc1,nr,nr1;
  struct point **neighbor;

  nc = ncol;
  nr = nrow;
  nr1 = nrow + 1;
  nc1 = ncol + 1;
  
  neighbor = (struct point **) malloc(sizeof(struct point*)*(nr+2));
  if(neighbor==NULL){
    fprintf(stderr,"NeighborSet: error in memory allocation\n");
      exit(1);
  }
  
  neighbor[0] = (struct point*) malloc(sizeof(struct point)*(nr+2)*(nc+2));
  if(neighbor[0]==NULL){
    fprintf(stderr,"NeighborSet: error in memory allocation\n");
    exit(1);
  }
    
  for(i=1,j=nr+2;i<j;i++)
    neighbor[i] = neighbor[i-1] + nc + 2;

  for(i=1;i<=nr;i++)
    for(j=1;j<=nc;j++){
      neighbor[i][j].row = i;
      neighbor[i][j].col = j;
    }

  if (boundary == WRAP) {
    for (i=1; i <= nr; i++) {
      neighbor[i][0]=neighbor[i][nc];
      neighbor[i][nc1]=neighbor[i][1];
    }
    for (j=0; j <= nc1; j++) {
      neighbor[0][j]=neighbor[nr][j];
      neighbor[nr1][j]=neighbor[1][j];
    }
  }else if (boundary == FIXED) {
    for (i=1; i <= nr; i++){
      neighbor[i][0].row = i;
      neighbor[i][0].col = 0;
      neighbor[i][nc1].row = i;
      neighbor[i][nc1].col = nc1;
    }
    for (j=0; j <= nc1; j++) {
      neighbor[0][j].row = 0;
      neighbor[0][j].col = j;
      neighbor[nr1][j].row = nr1;
      neighbor[nr1][j].col = j;
    }
  }else if (boundary == ECHO) {
    for (i=1; i <= nr; i++) {
      neighbor[i][0]=neighbor[i][1];
      neighbor[i][nc1]=neighbor[i][nc];
    }
    for (j=0; j <= nc1; j++) {
      neighbor[0][j]=neighbor[1][j];
      neighbor[nr1][j]=neighbor[nr][j];
    }
  }

  /* Here we shift the plane */
  if(direction==NORTH||direction==NORTHEAST||direction==NORTHWEST){
    for(i=nr;i>=1;i--)
      neighbor[i]=neighbor[i-1];
  }
  if(direction==WEST||direction==NORTHWEST||direction==SOUTHWEST){
    for(i=1;i<=nr;i++)
      neighbor[i]=neighbor[i]-1;
  }
  if(direction==EAST||direction==NORTHEAST||direction==SOUTHEAST){
    for(i=1;i<=nr;i++)
      neighbor[i]=neighbor[i]+1;
  }
  if(direction==SOUTH||direction==SOUTHEAST||direction==SOUTHWEST){
    for(i=1;i<=nr;i++)
      neighbor[i]=neighbor[i+1];
  }

  return(neighbor);
}

void NeighFree(struct point** neighbor)
{
  free(neighbor[0]);
  free(neighbor);
}

/******************************************margolus*/

/* MargolusNeigh(void)

   You do "(MARGOLUS ***)margolus = MargolusNeigh(), then you get
   the Margolus neighborfood in margolus. Look at the definition
   of MARGOLUS in cash2.h. If the user wants to get the
   coordinate of CCW neighbor (counter-clock wise neighbor) of
   "plane[i][j]", then  the user can get the coordinate with
   "margolus[phase][i][j].m[CCW].row" and
   "margolus[phase][i][j].m[CCW].col". 
   The idea is similar to NeighSet(). 

   For the index of "m[]", you can use CW (clock-wise), OPP
   (opposite) and CCW (all are defined in cash2.h). "phase" can
   take value 0 or 1. If it is 0, it is even phase (i.e., [1][1]
   will be the upper left corner of the 2x2 block); otherwise odd
   phase (i.e., [1][1] will be the lower right corner of the 2x2
   block). For details, look at cash2.h.

   CAUTION: ECHO boundary is not supported. 
 */
MARGOLUS *** MargolusNeigh(void)
{
  MARGOLUS ***margolus;
  int nr,nc;
  int i,j;
  /* vector: [0]:->, [1]:down, [2]:<-, [3]:up . We will
     succesively add the vectors to the current coordinate. This
     means that we rotate the current coordinate in
     clock-wise. However, the order of addition has to be taken
     care of, depending on which position the current coordinate
     is lieing on in the Margolus neighborfood. If the current
     cell is in UL, then add vector[0] to [3]. If you add
     vector[0] to the coordinate of the current cell, then you
     get the coorinate of the CW of UL (i.e. UR). If you add [0]
     and [1], then you get OOP of UL (i.e. LR). If you add [0],
     [1], and [2], you get CCW of UL (i.e. LL). If the current
     cell is in LR, then you've got to add the vector in the
     following order: first [2] so you get LL, second [3] so you
     get UL, third [4] so you get UR. */
  int rvec[4]={0,1,0,-1};
  int cvec[4]={1,0,-1,0};
  int where;  /* UL=0, UR=1, LR=2, LL=3 */
  int position[4] = {0,1,3,2};/* UL=0, UR=1, LR=2, LL=3 */
  int phase;  /* if even, 0. if odd, 1. */
  
  /* MargolusNeigh() can be used only if ... */
  if (nrow & 1 || ncol & 1) {
    fprintf(stderr,"MargolusNeigh: Please use even values for nrow and ncol\n");
    exit(1);
  }

  if (boundary == ECHO){
    fprintf(stderr,"MargolusNeigh: ECHO boundary is not supported. Sorry\n");
    exit(1);
  }

  nr = nrow;
  nc = ncol;

  /* margolus[0] or [1] are for even/odd phase. */
  margolus = (MARGOLUS***) malloc(sizeof(MARGOLUS**)*2);
  if(margolus==NULL){
    fprintf(stderr,"MargolusNeigh: error in memory allocation\n");
      exit(1);
  }

  for(phase=0;phase<2;phase++){
    margolus[phase] = (MARGOLUS**) malloc(sizeof(MARGOLUS*)*(nr+2));
    if(margolus[phase]==NULL){
      fprintf(stderr,"MargolusNeigh: error in memory allocation\n");
      exit(1);
    }
    margolus[phase][0] = (MARGOLUS*) malloc(sizeof(MARGOLUS)*(nr+2)*(nc+2));
    if(margolus[phase][0]==NULL){
      fprintf(stderr,"MargolusNeigh: error in memory allocation\n");
      exit(1);
    }
    for(i=1,j=nr+2;i<j;i++)
      margolus[phase][i] = margolus[phase][i-1] + nc + 2;
  }

  if(boundary==WRAP){
    for(phase=0;phase<2;phase++)
      for(i=1;i<=nr;i++)
	for(j=1;j<=nc;j++){
	  /* "where" represents the place in the Margolus
	     neighborhood (UL:0 UR:1 LR:2 LL:3).
	     if phase=0 (even phase), (1,1)=0 (1,2)=1 (2,2)=2 (2,1)=3 etc
	     if phase=1 (odd phase),  (1,1)=3 (1,2)=2 (2,2)=1 (2,1)=0 etc 
	     Note that "x=2*((i+1+phase)%2) + (j+1+phase)%2" gives
	     us the following:
	     even phase (phase=0) : (i,j)=x,
	     (1,1)=0; (1,2)=1;
	     (2,1)=2; (2,2)=3.
	     odd phase (phase=1) : (i,j)=x,
	     (1,1)=3; (1,2)=2;
	     (2,1)=1; (2,2)=0.
	     
	     Note that 'x' does not give us the right number for
	     the place in the Margolus neighborhood (the right
	     value of "where"). But, "position[x]" will gives us
	     the right number with respect to the position in the
	     Margolus neighborhood such that

	     in the even phase:
	     position[0]=0(UL); position[1]=1(UR);
	     position[2]=3(LL); position[3]=2(LR); 

	     in the odd phase:
	     position[3]=2(LR); position[2]=3(LL);
	     position[1]=1(UR); position[0]=0(UL);.
	  */
	  where = position[2*((i+1+phase)%2) + (j+1+phase)%2];

	  /* 1st click clock wise. Turn the (i,j) in clock-wise
	     once, then that coordinate is that of the clock-wise
	     neighbor of (i,j) */
	  /* Note the following. Given we have the following array.
	     1  2 3 .. nr-1  nr (original row)
	     Now we want to shift this to , for instance, the left, so that we have,
	     nr 1 2 .. nr-2  nr-1 (original row + (-1)-vector)
	     We have to do this if boundary==WRAP.

	     We can do this as follows:
	     1    2    3     ...  nr-1    nr   (original)
	     0    1    2     ...  nr-2    nr-1 (-1; shift left, we do this because we want an array starting from 0)
	     nr   nr+1 nr+2  ...  nr-2+nr nr-1+nr (+nr, we do this because the left most one shold not get negative number)
	     nr-1 nr   nr+1  ...  nr-3+nr nr-2+nr (-1;shift left; this is is the real shifting. It can be anything)
	     nr-1 0    1     ...  nr-3    nr-2    (%nr)
	     nr   1    2     ...  nr-2    nr-1    (+1;shift right, this will cancell out the first shift)

	     It's complicating, but we can't use shifting as we
	     did in NeighSet().
	  */
	  margolus[phase][i][j].m[CW].row=(nr+(i-1)+rvec[where])%nr + 1;
	  margolus[phase][i][j].m[CW].col=(nc+(j-1)+cvec[where])%nc + 1;

	  /* 2nd click clock wise. Turn the (i,j) in clock-wise
	     twice, then the coorinate will be of that of the
	     oppsite of the (i,j)*/
	  margolus[phase][i][j].m[OPP].row=(nr+(i-1)+rvec[where]+rvec[(where+1)%4])%nr + 1;
	  margolus[phase][i][j].m[OPP].col=(nc+(j-1)+cvec[where]+cvec[(where+1)%4])%nc + 1;;

	  /* 3rd click clock wise. Turn the (i,j) in clock-wise
	     three times, then the coodinate will be of that of
	     the counter-clock of the (i,j) */
	  margolus[phase][i][j].m[CCW].row=(nr+(i-1)+rvec[where]+rvec[(where+1)%4]+rvec[(where+2)%4])%nr + 1;
	  margolus[phase][i][j].m[CCW].col=(nc+(j-1)+cvec[where]+cvec[(where+1)%4]+cvec[(where+2)%4])%nc + 1;
	}
  }
  /* The following will set Margolus neighbors for
     boundary==FIXED. */
  else {
    for(phase=0;phase<2;phase++)
      for(i=1;i<=nr;i++)
	for(j=1;j<=nc;j++){
	  where = position[2*((i+1+phase)%2) + (j+1+phase)%2];

	  /* 1st click clock wise. Turn the (i,j) in clock-wise
	     once, then that coordinate is the clock-wise neighbor
	     of (i,j) */
	  margolus[phase][i][j].m[CW].row=i+rvec[where];
	  margolus[phase][i][j].m[CW].col=j+cvec[where];

	  /* 2nd click clock wise. Turn the (i,j) in clock-wise
	     twice, then the coorinate will be in oppsite of the
	     (i,j)*/
	  margolus[phase][i][j].m[OPP].row=i+rvec[where]+rvec[(where+1)%4];
	  margolus[phase][i][j].m[OPP].col=j+cvec[where]+cvec[(where+1)%4];

	  /* 3rd click clock wise. Turn the (i,j) in clock-wise
	     three times, then the coodinate will be in
	     counter-clock of the (i,j) */
	  margolus[phase][i][j].m[CCW].row=i+rvec[where]+rvec[(where+1)%4]+rvec[(where+2)%4];
	  margolus[phase][i][j].m[CCW].col=j+cvec[where]+cvec[(where+1)%4]+cvec[(where+2)%4];
	}
  }
  
  
  /* debug. To see whether MargolusNeigh() sets the Margolus
     neighborhood correctly. */
  /* {
    int j;
    int pos;
    int phase;
    char whe[3][4] = {"CW","OPP","CCW"};
    for(phase=0;phase<2;phase++){
      printf("phase %d\n",phase);
      for(pos=0;pos<3;pos++){
	printf("%s\n",whe[pos]);
	for(i=1;i<=nrow;i++){
	  for(j=1;j<=ncol;j++){
	    printf("(%d,%d) ,",margolus[phase][i][j].m[pos].row,margolus[phase][i][j].m[pos].col);
	  }
	  printf("\n");
	}
	printf("\n");
      }
      printf("\n");
    }
    }*/
  return margolus;
}



void MargolusFree(MARGOLUS*** margolus)
{
  free(margolus[0][0]);
  free(margolus[0]);
  free(margolus);
}


/*
  Little macro for MargolusDiffusion
  CC=current cell; CW=clock
  wise; OPP=oppsite; CCW=counter clock wise;
  TEMP=temporary. Arguments have to be "TYPE2"
*/
#define ROTCW(CC,CW,OPP,CCW,TEMP) {TEMP = CC; CC = CCW; CCW = OPP; OPP = CW; CW = TEMP;}
#define ROTCCW(CC,CW,OPP,CCW,TEMP) {TEMP = CC; CC = CW; CW = OPP; OPP = CCW; CCW = TEMP;}


/* MargolusDiffusion()

   Give your plane and Margolus neighborhood plane, it will do
   one step of Margolus diffusion. 'boundary' must be either WRAP
   or FIXED. If the boundary is FIXED, Margolus diffusion will be
   suppressed in the boundary in the odd phase so that
   "boundaryvalue2" will not flow into the plane. If the user wants
   that the entities of the model diffuse away from the plan by
   the boundary, make the user's own Margolus diffusion after
   reading the following code. With the use of MargolusNeigh(),
   it must be very simple to impliment.

   This algorithm was 30% faster (though more memory comsuming)
   than the original implimentation of Margolus diffusion in
   CASH2003 with TYPE=int and TYPE2=int with the optimization
   option -O3. If the size of types are larger, I expect that the
   difference between these algorithm  will be even larger.
*/
void MargolusDiffusion(TYPE2** world,MARGOLUS*** margolus)
{
  int nr1,nc1;
  int i,j;
  int cc_r,cc_c,cw_r,cw_c,opp_r,opp_c,ccw_r,ccw_c;
  static int phase=0;
  TYPE2 temp;

  nr1 = nrow-1;
  nc1 = ncol-1;

  /* If boundary==WRAP, do Margolus diffusion in the entire
     space. If not, but the phase=0, then do the diffusion in the
     entire space too. The reason why we do diffusion if phase=0
     even if the boundary condition is not WRAP is that any
     Margolus blocks do not go over the boundary in even phase.*/
  if(phase==0 || boundary==WRAP){
    /* (i,j) goes through each Margolus-block once except for
       boundary. In even-phase (phase=0), (i,j) will be "UL" of
       Margolus-block. In odd-phase (phase=1), (i,j) will be "LR"
       of Margolus-block. */
    for(i=1;i<=nr1;i+=2)
      for(j=1;j<=nc1;j+=2){
	/* (i,j)=HERE */
	cc_r = i;
	cc_c = j;
	/* (cw_r,cw_c)=CW_of_(i,j) */
	cw_r = margolus[phase][i][j].m[CW].row;
	cw_c = margolus[phase][i][j].m[CW].col;
	/* (opp_r,opp_c)=OPP_of_(i,j) */
	opp_r = margolus[phase][i][j].m[OPP].row;
	opp_c = margolus[phase][i][j].m[OPP].col;
	/* (ccw_r,ccw_c)=CCW_f_(i,j) */
	ccw_r = margolus[phase][i][j].m[CCW].row;
	ccw_c = margolus[phase][i][j].m[CCW].col;

	if(0.5 < genrand_real1()) 
	  ROTCW(world[cc_r][cc_c],world[cw_r][cw_c],world[opp_r][opp_c],world[ccw_r][ccw_c],temp)
	else 
	  ROTCCW(world[cc_r][cc_c],world[cw_r][cw_c],world[opp_r][opp_c],world[ccw_r][ccw_c],temp)
      }
  }

  /* We do not do Margolus diffusion on the boundary if
     boundary!=WRAP and if phase==odd: We suppress the diffusion
     in the edge in this case. */
  else{
    for(i=3;i<=nr1;i+=2)
      for(j=3;j<=nc1;j+=2){
	/* (i,j)=HERE */
	cc_r = i;
	cc_c = j;
	/* (cw_r,cw_c)=CW_of_(i,j) */
	cw_r = margolus[phase][i][j].m[CW].row;
	cw_c = margolus[phase][i][j].m[CW].col;
	/* (opp_r,opp_c)=OPP_of_(i,j) */
	opp_r = margolus[phase][i][j].m[OPP].row;
	opp_c = margolus[phase][i][j].m[OPP].col;
	/* (ccw_r,ccw_c)=CCW_f_(i,j) */
	ccw_r = margolus[phase][i][j].m[CCW].row;
	ccw_c = margolus[phase][i][j].m[CCW].col;
	
	if(0.5 < genrand_real1()) 
	  ROTCW(world[cc_r][cc_c],world[cw_r][cw_c],world[opp_r][opp_c],world[ccw_r][ccw_c],temp)
	else 
	  ROTCCW(world[cc_r][cc_c],world[cw_r][cw_c],world[opp_r][opp_c],world[ccw_r][ccw_c],temp)
      }
  }

  /* Change the phase for next time */
  phase = (phase+1)%2;
}

