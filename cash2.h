#include "constants.h"

#ifndef CASH2
#define CASH2

/* CA's cell type. You define your own type here. structure is
   allowed */
typedef struct __type2 { 
  int val;
  int seq[MAXLEN]; //do char[], anything else is not worth it
  //char seq[MAXLEN];
  //char end5[MAXLEN], end3[MAXLEN];	//dangling ends
  double krec;						//intrinsic recognition rate
  double kcat;						//intrinsic catalytic activity
  double kmut;						//intrinsic mutation rate
  double kdelta;					//intrinsic mutational step
  
  //double ckrec;						//complementary of krec
  //double ckcat;						//complementary of kcat
  //double ckmut;						//complementary of mutation rate
  //double ckdelta;					//complementary of mutational step
  
  char sign;						//arbitrary sign to track =/- strands
  double ktime;						//time needed to replicate a molecule
  
  //Variable over the lifetime of the molecule
  int bonr;				//row of whom I'm in complex with, if no one, then =-1
  int bonc;				//col as above
  int term;				//side of the complex I'm in: 5' is replicase, 3' is template
  double kdiss;			//dissociation rate = 1- complex form. probability
  
  double trepl;		//if being replicated, time until end of replication
  int to_be_zero;	//if trepl<1. activate this to zero time the next time step
} TYPE2;


/* This structure is used in asynchronous updating to determine
   the order of the cells to update. Look at the functions in
   chash2.c which use "struct updorder" for details and the usage
   of this structure. Look at UpdOrdReset() in cash2.c  */
struct updorder { 
  int row;
  int col;
};

/* This structure is used in neighborhoods retrieving. Look
   at NeighSet() in cash2.c */
struct point { 
  int row;
  int col;
};

/* This structure is used in Margolus neighborhood
   retrieving. Read MarGolusNeigh(). If you say "MARGOLUS
   something_you_defined[*][i][j]" where [*] denote even-or-odd
   phase, then "----[*][i][j].m[CW].row" and
   "----[*][i][j].m[CW].col" give you the coordinate of the CW
   neighbor (see below) of Margolus neighborhood of the [i][j]
   cell in the phase [*] (even/odd, 0 or 1).

   -----------
   |HERE| CW |
   -----------
   | CCW| OPP|
   -----------

   Note that the exact place of "HERE" does not matter in that
   "CW", "OPP" and "CCW" notations are invariant for the position
   of "HERE".

   What is even and what is odd? Let "SYD" denote
   "something_you_defined". Then, SYD[0] and SYD[1] gives you tow
   alternative Margolus partitioning of the plane. Then, SYD[0]'s
   partition will be such that [0][0] cell of your plane will be
   the upper-left corner (UL) of the 2x2 square. SYD[1]'s
   partition will be such that [0][0] cell  of your plane will be
   the lower-right corner (LR) of the 2x2 square. See the
   followings too. 

   UL,UR,LL,LR are like,
   -----------
   | UL | UR |
   -----------
   | LL | LR |
   -----------


   in even phase (SYD[0]), the plane is divided as follows:
   -----------------
   | 1,1 | 1,2 | etc.
   -------------
   | 2,1 | 2,2 | etc.
   -------------
   | etc. etc.   etc.

   in odd phase (SYD[1]) and boundary=WRAP, the plane is divided as follows:
   ------------------
   |nr,nc| nr,1| etc.
   -------------
   | 1,nc|  1,1| etc.
   -------------
   | etc. etc.   etc.

   where "nr" is the number of row, and "nc" is the number of
   column. If boundary==FIXED, then read the avobe table such
   that "nr=nc=0".
*/
#define CW (0)
#define OPP (1)
#define CCW (2)
typedef struct __margolus{
  struct point m[3];/* 0:clockwise 1:opposite 2:counter-clockwise */
} MARGOLUS;


/*********************************************basic*/
TYPE2 **NewP2(void);
TYPE2 **New2(void);
int PlaneFree2(TYPE2**);
void UpdOrdReset(struct updorder**);
TYPE2 **Copy2(TYPE2 **,TYPE2 **);
TYPE2 **Fill2(TYPE2 **,TYPE2);

/*********************************************shift*/
TYPE2 **Boundaries2(TYPE2**);

/*****************************************neighbors*/
struct point **NeighSet(int); 
void NeighFree(struct point**);

/******************************************margolus*/
MARGOLUS ***MargolusNeigh(void);
void MargolusFree(MARGOLUS***);
void MargolusDiffusion(TYPE2**,MARGOLUS***);

#endif //CASH2
