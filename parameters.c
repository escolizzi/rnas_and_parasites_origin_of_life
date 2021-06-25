#include <limits.h>
#include "parameters.h"

//SIMULATION PARAMETERS
int par_inittime=0; // initial time
int par_maxtime=INT_MAX; //for how time steps the simulation goes on
int par_time_savemovie=1000; //how often we save figure for movie
int par_time_savedata=50000; //how often we save data
int par_time_backup=500000; //how often we backup data
long int par_ulseed=334608376; //seed for the random numb. gen. (long int!)
double global_scale=0.2; //scaling constant in interval [0,1) to further slow down the dynamics (esp for complex formation)

//CHEMISTRY PARAMETERS
double par_factor_asso = -0.05; //energy scaling factor compl form
double par_factor_diss = -0.05; //compl diss
double par_kdec = 0.03; //decay rate (= decay probability)
double par_kdiff = 0.1; //basic diffusion rate
double par_beta=1.2;  //USED ALSO FOR INITIALISE PAR BETA
double par_kdiss=0.1;

//initialisation catalyst
int init_val=1;
double init_krec=0.5;
double init_kcat=1.;
double init_kmut=0.005;
double init_kdelta=0.05;
char init_sign='+';
double init_ktime=0.;

/*
double init_ckrec=0.3;
double init_ckcat=1.;
double init_ckmut=0.005;
double init_ckdelta=0.05;
*/

//initialisation parasite
//double init_beta = 1.4; //parasite advantage, it is saved in krec!

//               0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 
int init_seq[MAXLEN]={1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};
char par_init_end5[MAXLEN]="AAAAUUUUCCCCGGG";
char par_init_end3[MAXLEN]="AAUUAAUUCCGGCCG";

//OUTPUT
int par_nplane=2;
char par_name_out_file[MAXLEN]="dataMut.txt";
char par_movie_directory_name[MAXLEN]="movieMut";
char par_backup_directory_name[MAXLEN]="backupMut";
char par_backup_file_name[BUFSIZE]="backupMut";
char par_input_file_name[BUFSIZE]="";

