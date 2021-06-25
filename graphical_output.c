#include <stdio.h>
//#include "cash2.h"
#include "graphical_output.h"
#include "parameters.h"
#include "my.h"

//  initialise RGB color map
//  index 0 must be black
//  index 1 must be white
// MAX NUMBER OF COLOUR THAT CAN BE DEFINED IS 256 (see color.c)
void InitColormap(void)
{
  int i, r,g,b;
  int max_steps, how_many=250;// +2: black and white
  double step_size;
  //for now we use a gradient: we go blue->magenta->red->yellow->green.
  //given that from blue to magenta there are 256 steps [0->255], and
  //given that we make 4 transitions to colours, max_steps=4*256.
  //if we want 'how_many' colours, then our step_size will be max_steps/how_many
  //and we increase/decrease red,green,blue by that much each step
  ColorRGB(0,0,0,0);
  ColorRGB(1,255,255,255);
  
  r=0;g=0;b=255;	//start from blue
  max_steps=4*256;
  step_size=(int)((float)max_steps/(float)how_many);
  for(i=0;i<how_many;i++){
    ColorRGB(i+2,r,g,b);
    
    if(i*step_size< 255 - step_size) r += step_size;	//increase red, we go to magenta, and so on...
    else if(i*step_size< 2*(255 - step_size)) b -= step_size; //decrease blue, go to red
    else if(i*step_size< 3*(255 - step_size)) g += step_size; //increase green, go to yellow
    else if(i*step_size< 4*(255 - step_size)) r -= step_size; //decrease red, go to green
        
  }
  
  ColorRGB(0,0,0,0);
  ColorRGB(1,255,255,255);	//yeah yeah, we are repeating ourselves, but it's just to make sure
  //ColorRGB(1,255,0,0);
  //ColorRGB(2,0,255,0);
  //ColorRGB(3,255,255,0);
  //ColorRGB(4,255,0,255);
}


int GetColorIndexFromCat(int val,double kcat)
{
  int color;
  
  if(val==1){
    //we have 250 colors, scale kcat to be in ]2,252]
    color=2 + (int)( 0.5 + 249.*kcat/2.);
  }else if(val==2){
    color=1;
//     //krec for par is in interval (0,inf.)
//     //let's assume that the value will be bound in (0.,5.)
//     if(kcat>5.) color=1;
//     else color=2+ (int) ( 0.5 + kcat*250./5. );
  }
  
  return color;
  
  
}

int GetColorIndexFromPar(int val,double kcat)
{
  int color;
  
  if(val==1){
    //we have 250 colors, scale kcat to be in ]2,252]
    color=1;
  }else if(val==2){
     //krec for par is in interval (0,inf.)
     //let's assume that the value will be bound in (0.,2.)
     if(kcat>2.) color=1;
     else color=2+ (int) ( 0.5 + kcat*249./2. );
  }
  
  return color;
  
  
}

void SaveMovie(TYPE2** world, int Time)
{
  int nr,nc,i,j,color_index_cat,color_index_par;
  int **tomovie;
  
  nr=nrow; nc=ncol;
  if(  NULL == ( tomovie=(int **)malloc( (nr+1) *sizeof(int*)) )  ){
   fprintf(stderr,"SaveMovie(): Error. Memory allocation unsuccessful.\n");
   return;
  }
  for(i=0; i<=nr;i++){
    if(NULL==(tomovie[i]=(int *)malloc( par_nplane*(nc+1) *sizeof(int)))) {
      fprintf(stderr,"SaveMovie(): Error. Memory allocation unsuccessful.\n");
      return;
    }
  }
  //here you'll save data
  for(i=0;i<nr;i++)for(j=0;j<nc;j++){ 
    if(world[i+1][j+1].val==0){
      color_index_cat=0;
      color_index_par=0;
    }else{
      color_index_cat=GetColorIndexFromCat(world[i+1][j+1].val,world[i+1][j+1].krec); //this is function that makes colours from some features...
      color_index_par=GetColorIndexFromPar(world[i+1][j+1].val,world[i+1][j+1].krec); //this is function that makes colours from some features...
    }
    tomovie[i+1][j+1]=color_index_cat;
    tomovie[i+1][j+nc+1]=color_index_par;
  }
  //printf("Hello from save data\n");
  //For pretty movies you can include margin... for now this is ok
  
  // 0 because we overwrite the entire old colormap, if you didn't, and only added colours after the pre-defined ones, then you would put 16... I think
  PlanePNG(tomovie,0);
  
  for(i=0; i<nr;i++) free(tomovie[i]);
  free(tomovie);
}

int SaveDataRandomOrder(TYPE2 **world, int Time,struct updorder* updorder_p)
{
  int i,length,irow,icol,control=0;
  //int k;
  TYPE2 *icel;
  FILE *fp;
  
  fp=fopen(par_name_out_file,"a");
  if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
  
  
  //for(i=1;i<=nrow;i++) for(j=1;j<=ncol;j++){
  
  length = nrow*ncol;
  for(i=0;i<length;i++){
    irow=(updorder_p+i)->row;
    icol=(updorder_p+i)->col;
    
    icel=&world[irow][icol];
    if(icel->val>0){
      control++;
      //fprintf(fp,"%d %c %f %f %f %s %s\n", Time, icel->sign, icel->kcat, icel->kmut, icel->kdelta, icel->end5,icel->end3);		//OLD MODEL
      //fprintf(fp,"%d %d %c %f %f %f %f %f ", Time, icel->val, icel->sign, icel->krec, icel->kcat,icel->ckcat, icel->kmut, icel->kdelta);
      fprintf(fp,"%d %d %c %f %f %f %f ", Time, icel->val, icel->sign, icel->krec, icel->kcat, icel->kmut, icel->kdelta);
      //for(k=0;k<MAXLEN;k++) fprintf(fp,"%d",icel->seq[k]);
      fprintf(fp,"\n");
    }
  }
  
  
  fclose(fp);
  //fprintf(stderr,"Pop size: %d\n",count);
  return control;
}

/* This function will call all the other outputting
   functions. This also checks if the population is extinct; if
   so, it will set stopG to be 1. */
void OutputBackup(TYPE2 **world,int Time)
{
  char fname[BUFSIZE];
  char command[512];
    
  sprintf(command,"%s %s","mkdir -p",par_backup_directory_name);
  if(system(command)==-1){
    fprintf(stderr,"OutputBackup: Failed to mkdir %s. Save here\n",par_backup_directory_name);
    sprintf(fname,"%s_t%d",par_backup_file_name,Time);
  }
  else
    sprintf(fname,"%s/%s_t%d",par_backup_directory_name,par_backup_file_name,Time);
  
  SaveData(world,fname,Time);
}

// SaveData writes backup files, almost all information is written, 
// including with whom icel is in complex with (coordinates) and incel->kdiss
void SaveData(TYPE2 **world,char *fname, int Time)
{
  int i,j;
  //int k;
  int nr=nrow,nc=ncol;
  TYPE2 *icel;
  FILE *fp;
  
  fp=fopen(fname,"a");
  if(fp==NULL) fprintf(stderr,"SaveData(): Warning. File not opened.\n");
  
  /* The size of field */
  fprintf(fp,"#%d %d X %d val sign bonr bonc term kdiss krec kcat kmut kdelta ckrec ckcat ckmut ckdelta\n",Time,nr,nc);
  
  for(i=1;i<=nrow;i++) for(j=1;j<=ncol;j++){
    if(world[i][j].val>0){
      icel=&world[i][j];
      //fprintf(fp,"%d %c %f %f %f %s %s\n", Time, icel->sign, icel->kcat, icel->kmut, icel->kdelta, icel->end5,icel->end3);		//OLD MODEL
      fprintf(fp,"%d %c %d %d %d %f %f %f %f %f ", icel->val, icel->sign, icel->bonr, icel->bonc, icel->term, icel->kdiss, icel->krec, icel->kcat, icel->kmut, icel->kdelta);
      //for(k=0;k<MAXLEN;k++) fprintf(fp,"%d",icel->seq[k]);
      fprintf(fp,"\n");
    }else fprintf(fp,"N\n");
  }
  fclose(fp);
  //fprintf(stderr,"Pop size: %d\n",count);
}

