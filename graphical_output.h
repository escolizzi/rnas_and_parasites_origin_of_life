#include "cash.h"
#include "cash2.h"

void InitColormap(void);
int GetColorIndexFromCat(int val,double kcat);
int GetColorIndexFromPar(int val,double kcat);
void SaveMovie(TYPE2 **world,int Time);
int SaveDataRandomOrder(TYPE2 **world, int Time,struct updorder* updorder_p);
void OutputBackup(TYPE2 **world,int Time);
void SaveData(TYPE2 **world,char *fname,int Time);

