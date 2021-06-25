//External parameters, regardless the model
extern int nrow,ncol,scale;
extern int boundary,boundaryvalue;
extern TYPE2 boundaryvalue2;


void Initial(TYPE2 **);
void Meteorites(TYPE2 **world);
void AllThatCanHappen(TYPE2 **world, struct updorder* updorder_p,struct point** neighbor[]);
void Complex_Formation(double kdiss, TYPE2 *cat,int irow, int icol,TYPE2 *tmpl,int ineirow,int ineicol);
void Complex_Dissociation(TYPE2 *x,TYPE2 *xn);
void ToReplication(TYPE2 *icel,TYPE2 *inei,TYPE2 *neinei);
void MolDecay(int howmany, ...);
void SwapCells(int order, ...);

void Reactions(TYPE2 **world, struct updorder* updorder_p,struct point** neighbor[]);
void Decay(TYPE2 **world);
void MyDiffusion(TYPE2 **a, struct updorder* updorder_p,struct point** neighbor[], double pdiff);

void Finish(TYPE2 **world,struct point **neighbor[],struct updorder *updorder_p);
