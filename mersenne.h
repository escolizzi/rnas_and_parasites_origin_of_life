/********************************************mersenne*/
#ifndef MERSENNE
#define MERSENNE
/* Meresenne Twister, call function. See 
http://www.math.keio.ac.jp/~matumoto/mt.html for details. */
void init_genrand(unsigned long ); /* Set a random seed */
void init_by_array(unsigned long [], unsigned long); /* Set a random seed */
//static void next_state(void);
unsigned long genrand_int32(void); /* unsigned 32-bit integers */
long genrand_int31(void); /* unsigned 31-bit integers */
double genrand_real1(void); /* uniform real in [0,1] (32-bit resolution) */
double genrand_real2(void); /* uniform real in [0,1) (32-bit resolution) */
double genrand_real3(void); /* uniform real in (0,1) (32-bit resolution) */
double genrand_res53(void); /* uniform real in [0,1) with 53-bit
			       resolution */

/* Distribution */
double bnldev(double pp, int n);

/* Other functions */
double gammln(double); /* ln(GammaFunction) */
double factln(int); /* ln(N!). N is an argument */
double combinat(int, int); /* Combinatorial. Coefficient of
			      binomial distribution */
int genrand_int(int min,int max); //an random int from [min,max]

/* Chose "nchoice" numbers of element from the length of "length"
   sequence. */
void rand_choice(int length, int nchoice, int *result);
/* Chose one of the probability element of array, of which size
   is "nitem"  */
int choose_claim(int nitem,double *array);
#endif // MERSENNE
