
void Replication(TYPE2 *child,TYPE2 *parent,TYPE2 *cat);
void PerfectReplication(TYPE2 *child,TYPE2 *parent);
void Mutations(TYPE2 *child,TYPE2 *cat);
void IntSeqSubstitutions(int *seq, double kmut);
void SeqSubstitutions(char *seq, double kmut);
char NucleotideMutation(char);
void IntSeqRepl(int *trgt, int *src);
void SeqRepl(char *trgt, char *src);
char RNAc2cc(char rna);
