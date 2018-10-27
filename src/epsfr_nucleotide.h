#ifndef __NUCLEOTIDE__
#define __NUCLEOTIDE__

// Basic functionality about the nucleotide naming, pairing, numbering, etc.
bool check_rna_pairing(char nuc1, char nuc2);
char num_2_nuc(int num);
int nuc_2_num(char nuc);
bool is_valid_nuc_char(char nuc);

#endif // __NUCLEOTIDE__