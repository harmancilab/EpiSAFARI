#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "epsfr_nucleotide.h"
#include <string.h>
#include "epsfr_ansi_string.h"

int nuc_2_num(char nuc)
{
	if(toupper(nuc) == 'A')
	{
		return(0);
	}
	else if(toupper(nuc) == 'C')
	{
		return(1);
	}
	else if(toupper(nuc) == 'G')
	{
		return(2);
	}
	else if(toupper(nuc) == 'U' ||
		toupper(nuc) == 'T')
	{
		return(3);
	}
	else
	{
		return(4);
		//fprintf(stderr, "Cannot convert %c to number.\n", nuc);
		//exit(0);
	}
}

char num_2_nuc(int num)
{
	char nucs[] = "ACGT";

	char nuc = nucs[num];

	return(nuc);
}
