/* This is a plugin for surge that counts how many hydrogen atoms
   are attached to carbon atoms. */

#define HELPTEXT2 \
" This version counts the hydrogen atoms attached to carbon atoms.\n"

static int carbonindex = -1;   /* index into element table */
static long long CHcount[5*MAXN+1];

#define SURGEPLUGIN_STEP3 \
 if (carbonindex < 0) carbonindex = elementindex("C"); \
 { int ii,CHval; CHval=0; for (ii = 0; ii < n; ++ii) \
      if (vcol[ii] == carbonindex) CHval += hyd[ii]; \
   ++CHcount[CHval]; } 

#define SURGEPLUGIN_SUMMARY \
 fprintf(stderr,"Counts by the number of hydrogens attached to carbons:\n"); \
 { int ii; for (ii = 0; ii <= 5*MAXN; ++ii) \
     if (CHcount[ii] > 0) fprintf(stderr," %2d : %lld\n",ii,CHcount[ii]); }
