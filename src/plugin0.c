/* This is a plugin for surge that removes molecules with
   two atoms having three or more common neighbours. 
   Also, arsenic at valence 3 and 5 is added. */

#define HELPTEXT2 " This version removes molecules with K(2,3).\n"

#define SURGEPLUGIN_INIT \
   addelement("As","As",3,3); addelement("Az","As",5,5);

#define SURGEPLUGIN_STEP0 \
 { int ii,jj; \
   for (jj = n; --jj >= 1; ) \
   for (ii = jj; --ii >= 0; ) \
       if (POPCOUNT(g[ii] & g[jj]) >= 3) return 1; }
