/* This is a plugin for surge that implements an extra option
   -F# or -F#:# for the number of atoms with exactly 4 distinct
   non-H neighbours. */

#define HELPTEXT2 \
" -F# -F#:# Specify number of atoms with exactly 4 non-H neighbours\n"

static boolean Fswitch = FALSE;
static long Fmin,Fmax;

#define SURGEPLUGIN_STEP1 \
 { int ii,Fval; Fval=0; \
   for (ii = 0; ii < n; ++ii) if (deg[ii] == 4) ++Fval; \
   if (Fswitch && (Fval < Fmin || Fval > Fmax)) return; }

#define SURGEPLUGIN_SWITCHES \
  SWRANGE('F',":-",Fswitch,Fmin,Fmax,"surge -F")
