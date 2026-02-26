/* This is a plugin for surge that implements an extra option
   -V# or -V#:# for the number of atoms with exactly 4 distinct
   non-H neighbours. */

#define HELPTEXT2 \
" -V# -V#:# Specify number of atoms with exactly 4 non-H neighbours\n"

static boolean Vswitch = FALSE;
static long Vmin,Vmax;

#define SURGEPLUGIN_STEP1 \
 { int ii,Vval; Vval=0; \
   for (ii = 0; ii < n; ++ii) if (deg[ii] == 4) ++Vval; \
   if (Vswitch && (Vval < Vmin || Vval > Vmax)) return; }

#define SURGEPLUGIN_SWITCHES \
  SWRANGE('V',":-",Vswitch,Vmin,Vmax,"surge -V")
