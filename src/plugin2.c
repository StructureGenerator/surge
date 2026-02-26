/* This is a plugin for surge that optionally forbids
   adjacent oxygen atoms. */

#define HELPTEXT2 \
" This version forbids adjacent oxygen atoms if -Y is given.\n"

static boolean Yswitch = FALSE;
#define SURGEPLUGIN_SWITCHES  SWBOOLEAN('Y',Yswitch)

static int oxygenindex = -1;
#define SURGEPLUGIN_STEP2 \
 if (oxygenindex < 0) oxygenindex = elementindex("O"); \
 if (Yswitch) { int ii; for (ii = 0; ii < ne; ++ii) \
   if (vcol[edge[ii].x] == oxygenindex \
    && vcol[edge[ii].y] == oxygenindex) return; }
