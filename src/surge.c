/* This is a molecule generator based on geng.
   Version 1.0, November 11, 2021.

   Unix-style compilation command would be:

     gcc -o surge -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
         -march=native -mtune=native -DPREPRUNE=surgepreprune \
         -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
         surge.c geng.c planarity.c nautyW1.a

   You can build-in gzip output using the zlib library (https://zlib.net).
   Add -DZLIB to the compilation, and link with the zlib library either
   by adding -lz or libz.a . This will activate the -z command to gzip
   the output.

     gcc -o surge -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
         -march=native -mtune=native -DPREPRUNE=surgepreprune -DZLIB \
         -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
         surge.c geng.c planarity.c nautyW1.a -lz

   There is a makefile in the package; edit the first few lines.

   This version works best with geng version 3.2 or later. To use
   with an earlier version, add -DOLDGENG to the compilation command.

***********************************************************************/

#ifdef ZLIB
#define USAGE \
  "[-oFILE] [-z] [-u|-A|-S] [-T] [-e#|-e#:#] [-d#] [-c#] [-m#/#] formula"
#else
#define USAGE \
  "[-oFILE] [-u|-A|-S] [-T] [-e#|-e#:#] [-d#] [-c#] [-m#/#] formula"
#endif

#define HELPUSECMD

#define HELPTEXT1 \
"Make chemical graphs from a formula. Version 1.0.\n" \
"  Known elements are C,B,N,P,O,S,H,Cl,F,Br,I at their lowest valences.\n" \
"  Higher valences can be selected using Nx (Nitrogen/5), Sx,Sy (Sulfur 4/6)\n" \
"   Px (Phosphorus/5).\n" \
"\n" \
"  formula = a formula like C8H6N2\n" \
"\n" \
"  -u    Just count, don't write\n" \
"  -S    Output in SMILES format\n" \
"  -A    Output in alphabetical format\n" \
"  -O#   Output stage: 1 after geng, 2 after vcolg, 3 after multig\n" \
"        Default is to write SDfile format\n" \
"  -e# -e#:#  Restrict to given range of distinct non-H bonds\n" \
"  -t# -t#:#  Limit number of rings of length 3\n" \
"  -f# -f#:#  Limit number of cycles of length 4\n" \
"  -p# -p#:#  Limit number of cycles of length 5\n" \
"  -b    Only rings of even length (same as only cycles of even length)\n" \
"  -T    Disallow triple bonds\n" \
"  -P    Require planarity\n" \
"  -d#   Maximum degree not counting bond multiplicity or hydrogens (default 4)\n" \
"  -c#   Maximum coordination number (default 4). This is the maximum number\n" \
"         of distinct atoms (including H) that an atom can be bonded to\n" \
"         Coordination number > 4 is only allowed if no neighbours are H\n" \
"  -B#,...,# Specify sets of substructures to avoid (details in manual)\n" \
"     1 = no triple bonds in rings up to length 7\n" \
"     2 = Bredt's rule for two rings ij with one bond in\n" \
"           common (33, 34, 35, 36, 44, 45)\n" \
"     3 = Bredt's rule for two rings ij with two bonds in\n" \
"           common (i,j up to 56)\n" \
"     4 = Bredt's rule for two rings of length 6 sharing three bonds\n" \
"     5 = no substructures A=A=A (in ring or not)\n" \
"     6 = no substructures A=A=A in rings up to length 8\n" \
"        For -B5 and -B6, the central atom only has 2 non-H neighbours\n" \
"     7 = no K_33 or K_24 structure\n" \
"     8 = none of cone of P4 or K4 with 3-ear\n" \
"     9 = no atom in more than one ring of length 3 or 4\n" \
"  -v     Write more information to stderr\n" \
"  -m#/#  Do only a part. The two numbers are res/mod where 0<=res<mod.\n" \
"  -oFILE Write the output to the given file rather than to stdout.\n" \
"  -E..  Define a new element (see the manual)\n" \
"  -z     Write output in gzip format (only if compiled with zlib)\n"

#define EXTRAUSAGE ""
#define HELPTEXT2 ""

/* Undocumented options:
  -G...  Anything to the end of the parameter is passed to geng
  -x     Used for development purposes; not useful for users
*/

#define MAXN WORDSIZE    /* Not bigger than WORDSIZE, which can be 32 or 64 */
#define MAXNE (2*MAXN)
#include "gtools.h"
#include "naugroup.h"
#include "planarity.h"
#include <ctype.h>
#ifdef ZLIB
#include "zlib.h"
#endif

#define DEBUG 0

static struct smilesstruct
{
    int item;
    int x,y,r;
} smilesskeleton[4*MAXN+6*MAXNE];
/* Values for the item field */
#define SM_ATOM  1  /* Atom number x */
#define SM_BOND  2  /* Bond x,y */
#define SM_OPEN  3  /* Open ( */
#define SM_CLOSE 4  /* Close ) */
#define SM_RING0 5  /* Broken bond x,y for new ring r */
#define SM_RING1 6  /* End of ring r */
static int smileslen;

typedef unsigned long long counter;  /* For counters that might overflow */

static int hydrogens;
static int nv;  /* number of atoms except H */
static int numtypes;  /* different elements except H */
#define FORMULALEN 50
static int elementtype[FORMULALEN],elementcount[FORMULALEN];
                        /* Not H, decreasing valence */
static char canonform[2*FORMULALEN]; /* The formula in canonical order */
static int valencesum;  /* Sum of valences for non-H */
static int maxtype[8]; /* maxtype[d] is maximum index into
                    elementtype/elementcount for vertices of degree d */

/* Used with -A */
static boolean alphabetic;
static int newlabel[MAXN];  /* New number of atom i */

#define BADLISTS 9     /* Number of defined bad lists */
static boolean bad1;    /* Avoid triple edges in rings up to length 7 */
      /* bad1 is turned off if -T is given */
static boolean bad2;    /* Bredt's rule for one common bond */
static boolean bad3;    /* Bredt's rule for two common bonds */
static boolean bad4;    /* Bredt's rule for three common bonds */
static boolean bad5;    /* Avoid =A= even if not in a ring */
static boolean bad6;    /* Avoid =A= in rings up to length 8 */
      /* Note that bad6 is turned off if bad5 is set */
static boolean bad7;    /* Avoid K_{2,4} and K_{3,3} */
static boolean bad8;    /* Avoid cone(P4) and K4 with 3-ear */
      /* bad8 is turned off if -t and -f options make it impossible */
static boolean bad9;    /* No atom on two rings of length 3 or 4 */
      /* bad9 is turned off if -t and -f options make it impossible */

static boolean needcoordtest;

static boolean needcycles;
static int maxcycles=0;
#define MAXCYCLES 200
static setword inducedcycle[MAXCYCLES];
static int cyclecount;

int GENG_MAIN(int argc, char *argv[]);  /* geng main() */
static int min1,min12,max34,max4; /* bounds on degree counts on geng output */

static counter gengout=0, genggood=0;
static counter vcolgnontriv=0,vcolgout=0;
static counter multignontriv=0,multigout=0;
static long maxvgroup,maxegroup;

static boolean uswitch;  /* suppress output */
static boolean verbose;  /* print more information to stderr */
static int outlevel;  /* 1 = geng only, 2 = geng+vcolg,
                       3 = geng+vcolg+multig, 4 = everything */
static boolean smiles;  /* output in SMILES format */
static int maxbond;  /* maximum mult -1 of bonds (1 if -t, else 2) */

static boolean planar;  /* Molecules must be planar */
static int maxcoord;  /* Maximum coordination number allowed */        //UNUSED
static boolean xswitch;  /* Undocumented, used for development */

/* In the following, the counts are only meaningful if the
   corresponding boolean is true. */
static boolean tswitch;
static long min3rings,max3rings;  /* number of rings of length 3 */
static int count3ring[MAXN+1];
static boolean fswitch;
static long min4rings,max4rings;  /* number of rings of length 4 */
static int count4ring[MAXN+1];
static boolean pswitch;
static long min5rings,max5rings;  /* number of rings of length 5 */
static int count5ring[MAXN+1];
static boolean bipartite;

/* The following is only used if bad9 is selected */
static setword cycle34verts[MAXN+1];  /* set of vertices on rings
                                  of length 3 or 4 */

static long vgroupsize;  /* vcolg group size */
static size_t vgroupalloc=0;  /* Space allocated for vertex groups */
static int *vgroup=NULL;  /* Store of vertex group */
static long vgroupcount;

static long egroupsize;  /* multig group size */
static size_t egroupalloc=0;  /* Space allocated for edge groups */
static int *egroup=NULL;  /* Store of edge group */

typedef struct edge
{
    int x,y;         /* The atoms with this bond */
    int maxmult;     /* This edge can have multiplicity 1..maxmult+1 */
    int allenemate1,allenemate2;
        /* Previous bond that cannot be multiple at the same time as this */
    setword xy;
} edgetype;
static int numedges;
static edgetype edge[MAXNE];
static int edgenumber[MAXN][MAXN]; 
static int deg[MAXN];  /* Simple graph degree */

static FILE *outfile;
#ifdef ZLIB
static gzFile gzoutfile;
#endif
static boolean gzip;

/* Macros for appending to a string using pointer p */
#define PUTINT(xx) { unsigned long ul = (xx); char *sp,s[15]; \
 if (ul == 0) *(p++) = '0'; \
 else { sp = s; while (ul) { *(sp++) = (ul % 10) + '0'; ul /= 10; } \
        while (sp > s) { *(p++) = *(--sp); } }}
#define SPC *(p++) = ' '
#define PUTSTR(xx) { char *sp = (xx); \
   while (*sp != '\0') *(p++) = *(sp++); }
#define PUTBND(xx) { int bnd = (xx); if (bnd == 0) *(p++) = '-'; \
   else if (bnd == 1) *(p++) = '='; else *(p++) = '#'; }

/******************************************************************/

#define MAXELEMENTS 30
static struct elementstruct
{
    char *inputname,*name;
    boolean organic;  /* Belongs to organic subset */
    int valence;
    int lowervalence;  /* Next lower valence, or 0 if none */
    int maxcoord; /* Maximum number of distinct neighbours including H */
    int index;   /* Used in -T style outputs */
} element[MAXELEMENTS] =
{   /* The order of listing does not matter.
       Other elements can be added to the table by following the same
       pattern, up to any limit. All extra elements must be marked
       as non-organic.  The inputname field must have the form X or Xx
       and must be unique. */
  { "C", "C", TRUE,  4,0,4,  0 },
  { "N", "N", TRUE,  3,0,3,  1 },
  { "Nx","N", TRUE,  5,3,4, 10 },
  { "P", "P", TRUE,  3,0,3,  3 },
  { "Px","P", TRUE,  5,3,5, 13 },
  { "B", "B", TRUE,  3,0,3,  9 },
  { "O", "O", TRUE,  2,0,2,  2 },
  { "S", "S", TRUE,  2,0,2,  4 },
  { "Sx","S", TRUE,  4,2,4, 11 },
  { "Sy","S", TRUE,  6,4,6, 12 },
  { "H", "H", FALSE, 1,0,1, 99 },
  { "F", "F", TRUE,  1,0,1,  5 },
  { "Cl","Cl",TRUE,  1,0,1,  6 },
  { "Br","Br",TRUE,  1,0,1,  7 },
  { "I", "I", TRUE,  1,0,1,  8 },
  { "Si","Si",FALSE, 4,0,4, 14 },
  { NULL,NULL,FALSE, 0,0,0, 0 }
};
static int numelements;  /* Actual number in the table */
static int maxindex; /* Max value of the index field except 99 */
#define ISHYDROGEN(i) (strcmp(element[i].name,"H") == 0)

/******************************************************************/

static int
elementindex(char *inputname)
/* Index into element[] of element with this element input name,
   or -1 if it isn't present. */
{
    int i;

    for (i = 0; i < numelements; ++i)
        if (strcmp(element[i].inputname,inputname) == 0)
            break;

    if (i < numelements) return i;
    else
    {
        fprintf(stderr,">E Unknown element %s\n",inputname);
        exit(1);
    }
}

static void
addelement(char *inputname, char *name, int valence, int maxcoord)
/* Add an element to element[]. inputname and name must be of the
form A or Aa. If name==NULL then name is the same as inputname.  If
maxcoord==0 or maxcoord > valence then maxcoord is the same as valence. */
{
    int i;

    if (numelements == MAXELEMENTS)
        gt_abort(">E increase MAXELEMENTS\n");

    if (name == NULL) name = inputname;
    if (maxcoord == 0 || maxcoord > valence) maxcoord = valence;

    if (inputname == NULL || valence <= 0)
        gt_abort(">E invalid parameters in addelement()\n");

    for (i = 0; i < numelements; ++i)
    {
        if (strcmp(element[i].inputname,inputname) == 0)
        {
            fprintf(stderr,">E element %s is already present\n",name);
            exit(1);
        }
    }

    if (!isupper(inputname[0]) || (inputname[1] != '\0' 
             && (!islower(inputname[1]) || inputname[2] != '\0')))
        gt_abort(">E element names must have form A or Aa\n");
    if (!isupper(name[0]) || (name[1] != '\0' 
             && (!islower(name[1]) || name[2] != '\0')))
        gt_abort(">E element names must have form A or Aa\n");

    if (3*maxcoord < valence || valence < 0)
        gt_abort(">E impossible maxcoord/valence\n");

    element[numelements].inputname = strdup(inputname);
    element[numelements].name = strdup(name);
    element[numelements].valence = valence;
        /* lowervalence is only used for SMILES output of
           atoms in the organic subset and since those are
           in the table already we don't need a value to
           be set here. */
    element[numelements].lowervalence = 0;
    element[numelements].maxcoord = maxcoord;
    element[numelements].organic = FALSE;
    element[numelements].index = ++maxindex;

    if (verbose)
	fprintf(stderr,"Added element %s input %s, valence=%d maxcoord=%d\n",
		element[numelements].name,element[numelements].inputname,
                element[numelements].valence,element[numelements].maxcoord);
    ++numelements;
}

#undef HELPTEXT2
#ifdef SURGEPLUGIN
/* Note that the value of SURGEPLUGIN must be a filename in quotes, so on
   the compilation command you will probably need something like
   -DSURGEPLUGIN='"myplugin.h"' . */
#include SURGEPLUGIN
#endif

#ifndef HELPTEXT2
#define HELPTEXT2 ""
#endif

/******************************************************************/

#if 0
static unsigned long long
molcode(int *vcol, int *mult, int n, int ne)
/* This makes an isomorph-invariant code for a molecule.  It is intended
for debugging and can't be relied on to distinguish between different
molecules. It only works up to MAXN/2 atoms (but you can compile the
program with WORDSIZE=64 if you really need to). */
{
    graph g[MAXN],h[MAXN];
    int lab[MAXN],ptn[MAXN],orbits[MAXN],weight[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[2*MAXN];
    int i,x,y,nv;
    unsigned long long ans1,ans2;

    nv = 2*n;
    if (nv >= MAXN) gt_abort(">E surge : too big for longcode\n");

    for (i = 0; i < n; ++i) { weight[i] = vcol[i]; weight[n+i] = n+vcol[i]; }
    setlabptn(weight,lab,ptn,nv);

    for (i = 0; i < nv; ++i) g[i] = 0;

    for (i = 0; i < n; ++i) { g[i] |= bit[n+i]; g[n+i] |= bit[i]; }

    for (i = 0; i < ne; ++i)
    {
        x = edge[i].x;
        y = edge[i].y;
        g[x] |= bit[y]; g[y] |= bit[x];
        if (mult[i] > 0) { g[n+x] |= bit[n+y]; g[n+y] |= bit[n+x]; }
        if (mult[i] > 1) 
        {
            g[x] |= bit[y+n]; g[y+n] |= bit[x];
            g[x+n] |= bit[y]; g[y] |= bit[x+n];
        }
    }

    options.defaultptn = FALSE;
    options.getcanon = TRUE;

    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,2*MAXN,1,nv,h);
 
    ans1 = n;
    ans2 = ne;
    for (i = 0; i < nv; ++i)
    {
        ans1 = 177*ans1 + (h[i] >> (WORDSIZE-nv));
        ans2 = 1237*ans2 + (h[i] >> (WORDSIZE-nv));
    }

    return ans1^ans2;
}
#endif

/******************************************************************/

void dummy()
{
}

/******************************************************************/

static boolean
isplanar(graph *g, int n)
/* Check if g is planar, assuming g is connected */
{
    t_ver_sparse_rep V[MAXN];
    t_adjl_sparse_rep A[2*MAXNE+1];
    t_dlcl **dfs_tree,**back_edges,**mult_edges;
    t_ver_edge *embed_graph;
    int i,j,k,pop,nv,ne,c;
    int edge_pos,v,w;
    setword ww;
    boolean ans;
    graph h[MAXN];
    int newlab[MAXN];
    setword queue;

    queue = 0;
    ne = 0;
    for (i = 0; i < n; ++i)
    {
        h[i] = g[i];
        pop = POPCOUNT(h[i]);
        ne += pop;
        if (pop <= 2) queue |= bit[i];
    }
    ne /= 2;
    nv = n;

    while (queue && ne >= nv+3)
    {
        TAKEBIT(i,queue);
        pop = POPCOUNT(h[i]);
        if (pop == 1)  /* i--j with deg(i)=1 */
        {
            j = FIRSTBITNZ(h[i]);
            h[i] = 0;
            h[j] &= ~bit[i];
            --nv;
            --ne;
            queue |= bit[j];
        }
        else if (pop == 2)  /* j--i--k with deg(i)=2 */
        {
            j = FIRSTBITNZ(h[i]);
            k = FIRSTBITNZ(h[i] & ~bit[j]);
            h[i] = 0;
            h[j] &= ~bit[i];
            h[k] &= ~bit[i];
            --nv;
            if ((h[j] & bit[k]))
            {
                ne -= 2;
                queue |= (bit[j] | bit[k]);
            }
            else
            {
                --ne;
                h[j] |= bit[k];
                h[k] |= bit[j];
            }
        }
    }

    if (ne <= nv + 2) return TRUE;
    if (nv == 5 && ne <= 9) return TRUE;

    nv = 0;
    for (i = 0; i < n; ++i) if (h[i] != 0) newlab[i] = nv++;

    k = 0;
    for (i = 0; i < n; ++i)
    if (h[i] != 0)
    {
        V[newlab[i]].first_edge = k;
        ww = h[i];
        while (ww)
        {
            TAKEBIT(j,ww);
            A[k].end_vertex = newlab[j];
            A[k].next = k+1;
            ++k;
        }
        A[k-1].next = NIL;
    }

    ne = k/2;

    ans = sparseg_adjl_is_planar(V,nv,A,&c,&dfs_tree,&back_edges,
             &mult_edges,&embed_graph,&edge_pos,&v,&w);

    sparseg_dlcl_delete(dfs_tree,nv);
    sparseg_dlcl_delete(back_edges,nv);
    sparseg_dlcl_delete(mult_edges,nv);
    embedg_VES_delete(embed_graph,nv);

    return ans;
} 

/******************************************************************/

static void
SMILESoutput(int *vcol, int n, int *hyd, int *mult, int ne)
/* Write molecules in SMILES format */
{
    char *p,line[20*MAXNE];
    int i,x,y,r,m;
    const struct elementstruct *thiselement;

    p = line;

    for (i = 0; i < smileslen; ++i)
    {
        x = smilesskeleton[i].x;
        y = smilesskeleton[i].y;
        switch(smilesskeleton[i].item)
        {
         case SM_ATOM :
            thiselement = &element[vcol[x]];
#if 0
 /* This version seems to meet the OpenSmiles standard, but some
    problems with obabel have not been tracked down. */
            if (!thiselement->organic ||
                    thiselement->valence - hyd[x] <= thiselement->lowervalence)
#else
 /* This version gives explicit H for atoms in higher valences even
    if they are in the organic subset. */
            if (!thiselement->organic ||
                 (thiselement->lowervalence > 0 && hyd[x] > 0))
#endif
            {
                *(p++) = '[';
                PUTSTR(thiselement->name);
                if (hyd[x] > 0) *(p++) = 'H';
                if (hyd[x] > 1) PUTINT(hyd[x])
                *(p++) = ']';
            }
            else 
                PUTSTR(thiselement->name);
            break;
         case SM_BOND :
            m = mult[edgenumber[x][y]];
            if      (m == 1) *(p++) = '=';
            else if (m == 2) *(p++) = '#';
            break;
         case SM_OPEN :
            *(p++) = '(';
            break;
         case SM_CLOSE :
            *(p++) = ')';
            break;
         case SM_RING0 :
            m = mult[edgenumber[x][y]];
            if      (m == 1) *(p++) = '=';
            else if (m == 2) *(p++) = '#';
            r = smilesskeleton[i].r;
            if (r < 10)
                *(p++) = '0' + r;
            else
            {
                *(p++) = '%';
                *(p++) = '0' + r/10;
                *(p++) = '0' + r%10;
            }
            break;
         case SM_RING1 :
            r = smilesskeleton[i].r;
            if (r < 10)
                *(p++) = '0' + r;
            else
            {
                *(p++) = '%';
                *(p++) = '0' + r/10;
                *(p++) = '0' + r%10;
            }
            break;
        }
    }

    *(p++) = '\n';
    *p = '\0';

#ifdef ZLIB
    if (gzip)
    {
        if (gzputs(gzoutfile,line) < 0)
            gt_abort(">E surge : zlib output error\n");
        return;
    }
#endif

    if (fputs(line,outfile) == EOF) gt_abort(">E surge : output error\n");
}

/******************************************************************/

static void
SDFformat(int *vcol, int n, int *hyd, int *mult, int ne)
/* Write molecules in SDF format */
{
    int i,j;
    
#ifdef ZLIB
    if (gzip)
    {
        gzprintf(gzoutfile,"\nSurge 1.0\n\n");
        gzprintf(gzoutfile,"%3d%3d  0  0  0  0            999 V2000\n",n,ne);

        for (i = 0; i < n; ++i)
            gzprintf(gzoutfile,"    0.0000    0.0000    0.0000 %-2s"
                   "  0  0  0  0  0%3d  0  0  0  0  0  0\n",
                   element[vcol[i]].name,element[vcol[i]].valence);
    
        for (i = 0; i < ne; ++i)
            gzprintf(gzoutfile,"%3d%3d%3d  0  0  0  0\n",
                edge[i].x+1,edge[i].y+1,mult[i]+1);

        gzprintf(gzoutfile,"M  END\n$$$$\n");

        return;
    }
#endif

    fprintf(outfile,"\nSurge 1.0\n\n");
    fprintf(outfile,"%3d%3d  0  0  0  0            999 V2000\n",n,ne);

    for (i = 0; i < n; ++i)
        fprintf(outfile,"    0.0000    0.0000    0.0000 %-2s"
               "  0  0  0  0  0%3d  0  0  0  0  0  0\n",
               element[vcol[i]].name,element[vcol[i]].valence);

    for (i = 0; i < ne; ++i)
        fprintf(outfile,"%3d%3d%3d  0  0  0  0\n",
            edge[i].x+1,edge[i].y+1,mult[i]+1);

    fprintf(outfile,"M  END\n$$$$\n");
}

/****************************************************************/

static void
multigoutput(int *vcol, int n, int *mult, int ne)
/* Write output equal to the multig -T format. */
{
    char line[10+60*MAXN],*p;
    int i;

    p = line;
    PUTINT(n); SPC; PUTINT(ne);
    for (i = 0; i < n; ++i)
    {
        SPC;
        PUTINT(element[vcol[i]].index);
    }
    SPC;
    for (i = 0; i < ne; ++i)
    {
        SPC; PUTINT(edge[i].x);
        SPC; PUTINT(edge[i].y);
        SPC; PUTINT(mult[i]);
    }
    *(p++) = '\n';
    *p = '\0';

#ifdef ZLIB
    if (gzip)
    {
        if (gzputs(gzoutfile,line) < 0)
            gt_abort(">E surge : zlib output error\n");
        return;
    }
#endif

    if (fputs(line,outfile) == EOF) gt_abort(">E surge : output error\n");
}

/****************************************************************/

static void
alphabeticoutput(int *vcol, int n, int *mult, int ne)
/* Write alphabetic output */
{
    char line[10+60*MAXN],*p;
    int i,xx,yy;

    p = line;
    PUTINT(n);
    SPC;
    PUTINT(ne);
    SPC;
    PUTSTR(canonform);

    for (i = 0; i < ne; ++i)
    {
        SPC;
        xx = newlabel[edge[i].x];
        yy = newlabel[edge[i].y];
        if (xx < yy)
        {
            PUTINT(xx); PUTBND(mult[i]); PUTINT(yy);
        }
        else
        {
            PUTINT(yy); PUTBND(mult[i]); PUTINT(xx);
        }
    }
    *(p++) = '\n';
    *p = '\0';

#ifdef ZLIB
    if (gzip)
    {
        if (gzputs(gzoutfile,line) < 0)
            gt_abort(">E surge : zlib output error\n");
        return;
    }
#endif

    if (fputs(line,outfile) == EOF) gt_abort(">E surge : output error\n");
}

/******************************************************************/

static void
gotone(int *vcol, int n, int *hyd, int *mult, int ne, int level)
/* Now we have a completed molecule.
   deg[0..n-1] is the simple graph degrees
   hyd[0..n-1] is the number of implicit hydrogens
*/
{
    int i,val;

    for (i = level; i < ne; ++i) mult[i] = 0;

    if (needcoordtest)
        for (i = 0; i < n; ++i)
        {
            if (deg[i] + hyd[i] > element[vcol[i]].maxcoord) return;
            if (deg[i] + hyd[i] > 4 && hyd[i] > 0) return;
        }

#ifdef SURGEPLUGIN_STEP3
    SURGEPLUGIN_STEP3
#endif

    ++multigout;
    if (uswitch) return;

    if (outlevel == 3)
        multigoutput(vcol,n,mult,ne);
    else if (alphabetic)
        alphabeticoutput(vcol,n,mult,ne);
    else if (smiles)
        SMILESoutput(vcol,n,hyd,mult,ne);
    else
        SDFformat(vcol,n,hyd,mult,ne);
}

/******************************************************************/

static boolean
testemax(int *mult, int ne, int level)
/* Test if edge colouring is maximum wrt group. */
{
    int *gp,i,j;
    long kgp;

    for (i = level; i < ne; ++i) mult[i] = 0;

    /* kgp really starts at 1 on the next line */
    for (kgp = 1, gp = egroup; kgp < egroupsize; ++kgp, gp += ne)
    {
        for (i = 0; i < ne; ++i)
        {
            j = gp[i];
            if      (mult[j] > mult[i]) return FALSE;
            else if (mult[j] < mult[i]) break;
        }
    }

    return TRUE;
}

/***************************************************************************/

static void
/* Recursive scan for multiplying edges 
   Version which checks allenemate for -B3 */
escan2(int level, int needed, int *vcol, int *hyd, int *prev, int n, int *mult, int ne)
{
    int lev,maxlev,k,max,x,y;

    if (needed == 0)
    {
        if (egroupsize > 1 && !testemax(mult,ne,level))
            return;
        gotone(vcol,n,hyd,mult,ne,level);
        return;
    }
    else
    {
        maxlev = ne + 1 - (needed+maxbond-1)/maxbond;
        for (lev = level; lev < maxlev; ++lev)  
        {
            x = edge[lev].x;
            y = edge[lev].y;
            max = edge[lev].maxmult;

            if (needed < max) max = needed;
            if (hyd[x] < max) max = hyd[x];
            if (hyd[y] < max) max = hyd[y];
            if (prev[lev] >= 0 && mult[prev[lev]] < max) max = mult[prev[lev]];
            if (edge[lev].allenemate1 >= 0 && mult[edge[lev].allenemate1] >= 1)
                 max = 0;
            if (edge[lev].allenemate2 >= 0 && mult[edge[lev].allenemate2] >= 1)
                 max = 0;

            for (k = 1; k <= max; ++k)
            {
                mult[lev] = k;
                hyd[x] -= k;
                hyd[y] -= k;
                escan2(lev+1,needed-k,vcol,hyd,prev,n,mult,ne);
                hyd[x] += k;
                hyd[y] += k;
            }

            mult[lev] = 0;
        }
    }

    return;
}

/***************************************************************************/

static void
/* Recursive scan for multiplying edges */
escan(int level, int needed, int *vcol, int *hyd,
      int *prev, int n, int *mult, int ne)
{
    int lev,maxlev,k,max,x,y;

    if (needed == 0)
    {
        if (egroupsize > 1 && !testemax(mult,ne,level))
            return;
        gotone(vcol,n,hyd,mult,ne,level);
        return;
    }
    else
    {
        maxlev = ne + 1 - (needed+maxbond-1)/maxbond;
        for (lev = level; lev < maxlev; ++lev)  
        {
            x = edge[lev].x;
            y = edge[lev].y;
            max = edge[lev].maxmult;

            if (needed < max) max = needed;
            if (hyd[x] < max) max = hyd[x];
            if (hyd[y] < max) max = hyd[y];
            if (prev[lev] >= 0 && mult[prev[lev]] < max)
                max = mult[prev[lev]];

            for (k = 1; k <= max; ++k)
            {
                mult[lev] = k;
                hyd[x] -= k;
                hyd[y] -= k;
                escan(lev+1,needed-k,vcol,hyd,prev,n,mult,ne);
                hyd[x] += k;
                hyd[y] += k;
            }

            mult[lev] = 0;
        }
    }

    return;
}

/******************************************************************/

static void
findegroup(int *vcol, int n, int ne)
/* Set egroup to the set of vertex-colour-preserving vgroup elements */
{
    int *vgp,*egp,i,j,kgp;

    egp = egroup;

    /* kgp really starts at 1 on the next line */
    for (kgp = 1, vgp = vgroup; kgp < vgroupsize; ++kgp, vgp += n)
    {
        for (i = 0; i < n; ++i)
            if (vcol[vgp[i]] != vcol[i]) break;
        if (i == n)
            for (j = 0; j < ne; ++j)
                *(egp++) = edgenumber[vgp[edge[j].x]][vgp[edge[j].y]];
    }

    egroupsize = 1 + (egp - egroup) / ne;
}

/****************************************************************/

static void
colouredges(graph *g, int *vcolindex, int n)
/* This procedure receives graphs from the vcolg phase and
   colours the edges. */
{
    int hyd[MAXN];  /* Remaining degree of vertex */
    int i,j,k,ne;
    int mult[MAXNE];
    int prev[MAXNE]; /* If >= 0, earlier point that must have greater colour */
    int needed;  /* Extra edges needed */
    int iter[FORMULALEN];
    int vcol[MAXN];  /* index into element[] */

    ne = numedges;

    /* hyd[i] starts with the number of hydrogens needed if all bonds are
       single and is reduced as bonds are multiplied */

    needcoordtest = FALSE;
    for (i = 0; i < n; ++i)
    {
        vcol[i] = elementtype[vcolindex[i]];
        hyd[i] = element[vcol[i]].valence - deg[i];
        if (element[vcol[i]].valence > element[vcol[i]].maxcoord 
                    || (element[vcol[i]].valence > 4 && hyd[i] > 0))
            needcoordtest = TRUE;
    }

#ifdef SURGEPLUGIN_STEP2
    SURGEPLUGIN_STEP2
#endif

    needed = (valencesum - hydrogens)/2 - ne;  /* Extra edges needed */

    if (ne == 0 && needed > 0) return;

    if (alphabetic)
    {
        iter[0] = 0;
        for (i = 1; i < numtypes; ++i) iter[i] = iter[i-1] + elementcount[i-1];

        for (i = 0; i < n; ++i) newlabel[i] = iter[vcolindex[i]]++;
    }

    if (needed == 0)
    {
        gotone(vcol,n,hyd,mult,ne,0);
        return;
    }

    for (i = 0; i < ne; ++i) prev[i] = -1;

    for (i = 0; i < n; ++i)
    {
        if (deg[i] != 1) continue;
        /* Find most recent equivalent j */
        for (j = i; --j >= 0; )
            if (g[i] == g[j] && vcol[j] == vcol[i])
                break;

        if (j >= 0)
        {
            k = FIRSTBITNZ(g[i]);
            prev[edgenumber[i][k]] = edgenumber[j][k];
        }
    }

    if (vgroupsize == 1) 
        egroupsize = 1;
    else
    {
        if (egroupalloc < vgroupsize*ne)
        {
            if (egroup) free(egroup);
            if ((egroup = malloc((vgroupsize+48)*ne*sizeof(int))) == NULL)
                gt_abort(">E surge : Can't allocate space for egroup\n");
            egroupalloc = (vgroupsize+48) * ne;
        }
        findegroup(vcol,n,ne);
        if (vgroupsize % egroupsize != 0) gt_abort(">E egroup error\n");
    }

    if (egroupsize == 1 && needed == 1)
    {
        for (i = 0; i < ne; ++i) mult[i] = 0;
        for (i = 0; i < ne; ++i)
            if (prev[i] < 0 && edge[i].maxmult >= 1
                && hyd[edge[i].x] > 0 && hyd[edge[i].y] > 0)
            {
                mult[i] = 1;
                --hyd[edge[i].x];
                --hyd[edge[i].y];
                gotone(vcol,n,hyd,mult,ne,ne);
                ++hyd[edge[i].x];
                ++hyd[edge[i].y];
                mult[i] = 0;
            }
        return;
    }

    if (egroupsize != 1) ++multignontriv;
    if (egroupsize > maxegroup) maxegroup = egroupsize;

    if (bad5 || bad6) escan2(0,needed,vcol,hyd,prev,n,mult,ne);
    else              escan(0,needed,vcol,hyd,prev,n,mult,ne);
}

/******************************************************************/

static void
vcolgoutput(graph *g, int *vcolindex, int n)
/* Write output equal to the vcolg -T format. */
{
    char line[10+30*MAXN],*p;
    int i,j,ne;
    setword w;

    ne = 0; 
    for (i = 0; i < n; ++i) ne += POPCOUNT(g[i]);
    ne /= 2;

    p = line;
    PUTINT(n); SPC; PUTINT(ne);
    for (i = 0; i < n; ++i)
    {
        SPC;
        PUTINT(element[elementtype[vcolindex[i]]].index);
    }
    SPC;
    for (i = 0; i < n; ++i)
    {
        w = g[i] & BITMASK(i);
        while (w)
        {
            TAKEBIT(j,w);
            SPC; PUTINT(i); SPC; PUTINT(j);
        }
    }
    *(p++) = '\n';
    *p = '\0';

#ifdef ZLIB
    if (gzip)
    {
        if (gzputs(gzoutfile,line) < 0)
             gt_abort(">E surge : zlib output error\n");
        return;
    }
#endif

    if (fputs(line,outfile) == EOF) gt_abort(">E surge : output error\n");
}

/******************************************************************/

static boolean
testvmax(int *colindex, int n)
/* Test if vertex colouring is maximum wrt group. If so, return group.
   If not, return a safe level to return to. */
{
    int *gp,i,j;
    long kgp;

     /* kgp really starts at 1 on the next line */
    for (kgp = 1, gp = vgroup; kgp < vgroupsize; ++kgp, gp += n)
    {
        for (i = 0; i < n; ++i)
        {
            j = gp[i];
            if      (colindex[j] > colindex[i]) return FALSE;
            else if (colindex[j] < colindex[i]) break;
        }
    }

    return TRUE;
}
 
/**********************************************************************/

static void
vscan(int level, int *colindex, graph *g, int *prev,
                               int *maxcolindex, int *remain, int n)
/* Recursive vertex colour scan */
{
    int k,max;

    if (level == n)
    {
        if (vgroupsize == 1 || testvmax(colindex,n))
        {
            ++vcolgout;
            if (outlevel == 2)
            {
                if (!uswitch) vcolgoutput(g,colindex,n);
            }
            else
                colouredges(g,colindex,n);
        }
        return;
    }

    max = maxcolindex[level];
    if (prev[level] >= 0 && colindex[prev[level]] < max)
        max = colindex[prev[level]];

    for (k = 0; k <= max; ++k)   
    {
        if (remain[k] == 0) continue;
        colindex[level] = k;
        --remain[k];
        vscan(level+1,colindex,g,prev,maxcolindex,remain,n);
        ++remain[k];
    }
}

/**********************************************************************/

static void
storevgroup(int *p, int n)
/* Called by allgroup; store full group at vcolg phase */
{
    int *gp,i;

    if (vgroupcount == 0)
    {
        vgroupcount = 1;   /* Don't store identity */
        return;
    }

    gp = vgroup + n * (vgroupcount-1);
    for (i = 0; i < n; ++i) gp[i] = p[i];

    ++vgroupcount;
}

/**********************************************************************/

static void
colourvertices(graph *g, int n)
/* Main procedure for vcolg phase */
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[2*MAXN];
    grouprec *group;
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int prev[MAXN]; /* If >= 0, earlier point that must have greater colour */
    int weight[MAXN];
    int maxcolindex[MAXN];  /* Max colour index for each vertex */
    int remain[FORMULALEN];  /* Remaining number of colours for each vertex */
    int vcolindex[MAXN];  /* Index into elementtype[] */
    int i,j;

#ifdef SURGEPLUGIN_STEP1
    SURGEPLUGIN_STEP1
#endif

    for (i = 0; i < n; ++i)
    {
        prev[i] = -1;
        weight[i] = n*POPCOUNT(g[i]);
    }

    for (i = 0; i < n; ++i)
    {
        if (POPCOUNT(g[i]) != 1) continue;
        /* Find most recent equivalent j */
        for (j = i; --j >= 0; )
            if (g[j] == g[i]) break;

        if (j >= 0)
        {
            prev[i] = j;
            weight[i] = weight[j] + 1;
        }
    }

    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;
    options.defaultptn = FALSE;

    setlabptn(weight,lab,ptn,n);

    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,2*MAXN,1,n,NULL);

    if (stats.grpsize2 > 0 || stats.grpsize1 > 1e7)
        gt_abort(">E surge : vgroup size greater than 10^7 encountered\n");
    vgroupsize = stats.grpsize1 + 0.01;

    for (i = 0; i < numtypes; ++i) remain[i] = elementcount[i];
    for (i = 0; i < n; ++i) maxcolindex[i] = maxtype[deg[i]];

    if (vgroupsize == 1)  /* Trivial group */
    {
        vscan(0,vcolindex,g,prev,maxcolindex,remain,n);
        return;
    }

    ++vcolgnontriv;
    if (vgroupsize > maxvgroup) maxvgroup = vgroupsize;

    group = groupptr(FALSE);
    makecosetreps(group);

    if (vgroupalloc < (vgroupsize-1) * n)
    {
        if (vgroup) free(vgroup);
        if ((vgroup = malloc((vgroupsize+48)*n*sizeof(int))) == NULL)
            gt_abort(">E surge : Can't allocate space for vgroup\n");
        vgroupalloc = (vgroupsize+48) * n;
    }

    vgroupcount = 0;
    allgroup(group,storevgroup);

    if (vgroupcount != vgroupsize) gt_abort(">E surge : vgroup error\n");

   /* Check the logic of the next section. What about maxdeg? */
    if (numtypes == 1)
    {
        for (i = 0; i < n; ++i) vcolindex[i] = 0;
        ++vcolgout;
        colouredges(g,vcolindex,n);
        return;
    }

    j = n;      /* Can choose a better orbit? */
    for (i = 0; i < n; ++i)
        if (orbits[i] < i && orbits[i] < j) j = orbits[i];
    for (i = j + 1; i < n; ++i)
        if (orbits[i] == j) prev[i] = j;

    vscan(0,vcolindex,g,prev,maxcolindex,remain,n);
}

/******************************************************************/

extern int geng_maxe;

int
surgepreprune(graph *g, int n, int maxn)
/* This function is called by the PREPRUNE service of geng.
 It does -B9 which is more efficient here. If OLDGENG is 
 defined, it also speeds up the generation of connected
 graphs.  (For new geng, that function is built in.) */
{
    setword notvisited,queue;
    setword c34,w,ww,cyc;
    int ne,nc,i;

    if (bad9)    /* cycle34verts */
    {
        if (n <= 2)
            cycle34verts[n] = 0;
        else
        {
            c34 = cycle34verts[n-1];
            if (!tswitch)
            {
                w = g[n-1];
                while (w)
                {
                    TAKEBIT(i,w);
                    ww = g[i] & w;
                    if (POPCOUNT(ww) > 1) return 1;
                    if ((ww))
                    {
                        cyc = bit[n-1] | bit[i] | ww;
                        if ((c34 & cyc)) return 1;
                        c34 |= cyc;
                    }
                }
            }
            if (!fswitch)
            {
                for (i = n-1; --i >= 0;)
                {
                    w = g[i] & g[n-1];
                    if (POPCOUNT(w) > 2) return 1;
                    if (POPCOUNT(w) == 2)
                    {
                        cyc = bit[n-1] | bit[i] | w;
                        if ((c34 & cyc)) return 1;
                        c34 |= cyc;
                    }
                }
            }
            cycle34verts[n] = c34;
        }
    }

#ifdef OLDGENG
    if (n == maxn || geng_maxe - maxn >= 5) return 0;

    ne = 0;
    for (i = 0; i < n; ++i) ne += POPCOUNT(g[i]);
    ne /= 2;

    nc = 0;
    notvisited = ALLMASK(n);

    while (notvisited)
    {
        ++nc;
        queue = SWHIBIT(notvisited);
        notvisited &= ~queue;
        while (queue)
        {
            TAKEBIT(i,queue);
            notvisited &= ~bit[i];
            queue |= g[i] & notvisited;
        }
    }

    if (ne - n + nc > geng_maxe - maxn + 1) return 1;
#endif

    return 0;
}

int
surgeprune(graph *g, int n, int nmax)
/* This is a procedure that geng will call at each level
using the PRUNE service.
The options -t, -f, -p, -B7,8 are implemented here by
incrementally updating the required counts. */
{
    setword w,ax,ay,gx,gy,gxy,gxya,gi,gj;
    int i,j,x,y,k,a,b,extra;
    int i1,i2,i3,i4,d1,d2,d3,d4;
    int v[MAXN];

    if (tswitch)
    {
        if (n <= 2)
            count3ring[n] = 0;
        else
        {
            extra = 0;
            w = g[n-1];
            while (w)
            {
                TAKEBIT(i,w);
                extra += POPCOUNT(g[i]&w);
            }
            count3ring[n] = count3ring[n-1] + extra;
            if (count3ring[n] > max3rings) return 1;
        }
        if (n == nmax && count3ring[n] < min3rings)
            return 1;
    }

    if (fswitch)
    {
        if (n <= 3)
            count4ring[n] = 0;
        else
        {
            extra = 0;
            for (i = 0; i < n-1; ++i)
            {
                k = POPCOUNT(g[i]&g[n-1]);
                extra += k*k - k;
            }
            count4ring[n] = count4ring[n-1] + extra/2;
            if (count4ring[n] > max4rings) return 1;
        }
        if (n == nmax && count4ring[n] < min4rings)
            return 1;
    }

    if (pswitch)
    {
        if (n <= 4)
            count5ring[n] = 0;
        else
        {
            extra = 0;
            for (y = 1; y < n-1; ++y)
            {
                w = g[y] & ~BITMASK(y);
                while (w)
                {
                    TAKEBIT(x,w);
                    ax = (g[x] & g[n-1]) & ~bit[y];
                    ay = (g[y] & g[n-1]) & ~bit[x];
                    extra += POPCOUNT(ax)*POPCOUNT(ay) - POPCOUNT(ax&ay);
                }
            }
            count5ring[n] = count5ring[n-1] + extra;
            if (count5ring[n] > max5rings) return 1;
        }
        if (n == nmax && count5ring[n] < min5rings)
            return 1;
    }

    if (bad7 && n >= 6)
    {
          /* For K_33 we can assume that vertex n-1 is included */
        for (x = n-1; --x >= 1;)
        if (POPCOUNT(g[x]&g[n-1]) >= 3)
        {
             for (y = x; --y >= 0;)
                if (POPCOUNT(g[y]&g[x]&g[n-1]) >= 3) return 1;
        }
        
          /* But for K_24 we can't */
        for (x = n; --x >= 1;)
        if (POPCOUNT(g[x]) >= 4)
        {
            for (y = x; --y >= 0;)
                if (POPCOUNT(g[x]&g[y]) >= 4) return 1;
        }
    }

    if (bad8)
    {
       /* cone of P4 */
        if (n >= 5)
        {
            for (x = n; --x >= 0;)
            if (POPCOUNT(g[x]) == 4)
            {
              /* Think of a better way to do this */
                w = gx = g[x];
                TAKEBIT(i1,w);
                TAKEBIT(i2,w);
                TAKEBIT(i3,w);
                i4 = FIRSTBITNZ(w);
                d1 = POPCOUNT(g[i1]&gx);
                d2 = POPCOUNT(g[i2]&gx);
                d3 = POPCOUNT(g[i3]&gx);
                d4 = POPCOUNT(g[i4]&gx);
                if (d1 > 0 && d2 > 0 && d3 > 0 && d4 > 0
                   && ((d1 >= 2)+(d2 >= 2)+(d3 >= 2)+(d4 >= 2)) >= 2)
                      return 1;
            }
            else if (POPCOUNT(g[x]) > 4)
            {
                w = gx = g[x];
                i = 0;
                while (w)
                {
                    TAKEBIT(y,w);
                    if (POPCOUNT(g[y]&gx) >= 2) v[i++] = y;
                }
                for (--i; i >= 1; --i)
                for (j = 0; j < i; ++j)
                    if ((g[v[i]] & bit[v[j]]) &&
                              POPCOUNT((g[v[i]]|g[v[j]])&gx) >= 4)
                        return 1;
            }
        }

        /* K4 with a path of 3 edges between two of its vertices */

        if (n >= 6)
        {
            for (x = n; --x >= 1;)
            if (POPCOUNT(g[x]) >= 4)    
            {
                for (y = x; --y >= 0;)
                if (POPCOUNT(g[y]) >= 4 && (g[y]&bit[x]))
                {
                    gxy = g[x] & g[y];
                    while (gxy)
                    {
                        TAKEBIT(a,gxy);
                        gxya = gxy & g[a];
                        while (gxya)
                        {
                            TAKEBIT(b,gxya);
                            w = bit[x] | bit[y] | bit[a] | bit[b];
                            gi = g[x] & ~w;
                            gj = g[y] & ~w;
                            while (gi)
                            {
                                TAKEBIT(i,gi);
                                if ((g[i] & gj)) return 1;
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef SURGEPLUGIN_STEP0
    SURGEPLUGIN_STEP0
#endif

    return 0;
}

/******************************************************************/

static void
smilesdfs(graph *g, setword *seen, int v, int par, graph *back,
          struct smilesstruct *smilestemp, int *len)
/* Recursive DFS to collect SMILES information */
{
    setword gv,w;
    int k;
    boolean first;

    gv = g[v];
    first = TRUE;
    *seen |= bit[v];

    while (gv)
    {
        TAKEBIT(k,gv);
        if ((*seen & bit[k]))
        {
            if (k != par)
            {
                back[v] |= bit[k];
                back[k] |= bit[v];
            }
        }
        else
        {
            if (first)
            {
                if (POPCOUNT(g[k]) == 1)
                {
                    if ((w = (gv & ~*seen)))   /* really = */
                    {
                        gv |= bit[k];
                        k = FIRSTBITNZ(w);
                        gv &= ~bit[k];
                    }
                }
                smilesdfs(g,seen,k,v,back,smilestemp,len);
                first = FALSE;
            }
            else
            {
                smilestemp[(*len)++].item = SM_CLOSE;
                smilesdfs(g,seen,k,v,back,smilestemp,len);
                smilestemp[(*len)++].item = SM_OPEN;
            }
        }
    }

    smilestemp[*len].item = SM_ATOM;
    smilestemp[*len].x = v;
    ++*len;
    if (par >= 0)
    {
        smilestemp[*len].item = SM_BOND;
        smilestemp[*len].x = par;
        smilestemp[*len].y = v;
        ++*len;
    }
}

static void
makesmilesskeleton(graph *g, int n)
/* Make a skeleton SMILES structure for use in SMILESoutput */
{
    struct smilesstruct smilestemp[4*MAXN+6*MAXNE];
    graph back[MAXN],ring[MAXN];
    setword w,seen;
    int len,ringnumber;
    int i,j,v;

    for (v = n; --v >= 0; ) if (POPCOUNT(g[v]) == 1) break;
    if (v < 0) v = n-1;

    len = 0;
    seen = 0;
    for (i = 0; i < n; ++i) back[i] = 0;
    for (i = 0; i < n; ++i) ring[i] = 0;

    smilesdfs(g,&seen,v,-1,back,smilestemp,&len);

    smileslen = 0;
    ringnumber = 0;
    for (i = len; --i >= 0; )
    {
        smilesskeleton[smileslen++] = smilestemp[i];
        if (smilestemp[i].item == SM_ATOM)
        {
            v = smilestemp[i].x;
            w = ring[v];
            while (w)
            {
                TAKEBIT(j,w);
                smilesskeleton[smileslen].item = SM_RING1; 
                smilesskeleton[smileslen].r = j+1;
                ++smileslen;
            }
            w = back[v];
            while (w)
            {
                TAKEBIT(j,w);
                ++ringnumber;
                smilesskeleton[smileslen].item = SM_RING0; 
                smilesskeleton[smileslen].x = v;
                smilesskeleton[smileslen].y = j;
                smilesskeleton[smileslen].r = ringnumber;
                ++smileslen;
                ring[j] |= bit[ringnumber-1];
                back[j] &= ~bit[v];
            }
        }
    }

#ifdef SMILESSKELETON
  /* This will print the SMILES skeleton to stdout */

    fprintf(stdout,"len=%d",smileslen);
    for (i = 0; i < smileslen; ++i)
    {
        switch (smilesskeleton[i].item)
        {
         case SM_ATOM :
            fprintf(stdout," %d",smilesskeleton[i].x);
            break;
         case SM_BOND : 
            fprintf(stdout," %d-%d",smilesskeleton[i].x,smilesskeleton[i].y);
            break;
         case SM_OPEN :
            fprintf(stdout," (");
            break;
         case SM_CLOSE :
            fprintf(stdout," )");
            break;
         case SM_RING0 :
            fprintf(stdout," R%d:%d-%d",smilesskeleton[i].r,
                           smilesskeleton[i].x,smilesskeleton[i].y);
            break;
         case SM_RING1 :
            fprintf(stdout," R%d",smilesskeleton[i].r);
        }
    }
    fprintf(stdout,"\n");
#endif
}

/******************************************************************/

static void
inducedpaths(graph *g, int origin, int start, setword body,
                                        setword last, setword path)
/* Number of induced paths in g starting at start, extravertices within
 * body and ending in last.
 * {start}, body and last should be disjoint. */
{
    setword gs,w;
    int i;

    gs = g[start];

    w = gs & last;
    while (w)
    {
        TAKEBIT(i,w);
        inducedcycle[cyclecount++]
           = path | bit[edgenumber[start][i]] | bit[edgenumber[origin][i]];
    }

    w = gs & body;
    while (w)
    {
        TAKEBIT(i,w);
        inducedpaths(g,origin,i,body&~gs,last&~bit[i]&~gs,
                                         path|bit[edgenumber[start][i]]);
    }
}

static void
findinducedcycles(graph *g, int n)
/* Find all the induced cycles */
{
    setword body,last,cni;
    int i,j;

#if 0
    body = 0;
    for (i = 0; i < n; ++i)
        if (POPCOUNT(g[i]) > 1) body |= bit[i];
#else
    body = ALLMASK(n);
#endif

    cyclecount = 0;

    for (i = 0; i < n-2; ++i)
    {
        body &= ~bit[i];
        last = g[i] & body;
        cni = g[i] | bit[i];
        while (last)
        {
            TAKEBIT(j,last);
            inducedpaths(g,i,j,body&~cni,last,bit[edgenumber[i][j]]);
        }
        if (cyclecount > 3*MAXCYCLES/4) gt_abort(">E increase MAXCYCLES\n");
    }

    if (cyclecount > maxcycles) maxcycles = cyclecount;

#if 0
    printf("cyclecount=%d\n",cyclecount);

    {
        setword cyc; int i,j;

        for (i = 0; i < cyclecount; ++i)
        {
            cyc = inducedcycle[i];
            while (cyc)
            {
                TAKEBIT(j,cyc);
                printf(" %d-%d",edge[j].x,edge[j].y);
            }
            printf("\n");
        }
    }
#endif
}

/******************************************************************/

void
surgeproc(FILE *outfile, graph *gin, int n)
/* This is called by geng for each graph. */
{
    int i,j,k,d,n1,n12,n34,n4,ne;
    int isize,jsize;
    graph g[MAXN];
    setword w,wxy,ww,pw,cyc,cycle8;
    int x,y,e1,e2,e3;

    n1 = n12 = n34 = n4 = 0;

    ++gengout;

    for (i = 0; i < n; ++i)
    {
        d = POPCOUNT(gin[i]);   /* deg[i] is not defined yet */
        if      (d == 1) { ++n1; ++n12; }
        else if (d == 2) ++n12;
        else if (d == 3) ++n34;
        else             { ++n34; ++n4; }
    }

    if (n > 1 && (n1 < min1 || n12 < min12 || n34 > max34 || n4 > max4))
        return;

    if (planar && !isplanar(gin,n)) return;    /* Try later */

    ++genggood;

   /* Reverse to put higher degrees first */

    for (i = 0; i < n; ++i)
    {
        w = gin[n-i-1];
        pw = 0;
        while (w)
        {
            TAKEBIT(j,w);
            pw |= bit[n-j-1];
        }
        g[i] = pw;
    }

   /* Make the edge list with default parameters */

    ne = 0;
    for (i = 0; i < n; ++i)
    {
        deg[i] = POPCOUNT(g[i]);
        w = g[i] & BITMASK(i);
        while (w)
        {
            TAKEBIT(j,w);
            edge[ne].x = i;
            edge[ne].y = j;
            edge[ne].xy = bit[i] | bit[j];
            edge[ne].maxmult = maxbond;
            edge[ne].allenemate1 = edge[ne].allenemate2 = -1;
            edgenumber[i][j] = edgenumber[j][i] = ne;
            ++ne;
        }
    }
    numedges = ne;

    if (needcycles)
    {
        if (ne > WORDSIZE)
            gt_abort(">E surge : too many edges for badlists\n");
        findinducedcycles(g,n);
    }

    if (bad1)  /* no triple bonds in rings smaller than 7 */
    {
        for (i = 0; i < cyclecount; ++i)
        {
            cyc = inducedcycle[i];
            if (POPCOUNT(cyc) <= 7)
                while (cyc)
                {
                    TAKEBIT(j,cyc);
                    if (edge[j].maxmult == 2) edge[j].maxmult = 1;
                }
        }
    }

    if (bad2)  /* Bredt's rule for one common bond */
    {
        for (i = 0; i < cyclecount-1; ++i)
        {
            isize = POPCOUNT(inducedcycle[i]);
            if (isize > 6) continue;
            for (j = i+1; j < cyclecount; ++j)
            {
                jsize = POPCOUNT(inducedcycle[j]);
                if (jsize > 6) continue;

                w = inducedcycle[i] & inducedcycle[j];
                if (POPCOUNT(w) != 1) continue;

                if (isize*jsize <= 15)
                    edge[FIRSTBITNZ(w)].maxmult = 0;

                if (isize+jsize <= 9)
                {
                    wxy = edge[FIRSTBITNZ(w)].xy;
                    ww = (inducedcycle[i] | inducedcycle[j]) & ~w;
                    while (ww)
                    {
                        TAKEBIT(k,ww);
                        if ((edge[k].xy & wxy)) edge[k].maxmult = 0;
                    }
                }
            }   
        }
    }

    if (bad3)  /* Bredt's rule for two common bonds */
    {
        for (i = 0; i < cyclecount-1; ++i)
        {
            isize = POPCOUNT(inducedcycle[i]);
            if (isize == 3 || isize > 6) continue;

            for (j = i+1; j < cyclecount; ++j)
            {
                jsize = POPCOUNT(inducedcycle[j]);
                if (jsize == 3 || jsize > 6 || isize+jsize == 12) continue;

                w = inducedcycle[i] & inducedcycle[j];
                if (POPCOUNT(w) != 2) continue;

                ww = w;
                TAKEBIT(k,ww);
                edge[k].maxmult = 0;
                edge[FIRSTBITNZ(ww)].maxmult = 0;
                wxy = edge[k].xy ^ edge[FIRSTBITNZ(ww)].xy;
                
                ww = (inducedcycle[i] | inducedcycle[j]) & ~w;
                while (ww)
                {
                    TAKEBIT(k,ww);
                    if ((edge[k].xy & wxy)) edge[k].maxmult = 0;
                }
            }   
        }
    }

    if (bad4) /* Bredt's rule for two hexagons with 3 bonds in common */
    {
        for (i = 0; i < cyclecount-1; ++i)
        {
            isize = POPCOUNT(inducedcycle[i]);
            if (isize != 6) continue;

            for (j = i+1; j < cyclecount; ++j)
            {
                jsize = POPCOUNT(inducedcycle[j]);
                if (jsize != 6) continue;

                w = inducedcycle[i] & inducedcycle[j];
                if (POPCOUNT(w) != 3) continue;

                ww = inducedcycle[i] | inducedcycle[j];

                TAKEBIT(e1,w);
                TAKEBIT(e2,w);
                e3 = FIRSTBITNZ(w);

                wxy = edge[e1].xy ^ edge[e2].xy ^ edge[e3].xy;
                while (ww)
                {
                    TAKEBIT(k,ww);
                    if ((edge[k].xy & wxy)) edge[k].maxmult = 0;
                }
            }
        }
    }

    if (bad5)  /* No A=A=A, whether in ring or not */
    {
        for (i = 0; i < n; ++i)
        if (deg[i] == 2)
        {
            x = FIRSTBITNZ(g[i]);
            y = FIRSTBITNZ(g[i]&~bit[x]);
            e1 = edgenumber[i][x];
            e2 = edgenumber[i][y];
            if (edge[e2].allenemate1 < 0)
            {
                 if (e1 < e2) edge[e2].allenemate1 = e1;
                 else         edge[e1].allenemate1 = e2;
            }
            else
            {
                 if (e1 < e2) edge[e2].allenemate2 = e1;
                 else         edge[e1].allenemate2 = e2;
            }
        }
    }

    if (bad6)  /* No A=A=A in rings up to length 8 */
    {                           
        cycle8 = 0; 
        for (i = 0; i < cyclecount; ++i)
            if (POPCOUNT(inducedcycle[i]) <= 8) cycle8 |= inducedcycle[i];

        for (i = 0; i < n; ++i)
        if (deg[i] == 2)
        {
            x = FIRSTBITNZ(g[i]);
            y = FIRSTBITNZ(g[i]&~bit[x]);
            e1 = edgenumber[i][x];
            if (!(bit[e1] & cycle8)) continue;
            e2 = edgenumber[i][y];
            if (edge[e2].allenemate1 < 0)
            {
                 if (e1 < e2) edge[e2].allenemate1 = e1;
                 else         edge[e1].allenemate1 = e2;
            }
            else 
            {
                 if (e1 < e2) edge[e2].allenemate2 = e1;
                 else         edge[e1].allenemate2 = e2;
            }
        }
    }

    if (outlevel == 1)
    {
        if (!uswitch) writeg6(outfile,g,1,n);
        return;
    }

   /* Make a SMILES skeleton structure for later use */
    if (smiles) makesmilesskeleton(g,n);

    colourvertices(g,n);
}

/****************************************************************/

static void
decode_formula(char *formula, int *nv,
               int *mine, int *maxe, int *maxd, int *maxc)
/* Parse the input formula. The number of hydrogens goes to hydrogens.
   The other distinct elements go to elementlist[0..numtypes-1] and 
   elementcount[0..numtypes-1].
   *mine and *maxe have an edge range from -e and are updated.
   *mind and *maxc come from -d and -c and are updated.
   *nv gets the number of non-H atoms.
*/
{
    int i,j,d,mult,val,cnt,totval,dbe,forced;
    int maxvcoord,localmine,localmaxe,localmaxd,xi,yi;
    char *s1,*s2,*p;
    int count[FORMULALEN];

    if (numelements > FORMULALEN)
        gt_abort(">E surge : increase FORMULALEN\n");

    /* First we fill in count[*], which is parallel to element[*] */

    for (i = 0; i < numelements; ++i) count[i] = 0;

    for (s1 = formula; *s1 != '\0'; s1 = s2)
    {
        if (!isupper(*s1)) gt_abort(">E surge : unknown element name\n");
        for (s2 = s1+1; islower(*s2); ++s2) {}
        for (i = 0; i < numelements; ++i)
        {
            for (j = 0; element[i].inputname[j] != '\0'
                  && s1+j != s2 && element[i].inputname[j] == s1[j]; ++j) {}
            if (element[i].inputname[j] == '\0' && s1+j == s2) break;
        }
        if (i == numelements) gt_abort(">E surge : unknown element name\n");
        s1 = s2;
        if (!isdigit(*s2))
            ++count[i];
        else
        {
            mult = *s2 - '0';
            for (s2 = s1+1; isdigit(*s2); ++s2) mult = 10*mult+(*s2-'0');
            count[i] += mult;
        }
    }

    /* Next we collect elements actually used into elementtype[0..numtypes-1]
       and elementcount[0..numtypes-1], except for H which we just count. */

    numtypes = hydrogens = 0;
    for (i = 0; i < numelements; ++i)
    {
        cnt = count[i];
        if (cnt > 0)
        {
            if (ISHYDROGEN(i))
                hydrogens = cnt;
            else
            {
                elementtype[numtypes] = i;
                elementcount[numtypes] = cnt;
                ++numtypes;
            }
        }
    }

    /* Next we adjust *maxd and *maxc, as well as the maxcoord
       fields of elements */

    maxvcoord = 0;
    for (i = 0; i < numtypes; ++i)
        if (element[elementtype[i]].maxcoord > maxvcoord)
            maxvcoord = element[elementtype[i]].maxcoord;
    if (maxvcoord < *maxc)
        *maxc = maxvcoord;
    else if (maxvcoord > *maxc)
        for (i = 0; i < numtypes; ++i)
            if (element[elementtype[i]].maxcoord > *maxc)
                element[elementtype[i]].maxcoord = *maxc;
    if (*maxd > *maxc) *maxd = *maxc;

    /* Next we find some bounds on the number of vertices
       with various simple degrees. */

    min1 = min12 = max34 = max4 = 0;
    for (i = 0; i < numelements; ++i)
    {
        if (!ISHYDROGEN(i))
        {
            cnt = count[i];
            val = element[i].maxcoord;
            if (val <= 1) min1 += cnt;   // Check logic
            if (val <= 2) min12 += cnt;
            if (val >= 3) max34 += cnt;
            if (val >= 4) max4 += cnt;
                // Could add max5, could use both bounds everywhere
        }
    }

    /* Now sort by decreasing maximum coordination number */

    for (i = 1; i < numtypes; ++i)  /* really 1 */
    {
        xi = elementtype[i];
        yi = elementcount[i];
        for (j = i; element[elementtype[j-1]].maxcoord < element[xi].maxcoord; )
        {
            elementtype[j] = elementtype[j-1];
            elementcount[j] = elementcount[j-1];
            if (--j < 1) break;
        }
        elementtype[j] = xi;
        elementcount[j] = yi;
    }

    /* Make "canonical" molecule name (used for -A output) */

    p = canonform;
    for (i = 0; i < numtypes; ++i)
    {
        PUTSTR(element[elementtype[i]].inputname);
        if (elementcount[i] > 1) PUTINT(elementcount[i]);
    }
    if (hydrogens > 0) PUTSTR("H");
    if (hydrogens > 1) PUTINT(hydrogens);
    *p = '\0';

    /* Calculate *nv, totalval, forced which all exclude H */

    *nv = forced = 0;
    totval = hydrogens;
    for (i = 0; i < numtypes; ++i)
    {
        j = elementtype[i];
        cnt = elementcount[i];
        *nv += cnt;
        totval += cnt * element[j].valence;
        if (element[j].valence > element[j].maxcoord)
            forced += element[j].valence - element[j].maxcoord;
    }
    forced = (forced+1) / 2;

    if ((totval & 1)) gt_abort(">E surge : impossible parity\n");
    dbe = totval / 2 - (*nv + hydrogens - 1) - forced;
    if (dbe < 0) gt_abort(">E surge : negative DBE\n");
    if (*nv > MAXN) gt_abort(">E surge : too many non-hydrogen atoms\n");
    if (*nv == 0) gt_abort(">E surge : only hydrogen\n");
    valencesum = totval - hydrogens;

    localmine = *nv - 1;
    localmaxe = localmine + dbe;
    if (localmaxe > *nv * (*nv-1) / 2) localmaxe = *nv * (*nv-1) / 2;

    if (localmine > *mine) *mine = localmine;
    if (localmaxe < *maxe) *maxe = localmaxe;
    if (*mine > *maxe) gt_abort(">E surge : edge range is empty\n");

    fprintf(stderr,"%s  ",canonform);
    fprintf(stderr,"H=%d",hydrogens);
    for (i = 0; i < numtypes; ++i)
        fprintf(stderr," %s=%d",element[elementtype[i]].inputname,elementcount[i]);

    fprintf(stderr,"  nv=%d edges=%d-%d DBE=%d maxd=%d maxc=%d\n",
            *nv,*mine,*maxe,dbe,*maxd,*maxc);
    if (*maxe > MAXNE) gt_abort(">E surge : too many edges\n");

    for (d = 1; d <= *maxd; ++d)
    {
        for (i = 0; i < numtypes; ++i)
        {
            val = element[elementtype[i]].maxcoord;
            if (d <= val) maxtype[d] = i;
        }
    }
}

/****************************************************************/

static void
start_geng(int n, int maxd, int maxc,
           int mine, int maxe, char *extra1, long res, long mod)
/* start geng with arguments, extra1 before n, extra2 after n */
{
    int i,geng_argc,mind;
    char *geng_argv[20];
    char arga[30],argb[30];
    char resmod[40];
    char gengargs[80];
    char edgecount[40];

    mind = 1;
    if (hydrogens == 0)
    {
        for (i = 0; i < numtypes; ++i)
            if (element[elementtype[i]].valence < 4) break;
        if (i == numtypes) mind = 2;
    }

    if (n == 1) mind = 0;

    sprintf(arga,"-qcd%dD%d",mind,maxd);
    sprintf(argb,"%d",n);
    sprintf(edgecount,"%d:%d",mine,maxe);

    geng_argv[0] = "geng_surge";
    geng_argv[1] = arga;
    geng_argc = 2;

    if (tswitch && max3rings == 0)
    {
        geng_argv[geng_argc++] = "-t";
        tswitch = FALSE;
    }

    if (fswitch && max4rings == 0)
    {
        geng_argv[geng_argc++] = "-f";
        fswitch = FALSE;
    }

    if (bipartite) geng_argv[geng_argc++] = "-b";

    if (extra1)
    {
        snprintf(gengargs,78,"-%s",extra1);
        geng_argv[geng_argc++] = gengargs;
    }
    geng_argv[geng_argc++] = argb;
    geng_argv[geng_argc++] = edgecount;
    if (mod > 1)
    {
        sprintf(resmod,"%ld/%ld",res,mod);
        geng_argv[geng_argc++] = resmod;
    }
    geng_argv[geng_argc] = NULL;

    if (verbose)
    {
        fprintf(stderr,">geng");
        for (i = 1; geng_argv[i] != NULL; ++i)
            fprintf(stderr," %s",geng_argv[i]);
        fprintf(stderr,"\n");
    }

    geng_main(geng_argc,geng_argv);
}

/*******************************************************************/

static void
processEswitch(char **ps, char *id)
/* Process -E starting at *ps = the character after E, and update *ps.
   The value has the form <inputname>[,<name>]<valence>[<maxcoord>],
   where <inputname> and <name> are either <upper> or <upper><lower> and
   <valence> and <maxcoord> are single digits. Inputname comes first. */
{
    char inputname[3],name[3];
    int valence,maxcoord;
    char *s;
    int state,err;

    s = *ps;
    state = 0;

    while (state < 5)
    {
        switch (state) 
        {
           case 0:
                if (!isupper(*s))
                {
                    state = 6;
                    break;
                }
                inputname[0] = *s++;
                if (islower(*s))
                {
                    inputname[1] = *s++;
                    inputname[2] = '\0';
                }
                else
                    inputname[1] = '\0';
                state = 1;
                break;
           case 1:
                if (isupper(*s))
                    state = 2;
                else
                {
                    name[0] = inputname[0];
                    name[1] = inputname[1];
                    name[2] = inputname[2];
                    state = 3;
                }
                break;
            case 2:
                name[0] = *s++;
                if (islower(*s))
                {
                    name[1] = *s++;
                    name[2] = '\0';
                }
                else
                    name[1] = '\0';
                state = 3;
                break;
            case 3:
                if (!isdigit(*s))
                {
                    state = 7;
                    break;
                }
                valence = *s++ - '0';
                if (!isdigit(*s))
                    maxcoord = valence;
                else
                    maxcoord = *s++ - '0';
                state = 5;
                break;
        }
    }

    if (state == 6)
    {
        fprintf(stderr,">E %s : bad element name\n",id);
        exit(1);
    }
    if (state == 7)
    {
        fprintf(stderr,">E %s : bad valence or maxcoord\n",id);
        exit(1);
    }

    *ps = s;

    addelement(inputname,name,valence,maxcoord);
}

#define SWELEMENT(c,id) if (sw==c) {processEswitch(&arg,id);}

/*******************************************************************/

int
main(int argc, char *argv[])
{
    int argnum,i,j;
    boolean badargs,Gswitch,mswitch,Oswitch,eswitch,notriples;
    boolean oswitch,Bswitch,cswitch,Dswitch;
    char *extra1,*extra2,*formula,*arg,sw,*outfilename;
    long res,mod;
    int mine,maxe,maxd,maxc;
    long eminval,emaxval;
    double t1,t2;
    long badlist[BADLISTS];
    int badlen;

    HELP;

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    maxindex = 0;
    for (i = 0; i < MAXELEMENTS; ++i)
        if (element[i].inputname == NULL) break;
        else if (element[i].index < maxindex && !ISHYDROGEN(i))
            maxindex = element[i].index;
    numelements = i;

#ifdef SURGEPLUGIN_INIT
    SURGEPLUGIN_INIT;
#endif

    argnum = 0;
    badargs = verbose = Gswitch = mswitch = FALSE;
    uswitch = eswitch = notriples = smiles = FALSE;
    oswitch = gzip = alphabetic = Bswitch = FALSE;
    tswitch = fswitch = pswitch = bipartite = FALSE;
    cswitch = planar = xswitch = Dswitch = FALSE;
    Oswitch = FALSE; outlevel = 4;
    extra1 = extra2 = formula = NULL;
    bad1 = bad2 = bad3 = bad4 = bad5 = bad6 = bad7 = bad8 = bad9 = FALSE;

    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                if (sw == 'G')
                {
                    if (Gswitch)
                        gt_abort(">E surge: -G is only allowed once\n");
                    Gswitch = TRUE;
                    extra1 = arg;
                    break;
                }
                else if (sw == 'o')
                {
                    if (oswitch)
                        gt_abort(">E surge : -o is only allowed once\n");
                    oswitch = TRUE;
                    outfilename = arg;
                    break;
                }
                else SWRANGE('m',"/",mswitch,res,mod,"surge -m")
                else SWINT('O',Oswitch,outlevel,"surge -O")
                else SWINT('c',cswitch,maxc,"surge -c")
                else SWINT('d',Dswitch,maxd,"surge -d")
                else SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('v',verbose)
                else SWBOOLEAN('T',notriples)
                else SWBOOLEAN('S',smiles)
                else SWBOOLEAN('z',gzip)
                else SWBOOLEAN('A',alphabetic)
                else SWBOOLEAN('b',bipartite)
                else SWBOOLEAN('P',planar)
                else SWBOOLEAN('x',xswitch)
                else SWSEQUENCEMIN('B',",",Bswitch,badlist,1,BADLISTS,badlen,"surge -B")
                else SWRANGE('e',":-",eswitch,eminval,emaxval,"surge -e")
                else SWRANGE('t',":-",tswitch,min3rings,max3rings,"surge -t")
                else SWRANGE('f',":-",fswitch,min4rings,max4rings,"surge -f")
                else SWRANGE('p',":-",pswitch,min5rings,max5rings,"surge -p")
                else SWELEMENT('E',"surge -E")
#ifdef SURGEPLUGIN_SWITCHES
                else SURGEPLUGIN_SWITCHES
#endif
                else badargs = TRUE;

                if (Bswitch)
                {
                    for (i = 0; i < badlen; ++i)
                    {
                        if (badlist[i] < 1 || badlist[i] > BADLISTS)
                        gt_abort(">E surge : invalid bad list number\n");
                        if      (badlist[i] == 1) bad1 = TRUE;
                        else if (badlist[i] == 2) bad2 = TRUE;
                        else if (badlist[i] == 3) bad3 = TRUE;
                        else if (badlist[i] == 4) bad4 = TRUE;
                        else if (badlist[i] == 5) bad5 = TRUE;
                        else if (badlist[i] == 6) bad6 = TRUE;
                        else if (badlist[i] == 7) bad7 = TRUE;
                        else if (badlist[i] == 8) bad8 = TRUE;
                        else if (badlist[i] == 9) bad9 = TRUE;
                          /* Don't forget initialization if you add more */
                    }
                    Bswitch = FALSE;
                }
            }
        }
        else
        {
            ++argnum;
            if      (argnum == 1) formula = arg;
            else badargs = TRUE;
        }
    }

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE EXTRAUSAGE);
        GETHELP;
        exit(1);
    }

    if (Oswitch && (outlevel <= 0 || outlevel >= 5))
        gt_abort(">E surge : unknown value for -O\n");

    if (!cswitch) maxc = 4;
    if (!Dswitch) maxd = 4;

#ifndef ZLIB
    if (gzip)
        gt_abort(">E surge : -z is only allowed if zlib is compiled in\n");
#endif

    if (uswitch) gzip = oswitch = alphabetic = FALSE;

    if (alphabetic && smiles)
        gt_abort(">E surge : -A and -S are incompatible\n");

    if (!oswitch || (oswitch && strcmp(outfilename,"-") == 0))
        outfilename = "stdout";

    if (bad5) bad6 = FALSE;        /* bad6 is a subset of bad5 */
    if (notriples) bad1 = FALSE;
    if (tswitch && fswitch && max3rings+max4rings <= 1)
        bad9 = FALSE;

    needcycles = (bad1 || bad2 || bad3 || bad4 || bad6); 

    if (fswitch && max4rings < 6) bad7 = FALSE;

    if (tswitch && max3rings < 3) bad8 = FALSE;
    if (fswitch && max4rings < 2) bad8 = FALSE;
    if (pswitch && max5rings == 0) bad8 = FALSE;

    if (gzip)
    {
#ifdef ZLIB
        if (strcmp(outfilename,"stdout") == 0)
            gzoutfile = gzdopen(fileno(stdout),"wb");
        else
            gzoutfile = gzopen(outfilename,"wb");
        if (!gzoutfile)
            gt_abort(">E surge : unable to open compressed stream\n");
        gzbuffer(gzoutfile,1<<16);  /* Remove this line if gzbuffer()
           is not found; it means you have a old version of zlib. */
#endif
    }
    else
    {
        if (strcmp(outfilename,"stdout") == 0)
            outfile = stdout;
        else
            outfile = fopen(outfilename,"w");
        if (!outfile)
            gt_abort(">E surge : can't open output file\n");
    }

    maxbond = (notriples ? 1 : 2);

    if (!Oswitch) outlevel = 4;

    if (mswitch)
    {
        if (res < 0 || res >= mod)
            gt_abort(">E surge : -mres/mod needs 0 <= res < mod\n");
    }
    else
    {
        res = 0;
        mod = 1;
    }

    if (badargs || argnum != 1)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (eswitch)
    {
        mine = (int)eminval;
        maxe = (int)emaxval;
    }
    else
    {
        mine = 0;
        maxe = NOLIMIT;
    }

    decode_formula(formula,&nv,&mine,&maxe,&maxd,&maxc);

    t1 = CPUTIME;
    start_geng(nv,maxd,maxc,mine,maxe,extra1,res,mod);
#ifdef ZLIB
    if (gzip)
        if (gzclose(gzoutfile) != Z_OK)
            gt_abort(">E surge : error on closing compressed stream\n");
#endif
    t2 = CPUTIME;

    if (vgroup) free(vgroup);  /* Make valgrind happy */
    if (egroup) free(egroup);

    if (verbose)
    {
        fprintf(stderr,">G geng made %lld graphs, %lld accepted\n",
                       gengout,genggood);
        if (outlevel > 1)
            fprintf(stderr,">V vcolg %lld nontrivial groups, max size"
              " %ld, made %lld graphs\n",vcolgnontriv,maxvgroup,vcolgout);
        if (outlevel > 2)
            fprintf(stderr,">M multig %lld nontrivial groups, max size"
              " %ld, made %lld graphs\n",multignontriv,maxegroup,multigout);
    }

    if (needcycles) fprintf(stderr,"Max cycles = %d\n",maxcycles);

#ifdef SURGEPLUGIN_SUMMARY
    SURGEPLUGIN_SUMMARY
#endif

    fprintf(stderr,">Z %s %llu -> %llu -> %llu in %.2f sec\n",
        (uswitch ? "generated" : "wrote"),
        gengout,vcolgout,multigout,t2-t1);

    return 0;
}
