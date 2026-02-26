#define WORDSIZE 64
#define MAXN WORDSIZE
#include "gtools.h"

#define USAGE \
 "canonsdf [-A|-w|-u] [-p#|-p#:#] [-e#|-e#:#] [-d#] [-c#] [-a] [-Hcode]"

#define HELPTEXT \
" Process molecules in SDF format.\n" \
"\n" \
" -a  Allow more than 4 neighbours to include H\n" \
"      and N with more than 4 neighbours (default not allowed)\n" \
" -d# Maximum degree for simple graph without H (default 4)\n" \
" -c# Maximum coordination number without H (default 4)\n" \
" -e# -e#:# Number of bonds between non-H atoms\n" \
"\n" \
" -p# -p#:#  Only process one input or range of inputs (first is 1)\n" \
" -Hcode Output is restricted to the input with this hash code.\n" \
"   The code continues to the end of the argument. If this option\n" \
"   is given, all other constraints are ignored.\n" \
"\n" \
" -w  Write SDF to output identical to the input\n" \
" -A  Write an ASCII molecule decription\n" \
" -R  Use aromaticity to form canonical form\n" \
" -s  Form canonical form using only simple graph include hydrogens\n" \
" -S  Form canonical form using only simple graph exclude hydrogens\n" \
" -u  Just count\n" \
" The default is to write a 16-character hash code\n"


static char *atom[]
  = { "H","C","N","O","P","S","F","Cl","Br","I","B","Si" };
static int val1[]
  = {  1,  4,  3,  2,  3,  2,  1,  1,   1,   1,  1,  4 };
static int val2[]
  = {  0,  0,  5,  0,  5,  4,  0,  0,   0,   0,  0,  0 };
static int val3[]
  = {  0,  0,  0,  0,  0,  6,  0,  0,   0,   0,  0,  0 };
#define NUMATOMS (sizeof(atom)/sizeof(atom[0]))
#define HYDROGEN 0
#define CARBON 1
#define NITROGEN 2

static int fuzz1[] = {0x37541,0x61532,0x05257,0x26416};
static int fuzz2[] = {0x06532,0x70236,0x35523,0x62437};

#define FUZZ1(x) ((x) ^ fuzz1[(x)&3])
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3])
#define MV(x) ((x) >> (WORDSIZE-n))

typedef char card[84];
static card line[200];
static unsigned long long nin;

static int
atomcode(char *name)
{
    int i;

    for (i = 0; i < NUMATOMS; ++i)
        if (strcmp(atom[i],name) == 0) break;

    if (i == NUMATOMS)
    {
        fprintf(stderr,">E unknown element %s\n",name);
        exit(1);
    }

    return i;
}

static int
roundval(int index, int val)
/* Round val up to the next legal value */
{
    if (val1[index] >= val) return val1[index];
    if (val2[index] >= val) return val2[index];
    if (val3[index] >= val) return val3[index];

    fprintf(stderr,"index=%d val=%d\n",index,val);
    gt_abort(">E can't determine valence\n");
}

static int
readmolecule()
/* Read one SDF molecule into line[*], return number of lines (0 if none) */
{
    int numlines;

    if (fgets(line[0],82,stdin) == NULL) return 0;
    numlines = 1;
    while (fgets(line[numlines],82,stdin) != NULL)
    {
        if (strcmp(line[numlines],"$$$$\n") == 0
                 || strcmp(line[numlines],"$$$$\r\n") == 0
                 || strcmp(line[numlines],"$$$$\r\r\n") == 0)
        {
            if (numlines < 6) gt_abort(">E too few lines\n");
            return numlines+1;
        }
        ++numlines;
    }

    gt_abort(">E EOF encountered\n");
}

static void
writemolecule(int numlines)
{
    int i;

    for (i = 0; i < numlines; ++i) fputs(line[i],stdout);
}

static void
writeascii(int na, int *which, int *hyd,
                   int nb, int *u, int *v, int *mult)
/* Write molecule in ASCII */
{
    int i;

    for (i = 0; i < na; ++i)
    {
        printf("%d:%s",i,atom[which[i]]);
        if (hyd[i] == 1) printf("H ");
        else if (hyd[i] > 1) printf("H%d ",hyd[i]);
        else printf(" ");
    }

    for (i = 0; i < nb; ++i)
    {
        printf(" %d%c%d",u[i],(mult[i]==1?'-':mult[i]==2?'=':'#'),v[i]);
    }
    printf("\n");
}

static void
makecanoncode(int na, int *which, int *hyd,
                   int nb, int *u, int *v, int *mult, char *code)
/* Make a canonical hashcode for the molecule */
{
    int i,j,nhyd;
    int n,nh;
    graph g[MAXN],h[MAXN];
    int lab[MAXN],ptn[MAXN],orbits[MAXN],weight[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    set workspace[200];
    unsigned long long code1,code2;
    char *p;
    char *glyph = 
       "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789@%";

    if (strlen(glyph) != 64) gt_abort(">E glyph error\n");

    nhyd = 0;
    for (i = 0; i < na; ++i) nhyd += hyd[i];

    if (2*na + nhyd > MAXN) gt_abort(">E molecule too big\n");

    n = 2*na + nhyd;
    for (i = 0; i < na; ++i)
    {
        weight[i] = which[i];
        weight[na+i] = which[i] + 128;
    }
    for (i = 2*na; i < n; ++i) weight[i] = 0;
    setlabptn(weight,lab,ptn,n);

    EMPTYSET(g,n);
#define ADDE(ii,jj) { g[ii] |= bit[jj]; g[jj] |= bit[ii]; }

    nh = 2*na;
    for (i = 0; i < na; ++i)
    {
        ADDE(i,i+na);
        for (j = 0; j < hyd[i]; ++j)
        {
            ADDE(i,nh);
            ++nh;
        }
    }

    if (nh != n) gt_abort(">E huh?\n");

    for (i = 0; i < nb; ++i)
    {
        if (mult[i] == 1 || mult[i] == 3)
            ADDE(u[i],v[i]);
        if (mult[i] == 2 || mult[i] == 3)
            ADDE(na+u[i],na+v[i]);
    }

    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,200,1,n,h);

    code1 = 2 + 13*(nhyd+1);
    code2 = n + 17*nb;
    for (i = 0; i < n; ++i)
    {
        code1 = 1237LLU * code1 + FUZZ1(MV(h[i]));
        code2 = 1233457LLU * code2 + FUZZ2(MV(h[i]));
    }

    p = code;
    for (i = 0; i < 8; ++i)
    {
        *p++ = glyph[code1%64];
        code1 >>= 8;
    } 
    for (i = 0; i < 8; ++i)
    {
        *p++ = glyph[code2%64];
        code2 >>= 8;
    } 

    *p = '\0';
}

static void
makesimplecode(int na, int *which, int *hyd, boolean Sswitch,
                   int nb, int *u, int *v, char *code)
/* Make a canonical hashcode for the molecule ignoring bond multiplicity. */
{
    int i,j,nhyd;
    int n,nh;
    graph g[MAXN],h[MAXN];
    int lab[MAXN],ptn[MAXN],orbits[MAXN],weight[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    set workspace[200];
    unsigned long long code1,code2;
    char *p;
    char *glyph = 
       "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789@%";

    if (strlen(glyph) != 64) gt_abort(">E glyph error\n");

    nhyd = 0;

    if (!Sswitch)
    {
        for (i = 0; i < na; ++i) nhyd += hyd[i];
        if (na + nhyd > MAXN) gt_abort(">E molecule too big\n");
    }

    n = na + nhyd;
    for (i = 0; i < na; ++i) weight[i] = which[i];
    for (i = na; i < n; ++i) weight[i] = 0;
    setlabptn(weight,lab,ptn,n);

    EMPTYSET(g,n);
#define ADDE(ii,jj) { g[ii] |= bit[jj]; g[jj] |= bit[ii]; }

    nh = na;
    if (!Sswitch)
    {
        for (i = 0; i < na; ++i)
        {
            for (j = 0; j < hyd[i]; ++j)
            {
                ADDE(i,nh);
                ++nh;
            }
        }
    }

    if (nh != n) gt_abort(">E huh?\n");

    for (i = 0; i < nb; ++i) ADDE(u[i],v[i]);

//writeg6(stdout,g,1,n);

    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,200,1,n,h);

    code1 = 2 + 13*(nhyd+1);
    code2 = n + 17*nb;
    for (i = 0; i < n; ++i)
    {
        code1 = 1237LLU * code1 + FUZZ1(MV(h[i]));
        code2 = 1233457LLU * code2 + FUZZ2(MV(h[i]));
    }

    p = code;
    for (i = 0; i < 8; ++i)
    {
        *p++ = glyph[code1%64];
        code1 >>= 8;
    } 
    for (i = 0; i < 8; ++i)
    {
        *p++ = glyph[code2%64];
        code2 >>= 8;
    } 

    *p = '\0';
}

#define MAXAROMCODES 200
static struct { int mult[2*MAXN]; char code[20]; } aromqueue[MAXAROMCODES];
static int Cpos,aromhead,aromtail;
static struct twostep { int v1,v3; int e1,e2; setword v23; } twosteps[MAXN][4];
    /* Paths of three carbon atoms of the form C-C=C. */
static int numtwosteps[MAXN];

static void
maketwosteps(int na, int nb, int *which, int *u, int *v, int *mult)
/* Find all 2-paths C-C=C.  Could remove tree-like bits. */
{
    int i,j,xi,yi;
    int x[4],ex[4],y[4],ey[4],nx,ny;

    for (i = 0; i < na; ++i) numtwosteps[i] = 0;

    for (i = 0; i < na; ++i)  /* i is the central atom */
    {
        nx = ny = 0;
        if (which[i] == Cpos)
        {
            for (j = 0; j < nb; ++j)
            {
                if (u[j] == i && which[v[j]] == Cpos)
                {
                    if (mult[j] == 1) { ex[nx] = j; x[nx++] = v[j]; }
                    if (mult[j] == 2) { ey[ny] = j; y[ny++] = v[j]; }
                }
                else if (v[j] == i && which[u[j]] == Cpos)
                {
                    if (mult[j] == 1) { ex[nx] = j; x[nx++] = u[j]; }
                    if (mult[j] == 2) { ey[ny] = j; y[ny++] = u[j]; }
                }
            }
            for (xi = 0; xi < nx; ++xi)
            for (yi = 0; yi < ny; ++yi)
            {
                twosteps[x[xi]][numtwosteps[x[xi]]].e1 = ex[xi];
                twosteps[x[xi]][numtwosteps[x[xi]]].e2 = ey[yi];
                twosteps[x[xi]][numtwosteps[x[xi]]].v1 = x[xi];
                twosteps[x[xi]][numtwosteps[x[xi]]].v3 = y[yi];
                twosteps[x[xi]][numtwosteps[x[xi]]].v23 = bit[i] | bit[y[yi]];
                ++numtwosteps[x[xi]];
            }
        }
    }
}

static void
gotaromcycle(int level, struct twostep *path[],
    int na, int *which, int *hyd, int nb, int *u, int *v, int *mult)
{
    int newmult[2*MAXN];
    char code[20];
    int i,j;

    for (i = 0; i < nb; ++i) newmult[i] = mult[i];
    for (i = 0; i < level; ++i)
    {
        newmult[path[i]->e1] = 2;
        newmult[path[i]->e2] = 1;
    }
    makecanoncode(na,which,hyd,nb,u,v,newmult,code);

    for (j = 0; j < aromtail; ++j)
        if (strcmp(code,aromqueue[j].code) == 0) break;
    if (j == aromtail)
    {
        if (aromtail == MAXAROMCODES)
            gt_abort(">E increase MAXAROMCODES\n");
        strcpy(aromqueue[aromtail].code,code);
        for (i = 0; i < nb; ++i)
            aromqueue[aromtail].mult[i] = newmult[i];
        ++aromtail;
    }
}

static void
scancycles(int level, int first, struct twostep *path[], setword avail,
    int na, int *which, int *hyd, int nb, int *u, int *v, int *mult)
{
    int last,i;

    last = (level == 0 ? first : path[level-1]->v3);

    for (i = 0; i < numtwosteps[last]; ++i)
    {
        if (!(twosteps[last][i].v23 & ~avail))
        {
            path[level] = &twosteps[last][i];
            scancycles(level+1,first,path,
                    avail&~(twosteps[last][i].v23|bit[last]),
                    na,which,hyd,nb,u,v,mult);
        }
        else if (twosteps[last][i].v3 == first && level%2 == 0
                  && (twosteps[last][i].v23 & ~avail) == bit[first])
        {
            path[level] = &twosteps[last][i];
            gotaromcycle(level+1,path,na,which,hyd,nb,u,v,mult);
        }
    }
}
 
static void
makearomcode(int na, int *which, int *hyd,
                   int nb, int *u, int *v, int *mult, char *code)
/* Make a canonical hashcode for the aromatic type of the molecule */
{
    int i;
    struct twostep *path[2*MAXN];  /* Sequence of twosteps */

    Cpos = atomcode("C");
    for (i = 0; i < nb; ++i) aromqueue[0].mult[i] = mult[i];
    makecanoncode(na,which,hyd,nb,u,v,mult,aromqueue[0].code);
    aromhead = 0;
    aromtail = 1;

    while (aromhead < aromtail)
    {
        maketwosteps(na,nb,which,u,v,aromqueue[aromhead].mult);

        for (i = 0; i <= na-6; ++i)
            if (which[i] == Cpos)
                scancycles(0,i,path,ALLBITS>>(i+1),na,which,hyd,nb,
                                        u,v,aromqueue[aromhead].mult);
        ++aromhead;
    }

    strcpy(code,aromqueue[0].code);
    for (i = 1; i < aromtail; ++i)
        if (strcmp(code,aromqueue[i].code) > 0)
            strcpy(code,aromqueue[i].code);
}

static void
decode(int numlines, int *n, int *which, int *hyd,
       int *nb, int *u, int *v, int *mult)
/* which[0..*n-1] = non-H atoms
   hyd[0..*n-1] = number of hydrogens attached to each
   *nb = number of non-H bonds
   u[0..*nb-1]-v[0..*nb-1] = bonds
   mult[0..*nb-1] = multiplicities (1,2,3)
*/
{
    int i,j,nv,ne;
    int val[200],needval[200];
    double x,y,z;
    int j1,j2,j3,j4,j5,j6;
    char s[5];

   /* See SDFformats.txt for which formats are handled */

    if (sscanf(line[3],"%d%d",&nv,&ne) != 2)
        gt_abort(">E Can't read count line\n");

    *n = 0;

    for (i = 0; i < nv; ++i) 
    {
        if (sscanf(line[4+i],"%lf%lf%lf%s%d%d%d%d%d%d",&x,&y,&z,
                                    s,&j1,&j2,&j3,&j4,&j5,&j6) != 10)
            gt_abort(">E Can't read atom line\n");
        which[i] = atomcode(s);
        if (which[i] != HYDROGEN)
        {
            if (j4 > 0) hyd[i] = j4-1;
            else        hyd[i] = 0;
            ++*n;

            if (j6 > 0) val[i] = j6;
            else        val[i] = 0;

            needval[i] = (j4 == 0) && (j6 == 0);
        }
    }

    *nb = 0;

    for (i = 0; i < ne; ++i)
    {
        if (sscanf(line[4+nv+i],"%d%d%d",&j1,&j2,&j3) != 3)
            gt_abort(">E Can't read bond line\n");
        if (j3 < 1 || j3 > 4)
            gt_abort_2(">E %llu: irregular multiplicity\n%s\n",nin,line[4+nv+i]);
        if (which[j1-1] != HYDROGEN && which[j2-1] != HYDROGEN)
        {
            u[*nb] = j1-1;
            v[*nb] = j2-1;
            mult[*nb] = j3;
            ++*nb;
        }
        else if (which[j1-1] == HYDROGEN)
            ++hyd[j2-1];
        else if (which[j2-1] == HYDROGEN)
            ++hyd[j1-1];
        else
            gt_abort(">E the impossible happened\n");
    }

#if 0
{ int ii;
printf("----\nn=%d nb=%d nv=%d ne=%d\n",*n,*nb,nv,ne);
printf("atom="); for(ii=0;ii<*n;++ii)printf(" %d",which[ii]);printf("\n");
printf(" hyd="); for(ii=0;ii<*n;++ii)printf(" %d",hyd[ii]);printf("\n");
printf("need="); for(ii=0;ii<*n;++ii)printf(" %d",needval[ii]);printf("\n");
printf(" val="); for(ii=0;ii<*n;++ii)printf(" %d",val[ii]);printf("\n");
printf("edges="); for(ii=0;ii<*nb;++ii)printf(" %d%c%d",
  u[ii],(mult[ii]==1?'-':mult[ii]==2?'=':mult[ii]==3?'#':'*'),v[ii]);
  printf("\n");
}
#endif

    for (i = 0; i < *n; ++i)
    {
        boolean saw4;
        int q;
        /* If aromatic edges are incident, we take one of them
           to be double and the rest of them to be single. */

        q = 0;
        saw4 = FALSE;
        for (j = 0; j < *nb; ++j)
            if (u[j] == i || v[j] == i)
            {
                if (mult[j] == 4)
                {
                    if (!saw4)
                    {
                        q += 2;
                        saw4 = TRUE;
                    }
                    else
                        ++q;
                }
                else
                        q += mult[j];
            }

        if (needval[i])
        {
            val[i] = roundval(which[i],q+hyd[i]);
            hyd[i] = val[i] - q;
        }
        else if (val[i] > 0)
            hyd[i] = val[i] - q;
        else
            val[i] = roundval(which[i],q+hyd[i]);
    }

#if 0
{ int ii;
printf("----\nn=%d nb=%d nv=%d ne=%d\n",*n,*nb,nv,ne);
printf("atom="); for(ii=0;ii<*n;++ii)printf(" %d",which[ii]);printf("\n");
printf(" hyd="); for(ii=0;ii<*n;++ii)printf(" %d",hyd[ii]);printf("\n");
printf("need="); for(ii=0;ii<*n;++ii)printf(" %d",needval[ii]);printf("\n");
printf(" val="); for(ii=0;ii<*n;++ii)printf(" %d",val[ii]);printf("\n");
printf("edges="); for(ii=0;ii<*nb;++ii)printf(" %d%c%d",
  u[ii],(mult[ii]==1?'-':mult[ii]==2?'=':mult[ii]==3?'#':'*'),v[ii]);
  printf("\n");
}
#endif

  /* Try to find a plausible Kekule structure for aromatic edges.
   * This is a real kludge and frequently doesn't work. */

    { int k,laste,queue[2*MAXN],head,tail,ul,vl,ui,vi,du,dv;

        for (k = 0; k < *nb; ++k)
        {
            for (i = 0; i < *nb; ++i) if (mult[i] == 4) break;
            if (i == *nb) break;

            queue[0] = i;
            mult[i] = 1;
            head = 0;
            tail = 1;
            while (head < tail)
            {
                laste = queue[head++];
                ul = u[laste]; vl = v[laste];
                for (i = 0; i < *nb; ++i)
                {
                    ui = u[i]; vi = v[i];
                    if (mult[i] == 4 && (ui == ul || ui == vl ||
                                         vi == ul || vi == vl))
                    {
                        du = val[ui] - hyd[ui];
                        dv = val[vi] - hyd[vi];
                        for (j = 0; j < *nb; ++j)
                        {
                            if (j != i && (u[j] == ui || v[j] == ui))
                                du -= (mult[j] == 4 ? 1 : mult[j]); 
                            if (j != i && (u[j] == vi || v[j] == vi))
                                dv -= (mult[j] == 4 ? 1 : mult[j]); 
                        }
                        if (du <= 1 || dv <= 1) mult[i] = 1;
                        else                    mult[i] = 3 - mult[laste];
                        queue[tail++] = i;
                    }
                }
            }
        }
    }

#if 0
{ int ii;
printf("----\nn=%d nb=%d nv=%d ne=%d\n",*n,*nb,nv,ne);
printf("atom="); for(ii=0;ii<*n;++ii)printf(" %d",which[ii]);printf("\n");
printf(" hyd="); for(ii=0;ii<*n;++ii)printf(" %d",hyd[ii]);printf("\n");
printf("need="); for(ii=0;ii<*n;++ii)printf(" %d",needval[ii]);printf("\n");
printf(" val="); for(ii=0;ii<*n;++ii)printf(" %d",val[ii]);printf("\n");
printf("edges="); for(ii=0;ii<*nb;++ii)printf(" %d%c%d",
  u[ii],(mult[ii]==1?'-':mult[ii]==2?'=':mult[ii]==3?'#':'*'),v[ii]);
  printf("\n");
}
#endif

}

static boolean
isgood(int n, int *which, int *hyd, int nb, int *u, int *v, int *mult,
     boolean aswitch, boolean Rswitch, boolean sswitch, boolean Sswitch,
     int dmax, int cmax, long mine, long maxe, char *neededhash)
{
    int deg[MAXN];
    int i,coord;
    char hashcode[20];

    if (neededhash)
    {
        if (Rswitch)
            makearomcode(n,which,hyd,nb,u,v,mult,hashcode);
        else if (sswitch)
            makesimplecode(n,which,hyd,Sswitch,nb,u,v,hashcode);
        else
            makecanoncode(n,which,hyd,nb,u,v,mult,hashcode);
        return (strcmp(neededhash,hashcode) == 0);
    }

    if (nb < mine || nb > maxe) return FALSE;

    for (i = 0; i < n; ++i) deg[i] = 0;

    for (i = 0; i < nb; ++i)
    {
        ++deg[u[i]];
        ++deg[v[i]];
    }

    for (i = 0; i < n; ++i)
    {
        coord = deg[i] + hyd[i];
        if (!aswitch)
        {
            if (coord > 4 && hyd[i] > 0) return FALSE;
            if (which[i] == NITROGEN && coord > 4) return FALSE;
        }
        if (deg[i] > dmax) return FALSE;
        if (coord > cmax) return FALSE;
    }

    return TRUE;
}

int
main(int argc, char *argv[])
{
    unsigned long long nout;
    boolean badargs;
    boolean eswitch,pswitch,aswitch,dswitch,cswitch,Sswitch;
    boolean Hswitch,wswitch,uswitch,Aswitch,Rswitch,sswitch;
    int numlines,maxd,maxc;
    int j,n,nb;
    int u[2*MAXN],v[2*MAXN],mult[2*MAXN],which[MAXN],hyd[MAXN];
    char sw,*arg,*neededhash;
    long startindex,endindex,mine,maxe;
    char hashcode[20];

    HELP;

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    aswitch = dswitch = cswitch = wswitch = uswitch = Aswitch = FALSE;
    Rswitch = eswitch = pswitch = Hswitch = sswitch = badargs = FALSE;
    Sswitch = FALSE;

    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('w',wswitch)
                else SWBOOLEAN('A',Aswitch)
                else SWBOOLEAN('a',aswitch)
                else SWBOOLEAN('R',Rswitch)
                else SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('S',Sswitch)
                else SWINT('c',cswitch,maxc,"canonsdf -c")
                else SWINT('d',dswitch,maxd,"canonsdf -d")
                else SWRANGE('p',":-",pswitch,
                        startindex,endindex,"canonsdf -p")
                else SWRANGE('e',":-",eswitch,mine,maxe,"canonsdf -e")
                else if (sw == 'H')
                {
                    Hswitch = TRUE;
                    neededhash = arg;
                    break;
                }
                else badargs = TRUE;
            }
        }
        else
            badargs = TRUE;
    }

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!cswitch) maxc = 4;
    if (!dswitch) maxd = 4;
    if (pswitch && startindex < 0) startindex = 0;
    if (!Hswitch) neededhash = NULL;

    if ((uswitch!=0)+(Aswitch!=0)+(wswitch!=0) > 1)
        gt_abort(">E -v, -w and -u are incompatible\n");
    if ((Rswitch!=0)+(sswitch!=0)+(Sswitch!=0) > 1)
        gt_abort(">E -s, -S and -R are incompatible\n");
    if (Sswitch) sswitch = TRUE;

    if (!eswitch)
    {
        mine = 0;
        maxe = NAUTY_INFINITY;
    }

    nin = nout = 0;
    while ((numlines = readmolecule()) != 0)
    {
        ++nin;
        if (!pswitch ||
            (nin >= startindex && nin <= endindex))
        {
            decode(numlines,&n,which,hyd,&nb,u,v,mult);
            if (isgood(n,which,hyd,nb,u,v,mult,aswitch,Rswitch,sswitch,
                        Sswitch,maxd,maxc,mine,maxe,neededhash))
            {
                ++nout;
                if (!uswitch)
                {
                    if (wswitch) writemolecule(numlines);
                    else if (Aswitch) writeascii(n,which,hyd,nb,u,v,mult);
                    else if (Rswitch)
                    {
                        makearomcode(n,which,hyd,nb,u,v,mult,hashcode);
                        printf("%s\n",hashcode);
                    }
                    else if (sswitch)
                    {
                        makesimplecode(n,which,hyd,Sswitch,nb,u,v,hashcode);
                        printf("%s\n",hashcode);
                    }
                    else
                    {
                        makecanoncode(n,which,hyd,nb,u,v,mult,hashcode);
                        printf("%s\n",hashcode);
                    }
                }
            }
        }
    }

    if (uswitch)
        fprintf(stderr,
         ">Z read %lld molecules; %llu removed, %llu remaining\n",
         nin,nin-nout,nout);
    else
        fprintf(stderr,
         ">Z read %lld molecules; %llu removed, %llu written\n",
         nin,nin-nout,nout);

    exit(0);
}
