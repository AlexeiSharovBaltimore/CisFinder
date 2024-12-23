/*************************************

./patternFind -i ../../data/public-Myc_binding.fa -o ../output/1944.txt -FDR 0.05 -ratio 2 -n 100 -c ../../data/sharov-Smad1_flank.fa -rep ../../data/public-repeats_mouse -score 1

* patternFind is a part of CisFinder software
* http://lgsun.grc.nia.nih.gov/CisFinder
*
* Function: generates position frequency matrixes (PFM) for motifs over-represented in the
* test DNA sequence compared to control sequence
* 
* Syntax:
* patternFind -i inputFasta -o outputFile [-c controlFile, -pos positionFile,
* -ratio minEnrichmentRatio, -FDR maxFDR, -maxlen maxSequenceLength, -len motifLength,
* -strand strandOption, -score scoreOption, -userep, -getrep, -one, -cg, -brief,
* -n numberOfMotifs]
* 
* Comments:
* (a) If controlFile is not specified, then random sequence generated using 3rd order
* Markov chain is used as a control.
* (b) positionFile - stores position of each motif match in the test sequence file.
* It can be later used for motif clustering on the basis of their co-occurrence.
* (c) maxFDR = maximal False Discovery Rate (FDR) threshold. The program generates at least 
* 100 motifs even if they are not significant; additional motifs are included only 
* if they are significant (i.e. FDR < FDR threshold)
* (d) motifLength = 8 as a default. Possible values are 6, 8, and 10 only.
* (e) minEnrichmentRatio = minimum enrichment ratio, default = 1.5
* (f) strandOption: 0 = search both strands, 1 = search positive strand, 2 = search
* positive strand and use negative strand as a control.
* (g) scoreOption:
* 0 = use z-score (z) for motif over-representation for ordering motifs in the output file
* 1 = use z*(ratio-1) for sorting motifs, where ratio = over-representation ratio.
* 2 = use z*info for sorting motifs, where info = information content.
* 3 = use z*(ratio-1)*info
* 4 = use z*(1-selfsim) for sorting motifs, where selfsim = self-similarity of motif.
* 5 = use z*(ratio-1)*(1-selfsim)
* 6 = use z*info*(1-selfsim)
* 7 = use z*(ratio-1)*info*(1-selfsim)
* (e) userep: use repeats in sequence (lower-case in sequence)
* (f) getrep: generate repeat output file (motifs over-represented in repeats).
* (g) one: consider not more than 1 motif occurrence per sequence
* (h) cg: adjust motif abundance to C/G and CpG accurrence in test and control sequences
* (i) brief: do not generate PFM.
* (j) numberOfMotifs = maximum number of motifs to be generated (default = 500)
* See file formats in the readme.txt file
* 
* Author: Alexei Sharov   10/08/2008
* Lab. of Genetics, National Institute on Aging (NIA/NIH)
* Phone: 410-558-8556   Email: sharoval@mail.nih.gov
* 
* See disclaimer in the "readme.txt" file
* If you use this program in your research, please cite:
* Sharov, A.A. and Ko, M.S.H. 2008. Development of CisFinder - an express 
* algorithm to identify over-represented DNA sequence motifs and its application 
* to genome-wide chromatin immunoprecipitation assays for embryonic stem cells.
* BMC Bioinformatics (in review).
**************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define MAXSEGM 10
#define NGAPS 9
#define MAXSEQ 200000
#define LENGTH 10
#define MAXCG LENGTH+1
#define ITMAX 200
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);

typedef struct template_st{
	int n;
	int len[MAXSEGM];
	int start[MAXSEGM];
	int lenTotal;
}TEMPLATE;
typedef struct matrix_st{
	int gap;
	int *dir;
	int len;
	int offset1;		//offset forward
	int offset2;		//offset back
	long code;
	long freq;
	long freqAll;
	float ratio;
	float info;
	float selfSim;
	float score;
	float zvalue;
	float FDR;
	float p;
	int m[LENGTH][4];	//matrix
	char pattern[20];
	char patternRev[20];
	int nmatch;
	long *pos;
	long *match;
	int  *mmm1;
	int  *mmm2;
}MATRIX;
typedef struct param_st{
	char inputFile[300];
	char outputFile[300];
	char controlFile[300];
	char repeatFile[300];
	char positionFile[300];
	int N_matrix;
	int use_repeats;
	int get_repeat_file;
	int strand;
	int score_info;
	int score_ratio;
	int score_selfSim;
	int adjust_cg;
	int presence_option;
	int brief;
	int length;
	float min_ratio;
	float FDRthresh;
	long max_length;
}PARAM;
typedef struct repeat_st{
	float ratio;
	char *repeatType;
	int len;
	int *mmm1;
}REPEAT;
int lengthGlobal=8;
long maxintGlobal=1;

static char codes[16]={'O','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
static char revcodes[16]={'O','T','G','K','C','Y','S','B','A','W','R','D','M','H','V','N'};

/**********************************************/
double         normal_distribution     (float x)
/**********************************************
DESCRIPTION:
   returns an integral of normal probability distribution
   (y changes from 0 to 1; y = 0.5 for x = 0)
   Function translated from FORTRAN code (Numerical Recipes, Press W.H.
   et al. 1986)
*/
{
   static double
      a1=-1.26551223,
      a2= 1.00002368,
      a3= 0.37409196,
      a4= 0.09678418,
      a5=-0.18628806,
      a6= 0.27886807,
      a7=-1.13520398,
      a8= 1.48851587,
      a9=-0.82215223,
      a10=0.17087277;
   double z, t, y;

	z = fabs((double)x)/sqrt(2.);
   t = 1.0 / (1.0 + 0.5 * z);
   y = t*exp(-z * z + a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t *
     (a6 + t * (a7 + t * (a8 + t * (a9 + t * a10)))))))));
   if(x < 0.0) y = 2.0 - y;
   y = 1.0 - 0.5 * y;
   return(y);
}

/***********************************************/
float gammln(float xx)
/***********************************************
* Returns the value ln[gamma(xx)] for xx > 0. */
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) 
		ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

/***********************************************/
float gammq(float a, float x)
/***********************************************/
{
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0){ printf("Invalid arguments in routine gammq\n"); exit(0); }
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

/***********************************************/
void gser(float *gamser, float a, float x, float *gln)
/***********************************************/
{
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0){ printf("x less than 0 in routine gser\n"); exit(0); }
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf("a too large, ITMAX too small in routine gser\n");
		return;
	}
}

/***********************************************/
void gcf(float *gammcf, float a, float x, float *gln)
/***********************************************/
{
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) printf("a too large: %f, ITMAX too small in gcf\n",a);
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

/*************************************************************************/
void sortem   (long ie, float *a, int iperm, float *b,
		float *c, float *d, float *e, float *f, float *g, float *h)
/*************************************************************************/
{
	 long i, j, k, m, p, q, iring;
	 long lt[64], ut[64];
	 float ta, tb, tc, td, te, tf, tg, th, xa, xf, xg;
	 float xh, xe, xd, xc, xb;
	 int i__1;

	 /* Function Body */
	 j = ie-1;
	 m = 1;
	 i = 0;
	 iring = iperm + 1;
	 if (iperm > 7) {
		 iring = 1;
	 }

/* If this segment has more than two elements  we split it */
L10:	 if ((i__1 = j - i - 1) < 0) {
	goto L100;
	 } else if (i__1 == 0) {
	goto L90;
	 } else {
	goto L15;
	 }

/* p is the position of an arbitrary element in the segment we choose the 
* middle element. Under certain circumstances it may be advantageous 
* to choose p at random. */

L15:
	 p = (j + i) / 2;
	 ta = a[p];
	 a[p] = a[i];
	 switch (iring) {
	case 1:  goto L21;
	case 2:  goto L19;
	case 3:  goto L18;
	case 4:  goto L17;
	case 5:  goto L16;
	case 6:  goto L161;
	case 7:  goto L162;
	case 8:  goto L163;
	 }
L163:	 th = h[p];
	 h[p] = h[i];
L162:	 tg = g[p];
    g[p] = g[i];
L161:	 tf = f[p];
    f[p] = f[i];
L16:	 te = e[p];
	 e[p] = e[i];
L17:	 td = d[p];
    d[p] = d[i];
L18:	 tc = c[p];
    c[p] = c[i];
L19:	 tb = b[p];
	 b[p] = b[i];
L21: /* Start at the beginning of the segment, search for k such that a(k)>t */
    q = j;
    k = i;
L20:	 ++k;
	 if (k > q) {
	goto L60;
	 }
    if (a[k] <= ta) {
	goto L20;
    }
/* Such an element has now been found now search for a q such that a(q)<t 
* starting at the end of the segment. */
L30:  if (a[q] < ta) {
	goto L40;
	 }
    --q;
    if (q > k) {
	goto L30;
    }
    goto L50;

/* a(q) has now been found. we interchange a(q) and a(k) */

L40: xa = a[k];
    a[k] = a[q];
    a[q] = xa;
    switch (iring) {
	case 1:  goto L45;
	case 2:  goto L44;
	case 3:  goto L43;
	case 4:  goto L42;
	case 5:  goto L41;
	case 6:  goto L411;
	case 7:  goto L412;
	case 8:  goto L413;
    }
L413:     xh = h[k];
    h[k] = h[q];
	 h[q] = xh;
L412:	 xg = g[k];
    g[k] = g[q];
	 g[q] = xg;
L411:    xf = f[k];
	 f[k] = f[q];
    f[q] = xf;
L41:	 xe = e[k];
	 e[k] = e[q];
    e[q] = xe;
L42:	 xd = d[k];
    d[k] = d[q];
    d[q] = xd;
L43:    xc = c[k];
	 c[k] = c[q];
	 c[q] = xc;
L44:	 xb = b[k];
    b[k] = b[q];
	 b[q] = xb;
L45:
/* Update q and search for another pair to interchange: */
    --q;
   goto L20;
L50:    q = k - 1;
L60:
/* The upwards search has now met the downwards search: */
    a[i] = a[q];
    a[q] = ta;
   switch (iring) {
	case 1:  goto L65;
	case 2:  goto L64;
	case 3:  goto L63;
	case 4:  goto L62;
	case 5:  goto L61;
	case 6:  goto L611;
	case 7:  goto L612;
	case 8:  goto L613;
    }
L613:	 h[i] = h[q];
	 h[q] = th;
L612:    g[i] = g[q];
	 g[q] = tg;
L611:    f[i] = f[q];
    f[q] = tf;
L61:    e[i] = e[q];
	 e[q] = te;
L62:	 d[i] = d[q];
    d[q] = td;
L63:    c[i] = c[q];
    c[q] = tc;
L64:    b[i] = b[q];
    b[q] = tb;
L65:

/* The segment is now divided in three parts: (i,q-1),(q),(q+1,j) */
/* store the position of the largest segment in lt and ut */
    if (q << 1 <= i + j) {
	goto L70;
 }
 lt[m - 1] = i;
 ut[m - 1] = q - 1;
 i = q + 1;
 goto L80;
L70:	 lt[m - 1] = q + 1;
	 ut[m - 1] = j;
	 j = q - 1;
/* Update m and split the new smaller segment */
L80:	 ++m;
	 goto L10;

/* We arrive here if the segment has  two elements we test to see if */
/* the segment is properly ordered if not, we perform an interchange */
L90:
    if (a[i] <= a[j]) {
	goto L100;
    }
	 xa = a[i];
	 a[i] = a[j];
	 a[j] = xa;
	 switch (iring) {
	case 1:  goto L95;
	case 2:  goto L94;
	case 3:  goto L93;
	case 4:  goto L92;
	case 5:  goto L91;
	case 6:  goto L911;
	case 7:  goto L912;
	case 8:  goto L913;
    }
L913:	 xh = h[i];
	 h[i] = h[j];
	 h[j] = xh;
L912:    xg = g[i];
    g[i] = g[j];
	 g[j] = xg;
L911:	 xf = f[i];
    f[i] = f[j];
	 f[j] = xf;
L91:	 xe = e[i];
    e[i] = e[j];
    e[j] = xe;
L92:	 xd = d[i];
	 d[i] = d[j];
    d[j] = xd;
L93:	 xc = c[i];
	 c[i] = c[j];
	 c[j] = xc;
L94:    xb = b[i];
    b[i] = b[j];
	 b[j] = xb;
L95:

/* If lt and ut contain more segments to be sorted repeat process: */
L100:	 --m;
	 if (m <= 0) {
	goto L110;
	 }
	 i = lt[m - 1];
	 j = ut[m - 1];
	 goto L10;
L110:	 return;
} /* sortem_ */

/***********************************************/
void check (void *x)
/***********************************************/
{
	if (!x){
		printf("Out of memory\n");
		exit(0);
	}
}

/***********************************************/
void error_message (char *message)
/***********************************************/
{
	printf("ERROR: %s\n", message);
	exit(0);
}

/*************************************/
long   encode (int *ref)
/*************************************/
{
long k = 0, mult = 1, j;

for(j=0; j<lengthGlobal; ++j){
	k += ref[j]*mult;
	mult *= 4;
}
return(k);
}

/*************************************/
void   decode (long x, int *ref)
/*************************************/
{
long k = x, j;
for(j=0; j<lengthGlobal; ++j){
	ref[j] = k%4;
	k /= 4;
}
return;
}

/*************************************/
long   complementary (long x)
/*************************************/
{
int pos[LENGTH], pos1[LENGTH], j;
decode(x,pos);
for(j=0; j<lengthGlobal; ++j){
	pos1[lengthGlobal-j-1] = 3 - pos[j];
}
return(encode(pos1));
}

/*************************************/
float   get_ratio (int *ref)
/*************************************/
{
int i, x, max=-1000, min=10000000;
float ratio;
for(i=0; i<4; ++i){
	x = ref[i];
	if(min>x) min=x;
	if(max<x) max=x;
}
if(min<1) min = 1;
ratio = (float)max/min;
//printf("%.4f %d %d %d %d %d %d\n",ratio,min,max,ref[0],ref[1],ref[2],ref[3]);
return(ratio);
}

/*************************************/
void  get_self_similarity (MATRIX *Matrix)
/*************************************/
{
int i,j,k, off, nMatch, len;
float r, rmax, *M;

len = Matrix->len;
check(M = (float*)malloc(len*4*sizeof(float)));
for(j=0; j<len; ++j){
	float pp[4], sum;
	sum = 0;
	for(k=0; k<4; ++k){
		pp[k] = Matrix->mmm1[j*4+k];
		if(pp[k] < 0){ pp[k] = 0; }
		pp[k] = log(pp[k]+1);
		sum += pp[k];
	}
	sum /= 4;
	for(k=0; k<4; ++k){
		M[j*4+k] = pp[k]-sum;
	}
}
rmax = 0;
for(off=1; off<4; ++off){
	float r, sxx=0, syy=0, sxy=0, x, y;
	for(j=0; j<(len-off)*4; ++j){
		x = M[j];
		y = M[j+off*4];
		sxy += x*y;
		sxx += x*x;
		syy += y*y;
	}
	r = sxy;
	if(sxx && syy) r = sxy/sqrt(sxx*syy);
	if(rmax<r) rmax = r;
}
Matrix->selfSim = rmax;
free(M);
return;
}

/***********************************************/
int  *convert_to_numbers (char *x, int repeats, long *totalLen, long *length)
/***********************************************/
{
long len, i, pos;
int *y, xx, old=-1;
char ch;

len = strlen(x);
check(y = (int*)malloc(len*sizeof(int)));
*totalLen += len;
pos=0;
if(repeats==0){
	for(i=0;i<len;++i){
		ch = x[i];
		if(ch=='A'){ xx=0; }
		else if(ch=='C'){ xx=1; }
		else if(ch=='G'){ xx=2; }
		else if(ch=='T'){ xx=3; }
		else{ xx=-1; --*totalLen; }
		if(xx>=0 || old!=-1){ y[pos++]=xx; }
		old = xx;
	}
}
else if(repeats==1){
	for(i=0;i<len;++i){
		ch = x[i];
		if(ch=='A'||ch=='a'){ xx=0; }
		else if(ch=='C'||ch=='c'){ xx=1; }
		else if(ch=='G'||ch=='g'){ xx=2; }
		else if(ch=='T'||ch=='t'){ xx=3; }
		else{ xx=-1; --*totalLen; }
		if(xx>=0 || old!=-1){ y[pos++]=xx; }
		old = xx;
	}
}
else if(repeats==2){
	for(i=0;i<len;++i){
		ch = x[i];
		if(ch=='a'){ xx=0; }
		else if(ch=='c'){ xx=1; }
		else if(ch=='g'){ xx=2; }
		else if(ch=='t'){ xx=3; }
		else{ xx=-1; --*totalLen; }
		if(xx>=0 || old!=-1){ y[pos++]=xx; }
		old = xx;
	}
}
*length = pos;
return(y);
}

/*************************************/
void  count_CG  (int *seq, int *nCG, int *nCpG)
/*************************************/
{
int j,x=0,xold;

*nCG=0;
*nCpG = 0;
for(j=0; j<lengthGlobal; ++j){
	xold = x;
	x = seq[j];
	if(x==1 || x==2){ ++*nCG; }
	if(x==2 && xold==1){ ++*nCpG; }
}
return;
}

/*************************************/
void  get_frequency_CG (long *freq, float *frCG, float *frCpG)
/*************************************/
{
long i,j;
int nCG, nCpG, pos[LENGTH];
long cg[MAXCG], cpg[MAXCG], sum=0;

for(i=0; i<=lengthGlobal; ++i){
	cg[i]=0;
	cpg[i]=0;
}
for(i=0; i<maxintGlobal; ++i){
	decode(i,pos);
	count_CG (pos, &nCG, &nCpG);
	cg[nCG] += freq[i];
	cpg[nCG] += freq[i]*nCpG;
}
for(i=0; i<=lengthGlobal; ++i){
	sum += cg[i];
}
for(i=0; i<=lengthGlobal; ++i){
	frCG[i] = (float)cg[i]/sum;
	frCpG[i] = (float)cpg[i]/cg[i];
	//printf("%d %.4f %.4f\n",i,frCG[i],frCpG[i]);
}
return;
}

/***********************************************/
void  reverse_string (char *str) 
/***********************************************/
{
int i, len;
char *newstr;
len=strlen(str);
check(newstr = (char*)malloc((len+1)*sizeof(char)));
for(i=0; i<len; ++i){ newstr[len-i-1]=str[i]; }
memcpy(str,newstr,len);
free(newstr);
}

/***********************************************/
void  count_frequency (long nseq, long *lengthAll, int **seqAll, long *freq, TEMPLATE *template, int presence_option) 
/***********************************************/
{
int i, j, k, n, start, max, len, lenSum=0, *seq;
long iseq, mult[MAXSEGM], multgap[MAXSEGM], sum[MAXSEGM], sumAll, lastMissing, pos, length;
long *found;
int x;

long count1=0, count2=0;
int lenTotal = template->lenTotal;
n = template->n;
check(found = (long*)malloc(maxintGlobal*sizeof(long)));
for(iseq=0; iseq<nseq; ++iseq){
	seq = seqAll[iseq];
	length = lengthAll[iseq];
	lastMissing=-1;
	max = 0;
	lenSum = 0;
	memset(found,0,maxintGlobal*sizeof(long));
	for(i=0; i<n; ++i){
		len = template->len[i];
		mult[i] = 1;
		multgap[i] = 1;
		for(j=0; j<lenSum; ++j){
			multgap[i] *= 4;
		}
		lenSum += len;
		start = template->start[i];
		if(max < start+len){ max=start+len; }
		if(max >= length){ break; }
		k = 0;
		for(j=0; j<len; ++j){
			x = seq[j+start];
			if(x>0){
				k += x*mult[i];
			}else if(x<0){
				lastMissing = j+start;
			}
			mult[i] *= 4;
			//printf("%d %d %d %d %d %d %d %d %d\n",i,j,x,start,len,lenSum,k,mult[i],multgap[i]);
		}
		sum[i] = k;
	}
	if(max >= length){ continue; }
	for(pos=0; pos<length-max; ++pos){
		sumAll = 0; 
		for(i=0; i<n; ++i){
			sumAll += sum[i]*multgap[i];
			start = template->start[i];
			len = template->len[i];
			x = seq[pos+start+len];
			if(x>0){
				sum[i] += x*mult[i];
			}else if(x<0 && i==n-1 && lastMissing<pos+start+len){
				lastMissing = pos+start+len;
			}
			sum[i] /= 4;
			//printf("%d %d %d %d %d %d %d %d %d\n",pos,x,i,start,len,sum[i],sumAll,mult[i],multgap[i]);
		}
		//if(sumAll >= maxintGlobal || sumAll<0) printf("ERR %d %d\n",pos,sumAll);
		++count1;
		if(pos > lastMissing){
			if(presence_option){
				if(!found[sumAll]){
					found[sumAll]=1;
					freq[sumAll]++;
				}
			}else{
				//printf("%d %d\n",pos,found[sumAll]);
				if(pos >= found[sumAll]){
					++count2;
					freq[sumAll]++;
					found[sumAll] = pos+lenTotal;
				}
			}
		}
	}
	sumAll = 0; 
	for(i=0; i<n; ++i){
		sumAll += sum[i]*multgap[i];
	}
	if(length-max > lastMissing){
		if(presence_option){
			if(!found[sumAll]){
				found[sumAll]=1;
				freq[sumAll]++;
			}
		}else{
			if(pos >= found[sumAll]){
				freq[sumAll]++;
				found[sumAll] = pos+lenTotal;
			}
		}
	}
}
free(found);
return;
}

/***********************************************/
void   find_matches (int **sequence1, long *seqlen1, long nseq1, int **sequence2, long *seqlen2, long nseq2, 
	TEMPLATE *template, int *i2matrix, long *complement, MATRIX **matrix, int matStart, int matAdd, int strandOption, float *adjust, int presence_option)
/***********************************************/
{
long i, j, k, n, start, pos, max, len, iseq, code, code1, index;
long mult[MAXSEGM], multgap[MAXSEGM], sum[MAXSEGM];
int x;
long *posMatch, *codeMatch, *closeMatch, Nmatch;
MATRIX *pmat;

for(i=matStart; i<matStart+matAdd; ++i){
	j = matrix[i]->freqAll+1;
	check(matrix[i]->match = (long*)malloc(j*sizeof(long)));
	check(matrix[i]->pos = (long*)malloc(j*sizeof(long)));
	check(matrix[i]->dir = (int*)malloc(j*sizeof(int)));
	check(matrix[i]->mmm1 = (int*)calloc(30*4,sizeof(int)));
	check(matrix[i]->mmm2 = (int*)calloc(30*4,sizeof(int)));
	matrix[i]->nmatch = 0;
}
n = template->n;
int lenTotal = template->lenTotal;
int patternFullLen = lenTotal+4;
for(iseq=0; iseq<nseq1; ++iseq){
	int *seq;
	long length, lenSum=0, lastMissing=0;
	seq = sequence1[iseq];
	length = seqlen1[iseq];
	max = 0;
	for(i=0; i<n; ++i){
		mult[i] = 1;
		k = 0;
		multgap[i] = 1;
		for(j=0; j<lenSum; ++j){
			multgap[i] *= 4;
		}
		len = template->len[i];
		lenSum += len;
		start = template->start[i];
		if(max < start+len){ max=start+len; }
		if(max >= length) break;
		for(j=0; j<len; ++j){
			x = seq[j+start];
			if(x>0){
				k += x*mult[i];
			}else if(x<0){
				lastMissing = j+start;
			}
			mult[i] *= 4;
		}
		sum[i] = k;
	}
	if(max >= length) continue;
	for(pos=0; pos<length-max; ++pos){
		code = 0; 
		for(i=0; i<n; ++i){
			code += sum[i]*multgap[i];
			start = template->start[i];
			len = template->len[i];
			x = seq[pos+start+len];
			if(x>0){
				sum[i] += x*mult[i];
			}else if(x<0 && i==n-1 && lastMissing<pos+start+len){
				lastMissing = pos+start+len;
			}
			sum[i] /= 4;
		}
		if(pos <= lastMissing) continue;
		int found = 0;
		index = i2matrix[code];
		if(index>=0){
			found = 1;
			pmat = matrix[index];
			int jjj = pmat->nmatch;
			if(presence_option && jjj>0 && pmat->match[jjj-1]==iseq){ continue; }
			//printf("%d %d %d %d %d\n",iseq,pos,jjj,code,index);
			/* Remove matches that are too close to each other (e.g. in a simple repeat). */
			if(jjj>=1 && pmat->match[jjj-1]==iseq && pos - pmat->pos[jjj-1]<lenTotal){
			}else{
				pmat->dir[jjj] = 1;
				pmat->match[jjj] = iseq;
				pmat->pos[jjj] = pos;
				pmat->nmatch++;
			}
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-2+i;
				if(pos1<0 || pos1>=length) continue;
				x = seq[pos1];
				if(x<0) continue;
				pmat->mmm1[i*4+x]++;
			}
		}
		code1 = complement[code];
		index = i2matrix[code1];
		if(index>=0 && (!found || code1==code) && !strandOption){
			pmat = matrix[index];
			int jjj = pmat->nmatch;
			if(presence_option && jjj>0 && pmat->match[jjj-1]==iseq){ continue; }
			if(!found){
				if(jjj>=1 && pmat->match[jjj-1]==iseq && pos - pmat->pos[jjj-1]<lenTotal){
				}else{
					pmat->dir[jjj] = -1;
					pmat->match[jjj] = iseq;
					pmat->pos[jjj] = pos;
					pmat->nmatch++;
				}
			}
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-3+patternFullLen-i;
				if(pos1<0 || pos1>=length) continue;
				x = 3-seq[pos1];
				if(x<0) continue;
				pmat->mmm1[i*4+x]++;
			}

		}
	}
	code = 0; 
	for(i=0; i<n; ++i){
		code += sum[i]*multgap[i];
	}
	if(length-max > lastMissing){
		int found = 0;
		index = i2matrix[code];
		if(index>=0){
			found = 1;
			pmat = matrix[index];
			int jjj = pmat->nmatch;
			if(presence_option && jjj>0 && pmat->match[jjj-1]==iseq){ continue; }
			/* Remove matches that are too close to each other (e.g. in a simple repeat). */
			if(jjj>=1 && pmat->match[jjj-1]==iseq && pos - pmat->pos[jjj-1]<lenTotal){
			}else{
				pmat->dir[jjj] = 1;
				pmat->match[jjj] = iseq;
				pmat->pos[jjj] = pos;
				pmat->nmatch++;
			}
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-2+i;
				if(pos1<0 || pos1>=length) continue;
				x = seq[pos1];
				if(x<0) continue;
				pmat->mmm1[i*4+x]++;
			}
		}
		code1 = complement[code];
		index = i2matrix[code1];
		if(index>=0 && (!found || code1==code) && !strandOption){
			pmat = matrix[index];
			int jjj = pmat->nmatch;
			if(presence_option && jjj>0 && pmat->match[jjj-1]==iseq){ continue; }
			if(!found){
				if(jjj>=1 && pmat->match[jjj-1]==iseq && pos - pmat->pos[jjj-1]<lenTotal){
				}else{
					pmat->dir[jjj] = -1;
					pmat->match[jjj] = iseq;
					pmat->pos[jjj] = pos;
					pmat->nmatch++;
				}
			}
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-3+patternFullLen-i;
				if(pos1<0 || pos1>=length) continue;
				x = 3-seq[pos1];
				if(x<0) continue;
				pmat->mmm1[i*4+x]++;
			}
		}
	}
}
if(strandOption>=2){ return; }

/* Process control sequences */
check(posMatch = (long*)malloc(maxintGlobal*sizeof(long)));
check(codeMatch = (long*)malloc(maxintGlobal*sizeof(long)));
check(closeMatch = (long*)calloc(maxintGlobal,sizeof(long)));
for(iseq=0; iseq<nseq2; ++iseq){
	int *seq;
	long length, lenSum=0, lastMissing=0;
	seq = sequence2[iseq];
	length = seqlen2[iseq];
	max = 0;
	Nmatch = 0;
	for(i=0; i<n; ++i){
		mult[i] = 1;
		k = 0;
		multgap[i] = 1;
		for(j=0; j<lenSum; ++j){
			multgap[i] *= 4;
		}
		len = template->len[i];
		lenSum += len;
		start = template->start[i];
		if(max < start+len){ max=start+len; }
		if(max >= length) break;
		for(j=0; j<len; ++j){
			x = seq[j+start];
			if(x>0){
				k += x*mult[i];
			}else if(x<0){
				lastMissing = j+start;
			}
			mult[i] *= 4;
		}
		sum[i] = k;
	}
	if(max >= length) continue;
	for(pos=0; pos<length-max; ++pos){
		code = 0; 
		for(i=0; i<n; ++i){
			code += sum[i]*multgap[i];
			start = template->start[i];
			len = template->len[i];
			x = seq[pos+start+len];
			if(x>0){
				sum[i] += x*mult[i];
			}else if(x<0 && i==n-1 && lastMissing<pos+start+len){
				lastMissing = pos+start+len;
			}
			sum[i] /= 4;
		}
		if(pos <= lastMissing) continue;
		int found = 0, tooClose = 0;
		index = i2matrix[code];
		if(index>=0){
			found = 1;
			for(i=Nmatch-1; i>=0; --i){
				if(pos - posMatch[i] > lenTotal) break;
				if(index==codeMatch[i]) tooClose = 1;
			}
			if(tooClose) closeMatch[index]++;
			posMatch[Nmatch] = pos;
			codeMatch[Nmatch++] = index;
			pmat = matrix[index];
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-2+i;
				if(pos1<0 || pos1>=length) continue;
				x = seq[pos1];
				if(x<0) continue;
				pmat->mmm2[i*4+x]++;
			}
		}
		code1 = complement[code];
		index = i2matrix[code1];
		if(index>=0 && (!found || code1==code) && !strandOption){
			if(!found){
				for(i=Nmatch-1; i>=0; --i){
					if(pos - posMatch[i] > lenTotal) break;
					if(index==codeMatch[i]) tooClose = 1;
				}
				if(tooClose) closeMatch[index]++;
				posMatch[Nmatch] = pos;
				codeMatch[Nmatch++] = index;
			}
			pmat = matrix[index];
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-3+patternFullLen-i;
				if(pos1<0 || pos1>=length) continue;
				x = 3-seq[pos1];
				if(x<0) continue;
				pmat->mmm2[i*4+x]++;
			}
		}
	}
	code = 0; 
	for(i=0; i<n; ++i){
		code += sum[i]*multgap[i];
	}
	if(length-max > lastMissing){
		int found = 0, tooClose = 0;
		index = i2matrix[code];
		if(index>=0){
			found = 1;
			for(i=Nmatch-1; i>=0; --i){
				if(pos - posMatch[i] > lenTotal) break;
				if(index==codeMatch[i]) tooClose = 1;
			}
			if(tooClose) closeMatch[index]++;
			pmat = matrix[index];
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-2+i;
				if(pos1<0 || pos1>=length) continue;
				x = seq[pos1];
				pmat->mmm2[i*4+x]++;
			}
		}
		code1 = complement[code];
		index = i2matrix[code1];
		if(index>=0 && (!found || code1==code) && !strandOption){
			if(!found){
				for(i=Nmatch-1; i>=0; --i){
					if(pos - posMatch[i] > lenTotal) break;
					if(index==codeMatch[i]) tooClose = 1;
				}
				if(tooClose) closeMatch[index]++;
			}
			pmat = matrix[index];
			for(i=0; i<patternFullLen; ++i){
				int pos1 = pos-3+patternFullLen-i;
				if(pos1<0 || pos1>=length) continue;
				x = 3-seq[pos1];
				pmat->mmm2[i*4+x]++;
			}
		}
	}
}
for(i=matStart; i<matStart+matAdd; ++i){
	if(!closeMatch[i] || presence_option) continue;
	pmat = matrix[i];
	code = pmat->code;
	pmat->freq += closeMatch[i]*adjust[code];
	//printf("%d %d %d\n",i,code,closeMatch[i]);
}
free(posMatch);
free(codeMatch);
free(closeMatch);
return;
}

/***********************************************/
char  *truncate_filename(char *filename) 
/***********************************************/
{
char *ch;
ch = strrchr(filename,'/');
if(!ch) ch = filename;
else ++ch;
return(ch);
}

/***********************************************/
char *copy_string(char *str)
/***********************************************/
{
char *ch=NULL;
int len;
len = strlen(str);
if(len){
	check(ch = (char*)malloc((len+1)*sizeof(char)));
	strcpy(ch,str);
}
return(ch);
}

/***********************************************/
int  split_string (char *string, char *items[], int num)
/***********************************************/
{
char *ch;
int i=0;

ch = strchr(string,'\n');
if(ch) *ch = '\0';
ch = string;
while(1){
	items[i] = ch;
	ch = strchr(ch,' ');
	if(ch) *ch = '\0';
	else break;
	if(i>=num-1) break;
	ch++;
	i++;
}
return(i+1);
}

/***********************************************/
long  get_random_sequence(int **sequence1, long *seqlength1, long nseq1, int **sequence2, long *seqlength2)
/***********************************************/
{
long iseq, i, j, *sum, missing, lengthTotal=0;
float *transition;

check(transition = (float*)calloc(256,sizeof(float)));
check(sum = (long*)calloc(64,sizeof(long)));

for(iseq=0; iseq<nseq1; ++iseq){
	int *ref = sequence1[iseq];
	int num = 0;
	if(seqlength1[iseq] < 4) continue;
	missing = -10;
	for(i=0; i<3; ++i){
		int x = ref[i];
		if(x<0){ x=0; missing=i; }
		num = num*4+x;
	}
	//printf("%d  %d\n",iseq,num);
	for(i=3; i<seqlength1[iseq]; ++i){
		int x = ref[i];
		if(x<0){ x=0; missing=i; }
		if(i>missing+3){
			sum[num]++;
		}
		num = num*4+x;
		if(i>missing+3){
			transition[num] += 1;
		}
		x = ref[i-3];
		if(x<0) x=0;
		num -= x*64;
	}
}
for(i=0; i<64; ++i){
	//printf("%d",i);
	for(j=0; j<4; ++j){
		transition[i*4+j] /= sum[i];
		//printf("\t%.4f",transition[i*4+j]);
	}
	//printf("\n");
}
for(iseq=0; iseq<nseq1; ++iseq){
	seqlength2[iseq] = seqlength1[iseq];
	check(sequence2[iseq] = (int*)malloc(seqlength1[iseq]*sizeof(int)));
	int *ref = sequence2[iseq];
	int num = 0;
	for(i=0; i<3; ++i){
		ref[i] = rand()%4;
		num = num*4+ref[i];
	}
	for(i=3; i<seqlength2[iseq]; ++i){
		float p = (float)rand()/RAND_MAX;
		int k = 0;
		while(k<4 && p>=0){
			p -= transition[num*4+k];
			++k;
		}
		k--;
		if(k==4) k=3;
		ref[i] = k;
		num -= ref[i-3]*16;
		num = num*4+k;
		//printf("%d",k);
	}
	lengthTotal += seqlength2[iseq];
	//printf("\n\n");
}
free(transition);
free(sum);
return(lengthTotal);
}

/***********************************************/
void   read_array_int(char *string, int *x, int len)
/***********************************************/
{
char *ch=NULL;
char *items[80];
int i;

split_string (string, items, len);
for(i=0; i<len; ++i){
	x[i]=atoi(items[i]);
}
return;
}

/***********************************************/
void  read_repeats (char *filename, REPEAT *repeats) 
/***********************************************/
{
char *buffer, *ch, *repName;
int gap, len;
long ind, ind1;
long nRepeatNames=0, found;
float ratio;
FILE *fp;

check(buffer = (char*)malloc(3500*sizeof(char)));
check(repName = (char*)calloc(2000,sizeof(char)));
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("ERROR: Repeat file %s not found",ch); exit(0); }
int count=0;
while(fgets(buffer,3499,fp)){
	if(strstr(buffer,"Headers") || strstr(buffer,"Parameters")) continue;
	int num = sscanf(buffer,"%d%ld%f%d%s",&gap,&ind,&ratio,&len,repName);
	if(num<3){ continue; }
	//printf("%d %d %.4f %s\n",gap,ind,ratio,repName);
	ind1 = gap*maxintGlobal+ind;
	repeats[ind1].ratio = ratio;
	repeats[ind1].repeatType = copy_string(repName);
	repeats[ind1].len = len;
	if(len){
		check(repeats[ind1].mmm1 = (int*)malloc(len*4*sizeof(int)));
		fgets(buffer,3499,fp);
		read_array_int(buffer,repeats[ind1].mmm1,len*4);
	}
}
free(buffer);
free(repName);
return;
}

/*************************************/
int   compare_matrixes  (int* M1, int* M2, int N1, int N2, float *match, int *offset)
/*************************************/
{
int offset1, offset_start, offset_end, off, i, j, m, nMatch=0;
float match1, r;

*match=0;
*offset=0;
offset_start = -4;
offset_end = N1-N2+4;
if(N2>N1){
	offset_start = N1-4-N2;
	offset_end = 4;
}
for(off=offset_start; off<offset_end; ++off){
	int start=0, n=0;
	float sxy=0, sxx=0, syy=0;
	if(off>0){ start = off; }
	//printf("%d %d %d\n",N1,N2,off);
	for(i=start; i<N1 && i-off<N2; ++i){
		n++; 
		for(m=0; m<4; ++m){
			float x, y;
			//printf("%d %d\n",i,off);
			x = M1[i*4+m];
			y = M2[(i-off)*4+m];
			//printf("%d %d %.0f %.0f\n",i,m,x,y);
			sxy += x*y;
			sxx += x*x;
			syy += y*y;
		}
		if(sxy < 0 && n>=4){ break; }
	}

	if(sxy < 0 || n <6){ continue; }
	r = sxy/sqrt(sxx*syy);
	if(*match < r){
		*match=r;
		*offset=off;
		nMatch = n;
	}
}
return(nMatch);
}

/***********************************************/
long   read_file (char *filename, int rep, int **sequence, long *seqlength, long *nseq, char **seqNames, long max_length)
/***********************************************/
{
char *seq, *buffer,*seqname,*seqname1,*ch;
long seqlen=0, len, n1=0, parts=0;
FILE *fp;
long totalLen=0;

check(seq = (char*)malloc(MAXSEQ*sizeof(char)));
check(buffer = (char*)malloc(3300*sizeof(char)));
check(seqname = (char*)malloc(2000*sizeof(char)));
check(seqname1 = (char*)malloc(2000*sizeof(char)));
seq[0] = '\0';
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("ERROR: Input file %s not found",ch); exit(0); }
while(fgets(buffer,3299,fp)){
	len = strlen(buffer);
	if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; --len; }
	if(buffer[0]=='>'){
		if(seqlen > 0){
			if(seqNames){
				if(!parts) seqNames[n1] = copy_string(seqname);
				else{
					sprintf(seqname1,"%s-%d",seqname,++parts);
					seqNames[n1] = copy_string(seqname1);
				}
			}
			sequence[n1] = convert_to_numbers(seq, rep, &totalLen, &seqlen);
			if(seqlen<10){ free(sequence[n1]); }
			else{ seqlength[n1++] = seqlen; }
			if(n1>=MAXSEQ-1) break;
			if(totalLen >= max_length) break;
		}
		seq[0] = '\0';
		seqlen = 0;
		strcpy(seqname,&buffer[1]);
		parts = 0;
	}else{
		if(seqlen+len >= MAXSEQ){
			parts++;
			if(seqNames){
				sprintf(seqname1,"%s-%d",seqname,parts);
				seqNames[n1] = copy_string(seqname1);
			}
			sequence[n1] = convert_to_numbers(seq, rep, &totalLen, &seqlen);
			if(seqlen<10){ free(sequence[n1]); }
			else{ seqlength[n1++] = seqlen; }
			if(n1>=MAXSEQ-1) break;
			if(totalLen >= max_length) break;
			seq[0] = '\0';
			seqlen = 0;
		}
		seqlen += len;
		strcat(seq, buffer);
	}
}
if(totalLen < max_length && n1<MAXSEQ-1){
	sequence[n1] = convert_to_numbers(seq, rep, &totalLen, &seqlen);
	if(seqlen<10){ free(sequence[n1]); }
	else{ seqlength[n1++] = seqlen; }
}
fclose(fp);
if(n1>=MAXSEQ-1) printf("Warning: Too many sequences in input file. Loaded %d sequences\n", n1);
*nseq = n1;
printf("File %s loaded\n",ch);
free(seqname);
free(seqname1);
free(seq);
free(buffer);
return(totalLen);
}

/***********************************************/
void  get_templates (TEMPLATE *templates, int length)
/***********************************************/
{
int gap, segment;
TEMPLATE *t;
for(gap=0; gap<NGAPS; ++gap){
	t = &templates[gap];
	if(gap==0){
		t->n = 1;
		t->len[0]=length;
		t->start[0]=0;
		t->lenTotal = length;
	}else if(gap<6){
		segment = length/2;
		t->n = 2;
		t->len[0]=segment;
		t->start[0]=0;
		t->len[1]=segment;
		t->start[1]=segment+gap;
		t->lenTotal = gap+length;
	}else{
		segment = length/3;
		t->n = 3;
		t->len[0]=segment;
		t->start[0]=0;
		t->len[1]=length-2*segment;
		t->start[1]=segment+(gap-5);
		t->len[2]=segment;
		t->start[2]=length-segment+(gap-5)*2;
		t->lenTotal = length+(gap-5)*2;
	}
}
return;
}

/***********************************************/
PARAM *read_parameters (int nargs, char **argv)
/***********************************************/
{
PARAM *p;
int iarg=1, score=0, i;
static char *syntax = "patternFind -i inputFasta -o outputFile [-c controlFile, -pos positionFile, -ratio minEnrichmentRatio, -FDR maxFDR, -maxlen maxSequenceLength, -len motifLength, -strand strandOption, -score scoreOption, -userep, -getrep, -one, -cg, -brief, -n numberOfMotifs]\n";

if(nargs < 4){ printf("%s\n",syntax); exit(0); }
p = (PARAM*)calloc(1,sizeof(PARAM));
p->N_matrix = 500;
p->use_repeats = 0;
p->get_repeat_file = 0;
p->strand = 0;
p->score_info = 0;
p->score_ratio = 0;
p->score_selfSim = 0;
p->adjust_cg = 0;
p->min_ratio = 1.5;
p->presence_option = 0;
p->FDRthresh = 0.1;
p->max_length = 25000000;
p->length = 8;
while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1){ strcpy(p->inputFile,argv[++iarg]); }
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) strcpy(p->outputFile,argv[++iarg]);
	else if(!strcmp(argv[iarg],"-c") && iarg < nargs-1) strcpy(p->controlFile,argv[++iarg]);
	else if(!strcmp(argv[iarg],"-rep") && iarg < nargs-1) strcpy(p->repeatFile,argv[++iarg]);
	else if(!strcmp(argv[iarg],"-pos") && iarg < nargs-1) strcpy(p->positionFile,argv[++iarg]);
	else if(!strcmp(argv[iarg],"-ratio") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->min_ratio);
	else if(!strcmp(argv[iarg],"-FDR") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->FDRthresh);
	else if(!strcmp(argv[iarg],"-maxlen") && iarg < nargs-1) sscanf(argv[++iarg],"%ld",&p->max_length);
	else if(!strcmp(argv[iarg],"-len") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->length);
	else if(!strcmp(argv[iarg],"-userep")) p->use_repeats = 1;
	else if(!strcmp(argv[iarg],"-getrep")) p->get_repeat_file = 1;
	else if(!strcmp(argv[iarg],"-strand") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->strand);
	// strand=1 if search forward only; strand=2 if versus opposite strand
	else if(!strcmp(argv[iarg],"-one")) p->presence_option = 1;
	else if(!strcmp(argv[iarg],"-brief")) p->brief = 1;
	else if(!strcmp(argv[iarg],"-cg")) p->adjust_cg = 1;
	else if(!strcmp(argv[iarg],"-n") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->N_matrix);
	else if(!strcmp(argv[iarg],"-score") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&score);
	else{
		printf("ERROR: Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}
if(p->length!=6 && p->length!=8 && p->length!=10){ p->length=8; printf("Warning: Wrong motif length (6,8,or 10 are acceptable). Restore default = 8\n"); }
lengthGlobal = p->length;
for(i=0;i<p->length; ++i) maxintGlobal *= 4;

if(!p->inputFile[0] || !p->outputFile[0]){ printf("ERROR: %s\n",syntax); exit(0); }
if((score%2)) p->score_ratio=1; score /= 2;
if((score%2)) p->score_info=1; score /= 2;
if((score%2)) p->score_selfSim=1;
//printf("Score: ratio %d; info %d; selfSim %d.\n",p->score_ratio,p->score_info,p->score_selfSim);

if(p->N_matrix <1){ printf("WARNING: Wrong number of output matrixes. Using default N=500\n"); p->N_matrix=500; }
if(p->min_ratio <=1){ printf("WARNING: Minimum ratio should be >1. Using default=1.5\n"); p->min_ratio=1.5; }
return(p);
}

/***********************************************/
int main (int argc, char **argv) 
/***********************************************/
{
long *seqlen1, *seqlen2, *frequency1, *frequency2, *complement;
int **sequence1, **sequence2, gap, k, m, iii, pos[LENGTH], pos1[LENGTH], repeat_option=0;
long n1=0, n2=0, i, j, i1, num, ipos, count=0;
long *add, nadd=0;
int nMatrix=0, *total, *i2matrix, *difference, Npoints;
long length1=0, length2=0;
FILE *fp, *fp1;
TEMPLATE *templates;
float ratio, *ratioAll, *score, *index, *gapSort, maxEntropy, *pvalue;
float frCG1[MAXCG], frCG2[MAXCG], frCpG1[MAXCG], frCpG2[MAXCG];
float adjustLen, *adjust, x, y, log2;
MATRIX **matrix;
PARAM *p;
REPEAT *repeats;
char **seqNames, *ch, *ch1;

p = read_parameters(argc, argv);
check(frequency1 = (long*)malloc(maxintGlobal*sizeof(long)));
check(frequency2 = (long*)malloc(maxintGlobal*sizeof(long)));
check(complement = (long*)malloc(maxintGlobal*sizeof(long)));
check(difference = (int*)malloc(maxintGlobal*sizeof(int)));
check(total = (int*)malloc(maxintGlobal*sizeof(int)));
check(i2matrix = (int*)malloc(maxintGlobal*sizeof(int)));
check(seqlen1 = (long*)malloc(MAXSEQ*sizeof(long)));
check(seqlen2 = (long*)malloc(MAXSEQ*sizeof(long)));
check(sequence1 = (int**)malloc(MAXSEQ*sizeof(int*)));
check(sequence2 = (int**)malloc(MAXSEQ*sizeof(int*)));
check(ratioAll = (float*)malloc(maxintGlobal*sizeof(float)));
check(score = (float*)malloc(maxintGlobal*sizeof(float)));
check(pvalue = (float*)malloc(maxintGlobal*sizeof(float)));
check(index = (float*)calloc(maxintGlobal,sizeof(float)));
check(adjust = (float*)calloc(maxintGlobal,sizeof(float)));
check(matrix = (MATRIX**)malloc(maxintGlobal*sizeof(MATRIX*)));
check(add = (long*)malloc(100*sizeof(long)));
check(seqNames = (char**)malloc(MAXSEQ*sizeof(char*)));
check(templates = (TEMPLATE*)malloc(NGAPS*sizeof(TEMPLATE)));

if(p->get_repeat_file){ repeat_option=2; }
else if(p->use_repeats){ repeat_option=1; }
length1 = read_file(p->inputFile,repeat_option,sequence1,seqlen1,&n1,seqNames,p->max_length);
if(length1 < 10){ error_message("Test sequence is too short (<10)"); }
if(length1 < 5000){ printf("WARNING: Test sequence is too short (<5000)"); }
if(length1 > p->max_length){ printf("WARNING: Sequence is too long (>%d); truncated\n",p->max_length); }
if(p->presence_option){
	if(n1<2){ error_message("Presence option requires > one test sequence"); }
	if(n1<10){ printf("WARNING: Too few test sequences (<10)"); }
}
printf("Number of test sequences = %d. Total length = %d\n",n1,length1);
if(p->repeatFile[0] && !p->get_repeat_file){
	check(repeats = (REPEAT*)calloc(NGAPS*maxintGlobal,sizeof(REPEAT)));
	read_repeats(p->repeatFile,repeats);
}
if(p->controlFile[0] && p->strand<2){
	repeat_option=0;
	if(p->use_repeats && !p->get_repeat_file){ repeat_option=1; }
	length2 = read_file(p->controlFile,repeat_option,sequence2,seqlen2,&n2,NULL,p->max_length);
	if(p->presence_option){ // Adjust the length of control sequences
		int averTest = length1/n1;
		length2 = 0;
		for(i=0; i<n2; ++i){
			if(seqlen2[i] > averTest) seqlen2[i] = averTest;
			length2 += seqlen2[i];
		}
	}
}
else if(p->get_repeat_file){
	repeat_option=0;
	length2 = read_file(p->inputFile,repeat_option,sequence2,seqlen2,&n2,NULL,p->max_length);
}
if(n2){
	if(length2 < 10){ error_message("Control sequence is too short (<10)"); }
	if(length2 < 5000){ printf("WARNING: Control sequence is too short (<5000)"); }
	else if(length2 > p->max_length){ printf("WARNING: Control sequence is too long (>%d), truncated\n",p->max_length); }
	if(p->presence_option){
		if(n2<2){ error_message("Presence option requires > one control sequence"); }
		if(n2<10){ printf("WARNING: Too few control sequences (<10)"); }
	}
	printf("Number of control sequences = %d. Total length = %d\n",n2,length2);
}else{
	length2 = get_random_sequence(sequence1,seqlen1,n1,sequence2,seqlen2);
	n2 = n1;
}
maxEntropy = 2;
adjustLen = (float)length1/length2;
float min_control = (float)length2/maxintGlobal;
float min_total = 0.25*(1+p->min_ratio)*length1/maxintGlobal;

long N1hg = 0;
long N2hg = 0;
if(p->presence_option){
	adjustLen = (float)n1/n2;
	N1hg = n1;
	N2hg = n2;
	double pControl = 1.0 - exp(-1*(double)length2/n2/maxintGlobal);
	min_control = n2*pControl;
	double pTest = 1.0 - exp(-1*(double)length1/n1/maxintGlobal);
	min_total = 0.25*(1+p->min_ratio)*n1*pTest;
}
for(i=0; i<maxintGlobal; ++i){
	complement[i] = complementary(i);
}
log2 = log(2.0);

float *minfreq;
check(minfreq = (float*)malloc(30*sizeof(float)));
get_templates(templates, p->length);
for(gap=0; gap<NGAPS; gap++){
	int addMatrix=0, addMatrix1=0, jj;
	int lenTotal = templates[gap].lenTotal;
	for(i=0; i<maxintGlobal; ++i){
		frequency1[i] = 0;
		frequency2[i]=0;
		i2matrix[i] = -1;
	}
	count_frequency(n1,seqlen1,sequence1,frequency1,&templates[gap],p->presence_option);
	count_frequency(n2,seqlen2,sequence2,frequency2,&templates[gap],p->presence_option);
	if(!adjust[0]){
		if(p->adjust_cg){
			get_frequency_CG(frequency1,frCG1,frCpG1);
			get_frequency_CG(frequency2,frCG2,frCpG2);
		}
		for(i=0; i<maxintGlobal; ++i){
			int nCG, nCpG;
			adjust[i] = adjustLen;
			if(p->adjust_cg){
				decode(i, pos);
				count_CG(pos,&nCG,&nCpG);
				if(frCG2[nCG])
					adjust[i] *= frCG1[nCG]/frCG2[nCG];
				for(j=0; j<nCpG; ++j){
					if(frCpG2[nCpG])
						adjust[i] *= frCpG1[nCpG]/frCpG2[nCpG];
				}
				//printf("%d\t%.4f\n",i,adjust[i]);
			}
		}
	}
	for(i=0; i<maxintGlobal; ++i){
		j = complement[i];
		x = frequency1[i]+frequency1[j];
		y = frequency2[i]+frequency2[j];
		if(p->strand){
			x = frequency1[i];
			y = frequency2[i];
			if(p->strand >=2){
				y = frequency2[j];
			}
		}else if (i==j){
			x /= 2;
			y /= 2;
		}
		if(i!=j && !p->strand && y < 2*min_control) y = 2*min_control;
		else if(y < min_control) y = min_control;
		if(p->presence_option){
			N1hg = n1;
			N2hg = n2;
		}else{
			N1hg = length1/lenTotal;
			N2hg = length2/lenTotal;
		}
		if(x >= N1hg){ x = n1-1; }
		if(y >= N2hg){ y = n1-1; }
		double phg = (x+y)/(N1hg+N2hg);
		double qhg = x/N1hg;
		score[i] = 0;
		if(phg>0){
			score[i] = (qhg-phg)/sqrt(phg*(1-phg)*N2hg/(N1hg+N2hg-1)/N1hg);
		}
		y *= adjust[i];
		ratioAll[i] = (x+5)/(y+5);
		difference[i] = x-y;
		total[i] = x+y;
		//printf("%d\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\n",i,x,y/adjust[i],adjust[i],ratioAll[i],score[i]);
	}
	long nHypotheses = 0;
	for(i=0; i<maxintGlobal; ++i){
		float pp[LENGTH], entropy, info=0, maxAll[LENGTH], min_total1;
		j = complement[i];
		min_total1 = min_total;
		if(!p->strand && i!=j){ min_total1 *= 2; }
		if(total[i] < min_total1 || i2matrix[i] >= 0 || i2matrix[i] == -2) continue;
		if(!p->strand){ i2matrix[j] = -2; }
		nHypotheses++;
		if(ratioAll[i]<p->min_ratio) continue;
		if(score[i] < 1.643 || difference[i] < 8) continue;   /* To ensure p<0.05, single-tail test */

		jj = nMatrix+addMatrix;
		check(matrix[jj] = (MATRIX*)calloc(1,sizeof(MATRIX)));
		MATRIX *M = matrix[jj];
		decode(i, pos);
		for(k=0; k<lengthGlobal; ++k){
			pos1[k] = pos[k];
		}
		for(k=0; k<lengthGlobal; ++k){
			for(m=0; m<4; m++){
				pos1[k] = m;
				i1 = encode(pos1);
				M->m[k][m] = difference[i1];
				pos1[k] = pos[k];
			}
		}
		M->zvalue = score[i];
		M->gap = gap;
		M->code = i;
		M->freq = difference[i];
		M->freqAll = (total[i]+difference[i])/2;
		M->ratio = ratioAll[i];
		//printf("%d  %.4f\n",sum[i],coeff);
		i2matrix[i] = jj;
		if(!p->strand){ i2matrix[j] = jj; }
		++addMatrix;
	}
	find_matches(sequence1,seqlen1,n1,sequence2,seqlen2,n2,&templates[gap],i2matrix,complement,matrix,nMatrix,addMatrix,p->strand,adjust,p->presence_option);
	for(iii=nMatrix; iii<nMatrix+addMatrix; ++iii){
		float pp, entropy, info=0, sum1, sum2, sum3, sumLog1, sumLog2, sumLog3, max, Gtest=0, fr1, fr2, frSum;
		int *mmm1, *mmm2, code, df=0;
		MATRIX *M;
		M = matrix[iii];
		mmm1 = M->mmm1;
		mmm2 = M->mmm2;
		code = M->code;
		int n = templates[gap].n;
		int patternFullLen = templates[gap].start[n-1]+templates[gap].len[n-1]+4;
		int iblock=0;
		ipos = 0;
		for(j=0; j<patternFullLen; ++j){
			if(j<templates[gap].start[iblock]+2){
				sum1=0; sum2=0; sum3=0; sumLog1=0; sumLog2=0;
				for(k=0; k<4; ++k){
					fr1 = mmm1[j*4+k];
					fr2 = mmm2[j*4+k];
					sum1 += fr1;
					sum2 += fr2;
					frSum = fr1+fr2;
					sum3 += frSum;
					if(fr1) sumLog1 += fr1*log(fr1);
					if(fr2) sumLog1 += fr2*log(fr2);
					if(frSum) sumLog2 += frSum*log(frSum);
					fr1 -= fr2*adjust[code];
					if(fr1<0) fr1 = 0;
					mmm1[j*4+k] = fr1;
				}
				sumLog3 = sum1*log(sum1);
				if(sum2) sumLog3 += sum2*log(sum2);
				Gtest += 2*(sumLog1-sumLog2-sumLog3+sum3*log(sum3))/(1+(float)5/6/sum3);
				df += 3;
			}else if (j<templates[gap].start[iblock]+2+templates[gap].len[iblock]){
				for(k=0; k<4; ++k){
					fr1 = M->m[ipos][k];
					if(fr1<0) fr1 = 0;
					mmm1[j*4+k] = fr1;
				}
				++ipos;
			}else{
				if(iblock<n-1){
					++iblock;
					--j;
				}else{
					sum1=0; sum2=0; sum3=0; sumLog1=0; sumLog2=0;
					for(k=0; k<4; ++k){
						fr1 = mmm1[j*4+k]+1;
						fr2 = mmm2[j*4+k]+1;
						sum1 += fr1;
						sum2 += fr2;
						frSum = fr1+fr2;
						sum3 += frSum;
						if(fr1) sumLog1 += fr1*log(fr1);
						if(fr2) sumLog1 += fr2*log(fr2);
						if(frSum) sumLog2 += frSum*log(frSum);
						fr1 -= fr2*adjust[code];
						if(fr1<0) fr1 = 0;
						mmm1[j*4+k] = fr1;
					}
					sumLog3 = sum1*log(sum1);
					if(sum2) sumLog3 += sum2*log(sum2);
					Gtest += 2*(sumLog1-sumLog2-sumLog3+sum3*log(sum3))/(1+(float)5/6/sum3);
					df += 3;
				}
			}
		}
		if(Gtest<0) Gtest=0;
		M->p = 2*gammq(0.5*df, 0.5*Gtest);
		if(M->p > 1) M->p=1;
		//printf("%d\t%d\t%.3f\t%.6f\n",iii,df,Gtest,M->p);

		int start = 0, end = patternFullLen-1;
		while(start<2 && get_ratio(&mmm1[start*4]) < 3) ++start;
		while(end>patternFullLen-3 && get_ratio(&mmm1[end*4])<3) --end;
		if(end-start < 5){
			M->ratio = 0;
			continue;
		}
		M->offset1 = -2+start;
		M->offset2 = -2+(patternFullLen-end-1);
		patternFullLen = end-start+1;
		if(start > 0){
			memmove(mmm1,mmm1+start*4,patternFullLen*4*sizeof(int));
		}
		for(j=0; j<patternFullLen; ++j){
			int x1;
			minfreq[j] = 1000000;
			sum1 = 0;
			for(k=0; k<4; ++k){
				sum1 += mmm1[j*4+k];
			}
			for(k=0; k<4; ++k){
				x1 = mmm1[j*4+k];
				if(sum1==0){ x1=25; }
				else{ x1 = 100*(float)x1/sum1; }
				if(minfreq[j] > x1) minfreq[j]=x1;
				mmm1[j*4+k] = x1;
			}
		}
		// Contrasting of PFM using median of min frequency for each position
		sortem(patternFullLen, minfreq, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
		float median = minfreq[(patternFullLen+1)/2];
		for(j=0; j<patternFullLen; ++j){
			int x=0, y, x1, coef=1, i2, imax;
			max = -1;
			sum1 = 0;
			for(k=0; k<4; ++k){
				x1 = mmm1[j*4+k] - median;
				if(x1<0) x1=0;
				mmm1[j*4+k] = x1;
				if(max < x1) max = x1;
				sum1 += x1;
			}
			entropy = 0;
			for(k=0; k<4; ++k){
				x1 = mmm1[j*4+k];
				y = 0;
				if(x1 > max*0.5){
					y=1;
				}
				x += y*coef;
				coef *= 2;
				pp = (float)x1/sum1;
				if(pp>0) entropy -= pp*log(pp)/log2;
				mmm1[j*4+k] = 100*pp;
			}
			info += maxEntropy - entropy;
			M->pattern[j] = codes[x];
			M->patternRev[j] = revcodes[x];
		}
		M->pattern[patternFullLen]='\0';
		M->patternRev[patternFullLen]='\0';
		reverse_string(M->patternRev);
		M->len = patternFullLen;
		M->info = info;
		get_self_similarity(M);
		double pProduct = M->p * (1.0 - normal_distribution(M->zvalue));
		if(pProduct < 1.0E-30){ pProduct = 1.0E-30; }
		double chiSquare = -2.0*log(pProduct);
		df = 4;
		M->p = 2*gammq(0.5*df, 0.5*chiSquare);

		//if(M->ratio > 3) printf("%d  %.4f  %.4f  %.4f %e\n",iii,M->zvalue,Gtest,chiSquare,M->p);

		pvalue[iii-nMatrix] = M->p;
		index[iii-nMatrix] = iii;
		M->score = M->zvalue;
		if(p->score_ratio){ M->score *= (M->ratio-1); }
		if(p->score_info){ M->score *= info; }
		if(p->score_selfSim){ M->score *= (1.0 - M->selfSim); }
		//printf("%d\t%.4f\t%.4f\n",M->freqAll,M->ratio,M->score);
	}
	sortem(addMatrix, pvalue, 1, index, NULL, NULL, NULL, NULL, NULL, NULL);
	float FDR1 = 1;
	for(iii=addMatrix-1; iii>=0; --iii){
		MATRIX *M;
		M = matrix[(int)index[iii]];
		int rank = iii+1;
		float FDR = M->p/rank*nHypotheses;
		if(FDR >= FDR1) FDR = FDR1;
		else FDR1 = FDR;
		M->FDR = FDR;
		//printf("%d\t%.4f\n",M->freqAll,M->ratio);
		//printf("%d %.5f %.5f\n",iii,FDR);
	}
	//printf("%d %d %d\n",nMatrix,addMatrix,nHypotheses);
	nMatrix += addMatrix;
}
int nSignificant=0;
for(i=0; i<nMatrix; ++i){
	MATRIX *M;
	M = matrix[i];
	index[i] = i;
	score[i] = M->score;
	if(M->FDR <= p->FDRthresh && M->ratio >= p->min_ratio) ++nSignificant;
	//printf("%d %.0f %.4f %.4f %.4f %d - %.4f\n",i,index[i],score[i],M->info,M->ratio,M->freq,M->FDR);
}
sortem(nMatrix, score, 1, index, NULL, NULL, NULL, NULL, NULL, NULL);
fp = fopen(p->outputFile,"w");
if(!fp) error_message("Output file not opened");
ch = truncate_filename(p->inputFile);
if(p->positionFile[0]){
	fp1 = fopen(p->positionFile,"w");
	if(!fp1) error_message("Position file not opened");
	fprintf(fp1,"Parameters:\tmotifsGeneratedBy=patternScan(1.0)\tinput=%s\tNseq=%d",ch,n1);
	if(p->controlFile[0]){ ch = truncate_filename(p->controlFile); fprintf(fp1,"\tControl=%s",ch); }
	fprintf(fp1,"\tScanRepeats=%d\tStrand=%d\tScoreRatio=%d\tScoreInfo=%d\tScoreRepeat=%d\tAdjustCG=%d\tPresence=%d\tMinRatio=%.4f\tFDR=%.4f\n",
	   p->use_repeats,p->strand,p->score_ratio,p->score_info,p->score_selfSim,p->adjust_cg,p->presence_option,p->min_ratio,p->FDRthresh);
	fprintf(fp1,"Headers:\tMotifName\tSequenceNo\tSequenceName\tDir\tPosition\n");
}
if(!p->get_repeat_file){
	fprintf(fp,"Parameters:\tmotifsGeneratedBy=patternScan(1.0)\tinput=%s",ch);
	if(p->controlFile[0]){ ch = truncate_filename(p->controlFile); fprintf(fp,"\tControl=%s",ch); }
	if(p->repeatFile[0]){ ch = truncate_filename(p->repeatFile); fprintf(fp,"\tRepeats=%s",ch); }
	fprintf(fp,"\tScanRepeats=%d\tStrand=%d\tScoreRatio=%d\tScoreInfo=%d\tScoreRepeat=%d\tAdjustCG=%d\tPresence=%d\tMinRatio=%.4f\tFDR=%.4f\n",
	   p->use_repeats,p->strand,p->score_ratio,p->score_info,p->score_selfSim,p->adjust_cg,p->presence_option,p->min_ratio,p->FDRthresh);
	fprintf(fp,"Headers:\tName\tPattern\tPatternRev\tFreq\tRatio\tInfo\tScore\tp\tFDR");
	if(p->repeatFile[0]) fprintf(fp,"\tRepeat");
	fprintf(fp,"\n");
}
for(iii=nMatrix-1; iii>=0; --iii){
	char consensus[30], revconsensus[30];
	float info, score1, max, repeat1;
	int freq, *mmm1;
	MATRIX *M;
	i1 = index[iii];
	M = matrix[i1];
	mmm1 = M->mmm1;
	if(M->ratio < p->min_ratio) continue;
	if(M->FDR > p->FDRthresh && (nSignificant >= 100 || iii < nMatrix-100)) continue;
	score1 = score[iii];
	if(!score1) break;
	freq   = M->freqAll;
	ratio = M->ratio;
	info   = M->info;
	if(p->get_repeat_file){
		if(ratio<3) continue;
		long code = M->code;
		fprintf(fp,"%d\t%d",M->gap,code);
		fprintf(fp,"\t%.3f\t%d\n",ratio,M->len);
		for(j=0; j<M->len; ++j){
			for(k=0; k<4; ++k){
				fprintf(fp,"%d ",mmm1[j*4+k]);
			}
		}
		fprintf(fp,"\n");
		continue;
	}
	//fprintf(fp,"%d\t%d\n",M->code,M->gap);
	if(M->p>1.0e-10 && M->p<0.0002){
		fprintf(fp,">M%03d\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%.1e\t%.1e",++count,M->pattern,M->patternRev,freq,ratio,info,score1,M->p,M->FDR);
	}else{
		fprintf(fp,">M%03d\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f",++count,M->pattern,M->patternRev,freq,ratio,info,score1,M->p,M->FDR);
	}
	if(p->repeatFile[0]){
		long ind3 = M->gap*maxintGlobal+M->code;
		if(!repeats[ind3].ratio){
			fprintf(fp,"\t0\t");
		}else{
			float match=0; int offset;
			compare_matrixes(mmm1,repeats[ind3].mmm1,M->len,repeats[ind3].len,&match,&offset);
			if(match<0.75){
				fprintf(fp,"\t0\t");
			}else{
				fprintf(fp,"\t%.2f",repeats[ind3].ratio);
				if(repeats[ind3].repeatType){
					fprintf(fp,"\t%s",repeats[ind3].repeatType);
				}
			}
		}
	}
	fprintf(fp,"\n");
	if(!p->brief){
		for(j=0; j<M->len; ++j){
			fprintf(fp,"%d",j);
			for(k=0; k<4; ++k){
				fprintf(fp,"\t%d",mmm1[j*4+k]);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	if(p->positionFile[0]){
		for(i=0; i<M->nmatch; ++i){
			int pos0 = M->pos[i] + M->offset1;
			if(M->dir[i]<0) pos0 = M->pos[i] + M->offset2;
			fprintf(fp1,"M%03d\t%d\t%s\t%d\t%d\n",count,M->match[i],seqNames[M->match[i]],M->dir[i],pos0);
		}
	}
	if(count >= p->N_matrix){ break; }
}
fclose(fp);
if(p->positionFile[0]){
	fclose(fp1);
}
if(!count) error_message("No motifs found. Try using lower enrichment ratio.");
printf("patternFind done\n");
return(0);
}

