/*************************************
* patternCluster is a part of CisFinder software
* http://lgsun.grc.nia.nih.gov/CisFinder
*
* Function: clusters motifs with position frequency matrixes (PFM) 
* based on PFM similarity and/or motif co-occurrence in the test sequence.
* Similarity is measured by correlation of position-weight matrices (PWM) which
* are log-transformed PFMs.
* 
* Syntax:
* patternCluster -i inputMotifFile -o outputFile [-pos positionFile, 
* -match matchThreshold, -n numberOfMotifs, -repeat maxRepeatEnrichment, -posonly]
* 
* Comments:
* (a) positionFile - stores position of each motif match in the test sequence file.
* It is used for motif clustering on the basis of their co-occurrence.
* If positionFile is not specified, then clustering is done exclusively on the basis
* of similarity of motifs, otherwise co-occurrence method is used at least for
* clustering of motifs with high level of self-similarity (>0.5 on average). However
* if "posonly" option is used, then clustering is based exclusively on co-occurrence
* of motifs.
* (b) posonly = motifs are clustered exclusively based on co-occurrence.
* (c) matchThreshold = minimum similarity (correlation of PWMs) between motifs;
* default value = 0.75.
* (d) numberOfMotifs can be used to limit the number of input motifs
* (e) maxRepeatEnrichment (ratio threshold) can be used to filter out motifs with 
* high over-representation in repeats.
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
#define MAXLENGTH 150
#define MAXINT 65536
#define MAXDIST 50
#define NHEADERS 10

int DEBUG=0;

typedef struct matrix_st{
	int id;
	int num;
	int alloc;
	int len;
	int checked;
	int palindrome;
	int method;
	long freq;
	float p;
	float FDR;
	float ratio;
	float info;
	float score;
	float repeat;
	float selfSim;	//self-similarity
	int *memb;
	int *off;
	int *dir;
	float m[MAXLENGTH][4];
	char pattern[MAXLENGTH];
	char patternRev[MAXLENGTH];
	char *name;
	char *repeatName;
}MATRIX;
typedef struct pair_st{
	int ind1;
	int ind2;
	int off;
	int dir;
	int nMatch;
	int num;
	float match;
}PAIR;
typedef struct position_st{
	int nmatch;
	int alloc;
	int *dir;
	int *TFid;
	long *pos;
	long iseq;
}POSITION;
typedef struct param_st{
	char *inputFile;
	char *outputFile;
	char *positionFile;
	int N_matrix;
	float matchThresh;
	float max_repeat;
	int position_only;
	int headerIndex[NHEADERS];
	char *param;
}PARAM;
typedef struct cluster_st{
	int c1;
	int c2;
	float distance;
	long maxscore;
}CLUSTER;

static char *codes[16]={"O","A","C","M","G","R","S","V","T","W","Y","H","K","D","B","N"};
static char *revcodes[16]={"O","T","G","K","C","Y","S","B","A","W","R","D","M","H","V","N"};
static char *headerItems[NHEADERS]={"Name","Pattern","PatternRev","Freq","Ratio","Info","Score","p","FDR","Repeat"};

int add_to_cluster(MATRIX *cl, int ind1, int ind2, int dir, int off, MATRIX *Matrix, float threshold);
void   common_ref  (long pos1, int dir1, int len1, long pos2, int dir2, int len2, long *pos3, int *dir3);
void  hierarchical_clustering (float *input, int n, int *select2pattern, float *index1, FILE *fp);
char *copy_string(char *str);

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

/*************************************/
MATRIX *reverse_matrix  (MATRIX *M)
/*************************************/
{
int len, i, j;
MATRIX *M1;

check(M1 = (MATRIX*)malloc(sizeof(MATRIX)));
memcpy(M1, M, sizeof(MATRIX));
len = M->len;
memcpy(M1->pattern,M->patternRev,MAXLENGTH*sizeof(char));
memcpy(M1->patternRev,M->pattern,MAXLENGTH*sizeof(char));
for(i=len-1; i>=0; --i){
	for(j=3; j>=0; --j){
		M1->m[len-i-1][3-j] = M->m[i][j];
	}
}
return(M1);
}

/*************************************/
float   matrix_match  (MATRIX* M1, MATRIX* M2, int offset, int dir, int *nMatch)
/*************************************/
{
MATRIX *M3;
if(dir<0){
	M3 = reverse_matrix(M2);
}else{
	check(M3 = (MATRIX*)malloc(sizeof(MATRIX)));
	memcpy(M3, M2, sizeof(MATRIX));
}
int N1, N2, i, j, m, n=0;
float match, r, sxx=0, syy=0, sxy=0;

N1 = M1->len;
N2 = M2->len;
int start=0;
if(offset>0){ start = offset; }
for(i=start; i<N1 && i-offset<N2; ++i){
	n++; 
	for(m=0; m<4; ++m){
		float x, y;
		x = M1->m[i][m];
		y = M3->m[i-offset][m];
		sxy += x*y;
		sxx += x*x;
		syy += y*y;
	}
}
if(sxy <= 0 || !sxx || !syy){ free(M3); return(0); }
r = sxy/sqrt(sxx*syy);
*nMatch = n;
free(M3);
return(r);
}

/*************************************/
int   compare_matrixes  (MATRIX* M1, MATRIX* M2, float *match, int *offset, float thresh)
/*************************************/
{
int N1, N2, offset1, offset_start, offset_end, off, i, j, m, nMatch=0;
float match1, r;

N1 = M1->len;
N2 = M2->len;

*match=0;
*offset=0;
offset_start = -2;
offset_end = N1-N2+2;
if(N2>N1){
	offset_start = N1-2-N2;
	offset_end = 2;
}
for(off=offset_start; off<offset_end; ++off){
	int start=0, n=0, overhg=0, f1,f2,f3,f4;
	float sxy=0, sxx=0, syy=0;
	if(off>0){ start = off; }
	n = 0;
	for(i=start; i<N1 && i-off<N2; ++i){
		n++; 
		for(m=0; m<4; ++m){
			float x, y;
			//printf("%d %d\n",i,off);
			x = M1->m[i][m];
			y = M2->m[i-off][m];
			sxy += x*y;
			sxx += x*x;
			syy += y*y;
		}
		if(sxy < 0 && n>=4){ break; }
	}
	if(sxy < 0 || n <6){ continue; }
	r = sxy/sqrt(sxx*syy);
	if(*match<thresh && r>*match || r>thresh && r+(1.0-r)*0.5*(1.0-6.0/(float)n) > *match+(1.0- *match)*0.5*(1.0-6.0/(float)nMatch)){
		*match=r;
		*offset=off;
		nMatch = n;
	}
}
return(nMatch);
}

/*************************************/
int  align_matrixes_best (MATRIX* M1, MATRIX* M2, float *match, int *offset, int *dir, float thresh)
/*************************************/
{
float match1;
int offset1, nMatch, nMatch1;
MATRIX *M3;

nMatch = compare_matrixes(M1,M2,match,offset,thresh);
*dir = 1;
M3 = reverse_matrix(M2);
nMatch1 = compare_matrixes(M1,M3,&match1,&offset1,thresh);
free(M3);
if(*match<thresh && match1>*match || match1>thresh && match1+(1.0-match1)*0.5*(1.0-6.0/(float)nMatch1) > *match+(1.0 - *match)*0.5*(1.0-6.0/(float)nMatch)){
	*match = match1;
	*offset = offset1;
	*dir = -1;
	nMatch = nMatch1;
}
return(nMatch);
}

/*************************************/
void  get_self_similarity (MATRIX *Matrix, int nPattern, PARAM *p)
/*************************************/
{
int i, off, nMatch;
float r, rmax;

for(i=0; i<nPattern; ++i){
	rmax = 0;
	for(off=1; off<5; ++off){
		r = matrix_match(&Matrix[i],&Matrix[i],off,1,&nMatch);
		if(nMatch>=5 && rmax<r) rmax = r;
	}
	Matrix[i].selfSim = rmax;
}
return;
}

/**************************/
float   information (float *xr)
/**************************/
{
float sum=0, x, maxEntropy, pp[4], entropy=0;
int i;

maxEntropy = 2*log(2.0);
for(i=0; i<4; ++i){
	pp[i] = exp(xr[i]);
	sum += pp[i];
}
for(i=0; i<4; ++i){
	pp[i] /= sum;
	entropy -= pp[i]*log(pp[i]);
}
return(maxEntropy - entropy);
}

/*************************************/
void   update_matrix (MATRIX *cluster, MATRIX *Matrix, float matchThresh)
/*************************************/
{
int Npos=0, Nseq, **coordLink, *nLink, i, j, start, end, len, ind, off;
long freqTest=0, freqCont=0;
float Mnew[100][4], r, rmax;
MATRIX *M;

Nseq = cluster->num;
cluster->freq = 0;
cluster->repeat = 0;
for(i=0; i<Nseq; ++i){
	ind = cluster->memb[i];
	len = Matrix[ind].len;
	start = cluster->off[i];
	end = start+len;
	if(Npos<end){ Npos=end; }
	freqTest += Matrix[ind].freq;
	freqCont += Matrix[ind].freq/Matrix[ind].ratio;
	float repeatFold = Matrix[ind].repeat;
	if(repeatFold < 1) repeatFold=1;
	cluster->repeat += repeatFold * Matrix[ind].freq;
	cluster->freq += Matrix[ind].freq;
	if(cluster->score < Matrix[ind].score){
		cluster->score = Matrix[ind].score;
	}
}
if(Npos > MAXLENGTH){ Npos = MAXLENGTH; }

cluster->ratio = (float)freqTest/freqCont;
cluster->repeat /= cluster->freq;
if(cluster->repeat < 2) cluster->repeat = 0;
//printf("%.4f\n",cluster->repeat);


if(!Npos) error_message("In update_matrix");
cluster->len = Npos;
check(coordLink = (int**)malloc(Npos*sizeof(int*)));
check(nLink = (int*)calloc(Npos,sizeof(int)));
for(i=0; i<Npos; ++i){
	check(coordLink[i] = (int*)malloc(Nseq*sizeof(int)));
}
for(i=0; i<Nseq; ++i){
	ind = cluster->memb[i];
	len = Matrix[ind].len;
	start = cluster->off[i];
	end = start+len;
	if(end > MAXLENGTH){ end = MAXLENGTH; }
	for(j=start; j<end; ++j){
		coordLink[j][nLink[j]]=i;
		nLink[j] += 1;
	}
}
cluster->info = 0;
for(i=0; i<Npos; ++i){
	int k, dir, pos, name, Ninfo;
	float sumWgt = 0, newm[4]={0,0,0,0}, info, wgt, x;
	if(!nLink[i]){
		for(j=0; j<4; ++j){
			newm[j] += 1;
		}
	}
	for(k=0; k<nLink[i]; ++k){
		name = coordLink[i][k];
		ind = cluster->memb[name];
		M = &Matrix[ind];
		start = cluster->off[name];
		dir = cluster->dir[name];
		pos = i-start;
		if(dir<0){
			pos = M->len-pos-1;
		}
		info = 0;
		Ninfo = 0;
		for(j=pos-1; j<pos+2; ++j){
			if(j<0 || j==M->len) continue;
			info += information(&M->m[j][0]);
			Ninfo++;
		}
		wgt = M->freq * info/Ninfo;
		sumWgt += wgt;
		for(j=0; j<4; ++j){
			if(dir>0){ x = M->m[pos][j]; }
			else{ x = M->m[pos][3-j]; }
			newm[j] += x*wgt;
		}
	}
	for(j=0; j<4; ++j){
		if(sumWgt) newm[j] /= sumWgt;
		cluster->m[i][j] = newm[j];
	}
	cluster->info += information(&cluster->m[i][0]);
}
rmax = 0;
for(off=1; off<4; ++off){
	int nMatch=0;
	r = matrix_match(&Matrix[i],&Matrix[i],off,1,&nMatch);
	if(nMatch>=5 && rmax>r) rmax = r;
}
Matrix[i].selfSim = rmax;
for(i=0; i<Npos; ++i) free(coordLink[i]);
free(coordLink);
free(nLink);
return;
}

/***********************************************/
void  check_palindrome (MATRIX *M, float thresh)
/***********************************************/
{
int len, nMatch, offset;
float match;
MATRIX *M1;
M1 = reverse_matrix(M);
nMatch = compare_matrixes(M,M1,&match,&offset,thresh);
if(match > 0.7 && nMatch >=5) M->palindrome=1;
else M->palindrome=0;
free(M1);
return;
}

/**************************************************/
MATRIX *align_cluster_sim (int *members, int Nseq, int *remains, int *nRemains, MATRIX *Matrix, int nMatrix,
			PAIR *pairs, int nPairs, float matchThresh)
/**************************************************/
{
MATRIX *newCluster;
int *memberClust, i, j, ind1, ind2, dir, offset, nPairs1=0, count=0;
int *memberLinks, bestPair, found=0;
float match, *scorePair, *sorted;
PAIR *pairs1;

if(Nseq <= 1){
	return(NULL);
}
check(memberClust = (int*)calloc(nMatrix,sizeof(int)));
for(i=0; i<Nseq; ++i){
	memberClust[members[i]] = 1;
}
for(i=0; i<nPairs; ++i){
	ind1 = pairs[i].ind1;
	ind2 = pairs[i].ind2;
	if(memberClust[ind1] && memberClust[ind2]){
		found = 1;
		memberClust[ind1]++;
		memberClust[ind2]++;
	}
}
if(!found){
	free(memberClust);
	return(NULL);
}
check(pairs1 = (PAIR*)malloc(Nseq*sizeof(PAIR)));
check(newCluster = (MATRIX*)calloc(1,sizeof(MATRIX)));
check(newCluster->memb = (int*)malloc(Nseq*sizeof(int)));
check(newCluster->dir = (int*)malloc(Nseq*sizeof(int)));
check(newCluster->off = (int*)malloc(Nseq*sizeof(int)));
check(scorePair = (float*)malloc(nPairs*sizeof(float)));
check(sorted = (float*)malloc(nPairs*sizeof(float)));
for(i=0; i<nPairs; ++i){
	ind1 = pairs[i].ind1;
	ind2 = pairs[i].ind2;
	sorted[i] = i;
	if(memberClust[ind1] && memberClust[ind2]){
		scorePair[i] = (memberClust[ind1]+memberClust[ind2])*(pairs[i].match+(1.0-pairs[i].match)*0.5*(1.0-6.0/(float)pairs[i].nMatch));
	}else{
		scorePair[i] = 0;
	}
}
if(nPairs>1)
	sortem(nPairs, scorePair, 1, sorted, NULL, NULL, NULL, NULL, NULL, NULL);

/* put the best pair into the new cluster */
bestPair = sorted[nPairs-1];
newCluster->num = 2;
ind1 = pairs[bestPair].ind1;
ind2 = pairs[bestPair].ind2;
align_matrixes_best (&Matrix[ind1],&Matrix[ind2],&match,&offset,&dir,matchThresh);
newCluster->memb[0]=ind1;
newCluster->memb[1]=ind2;
newCluster->dir[0]=1;
newCluster->dir[1]=dir;
if(offset>=0){
	newCluster->off[0]=0;
	newCluster->off[1]=offset;
}else{
	newCluster->off[0]=-offset;
	newCluster->off[1]=0;
}
for(i=0; i<Nseq; ++i){
	memberClust[members[i]] = 1;
}
memberClust[pairs[bestPair].ind1] = 2;
memberClust[pairs[bestPair].ind2] = 2;
if(Nseq == 2){
	update_matrix(newCluster,Matrix,matchThresh);
	free(pairs1);
	free(memberClust);
	free(sorted);
	free(scorePair);
	return(newCluster);
}
check(sorted = (float*)realloc(sorted,Nseq*sizeof(float)));
check(scorePair = (float*)realloc(scorePair,Nseq*sizeof(float)));
while(1){
	int offset,dir, found=0, offsetAdjust=0, ind, num, k, nPairs1;
	update_matrix(newCluster, Matrix,matchThresh);
	nPairs1=0;
	for(i=0; i<Nseq; ++i){
		ind = members[i];
		if(memberClust[ind]==2 || Matrix[ind].len > newCluster->len+2){ continue; }
		int nnn = align_matrixes_best (newCluster, &Matrix[ind], &match, &offset, &dir,matchThresh);
		if(match < matchThresh){ continue; }
		//printf("%d %d %d %d %d %.4f\n",i,ind,offset,dir,nPairs1,match);
		pairs1[nPairs1].ind1=ind;
		pairs1[nPairs1].match=match;
		pairs1[nPairs1].dir=dir;
		pairs1[nPairs1].off=offset;
		sorted[nPairs1]=nPairs1;
		scorePair[nPairs1]=match+(1.0-match)*0.5*(1.0-6.0/(float)nnn);
		++nPairs1;
	}
	if(!nPairs1){
		break;
	}else if(nPairs1>1){
		sortem(nPairs1, scorePair, 1, sorted, NULL, NULL, NULL, NULL, NULL, NULL);
	}
	num = newCluster->num;
	for(j=nPairs1-1; j>=0 && j>nPairs1-1-newCluster->num*2; --j){
		k = sorted[j];
		ind = pairs1[k].ind1;
		offset = pairs1[k].off;
		dir = pairs1[k].dir;
		//printf("%d %d %d %d %.4f\n",j,ind,offset,dir,pairs1[k].match);
		memberClust[ind]=2;
		newCluster->memb[num] = ind;
		offset += offsetAdjust;
		newCluster->off[num] = offset;
		newCluster->dir[num] = dir;
		found = 1;
		if(offset<0){
			for(k=0; k<=num; ++k){
				newCluster->off[k] -= offset;
			}
			offsetAdjust -= offset;
		}
		num++;
	}
	newCluster->num = num;
	count++;
}
*nRemains=0;
for(i=0; i<Nseq; ++i){
	ind1 = members[i];
	if(memberClust[ind1]==1){
		remains[*nRemains]=ind1;
		*nRemains += 1;
	}
}
free(sorted);
free(scorePair);
free(pairs1);
free(memberClust);
return(newCluster);
}

/**************************************************/
MATRIX *align_cluster_pos (int *members, int Nseq, int *remains, int *nRemains, MATRIX *Matrix, int nPattern,
			PAIR *pairs, int nPairs, float matchThresh)
/**************************************************/
{
MATRIX *newCluster;
int i, j, ind1, ind2, dir, off, ipair;
int *memberLinks, *memberClust, found=0;
float match, *score, *sorted, threshold;

if(Nseq == 1){
	return(NULL);
}
check(memberClust = (int*)calloc(nPattern,sizeof(int)));
for(i=0; i<Nseq; ++i){
	memberClust[members[i]] = 1;
}
check(newCluster = (MATRIX*)calloc(1,sizeof(MATRIX)));
for(ipair=0; ipair<nPairs; ++ipair){
	ind1 = pairs[ipair].ind1;
	ind2 = pairs[ipair].ind2;
	off = pairs[ipair].off;
	dir = pairs[ipair].dir;
	if(!memberClust[ind1] || !memberClust[ind2]) continue;
	if(!newCluster->num){
		add_to_cluster(newCluster,-1,ind1,1,0,Matrix,matchThresh);
		add_to_cluster(newCluster,-1,ind2,dir,off,Matrix,matchThresh);
		update_matrix(newCluster, Matrix,matchThresh);
		memberClust[ind1] = 2;
		memberClust[ind2] = 2;
		found = 1;
		continue;
	}
	if(memberClust[ind1]==memberClust[ind2]) continue;
	if(memberClust[ind1]==1){
		if(dir>0) off = -off;
		else off = off+Matrix[ind2].len - Matrix[ind1].len;
		int swap=ind1; ind1=ind2; ind2=swap;
	}
	if(pairs[ipair].match < 0.05) continue;
	if(add_to_cluster(newCluster,ind1,ind2,dir,off,Matrix,matchThresh)){
		if(DEBUG) printf("m%d via m%d\n",ind2,ind1);
		update_matrix(newCluster, Matrix,matchThresh);
		memberClust[ind2] = 2;
		found = 1;
	}
}
if(!found){
	free(memberClust);
	free(newCluster);
	return(NULL);
}
*nRemains=0;
for(i=0; i<Nseq; ++i){
	ind1 = members[i];
	if(memberClust[ind1]==1){
		remains[*nRemains]=ind1;
		*nRemains += 1;
	}
}
free(memberClust);
return(newCluster);
}

/***********************************************/
void  reverse_string (char *str) 
/***********************************************/
{
int i, len;
char *newstr;
len=strlen(str);
check(newstr = (char*)calloc((len+1),sizeof(char)));
for(i=0; i<len; ++i){ newstr[len-i-1]=str[i]; }
if(len) memcpy(str,newstr,len*sizeof(char));
free(newstr);
}

/*************************************/
void   print_matrix  (MATRIX *M, FILE *fp)
/*************************************/
{
int i,j;
float pp[4], sum;
for(i=0; i<M->len; ++i){
	fprintf(fp,"%d",i);
	sum = 0;
	for(j=0; j<4; ++j){
		pp[j] = exp(M->m[i][j])-0.01;
		if(pp[j]<0) pp[j]=0;
		sum += pp[j];
	}
	for(j=0; j<4; ++j){
		pp[j] *= 100/sum;
		fprintf(fp,"\t%.0f",pp[j]);
	}
	fprintf(fp,"\n");
}
fprintf(fp,"\n");
return;
}

/*************************************/
void   print_cluster (int clusterNo, MATRIX **clusters, int icluster, FILE *fp, MATRIX *Matrix,
		PAIR *posPairs, int nPosPairs, int nPattern, PARAM *p)
/*************************************/
{
int len, clLen, coef,y,x, i,j,k,num, max_freq=0, *pattern2member,dir;
long off;
float max, x1, *index,*freq;
MATRIX *cluster;

cluster = clusters[icluster];
num = cluster->num;
clLen = cluster->len;
check(freq = (float*)calloc(clLen,sizeof(float)));
for(i=0;i<num;++i){
	off = cluster->off[i];
	j = cluster->memb[i];
	len = Matrix[j].len;
	for(k=0;k<len;++k){
		if(off+k<0 || off+k>=clLen) continue;
		float fr = Matrix[j].freq*Matrix[j].ratio;
		if(freq[off+k] < fr) freq[off+k] = fr;
	}
}
if(num > 1){
	for(i=0; i<clLen; ++i){
		max=-1;
		for(j=0; j<4; ++j){
			x1 = cluster->m[i][j];
			if(x1>max){ max=x1; }
		}
		x=0;
		coef=1;
		for(j=0; j<4; ++j){
			y = 1;
			if(cluster->m[i][j] < max-0.6931){ y=0; }
			x += y*coef;
			coef *= 2;
		}
		strcat(cluster->pattern,codes[x]);
		strcat(cluster->patternRev,revcodes[x]);
	}
	reverse_string(cluster->patternRev);
}
fprintf(fp,">C%03d\t%s\t%s",clusterNo,cluster->pattern,cluster->patternRev);
if(p->headerIndex[3]>=0) fprintf(fp,"\t%d",cluster->freq);
if(p->headerIndex[4]>=0) fprintf(fp,"\t%.4f",cluster->ratio);
fprintf(fp,"\t%.4f\t%.4f",cluster->info,cluster->score);
if(p->headerIndex[7]>=0) fprintf(fp,"\t%.4f",cluster->p);
if(p->headerIndex[8]>=0) fprintf(fp,"\t%.4f",cluster->FDR);
fprintf(fp,"\t%d\t%d\t%d",cluster->palindrome,num,cluster->method);
if(num==1) fprintf(fp,"\t%s",cluster->name);
else fprintf(fp,"\t-1");
if(p->headerIndex[9]>=0){ fprintf(fp,"\t%.2f\t",cluster->repeat);
	if(cluster->repeat && cluster->repeatName) fprintf(fp,"%s",cluster->repeatName);
}
fprintf(fp,"\n");
print_matrix(cluster, fp);
check(index = (float*)malloc(num*sizeof(float)));
for(i=0;i<num;++i){ index[i] = i; }
if(num > 2){
	float *input, match=0.5, *pairSim, *s, *ss;
	int ind1,ind2,dir1,dir2,off1,off2,nMatch=4,i1,i2,ipair,*pairMem1,*pairMem2,len1,len2;
	fprintf(fp,"Zvalue\n");
	for(i=0; i<clLen; ++i){
		fprintf(fp,"%.2f",sqrt(freq[i]));
		if(i<clLen-1) fprintf(fp,"\t");
	}
	fprintf(fp,"\n\n");
	check(input = (float*)calloc(num*num,sizeof(float)));
	for(i=0;i<num;++i){
		off1 = cluster->off[i];
		dir1 = cluster->dir[i];
		ind1 = cluster->memb[i];
		len1 = Matrix[ind1].len;
		for(j=i+1;j<num;++j){
			off2 = cluster->off[j];
			dir2 = cluster->dir[j];
			ind2 = cluster->memb[j];
			len2 = Matrix[ind2].len;
			common_ref(off1,dir1,Matrix[ind1].len,off2,dir2,Matrix[ind2].len,&off,&dir);
			if(off>len1-3 || off<3-len2){
				input[i*num+j] = 1;
				input[j*num+i] = 1;
			}else{
				match = matrix_match(&Matrix[ind1],&Matrix[ind2],off,dir,&nMatch);
				input[i*num+j] = 1/(match*nMatch+1);
				input[j*num+i] = input[i*num+j];
			}
		}
	}
	hierarchical_clustering(input,num,cluster->memb,index,NULL);
}
free(freq);
if(num > 1){
	char *pattern;
	fprintf(fp,"members\n");
	for(i=0; i<num; ++i){
		int imem, ind, start, dir;
		imem = index[i];
		ind = cluster->memb[imem];
		start = cluster->off[imem];
		dir = cluster->dir[imem];
		len = Matrix[ind].len;
		//printf("%d %d %d %d %d %d\n",icluster,start,len,clLen,imem,ind);
		fprintf(fp,"%s\t",Matrix[ind].name);
		for(j=0; j<start; ++j){
			fprintf(fp," ");
		}
		pattern = Matrix[ind].pattern;
		if(dir<0) pattern = Matrix[ind].patternRev;
		fprintf(fp,"%s\t",pattern);
		if(p->headerIndex[3]>=0) fprintf(fp,"%d",Matrix[ind].freq);
		fprintf(fp,"\t");
		if(p->headerIndex[4]>=0) fprintf(fp,"%.4f",Matrix[ind].ratio);
		fprintf(fp,"\t%d\t%d\t",dir,start);
		if(p->headerIndex[9]>=0) fprintf(fp,"%.2f\t",Matrix[ind].repeat);
		if(Matrix[ind].repeatName) fprintf(fp,"%s",Matrix[ind].repeatName);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}
free(index);
return;
}

/***********************************************/
void  print_intercluster(FILE *fp, MATRIX **clusters, float *sorted, int nClusters, MATRIX *Matrix, int nPattern,
		PAIR *posPairs, int nPosPairs, float matchThresh)
/***********************************************/
{
int i,j,i1,j1,k,dir,offset,offset1,ipair;
int found, found1=0, nMatch, len, len1;
float match;

found = 0;
for(ipair=0; ipair<nPosPairs; ++ipair){
	int ind1, ind2, cl1, cl2, overlay;
	ind1 = posPairs[ipair].ind1;
	ind2 = posPairs[ipair].ind2;
	offset = posPairs[ipair].off;
	nMatch = align_matrixes_best(&Matrix[ind1],&Matrix[ind2],&match,&offset1,&dir,matchThresh);
	if(offset == offset1) continue;
	overlay = Matrix[ind1].len-offset;
	if(posPairs[ipair].dir < 0) overlay = Matrix[ind2].len+offset;
	if(overlay >= 6) continue;
	if(!found){
		fprintf(fp,"Intercluster links\n");
		fprintf(fp,"Name1\tPattern1\tName2\tPattern2\tStrand\tOffset\tRelativeFreq\tFrequency\n");
		found = 1;
	}
	fprintf(fp,"%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%d\n",Matrix[ind1].name,Matrix[ind1].pattern,Matrix[ind2].name,Matrix[ind2].pattern,posPairs[ipair].dir,offset,posPairs[ipair].match,posPairs[ipair].num);
}
if(found) fprintf(fp,"\n");
return;
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
	ch = strchr(ch,'\t');
	if(ch) *ch = '\0';
	else break;
	if(i>=num-1) break;
	ch++;
	i++;
}
return(i+1);
}

/***********************************************/
void  parse_headers (char *headers, PARAM *p)
/***********************************************/
{
char *items[10], *ch;
int n, i, j;

n = split_string(headers,items,10);
for(j=0;j<NHEADERS;++j) p->headerIndex[j] = -1;
for(i=0;i<n;++i){
	for(j=0;j<NHEADERS;++j){
		if(!strcmp(items[i],headerItems[j])){
			p->headerIndex[j] = i;
		}
	}
}
return;
}

/***********************************************/
int  read_file (char *filename, MATRIX **Mp, PARAM *p)
/***********************************************/
{
char *buffer, *ch, *items[10];
FILE *fp;
MATRIX *M;
int N=0, i,pos,NNN;
float m[4];

check(buffer = (char*)malloc(3201*sizeof(char)));
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("Input file %s not found",ch); exit(0); }
while(fgets(buffer,3200,fp)){
	if(strstr(buffer,"Parameters")==buffer){
		p->param = copy_string(&buffer[12]);
		p->param[strlen(p->param-1)]='\0';
	}else if(strstr(buffer,"Headers")==buffer){
		parse_headers(&buffer[9],p);
	}
	if(buffer[0] == '>') N++;
}
rewind(fp);
if(!N) error_message("No input motifs found!");
NNN = N;
if(NNN > p->N_matrix){
	NNN = p->N_matrix;
	printf("Maximum number of motifs = %d; others are ignored\n",p->N_matrix);
}
check(*Mp = (MATRIX*)calloc(NNN,sizeof(MATRIX)));
M = *Mp;
N = 0;
float log2 = log(2.0);
fgets(buffer,3200,fp);
while(1){
	int done=0;
	float sd;
	while(buffer[0] != '>' && !done){
		if(!fgets(buffer,32000,fp)) done=1;
	}
	if(N>=NNN || done){ break; }

	int n = split_string(&buffer[1],items,10);
	if(p->headerIndex[9]>=0) sscanf(items[p->headerIndex[9]],"%f",&M[N].repeat);
	if(M[N].repeat <= p->max_repeat){
		if(p->headerIndex[0]>=0) M[N].name=copy_string(items[p->headerIndex[0]]);
		else M[N].name=copy_string(items[0]);
		if(p->headerIndex[3]>=0) sscanf(items[p->headerIndex[3]],"%ld",&M[N].freq);
		else M[N].freq=100;
		if(p->headerIndex[4]>=0) sscanf(items[p->headerIndex[4]],"%f",&M[N].ratio);
		else M[N].ratio=2;
		if(p->headerIndex[6]>=0) sscanf(items[p->headerIndex[6]],"%f",&M[N].score);
		if(p->headerIndex[7]>=0) sscanf(items[p->headerIndex[7]],"%f",&M[N].p);
		if(p->headerIndex[8]>=0) sscanf(items[p->headerIndex[8]],"%f",&M[N].FDR);
		if(M[N].repeat && p->headerIndex[9]<n-1)
			M[N].repeatName=copy_string(items[p->headerIndex[9]+1]);
	}
	pos=0;
	M[N].info=0;
	while(fgets(buffer,3200,fp)){
		float sum=0, max=0, sum1=0;
		if(strlen(buffer)<3 || buffer[0] == '>'){ break; }
		if(pos==149){ printf("WARNING: Motif %s is too long (>150), truncating.\n",M[N].name); ++pos; }
		if(pos>149 || M[N].repeat > p->max_repeat){ continue; }
		sscanf(buffer,"%d%f%f%f%f",&i,&m[0],&m[1],&m[2],&m[3]);
		if(i!=pos){ printf("Wrong file format in line: %s",buffer); exit(0); }
		for(i=0; i<4; ++i){
			sum += m[i];
		}
		float pseudocount=sum*0.01;
		if(!pseudocount || pseudocount>1) pseudocount=1;
		sum = 0;
		for(i=0; i<4; ++i){
			if(m[i]<0) error_message("Negative value in PFM");
			m[i] += pseudocount;
			if(max<m[i]) max=m[i];
			sum1 += m[i];
			M[N].m[pos][i] = log(m[i]); 
			sum += M[N].m[pos][i];
		}
		int x=0;
		int coef=1;
		float entropy=0;
		for(i=0; i<4; ++i){
			M[N].m[pos][i] -= sum/4;
			float pp = m[i]/sum1;
			entropy -= pp*log(pp)/log2;
			int y = 1;
			if(m[i] < max/2){ y=0; }
			x += y*coef;
			coef *= 2;
		}
		M[N].info += 2.0-entropy;
		strcat(M[N].pattern,codes[x]);
		strcat(M[N].patternRev,revcodes[x]);
		++pos;
	}
	if(M[N].repeat > p->max_repeat) continue;
	reverse_string(M[N].patternRev);
	M[N].len = pos;
	if(!M[N].score){
		M[N].score = M[N].info;
		if(M[N].ratio) M[N].score += (M[N].ratio - 1);
		if(M[N].freq > 0) M[N].score += sqrt(M[N].freq);
	}
	M[N].num=1;
	check(M[N].memb = (int*)malloc(sizeof(int)));
	check(M[N].off = (int*)calloc(1,sizeof(int)));
	check(M[N].dir = (int*)malloc(sizeof(int)));
	M[N].memb[0] = N;	
	M[N].dir[0] = 1;	
	N++;
}
fclose(fp);
printf("File %s loaded. Input patterns: %d\n",ch,N);
free(buffer);
return(N);
}

/***********************************************/
void   sort_positions(POSITION *Pp, int nSeq)
/***********************************************/
{
int iseq, i, j, nmatch, *dir, *TFid;
float *pos, *index;

check(pos = (float*)malloc(MAXINT*sizeof(float)));
check(index = (float*)malloc(MAXINT*sizeof(float)));
check(dir = (int*)malloc(MAXINT*sizeof(int)));
check(TFid = (int*)malloc(MAXINT*sizeof(int)));
for(iseq=0;iseq<nSeq;++iseq){
	POSITION *PP = &Pp[iseq];
	nmatch = PP->nmatch;
	memmove(dir,PP->dir,nmatch*sizeof(int));
	memmove(TFid,PP->TFid,nmatch*sizeof(int));
	for(i=0;i<nmatch; ++i){
		pos[i] = PP->pos[i];
		index[i] = i;
	}
	sortem(nmatch, pos, 1, index, NULL, NULL, NULL, NULL, NULL, NULL);
	for(i=0;i<nmatch; ++i){
		j = index[i];
		PP->pos[i] = pos[i];
		PP->dir[i] = dir[j];
		PP->TFid[i] = TFid[j];
		//printf("%d %d %d %d\n",iseq,PP->pos[i],PP->dir[i],PP->TFid[i]);
	}
}
free(pos);
free(index);
free(dir);
free(TFid);
return;
}

/***********************************************/
int  read_positions (char *filename, POSITION **Pp, MATRIX *Matrix, int nPattern)
/***********************************************/
{
char *buffer,*junk,*pos,*motifName,*oldName,*items[10],*ch;
FILE *fp;
POSITION *PP;
int N=0, i, j, NNN,*uploaded,num, ind=0;
float m[4];

check(buffer = (char*)calloc(3300,sizeof(char)));
check(motifName = (char*)calloc(500,sizeof(char)));
check(oldName = (char*)calloc(500,sizeof(char)));
check(junk = (char*)malloc(500*sizeof(char)));
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("ERROR: Input file %s not found",ch); exit(0); }
while(!strstr(buffer,"Parameters")){
	if(!fgets(buffer,3299,fp)) break;
}
int n = split_string(buffer,items,10);
for(i=0;i<n;++i){
	if(strstr(items[i],"Nseq")){
		sscanf(&items[i][5],"%d",&NNN);
		break;
	}
}
if(!NNN) error_message("Nseq not specified in position file");
check(*Pp = (POSITION*)calloc(NNN,sizeof(POSITION)));
PP = *Pp;
for(i=0;i<NNN;++i){
	int k = 200;
	PP[i].alloc = k;
	PP[i].nmatch = 0;
	check(PP[i].dir = (int*)malloc(k*sizeof(int)));
	check(PP[i].TFid = (int*)malloc(k*sizeof(int)));
	check(PP[i].pos = (long*)malloc(k*sizeof(long)));
}
while(fgets(buffer,3299,fp)){
	int TFid, SEQid, dir, k;
	long pos;
	if(strstr(buffer,"Headers") || strstr(buffer,"Parameters")) continue;
	sscanf(buffer,"%s%d%s%d%ld",motifName,&SEQid,junk,&dir,&pos);
	if(PP[SEQid].nmatch >= PP[SEQid].alloc){
		PP[SEQid].alloc += 300;
		k = PP[SEQid].alloc;
		check(PP[SEQid].dir = (int*)realloc(PP[SEQid].dir,k*sizeof(int)));
		check(PP[SEQid].TFid = (int*)realloc(PP[SEQid].TFid,k*sizeof(int)));
		check(PP[SEQid].pos = (long*)realloc(PP[SEQid].pos,k*sizeof(long)));
	}
	if(strcmp(motifName,oldName)){
		ind = -1;
		for(i=0;i<nPattern;++i){
			if(!strcmp(motifName,Matrix[i].name)){
				ind=i;
				break;
			}
		}
		strcpy(oldName,motifName);
	}
	if(ind == -1) continue;
	k = PP[SEQid].nmatch;
	PP[SEQid].dir[k] = dir;
	PP[SEQid].pos[k] = pos;
	PP[SEQid].TFid[k] = ind;
	++PP[SEQid].nmatch;
}
printf("File %s loaded. Input sequences: %d\n",ch,NNN);
free(buffer);
free(junk);
free(motifName);
free(oldName);
return(NNN);
}

/***********************************************/
void   common_ref  (long pos1, int dir1, int len1, long pos2, int dir2, int len2, long *pos3, int *dir3)
/***********************************************/
{
if(dir1>0){
	*dir3 = dir2;
	*pos3 = pos2-pos1;
}else{
	*dir3 = -dir2;
	*pos3 = (pos1+len1)-(pos2+len2);
}
return;
}

/***********************************************/
int  add_to_cluster  (MATRIX *cl, int ind1, int ind2, int dir, int off, MATRIX *Matrix, float threshold)
/***********************************************/
{
int i, off1, dir1, off2, dir2, index=0, nMatch;

while(index < cl->num && ind2 != cl->memb[index]) ++index;
if(index < cl->num) return(0);
if(cl->alloc < 0 || cl->num<0) error_message("Add_to_cluster");
if(cl->alloc <= cl->num){
	cl->alloc += 20;
	if(!cl->num){
		check(cl->memb = (int*)malloc(cl->alloc*sizeof(int)));
		check(cl->off = (int*)malloc(cl->alloc*sizeof(int)));
		check(cl->dir = (int*)malloc(cl->alloc*sizeof(int)));
	}else{
		check(cl->memb = (int*)realloc(cl->memb,cl->alloc*sizeof(int)));
		check(cl->off = (int*)realloc(cl->off,cl->alloc*sizeof(int)));
		check(cl->dir = (int*)realloc(cl->dir,cl->alloc*sizeof(int)));
	}
}
if(ind1 < 0){
	if(off<0){
		for(i=0; i<cl->num; ++i) cl->off[i] -= off;
		off = 0;
	}
	cl->memb[cl->num] = ind2;
	cl->dir[cl->num] = dir;
	cl->off[cl->num] = off;
	cl->num++;
	return(1);
}
float match;
index = 0;
while(index < cl->num && ind1 != cl->memb[index]) ++index;
if(index == cl->num) error_message("member not found");
off1 = cl->off[index];
dir1 = cl->dir[index];
dir2 = dir*dir1;
if(dir1>0){
	off2 = off+off1;
}else{
	off2 = off1+Matrix[ind1].len-Matrix[ind2].len-off;
	//printf("%d %d %d %d %d\n",off2,off1,Matrix[ind1].len,Matrix[ind2].len,off);
}
match = matrix_match(cl, &Matrix[ind2], off2, dir2, &nMatch);
if(DEBUG) printf("Match = %.4f %d %d %d %d %d %d\n",match,dir,off,dir2,off2,dir1,off1);
if(match >= threshold){
	if(off2<0){
		for(i=0; i<cl->num; ++i) cl->off[i] -= off2;
		off2 = 0;
	}
	cl->memb[cl->num] = ind2;
	cl->dir[cl->num] = dir2;
	cl->off[cl->num] = off2;
	cl->num++;
	return(1);
}
return(0);
}

/***********************************************/
PAIR  *count_TFpairs  (POSITION *Pp, int nSeq, MATRIX *Matrix, int nPattern, int maxDistance, int *nPosPairs)
/***********************************************/
{
int i, j, k, dir1, dir2, dir3, TFid1, TFid2, len1, len2, nPairs=0, **TFpair, **TFindex, ipair;
int *pair2TF1, *pair2TF2, n, **TFpairPos;
long iseq, pos1, pos2, pos3, nPairAlloc, nPairs1, nPairs2;
float *num, *match, *TF1, *TF2, *dir, *offset;
PAIR *pairs;

check(TFpair = (int**)malloc(nPattern*sizeof(int*)));
for(i=0; i<nPattern; ++i){
	check(TFpair[i] = (int*)calloc(nPattern,sizeof(int)));
}

for(iseq=0;iseq<nSeq;++iseq){
	POSITION *PP = &Pp[iseq];
	int nmatch = PP->nmatch;
	for(i=0;i<nmatch; ++i){
		pos1 = PP->pos[i];
		dir1 = PP->dir[i];
		TFid1 = PP->TFid[i];
		len1 = Matrix[TFid1].len;
		for(j=i+1;j<nmatch; ++j){
			pos2 = PP->pos[j];
			dir2 = PP->dir[j];
			TFid2 = PP->TFid[j];
			len2 = Matrix[TFid2].len;
			common_ref(pos1,dir1,len1,pos2,dir2,len2,&pos3,&dir3);
			//printf("%d - %d %d %d - %d %d %d - %d %d\n",iseq,pos1,dir1,len1,pos2,dir2,len2,pos3,dir3);
			if(pos3< -maxDistance-10 || pos3 > maxDistance+10) break;
			if(pos3< -maxDistance || pos3 > maxDistance) continue;
			TFpair[TFid1][TFid2]++;
			TFpair[TFid2][TFid1]++;
			//printf("%d %d %d %d %d %d %d\n",iseq,i,TFid1,TFid2,pos3,dir3,nPattern);
		}
	}
}

check(num = (float*)malloc(MAXINT*sizeof(float)));
check(TF1 = (float*)malloc(MAXINT*sizeof(float)));
check(TF2 = (float*)malloc(MAXINT*sizeof(float)));
for(i=0;i<nPattern; ++i){
	for(j=i;j<nPattern; ++j){
		int n = TFpair[i][j];
		if(i==j) n /= 2;
		if(n < 5) continue;
		num[nPairs] = n;
		TF1[nPairs] = i;
		TF2[nPairs] = j;
		++nPairs;
		if(nPairs>=MAXINT) break;
	}
	if(nPairs>=MAXINT) break;
}
sortem(nPairs, num, 2, TF1, TF2, NULL, NULL, NULL, NULL, NULL);

/* Fill up index table */
TFindex = TFpair;
for(i=0; i<nPattern; ++i){
	for(j=0; j<nPattern; ++j){
		TFindex[i][j] = 0;
	}
}
check(pair2TF1 = (int*)malloc(nPairs*sizeof(int)));
check(pair2TF2 = (int*)malloc(nPairs*sizeof(int)));
for(i=nPairs-1; i>=0; --i){
	j = TF1[i];
	k = TF2[i];
	TFindex[j][k] = nPairs-i;
	TFindex[k][j] = nPairs-i;
	pair2TF1[nPairs-i-1] = j;
	pair2TF2[nPairs-i-1] = k;
	//printf("%d %d %d %d %s %s\n",i,j,k,(int)num[i],Matrix[j].pattern,Matrix[k].pattern);
	//if(nPairs-i>50)exit(0);
}

/* Second search for exact position */
check(TFpairPos = (int**)malloc(nPairs*sizeof(int*)));
len1 = (maxDistance*2+1)*2;
for(i=0; i<nPairs; ++i){
	check(TFpairPos[i] = (int*)calloc(len1,sizeof(int)));
}

for(iseq=0;iseq<nSeq;++iseq){
	POSITION *PP = &Pp[iseq];
	int nmatch = PP->nmatch;
	for(i=0;i<nmatch; ++i){
		pos1 = PP->pos[i];
		dir1 = PP->dir[i];
		TFid1 = PP->TFid[i];
		len1 = Matrix[TFid1].len;
		for(j=i+1;j<nmatch; ++j){
			int ipair, iii;
			TFid2 = PP->TFid[j];
			if(!TFindex[TFid1][TFid2]) continue;
			ipair = TFindex[TFid1][TFid2]-1;
			pos2 = PP->pos[j];
			dir2 = PP->dir[j];
			len2 = Matrix[TFid2].len;
			common_ref(pos1,dir1,len1,pos2,dir2,len2,&pos3,&dir3);
			if(TFid2<TFid1) pos3 = -pos3;
			if(pos3< -maxDistance || pos3 > maxDistance) continue;
			if(pos3< -maxDistance-10 || pos3 > maxDistance+10) break;
			iii = ((dir3+1)/2)*(maxDistance*2+1)+pos3+maxDistance;
			//printf("%d %d %d %d %d %d %d %d - %d %d %d %d\n",iseq,i,j,k,dir3,pos3,ipair,iii, pos1,dir1,pos2,dir2);
			TFpairPos[ipair][iii]++;
		}
	}
}

nPairs1 = 0;
check(dir = (float*)malloc(MAXINT*sizeof(float)));
check(offset = (float*)malloc(MAXINT*sizeof(float)));
check(match = (float*)malloc(MAXINT*sizeof(float)));
for(ipair=0; ipair<nPairs; ++ipair){
	int dir0, off, iii, TFid1, TFid2, freqMin;
	TFid1 = pair2TF1[ipair];
	TFid2 = pair2TF2[ipair];
	int numMin = Matrix[TFid1].freq;
	if(numMin > Matrix[TFid2].freq) numMin = Matrix[TFid2].freq;
	for(dir0=0; dir0<2; ++dir0){
		for(off=-maxDistance;off<=maxDistance; ++off){
			int n0, n1;
			iii = dir0*(maxDistance*2+1)+off+maxDistance;
			n = TFpairPos[ipair][iii];
			if(n >= 5 && n >= 0.05*numMin){
				TF1[nPairs1] = TFid1;
				TF2[nPairs1] = TFid2;
				match[nPairs1] = (float)n/numMin;
				num[nPairs1] = n;
				dir[nPairs1] = dir0*2-1;
				offset[nPairs1] = off;
				++nPairs1;
				if(nPairs1>=MAXINT) break;
			}
			if(nPairs1>=MAXINT) break;
		}
		if(nPairs1>=MAXINT) break;
	}
	if(nPairs1>=MAXINT) break;
}

sortem(nPairs1, match, 5, TF1, TF2, dir, offset, num, NULL, NULL);
check(pairs = (PAIR*)malloc(nPairs1*sizeof(PAIR)));
for(ipair=0; ipair<nPairs1; ++ipair){
	i = nPairs1-1-ipair;
	if(num[i] < 0.05) break;
	pairs[ipair].match = match[i];
	pairs[ipair].num = num[i];
	pairs[ipair].ind1 = TF1[i];
	pairs[ipair].ind2 = TF2[i];
	pairs[ipair].dir = dir[i];
	pairs[ipair].off = offset[i];
}
nPairs2 = ipair;
check(pairs = (PAIR*)realloc(pairs,nPairs2*sizeof(PAIR)));
*nPosPairs=nPairs2;
for(i=0; i<nPairs; ++i){
	free(TFpairPos[i]);
}

free(TFpairPos);
for(i=0; i<nPattern; ++i){
	free(TFpair[i]);
}
free(TFpair);
free(match);
free(num);
free(TF1);
free(TF2);
free(dir);
free(offset);
free(pair2TF1);
free(pair2TF2);
return(pairs);
}

/******************************************/
void combine_clusters  (float *a, int *weight, int n, int c1, int c2)
/******************************************/
{
int i, j, i1, j1;
float w1, w2, w3, x1, x2;
float *a_new;

if(c1>=n || c2>=n || c1==c2 || n<2) error_message("Wrong cluster number");
a_new = (float*)calloc((n-1)*(n-1),sizeof(float));
for(i=0;i<n-1;++i){
	a_new[i*(n-1)+(n-2)]=0;
	a_new[(n-2)*(n-1)+i]=0;
}
w3 = weight[c1]+weight[c2];
w1 = (float)weight[c1]/w3;
w2 = 1.0 - w1;
for(i=0;i<n;++i){
	i1 = i;
	if(i>c1) i1--; 
	if(i>c2) i1--;
	if(i==c1 || i==c2) continue;
	weight[i1] = weight[i];
	for(j=0;j<n;++j){
		j1 = j;
		if(j>c1) j1--;
		if(j>c2) j1--;
		if(j==c1 || j==c2) continue;
		a_new[i1*(n-1)+j1] = a[i*n+j];
	}
	x1 = a[i*n+c1];
	x2 = a[i*n+c2];
	a_new[i1*(n-1)+(n-2)] = x1*w1 + x2*w2;
	a_new[(n-2)*(n-1)+i1] = a_new[i1*(n-1)+(n-2)];
}
weight[n-2] = w3;
for(i=0;i<(n-1)*(n-1);++i)
	a[i] = a_new[i];
free(a_new);
return;
}

/******************************************/
void  add_score(int n, int k, CLUSTER *c, float *score)
/******************************************/
{
long nstack, clust;
long *stack;
long c1,c2,add,score1,score2,swap;

stack = (long*)malloc(2*n*sizeof(long));
nstack=0;

c1 = c[k].c1;
score1 = 0;
if(c1 >= n){
	c1-=n;
	score1 = c[c1].maxscore;
}
score2 = 0;
c2 = c[k].c2;
if(c2 >= n){
	c2-=n;
	score2 = c[c2].maxscore;
}
if(score1 < score2){
	swap = c1;
	c1 = c2;
	c2 = swap;
	swap = score1;
	score1 = score2;
	score2 = swap;
}
add = score1+1;
if(score2==0)
	score[c2] += add;
else
	stack[nstack++] = c2;
c[k].maxscore = add + score2;
while(nstack){
	clust = stack[--nstack];
	c1 = c[clust].c1;
	c2 = c[clust].c2;
	if(c1 >= n){
		stack[nstack++] = c1-n;
	}else{
		score[c1] += add;
	}
	if(c2 >= n){
		stack[nstack++] = c2-n;
	}else{
		score[c2] += add;
	}
}
free(stack);
return;
}

/******************************************/
void  hierarchical_clustering (float *input, int n, int *select2pattern, float *index1, FILE *fp)
/******************************************/
{
long i, j, n1, i1, c1=0, c2=0;
float *b, distance, s=1, maxx=0;
int *weight;
float *score;
int *icluster;
CLUSTER *c;

b = (float*)calloc(n*n,sizeof(float));
weight = (int*)malloc(n*sizeof(int));
score = (float*)calloc(n, sizeof(float));
icluster = (int*)malloc(2*n*sizeof(int));
c = (CLUSTER*)calloc(n,sizeof(CLUSTER));
for(i=0;i<n*n;++i){
	if(maxx < input[i]) maxx = input[i];
}
for(i=0;i<n;++i){
	c[i].c1 = -1;
	weight[i] = 1;
	icluster[i] = i;
	icluster[i+n] = i+n;
	for(j=i+1;j<n;++j){
		b[i*n+j] = input[i*n+j]/maxx;
		b[j*n+i] = b[i*n+j];
	}
}
n1 = n;
while(n1 > 1){
	distance = 1000;
	for(i=0;i<n1;++i){
		for(j=i+1;j<n1;++j){
			if(distance > b[i*n1+j]){
				distance = b[i*n1+j];
				c1 = i;
				c2 = j;
			}
		}
	}
	c[n-n1].c1 = icluster[c1];
	c[n-n1].c2 = icluster[c2];
	c[n-n1].distance = distance;
	add_score(n, n-n1, c, score);
	combine_clusters(b, weight, n1, c1, c2);
	for(i=0;i<2*n;++i){
		i1 = i;
		if(i>c1) i1--; 
		if(i>c2) i1--;
		icluster[i1] = icluster[i];
	}
	--n1;
}
for(i=0;i<n;++i) index1[i] = i;
sortem(n, score, 1, index1, NULL, NULL, NULL, NULL, NULL, NULL);
for(i=0; i<n; ++i){
	j = index1[i];
	if(fp) fprintf(fp,"%d\t%d\n",(int)index1[i],select2pattern[j]);
	//printf("%d\t%d\n",(int)index1[i],select2pattern[j]);
}
for(i=0; c[i].c1 >= 0; ++i)
	if(fp) fprintf(fp,"%d\t%d\t%d\t%.4f\n", n+i, c[i].c1, c[i].c2, c[i].distance);
free (weight);
free (icluster);
free (score);
free (b);
free (c);
return;
}

/******************************************/
PAIR  *get_similarity_pairs (MATRIX *Matrix, int nPattern, int *nPairs1, float matchThresh)
/******************************************/
{
int i, j, ind, ind1, nPairAlloc = 1000;
int len, len1, nPairs=0, offset, dir, nMatch;
PAIR *pairs;
float infoRatio, match;

check(pairs = (PAIR*)malloc(nPairAlloc*sizeof(PAIR)));
for(ind=0; ind<nPattern; ++ind){
	for(ind1=ind+1; ind1<nPattern; ++ind1){
		nMatch = align_matrixes_best(&Matrix[ind],&Matrix[ind1],&match,&offset,&dir,matchThresh);
		if(match < matchThresh){ continue; }
		if(nPairs > nPairAlloc-2){
			nPairAlloc += 1000;
			check(pairs = (PAIR*)realloc(pairs,nPairAlloc*sizeof(PAIR)));
		}
		pairs[nPairs].ind1=ind;
		pairs[nPairs].ind2=ind1;
		pairs[nPairs].match=match;
		pairs[nPairs].dir=dir;
		pairs[nPairs].off=offset;
		pairs[nPairs].nMatch=nMatch;
		++nPairs;
		//printf("%d %d %d %d %.4f %s %s\n",ind,ind1,dir,offset,match,Matrix[ind].pattern,Matrix[ind1].pattern);
	}
}
*nPairs1 = nPairs;
return(pairs);
}

/******************************************/
int  **get_linkage (PAIR *pairs, int nPairs, int *nLinked, int nPattern)
/******************************************/
{
int i, j, ind, ind1, ind2;
int **linked, *nalloc;

check(linked = (int**)malloc(nPattern*sizeof(int*)));
check(nalloc = (int*)malloc(nPattern*sizeof(int)));
for(ind=0; ind<nPattern; ++ind){
	nalloc[ind] = 20;
	check(linked[ind] = (int*)malloc(nalloc[ind]*sizeof(int)));
	nLinked[ind]=0;
}
for(i=0; i<nPairs; ++i){
	ind = pairs[i].ind1;
	ind1 = pairs[i].ind2;
	linked[ind][nLinked[ind]++] = ind1;
	linked[ind1][nLinked[ind1]++] = ind;
	if(nLinked[ind]==nalloc[ind]){
		nalloc[ind] += 50;
		check(linked[ind] = (int*)realloc(linked[ind],nalloc[ind]*sizeof(int)));
	}
	if(nLinked[ind1]==nalloc[ind1]){
		nalloc[ind1] += 50;
		check(linked[ind1] = (int*)realloc(linked[ind1],nalloc[ind1]*sizeof(int)));
	}
}
free(nalloc);
return(linked);
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
PARAM *read_parameters (int nargs, char **argv)
/***********************************************/
{
PARAM *p;
int iarg=1, score=-1,i;

if(nargs < 4){ printf("patternCluster -i inputMotifFile -o outputFile [-pos positionFile, -match matchThreshold, -n numberOfMotifs, -repeat maxRepeatEnrichment, -posonly]\n"); exit(0); }
check(p = (PARAM*)calloc(1,sizeof(PARAM)));
p->N_matrix = 5000;
p->matchThresh=0.85;
p->position_only=0;
p->max_repeat = 1000;
while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1) p->inputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-pos") && iarg < nargs-1) p->positionFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-match") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->matchThresh);
	else if(!strcmp(argv[iarg],"-n") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->N_matrix);
	else if(!strcmp(argv[iarg],"-posonly")) p->position_only=1;
	else if(!strcmp(argv[iarg],"-repeat") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->max_repeat);
	else{
		printf("Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}
for(i=0; i<NHEADERS; ++i){ p->headerIndex[i]=-1; }
if(p->N_matrix <1) error_message("Wrong number of output matrixes shown after -n option");
if(p->matchThresh <0.5 || p->matchThresh>0.99) error_message("Match threshold should be from 0.5 to 0.99");
if(p->position_only && !p->positionFile) error_message("Position file needs to be specified");
return(p);
}

/***********************************************/
int main (int argc, char **argv) 
/***********************************************/
{
MATRIX *Matrix=NULL, **clusters;
PARAM *p;
PAIR *simPairs, *posPairs;
FILE *fp;
int nPattern, ind, ind1, ind2, nSimPairs=0, len, len1, i, j, nSeq=0, nPosPairs=0;
int icluster=1, nClusters=0, *matrixCluster, **linked, *nLinked, *remains, nRemains=0, *members, Nmembers=0;
int *remainsAdd, nRemainsAdd=0, ipair, *stack, nStack=0;
float infoRatio, match, *input;
int offset, dir, count;
float *sorted, *score, maxScore=0;
POSITION *Position;

p = read_parameters(argc, argv);
nPattern = read_file(p->inputFile,&Matrix,p);
if(nPattern<2) error_message("Input motif file is empty");
get_self_similarity(Matrix,nPattern,p);
simPairs = get_similarity_pairs(Matrix,nPattern,&nSimPairs,p->matchThresh);
if(p->positionFile){
	nSeq = read_positions (p->positionFile,&Position,Matrix,nPattern);
	sort_positions(Position,nSeq);
	posPairs = count_TFpairs(Position,nSeq,Matrix,nPattern,MAXDIST,&nPosPairs);
}

check(nLinked = (int*)calloc(nPattern,sizeof(int)));
if(p->position_only){
	linked = get_linkage(posPairs,nPosPairs,nLinked,nPattern);
}
else{
	linked = get_linkage(simPairs,nSimPairs,nLinked,nPattern);
}
check(matrixCluster = (int*)calloc(nPattern,sizeof(int)));
check(stack = (int*)calloc(MAXINT,sizeof(int)));
check(clusters = (MATRIX**)malloc(nPattern*sizeof(MATRIX*)));
check(remains = (int*)calloc(nPattern,sizeof(int)));
check(remainsAdd = (int*)calloc(nPattern,sizeof(int)));
check(members = (int*)malloc(nPattern*sizeof(int)));
for(ind=0; ind<nPattern; ++ind){
	check_palindrome(&Matrix[ind],p->matchThresh);
}
for(ind=0; ind<nPattern; ++ind){
	float avg_selfSim=0;
	Nmembers=0;
	if(matrixCluster[ind] || !nLinked[ind]){ continue; }
	stack[nStack++] = ind;
	while(nStack){
		ind1 = stack[--nStack];
		if(matrixCluster[ind1]){ continue; }
		avg_selfSim += Matrix[ind1].selfSim;
		members[Nmembers++]=ind1;
		matrixCluster[ind1] = icluster;
		for(i=0; i<nLinked[ind1]; ++i){
			ind2 = linked[ind1][i];
			if(!matrixCluster[ind2]){
				int i1, found=0;
				for(i1=0;i1<nStack && !found;++i1){
					if(stack[i1]==ind2) found=1;
				}
				if(!found){
					stack[nStack++] = ind2;
				}
			}
		}
	}
	nRemains=0;
	while(1){
		MATRIX *mp;
		int methodUsed=0;
		avg_selfSim /= Nmembers;
		if(p->positionFile && (p->position_only || avg_selfSim >=0.5)){
			mp = align_cluster_pos(members,Nmembers,remains,&nRemains,Matrix,nPattern,posPairs,nPosPairs,p->matchThresh);
			methodUsed=1;
		}else{
			mp = align_cluster_sim(members,Nmembers,remains,&nRemains,Matrix,nPattern,simPairs,nSimPairs,p->matchThresh);
		}
		if(mp){
			//printf("%d %d %d\n",icluster,mp->num,nRemains);
			check_palindrome(mp,p->matchThresh);
			mp->method = methodUsed;
			clusters[nClusters++] = mp;
		}else{
			for(i=0; i<Nmembers; ++i){
				clusters[nClusters++] = &(Matrix[members[i]]);
			}
			break;
		}
		avg_selfSim=0;
		if(nRemains && nRemains < Nmembers){
			for(i=0; i<nRemains; ++i){
				members[i] = remains[i];
				avg_selfSim += Matrix[members[i]].selfSim;
			}
			Nmembers = nRemains;
		}else{
			break;
		}
	}
	++icluster;
}

for(ind=0; ind<nPattern; ++ind){
	if(!matrixCluster[ind]){
		clusters[nClusters++] = &(Matrix[ind]);
	}
}
printf("Pairs done. N clusters = %d\n",nClusters);
check(sorted = (float*)malloc(nClusters*sizeof(float)));
check(score = (float*)malloc(nClusters*sizeof(float)));
for(i=0; i<nClusters; ++i){
	score[i] = clusters[i]->score*sqrt(clusters[i]->num);
	sorted[i] = i;
}
sortem(nClusters, score, 1, sorted, NULL, NULL, NULL, NULL, NULL, NULL);
fp = fopen(p->outputFile,"w");
if(!fp) error_message("Output file not opened");
fprintf(fp,"Parameters:\tClustersGeneratedBy=patternCluster(1.0)\t%s",p->param);
int method1=0;
if(p->positionFile) method1++;
if(p->position_only) method1++;
fprintf(fp,"\tmatchThresh=%.4f\tMethod=%d\tMaxRepeat=%.4f\n",p->matchThresh,method1,p->max_repeat);
fprintf(fp,"Headers:\tName\tPattern\tPatternRev");
if(p->headerIndex[3]>=0) fprintf(fp,"\tFreq");
if(p->headerIndex[4]>=0) fprintf(fp,"\tRatio");
fprintf(fp,"\tInfo\tScore");
if(p->headerIndex[7]>=0) fprintf(fp,"\tp");
if(p->headerIndex[8]>=0) fprintf(fp,"\tFDR");
fprintf(fp,"\tPalindrome\tNmembers\tMethod\tMotifName");
if(p->headerIndex[9]>=0) fprintf(fp,"\tRepeat\tRepeatType");
fprintf(fp,"\n");
count=0;
for(i=nClusters-1; i>=0; --i){
	j = sorted[i];
	clusters[j]->id = ++count;
	print_cluster(count,clusters,j,fp,Matrix,posPairs,nPosPairs,nPattern,p);
	//printf("%d %d %s %.4f\n",count,clusters[j]->num,clusters[j]->pattern,clusters[j]->score);
}
if(p->positionFile){
	print_intercluster(fp,clusters,sorted,nClusters,Matrix,nPattern,posPairs,nPosPairs,p->matchThresh);
}
fclose(fp);
return(0);
}



