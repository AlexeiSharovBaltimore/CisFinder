/*************************************
* patternTest is a part of CisFinder software
* http://lgsun.grc.nia.nih.gov/CisFinder
*
* Function: improve motifs using resampling method
* 
* Syntax:
* patternTest -i motifFile -f fastaFile -o outputFile [-c controlFile, 
* -prog progressFile, -maxlen maxSequenceLength, -strand strandOption, 
* -n numberOfMotifs, -method resampleMethod, -iter numIteartions, 
* -siter numSubiterations, -fp falsePositives, -userep, -one]
* 
* Comments:
* (a) If controlFile is not specified, then random sequence generated using 3rd order
* Markov chain is used as a control.
* (b) progressFile keeps records for each iteration and subiteration.
* (c) strandOption: 0 = search both strands, 1 = search positive strand, 2 = search
* positive strand and use negative strand as a control.
* (d) numberOfMotifs can be used to limit the number of input motifs (default=100)
* (e) resampleMethod: 1=regression method (default), 2=simple resampling, 3=difference method.
* (f) numIteartions = number of main iterations when a new set of motif matches 
* is generated (default = 3)
* (g) numSubiterations = number of sub-iterations which process the same set of motif
* matches, and only the match score changes (default = 3).
* (e) userep: use repeats in sequence (lower-case in sequence)
* (g) one: consider not more than 1 motif occurrence per sequence (with highest 
* match score)
* (h) falsePositives = number of expected false positives per 10000 bp in the
* control sequence (default = 5).
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
#define MAXSEQ 200000
#define MINDIST 15
#define SCORE 2
#define NHEADERS 9
#define NSTACK 60
#define METHOD_REGRESS 1
#define METHOD_RESAMPLE 2
#define METHOD_DIFF 3

int DEBUG=0;

typedef struct matrix_st{
	int id;
	int len;
	int freq;
	int offset;
	int num;  //number of sub-motifs
	float ratio;
	float info;
	float score;
	float repeat;
	float selfSim;	//self-similarity
	float scoreThresh;
	float scoreThreshCore;
	int threshDone;
	long threshMult;
	float FDR;
	float m[MAXLENGTH][4];
	char pattern[MAXLENGTH];
	char patternRev[MAXLENGTH];
	char name[MAXLENGTH];
}MATRIX;
typedef struct position_st{
	int num;
	int alloc;
	long *iseq;
	long *pos;
	int *dir;
	int *pattern;
	float *score;
}POSITION;
typedef struct lookup_st{
	int num;
	int alloc;
	int *imatr;
	int *dir;
	float *score;
}LOOKUP;
typedef struct param_st{
	char *inputFile;
	char *outputFile;
	char *testFile;
	char *controlFile;
	char *progressFile;
	int use_repeats;
	int strand;
	int presence_option;
	int N_matrix;
	int method;
	long max_length;
	int iterations;
	int subIterations;
	float falsePositives;	// Expected number of false positives per 10 Kb //
	int headerIndex[NHEADERS];
}PARAM;

static char *codes[16]={"O","A","C","M","G","R","S","V","T","W","Y","H","K","D","B","N"};
static char *headerItems[NHEADERS]={"Name","Pattern","Threshold","Freq","Ratio","Info","Score","FDR","Repeat"};

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

/*************************************/
long   complementary (long k)
/*************************************/
{
int pos[8], pos1[8], j, mult;

for(j=0; j<8; ++j){
	pos[j] = k%4;
	k /= 4;
}
for(j=0; j<8; ++j){
	pos1[8-j-1] = 3 - pos[j];
}
k = 0;
mult = 1;
for(j=0; j<8; ++j){
	k += pos1[j]*mult;
	mult *= 4;
}
return(k);
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
if(len) memcpy(str,newstr,len);
free(newstr);
}

/*************************************/
void  regression (float *yy, float *xx, int n, float *a, float *b)
/*************************************/
{
float sx=0,sy=0,syx=0,sxx=0,syy=0,x,y;
int i;

for(i=0; i<n; ++i){
	sx += xx[i];
	sy += yy[i];
}
sx /= n;
sy /= n;
for(i=0; i<n; ++i){
	x = xx[i]-sx;
	y = yy[i]-sy;
	syx += x*y;
	sxx += x*x;
	syy += y*y;
}
*b = syx/sxx;
*a = sy - (*b)*sx;
return;
}

/*************************************/
void   print_matrix  (MATRIX *M, FILE *fp)
/*************************************/
{
int i,j;
float pp[4], sum;
double log10 = log(10.0);
for(i=0; i<M->len; ++i){
	fprintf(fp,"%d",i);
	sum = 0;
	for(j=0; j<4; ++j){
		pp[j] = exp(M->m[i][j]*log10);
		sum += pp[j];
	}
	for(j=0; j<4; ++j){
		pp[j] *= 1000/sum;
		fprintf(fp,"\t%.0f",pp[j]);
	}
	fprintf(fp,"\n");
}
fprintf(fp,"\n");
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
int  read_input (char *filename, MATRIX **Mp, PARAM *p)
/***********************************************/
{
char *buffer, *ch, *items[10];
FILE *fp;
MATRIX *M;
int N=0, i,pos,NNN,len;
float m[4];

check(buffer = (char*)malloc(3500*sizeof(char)));
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("Input file %s not found",ch); exit(0); }
for(i=0; i<2; ++i){
	fgets(buffer,3499,fp);
	if(strstr(buffer,"Headers")==buffer){
		parse_headers(&buffer[9],p);
	}
}
rewind(fp);
while(fgets(buffer,3499,fp)){
	if(buffer[0] == '>') N++;
}
rewind(fp);
if(N > p->N_matrix){
	N = p->N_matrix;
	printf("Maximum number of motifs = %d; others are ignored\n",p->N_matrix);
}
NNN = N;
check(*Mp = (MATRIX*)malloc(NNN*sizeof(MATRIX)));
M = *Mp;
N = 0;
double log10 = log(10);
while(fgets(buffer,3499,fp)){
	int done=0;
	float sd, sumMax=0;
	while(buffer[0] != '>' && !done){
		if(!fgets(buffer,3499,fp)) done=1;
	}
	if(N>=NNN || done){ break; }
	int n = split_string(buffer,items,10);
	if(p->headerIndex[0]>=0) strcpy(M[N].name,items[p->headerIndex[0]]);
	else strcpy(M[N].name,items[0]);
	if(p->headerIndex[1]>=0) strcpy(M[N].pattern,items[p->headerIndex[1]]);
	if(p->headerIndex[2]>=0) strcpy(M[N].patternRev,items[p->headerIndex[2]]);
	if(p->headerIndex[3]>=0) sscanf(items[p->headerIndex[3]],"%d",&M[N].freq);
	if(p->headerIndex[4]>=0) sscanf(items[p->headerIndex[4]],"%f",&M[N].ratio);
	if(p->headerIndex[5]>=0) sscanf(items[p->headerIndex[5]],"%f",&M[N].info);
	if(p->headerIndex[6]>=0) sscanf(items[p->headerIndex[6]],"%f",&M[N].score);
	if(p->headerIndex[7]>=0) sscanf(items[p->headerIndex[7]],"%f",&M[N].FDR);
	if(p->headerIndex[8]>=0) sscanf(items[p->headerIndex[8]],"%f",&M[N].repeat);
	memmove(M[N].name,M[N].name+1,strlen(M[N].name)*sizeof(char));
	sd = sqrt((double)M[N].freq);
	pos=0;
	while(fgets(buffer,3499,fp)){
		float sum=0, max=0, x;
		if(strlen(buffer)<3){ break; }
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
			sum += m[i];
		}
		for(i=0; i<4; ++i){
			x = log(4*m[i]/sum)/log10;
			M[N].m[pos][i] = x;
			if(max < x) max = x;
		}
		sumMax += max;
		++pos;
	}
	M[N].len = pos;
	M[N].threshMult=1;
	M[N].scoreThresh = SCORE;
	if(M[N].scoreThresh > sumMax*0.75) M[N].scoreThresh = sumMax*0.75;
	if(M[N].scoreThresh < 1.301) M[N].scoreThresh = 1.301;
	N++;
}
fclose(fp);
ch = truncate_filename(filename);
printf("File %s loaded. Input patterns: %d\n",ch,N);
free(buffer);
return(NNN);
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

/***********************************************/
long read_file (char *filename, int rep, int **sequence, long *seqlength, long *nseq, char **seqNames)
/***********************************************/
{
char *seq, *buffer,*seqname,*seqname1, *ch;
long seqlen=0, len, n1=0, parts=0;
FILE *fp;
long totalLen=0;

check(seq = (char*)malloc(MAXSEQ*sizeof(char)));
check(buffer = (char*)malloc(3300*sizeof(char)));
check(seqname = (char*)malloc(1300*sizeof(char)));
check(seqname1 = (char*)malloc(1300*sizeof(char)));
seq[0] = '\0';
fp = fopen(filename,"r");
if(!fp){ printf("Input file %s not found",filename); exit(0); }
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
			if(n1>=MAXSEQ-1) error_message("Too many sequences in input file");
		}
		strcpy(seqname,&buffer[1]);
		seq[0] = '\0';
		seqlen = 0;
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
			if(n1>=MAXSEQ-1) error_message("Too many sequences in input file");
			seq[0] = '\0';
			seqlen = 0;
		}
		seqlen += len;
		strcat(seq, buffer);
	}
}
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
if(n1>=MAXSEQ-1) error_message("Too many sequences in input file");
fclose(fp);
*nseq = n1;
ch = truncate_filename(filename);
printf("File %s loaded\n",ch);
free(seq);
free(buffer);
free(seqname);
free(seqname1);
return(totalLen);
}

/***********************************************/
void  matrix_offset (MATRIX *Matrix, int nPattern, PARAM *p)
/***********************************************/
{
float maxEntropy = 2, pp[4], entropy, *info, infoCum, sum, sumMean, sumVar, sumMean1, sumVar1, sumMax;
int i, j, k;
MATRIX *M;

double log2 = log(2.0);
double log10 = log(10.0);
check(info = (float*)malloc(MAXLENGTH*sizeof(float)));
for(i=0; i<nPattern; ++i){
	float maxInfo=0, var=0;
	int imax=0;
	M = &Matrix[i];
	infoCum = 0;
	for(j=0; j<M->len; ++j){
		sum = 0;
		for(k=0; k<4; k++){
			pp[k] = exp(M->m[j][k]*log10);
			sum += pp[k];
		}
		entropy = 0;
		for(k=0; k<4; ++k){
			pp[k] /= sum;
			entropy -= pp[k]*log(pp[k])/log2;
		}
		info[j] = maxEntropy - entropy;
		infoCum += info[j];
		if(j>=8) infoCum -= info[j-8];
		//printf("%d %.4f %.4f\n",j,infoCum,info[j]);
		if(maxInfo < infoCum){
			maxInfo = infoCum;
			imax = j;
		}
	}
	M->offset = imax-7;
	if(M->offset<0) M->offset=0;
	sumMean=0;
	sumVar=0;
	sumMean1=0;
	sumVar1=0;
	sumMax=0;
	for(j=0; j<M->len; ++j){
		float sx=0, ssx=0, x, max=-1;
		for(k=0; k<4; k++){
			x = M->m[j][k];
			sx += x;
			ssx += x*x;
			if(max < x) max = x;
		}
		sumMax += max;
		ssx = (ssx - sx*sx/4)/3;
		sx /= 4;
		if(j>=M->offset && j<=imax){
			sumMean1 += sx;
			sumVar1 += ssx;
		}
		sumMean += sx;
		sumVar += ssx;
	}
	float logProp = log((float)(2+10.0/p->falsePositives)*p->falsePositives/10000);
	float z = 0.968 -0.3235*logProp -0.0037336*logProp*logProp;
	M->scoreThresh = sumMean + z*sqrt(sumVar);
	if(M->scoreThresh > sumMax*0.75) M->scoreThresh = sumMax*0.75;
	if(M->scoreThresh < 1.301) M->scoreThresh = 1.301;
	M->scoreThreshCore = M->scoreThresh - sumMean1 - 3*sqrt(sumVar1);
	//printf("%.4f %.4f %.4f %.4f %.4f\n",M->scoreThresh,sumMean,sumVar,M->scoreThreshCore,z);
}
free(info);
return;
}

/***********************************************/
void  make_lookup_table (LOOKUP *lookupTable, MATRIX *Matrix, int nPattern, long *complement, int strand_option)
/***********************************************/
{
long ind, ind1;
int i, j, pos[8], k, im, *seeds;
MATRIX *M;
LOOKUP *L,*L1;
float score, missing;

check(seeds = (int*)calloc(nPattern,sizeof(int)));
for(ind=0; ind<MAXINT; ++ind){
	L = &lookupTable[ind];
	L->num=0;
	if(!L->alloc){
		L->alloc = 5;
		check(L->score = (float*)malloc(L->alloc*sizeof(float)));
		check(L->imatr = (int*)malloc(L->alloc*sizeof(int)));
		check(L->dir = (int*)malloc(L->alloc*sizeof(int)));
	}
}
for(ind=0; ind<MAXINT; ++ind){
	ind1 = ind;
	for(i=0; i<8; ++i){
		k = ind1%4;
		pos[i] = k;
		ind1 = (ind1-k)/4;
	}
	L = &lookupTable[ind];
	ind1 = complement[ind];
	L1 = &lookupTable[ind1];
	for(im=0; im<nPattern; ++im){
		M = &Matrix[im];
		score = 0;
		for(i=0; i<8; ++i){
			j = i + M->offset;
			if(j<M->len) score += M->m[j][pos[i]];
		}
		if(score >= M->scoreThreshCore){
			//printf("A %d %.4f\n",ind,score);
			if(L->num==L->alloc){
				L->alloc += 20;
				check(L->score = (float*)realloc(L->score,L->alloc*sizeof(float)));
				check(L->imatr = (int*)realloc(L->imatr,L->alloc*sizeof(int)));
				check(L->dir = (int*)realloc(L->dir,L->alloc*sizeof(int)));
			}
			L->imatr[L->num] = im;
			L->score[L->num] = score;
			L->dir[L->num] = 1;
			L->num++;
			seeds[im]++;
			if(strand_option) continue;
			if(ind1==ind) continue;
			if(L1->num==L1->alloc){
				L1->alloc += 20;
				check(L1->score = (float*)realloc(L1->score,L1->alloc*sizeof(float)));
				check(L1->imatr = (int*)realloc(L1->imatr,L1->alloc*sizeof(int)));
				check(L1->dir = (int*)realloc(L1->dir,L1->alloc*sizeof(int)));
			}
			L1->imatr[L1->num] = im;
			L1->score[L1->num] = score;
			L1->dir[L1->num] = -1;
			L1->num++;
			seeds[im]++;
		}
	}
}
//for(im=0; im<nPattern; ++im) printf("%d ",seeds[im]);
//printf("\n");
//Remove redundancy in lookup table
for(ind=0; ind<MAXINT; ++ind){
	L = &lookupTable[ind];
	for(i=0; i<L->num; ++i){
		im = L->imatr[i];
		score = L->score[i];
		for(j=i+1; j<L->num; ++j){
			int im1 = L->imatr[j];
			if(im1 != im) continue;
			float score1 = L->score[j];
			if(score1 > score){
				L->score[i] = score1;
				score = score1;
				L->dir[i] = L->dir[j];
			}
			for(k=j; k<L->num-1; ++k){
				L->score[k] = L->score[k+1];
				L->imatr[k] = L->imatr[k+1];
				L->dir[k] = L->dir[k+1];
			}
			--j;
			--L->num;
		}
	}
}
free(seeds);
return;
}

/***********************************************/
int  remove_position (POSITION *P, int iseq, int pos, int len)
/***********************************************/
{
int i, j, k;

i=P->num-1;
while(i>=0 && P->pos[i] > pos && P->iseq[i]==iseq) --i; 
if(P->pos[i]==pos && P->iseq[i]==iseq){
	for(j=i; j<P->num-1; ++j){
		P->dir[j] = P->dir[j+1];
		P->score[j] = P->score[j+1];
		P->iseq[j] = P->iseq[j+1];
		P->pos[j] = P->pos[j+1];
		for(k=0; k<len; ++k){
			P->pattern[j*len+k] = P->pattern[(j+1)*len+k];
		}
	}
	--P->num;
}else{
	return(0);		
}
return(1);
}

/***********************************************/
int  find_motifs (int **sequence, long *seqlen, long nseq, LOOKUP *lookupTable,
		MATRIX *Matrix, int nPattern, POSITION *position, int adjustThresh, PARAM *p)
/***********************************************/
{
long iseq, pos, ind, lastMissing=0, totalLength=0;
int i, j, k, mult, *seq, dir, *pattern, N, x, threshChanged=0, im, done=0;
int nstack=0;
float score, score1, *stack;
MATRIX *M;
LOOKUP *L;
POSITION *P;

if(!adjustThresh) done=1;
check(stack = (float*)malloc(NSTACK*sizeof(float)));
check(pattern = (int*)malloc(MAXLENGTH*sizeof(int)));
for(im=0; im<nPattern; ++im){
	position[im].num = 0;
}
for(iseq=0; iseq<nseq; ++iseq){
	seq = sequence[iseq];
	ind = 0;
	mult = 1;
	lastMissing = -1;
	nstack=0;
	totalLength += seqlen[iseq]-8;
	for(pos=0; pos<seqlen[iseq]-8; ++pos){
		if(pos==0){
			for(j=0; j<8; ++j){
				x = seq[j];
				if(x>0) ind += x*mult;
				else if(x<0) lastMissing = j;
				mult *= 4;
			}
		}else{
			x = seq[pos+7];
			if(x>0) ind += x*mult;
			else if(x<0) lastMissing = pos+7;
			ind /= 4;
		}
		L = &lookupTable[ind];
		for(i=0; i<L->num; ++i){
			int start, pos1, im1, len1, missing;
			im = L->imatr[i];
			M = &Matrix[im];
			score = L->score[i];
			dir = L->dir[i];
			//printf("%d %d %d %d %.4f\n",iseq,pos,ind,dir,score);
			if(dir>0){
				start = pos - M->offset;
				if(start<=lastMissing || start + M->len >= seqlen[iseq]) continue;
			}else{
				start = pos+8+M->offset-1;
				if(start - (M->len-1)<=lastMissing || start >= seqlen[iseq]){ continue; }
			}
			missing = 0;
			for(j=0; j<M->len; ++j){
				if(dir>0){
					pattern[j] = seq[start+j];
					if(pattern[j] <0) missing=1;
				}else{
					pattern[j] = 3 - seq[start-j];
					if(pattern[j] >3) missing=1;
				}
				if(missing){ break; }
				if(j>=M->offset && j<M->offset+8) continue;
				score += M->m[j][pattern[j]];
			}
			if(score < M->scoreThresh || missing) continue;
			if(dir<0) start = start - M->len + 1;

			/* Check for overlap with other TFBS */
			int len = M->len;
			int redundant = 0;
			for(j=nstack-1; j>=0; j-=3){
				im1 = stack[j-2];
				pos1 = stack[j-1];
				score1 = stack[j];
				len1 = Matrix[im1].len;
				//printf("%d %.4f %d %.4f\n",pos1,score1,start,score);		
				if(!p->presence_option && pos-pos1 > MINDIST) break;
				if(im1!=im && start-pos1 >= len1*0.75 && pos1+len1-start < len*0.75) continue;
				if(score1 > score){
					redundant=1;
					break;
				}else{
					if(!remove_position(&position[im1],iseq,pos1,M->len)) error_message("Deletion Failed");
					memmove(&stack[j-2], &stack[j+1], (NSTACK-j-1)*sizeof(float));
					nstack -= 3;
				}				
			}
			if(redundant) continue;
			if(nstack >= NSTACK-3){
				memmove(stack, &stack[3], (NSTACK-3)*sizeof(float));
				nstack -= 3;
			}
			stack[nstack++] = im;
			stack[nstack++] = start;
			stack[nstack++] = score;

			P = &position[im];
			N = P->num;
			if(!P->alloc){
				P->alloc = 100;
				check(P->score = (float*)malloc(P->alloc*sizeof(float)));
				check(P->iseq = (long*)malloc(P->alloc*sizeof(long)));
				check(P->pos = (long*)malloc(P->alloc*sizeof(long)));
				check(P->dir = (int*)malloc(P->alloc*sizeof(int)));
				check(P->pattern = (int*)malloc(P->alloc*M->len*sizeof(int)));
			}else if(N == P->alloc){
				P->alloc += 100;
				check(P->score = (float*)realloc(P->score,P->alloc*sizeof(float)));
				check(P->iseq = (long*)realloc(P->iseq,P->alloc*sizeof(long)));
				check(P->pos = (long*)realloc(P->pos,P->alloc*sizeof(long)));
				check(P->dir = (int*)realloc(P->dir,P->alloc*sizeof(int)));
				check(P->pattern = (int*)realloc(P->pattern,P->alloc*M->len*sizeof(int)));
			}
			P->dir[N] = dir;
			P->score[N] = score;
			P->iseq[N] = iseq;
			P->pos[N] = start;
			for(j=0; j<M->len; ++j){
				P->pattern[N*M->len+j] = pattern[j];
			}
			P->num++;
		}
	}
	long testLength = 200*10000/p->falsePositives;
	if(!done && (totalLength > testLength || iseq==nseq-1)){
		for(im=0; im<nPattern; ++im){
			float logratio, expected;
			M = &Matrix[im];
			if(M->threshDone) continue;
			N = position[im].num;
			expected = (1+10.0/p->falsePositives)*totalLength*p->falsePositives/10000;
			if(N >= expected || M->scoreThresh < 1.301){
				float *scores;
				float x1 = M->scoreThresh;
				if(M->scoreThresh >= 1.301){
					check(scores = (float*)malloc(MAXINT*sizeof(float)));
					for(i=0;i<N;++i) scores[i] = position[im].score[i];
					sortem(N, scores, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
					M->scoreThresh = scores[N-(int)expected];
					free(scores);
				}
				if(M->scoreThresh < 1.301) M->scoreThresh = 1.301;
				M->scoreThreshCore += M->scoreThresh - x1;
				M->threshDone = 1;
				//printf("Done: %d %.4f\n",im,M->scoreThresh);
				continue;
			}
			expected = (2+10.0/p->falsePositives)*totalLength*p->falsePositives/10000;
			if(N<10) N=10;
			logratio = log((float)N/expected);
			//printf("%d - %d %.4f %.4f %d %.1f\n",im,N,Matrix[im].scoreThresh,logratio,totalLength,expected);
			M->scoreThresh += logratio*M->threshMult;
			M->scoreThreshCore += logratio*M->threshMult;
			M->threshMult *= 2;
			threshChanged++;
			//printf("%d - %d %.4f %.4f %d\n",im,N,M->scoreThresh,logratio,totalLength);
		}
		done = 1;
	}
	if(threshChanged) break;
}
free(pattern);
free(stack);
return(threshChanged);
}

/***********************************************/
void  adjust_scoreThresh_old(MATRIX *Matrix, int nPattern, POSITION *position1, POSITION *position2,
		float adjustLength)
/***********************************************/
{
long imatr, ipos, i, j, N1, N2, i1, i2, n1, n2;
float *score1, *score2, thresh, thresh1, maxDiff;

for(imatr=0; imatr<nPattern; ++imatr){
	N1 = position1[imatr].num;
	check(score1 = (float*)malloc(N1*sizeof(float)));
	memcpy(score1,position1[imatr].score,N1*sizeof(int));
	sortem(N1, score1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	N2 = position2[imatr].num;
	check(score2 = (float*)malloc(N2*sizeof(float)));
	memcpy(score2,position2[imatr].score,N2*sizeof(int));
	sortem(N2, score2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	i1=0;
	i2=0;
	thresh1 = 0;
	maxDiff = -1;
	while(1){
		float diff;
		thresh = score1[N1-1-i1];
		while(i1<N1 && score1[N1-1-i1] >= thresh) ++i1;
		while(i2<N2 && score2[N2-1-i2] >= thresh) ++i2;
		diff = i1-i2*adjustLength*1.5;
		//printf("%d\t%d\t%.4f\t%.1f\n",i1,i2,thresh,diff);
		if(maxDiff<diff && i1>30 || thresh > 8){
			maxDiff=diff;
			thresh1 = thresh;
			n1 = i1;
			n2 = i2*adjustLength;
		}
		if(i1>=N1-1) break;
		thresh = score1[N1-1-i1];
	}
	Matrix[imatr].scoreThresh = thresh1;
	Matrix[imatr].freq = n1-n2;
	Matrix[imatr].ratio = (float)n1/n2;
	Matrix[imatr].score = Matrix[imatr].freq/sqrt(n2);
	free(score1);
	free(score2);
}
return;
}

/***********************************************/
void  adjust_scoreThresh(MATRIX *Matrix, int nPattern, POSITION *position1, POSITION *position2,
	 long length1, long length2, PARAM *p)
/***********************************************/
{
long imatr, ipos, i, j, N1, N2, i1, i2, n1, n2;
float *score1, *score2, thresh, thresh1, maxDiff;

for(imatr=0; imatr<nPattern; ++imatr){
	N1 = position1[imatr].num;
	check(score1 = (float*)malloc(N1*sizeof(float)));
	memcpy(score1,position1[imatr].score,N1*sizeof(int));
	sortem(N1, score1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	N2 = position2[imatr].num;
	check(score2 = (float*)malloc(N2*sizeof(float)));
	memcpy(score2,position2[imatr].score,N2*sizeof(int));
	sortem(N2, score2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	i2 = length2*p->falsePositives/10000;
	if(i2>N2) i2 = N2-1;
	thresh = score2[N2-1-i2];
	i1=0;
	while(i1<N1 && score1[N1-1-i1] >= thresh) ++i1;
	n1 = i1;
	n2 = i2*length1/length2;
	Matrix[imatr].scoreThresh = thresh;
	Matrix[imatr].freq = n1-n2;
	Matrix[imatr].ratio = (float)n1/n2;
	Matrix[imatr].score = Matrix[imatr].freq/sqrt(n2);
	free(score1);
	free(score2);
}
return;
}

/***********************************************/
void  get_frequency(MATRIX *Matrix, int nPattern, POSITION *position, int **freqMatch)
/***********************************************/
{
int imatr, i, j, N, len;

for(imatr=0; imatr<nPattern; ++imatr){
	len = Matrix[imatr].len;
	for(i=0; i<len*4; ++i){
		freqMatch[imatr][i]=0;
	}
	N = position[imatr].num;
	for(i=0; i<N; ++i){
		if(position[imatr].score[i] < Matrix[imatr].scoreThresh) continue;
		for(j=0; j<len; ++j){
			int x = position[imatr].pattern[i*len+j];
			freqMatch[imatr][j*4+x]++;
		}
	}
}
return;
}

/***********************************************/
void  get_reverse_sequence (int **sequence1, long *seqlen1, long n1, int **sequence2)
/***********************************************/
{
long i, len, j;
int x;
for(i=0; i<n1; ++i){
	len = seqlen1[i];
	check(sequence2[i] = (int*)malloc(len*sizeof(int)));
	for(j=0; j<len; ++j){
		x = sequence1[i][len-j-1];
		if(x>=0) x = 3-x;
		sequence2[i][j] = x;
	}
}
return;
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
void  update_scores(MATRIX *Matrix, int nPattern, POSITION *position)
/***********************************************/
{
int imatr, i, j, N, len;
float score;
MATRIX *M;

for(imatr=0; imatr<nPattern; ++imatr){
	M = &Matrix[imatr];
	len = M->len;
	N = position[imatr].num;
	for(i=0; i<N; ++i){
		score = 0;
		for(j=0; j<len; ++j){
			int x = position[imatr].pattern[i*len+j];
			score += M->m[j][x];
		}
//printf("%.4f\t%.4f\n",position[imatr].score[i],score);
		position[imatr].score[i] = score;
	}
}
return;
}

/***********************************************/
void  print_progress (int iter, int iter1, int **freqMatch1, int **freqMatch2, POSITION *position1, POSITION *position2,
		MATRIX *Matrix, int nPattern, float adjust, FILE *fp)
/***********************************************/
{
int imatr, i, j, n1, n2, N1, N2, *fr1, *fr2;

if(!iter && !iter1){
	fprintf(fp,"Name\tIteration\tsubIteration\tN-TestGrThresh\tN-ContGrThresh\tN-TestAll\tN-ContAll\n");
}
for(imatr=0; imatr<nPattern; ++imatr){
	fr1 = freqMatch1[imatr];
	fr2 = freqMatch2[imatr];
	n1 = fr1[0]+fr1[1]+fr1[2]+fr1[3];
	n2 = (fr2[0]+fr2[1]+fr2[2]+fr2[3])*adjust;
	N1 = position1[imatr].num;
	N2 = position2[imatr].num*adjust;
	fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n",Matrix[imatr].name,iter+1,iter1+1,n1,n2,N1,N2,Matrix[imatr].scoreThresh);
	//printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n",Matrix[imatr].name,iter+1,iter1+1,n1,n2,N1,N2,Matrix[imatr].scoreThresh);
}
return;
}

/***********************************************/
void  adjust_matrix(MATRIX *Matrix, int nPattern, int **freqMatch1, int **freqMatch2, float adjust, int method)
/***********************************************/
{
int imatr, i, j, k, N1, N2, len, *fr1, *fr2;
float x, y, *test, *cont, a, b, pp[4], sum, residual;
MATRIX *M;

double log10 = log(10.0);
for(imatr=0; imatr<nPattern; ++imatr){
	len = Matrix[imatr].len;
	check(test = (float*)malloc(len*4*sizeof(float)));
	check(cont = (float*)malloc(len*4*sizeof(float)));
	fr1 = freqMatch1[imatr];
	fr2 = freqMatch2[imatr];
	N1 = fr1[0]+fr1[1]+fr1[2]+fr1[3];
	N2 = fr2[0]+fr2[1]+fr2[2]+fr2[3];
	M = &Matrix[imatr];
	if(method==METHOD_REGRESS){
		for(i=0; i<len; ++i){
			for(j=0; j<4; ++j){
				k = i*4+j;
				test[k] = log((float)fr1[k]/N1 + 0.01)/log10;
				cont[k] = log((float)fr2[k]/N2 + 0.01)/log10;
				//printf("%.4f\t%.4f\n",test[k],cont[k]);
			}
		}
		regression(test,cont,len*4,&a,&b);
	}
	for(i=0; i<len; ++i){
		sum = 0;
		for(j=0; j<4; ++j){
			k = i*4+j;
			if(method==METHOD_REGRESS){
				residual = test[k] - (a+b*cont[k]);
				M->m[i][j] += 0.5*residual;
				pp[j] = exp(M->m[i][j]*log10);
			}else if(method==METHOD_RESAMPLE){
				pp[j] = fr1[k];
				if(pp[j]<1) pp[j]=1;
			}else if(method==METHOD_DIFF){
				pp[j] = fr1[k] - adjust*fr2[k];
				if(pp[j]<1) pp[j]=1;
			}else{
				error_message("Method unknown");
			}
			sum += pp[j];
		}
		for(j=0; j<4; ++j){
			pp[j] /= sum;
			M->m[i][j] = log(pp[j]*4)/log10;
		}
	}
	free(test);
	free(cont);
}
return;
}

/***********************************************/
PARAM *read_parameters (int nargs, char **argv)
/***********************************************/
{
PARAM *p;
int iarg=1, score=-1;

if(nargs < 4){ printf("patternTest -i motifFile -f fastaFile -o outputFile [-c controlFile, -prog progressFile, -maxlen maxSequenceLength, -strand strandOption, -n numberOfMotifs, -method resampleMethod, -iter numIteartions, -siter numSubiterations, -fp falsePositives, -userep, -one]"); exit(0); }
check(p = (PARAM*)calloc(1,sizeof(PARAM)));
p->max_length = 25000000;
p->N_matrix = 100;
p->iterations = 3;
p->subIterations = 3;
p->method = 1;
p->falsePositives = 5;
while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1) p->inputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-f") && iarg < nargs-1) p->testFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-c") && iarg < nargs-1) p->controlFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-prog") && iarg < nargs-1) p->progressFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-maxlen") && iarg < nargs-1) sscanf(argv[++iarg],"%ld",&p->max_length);
	else if(!strcmp(argv[iarg],"-userep")) p->use_repeats = 1;
	else if(!strcmp(argv[iarg],"-strand") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->strand);
	else if(!strcmp(argv[iarg],"-method") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->method);
	else if(!strcmp(argv[iarg],"-fp") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->falsePositives);
	else if(!strcmp(argv[iarg],"-n") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->N_matrix);
	else if(!strcmp(argv[iarg],"-one")) p->presence_option = 1;
	else if(!strcmp(argv[iarg],"-iter") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->iterations);
	else if(!strcmp(argv[iarg],"-siter") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->subIterations);
	else{
		printf("Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}
if(!p->inputFile[0] || !p->outputFile[0]){ printf("ERROR: input or output file not specified\n"); exit(0); }
if(p->falsePositives < 0.01) error_message("False positives should be >0.01");
return(p);
}

/***********************************************/
int main (int argc, char **argv) 
/***********************************************/
{
MATRIX *Matrix=NULL;
PARAM *p;
FILE *fp, *fp1;
LOOKUP *lookupTable;
POSITION *position1, *position2;
long *seqlen1, *seqlen2, *complement, ind, ind1, ind2, n1, n2, iter, iter1;
int len, len1, i, j, im, nPattern, **freqMatch1, **freqMatch2, **sequence1, **sequence2;
long length1=0, length2=0;
float adjust, x, y, log2;
char **seqNames;

double log10 = log(10.0);
p = read_parameters(argc, argv);
nPattern = read_input (p->inputFile,&Matrix,p);
if(nPattern<=0) error_message("Input motif file is empty");
check(seqlen1 = (long*)malloc(MAXSEQ*sizeof(long)));
check(seqlen2 = (long*)malloc(MAXSEQ*sizeof(long)));
check(sequence1 = (int**)malloc(MAXSEQ*sizeof(int*)));
check(sequence2 = (int**)malloc(MAXSEQ*sizeof(int*)));
check(lookupTable = (LOOKUP*)calloc(MAXINT,sizeof(LOOKUP)));
check(complement = (long*)malloc(MAXINT*sizeof(long)));
check(position1 = (POSITION*)calloc(nPattern,sizeof(POSITION)));
check(position2 = (POSITION*)calloc(nPattern,sizeof(POSITION)));
check(freqMatch1 = (int**)malloc(nPattern*sizeof(int*)));
check(freqMatch2 = (int**)malloc(nPattern*sizeof(int*)));
check(seqNames = (char**)malloc(MAXSEQ*sizeof(char*)));
for(i=0; i<nPattern; ++i){
	check(freqMatch1[i] = (int*)malloc(Matrix[i].len*4*sizeof(int)));
	check(freqMatch2[i] = (int*)malloc(Matrix[i].len*4*sizeof(int)));
}
length1 = read_file(p->testFile,p->use_repeats,sequence1,seqlen1,&n1,seqNames);
if(length1 < 10){ error_message("Test sequence is too short (<10)"); }
if(length1 < 5000){ printf("WARNING: Test sequence is too short (<5000)"); }
if(length1 > p->max_length){ printf("WARNING: Sequence is too long (>%d); truncated\n",p->max_length); }
if(p->presence_option){
	if(n1<1){ error_message("Presence option requires > one test sequence"); }
	if(n1<10){ printf("WARNING: Too few test sequences (<10)"); }
}
printf("Number of test sequences = %d. Total length = %d\n",n1,length1);
if(p->strand>=2){
	get_reverse_sequence(sequence1,seqlen1,n1,sequence2);
	n2 = n1;
	length2 = length1;
	memcpy(seqlen2,seqlen1,n1*sizeof(long));
}else if(p->controlFile){
	length2 = read_file(p->controlFile,p->use_repeats,sequence2,seqlen2,&n2,NULL);
	if(p->presence_option){ // Adjust the length of control sequences
		int averTest = length1/n1;
		length2 = 0;
		for(i=0; i<n2; ++i){
			if(seqlen2[i] > averTest) seqlen2[i] = averTest;
			length2 += seqlen2[i];
		}
	}
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
for(i=0; i<MAXINT; ++i){
	complement[i] = complementary(i);
}
log2 = log(2.0);

fp = fopen(p->outputFile,"w");
if(!fp) error_message("Output file not opened");
if(p->progressFile){
	fp1 = fopen(p->progressFile,"w");
	if(!fp1) error_message("Progress file not opened");
}
adjust = (float)length1/length2;
for(iter=0; iter<=p->iterations; ++iter){
	int threshChanged = 1;
	//printf("Iter %d\n",iter+1);
	if(iter == p->iterations) p->use_repeats=1;
	matrix_offset(Matrix,nPattern,p);
	iter1=0;
	while(threshChanged && iter1<3){
		make_lookup_table(lookupTable,Matrix,nPattern,complement,p->strand);
		threshChanged = find_motifs(sequence2,seqlen2,n2,lookupTable,Matrix,nPattern,position2,1,p);
		++iter1;
	}
	if(threshChanged){
		find_motifs(sequence2,seqlen2,n2,lookupTable,Matrix,nPattern,position2,0,p);
	}
	find_motifs(sequence1,seqlen1,n1,lookupTable,Matrix,nPattern,position1,0,p);
	for(iter1=0; iter1<p->subIterations; ++iter1){
		//printf("Iter %d - %d %d %d\n",iter+1,iter1+1,position1[0].num,position2[0].num);
		if(iter1>0){
			update_scores(Matrix,nPattern,position1);
			update_scores(Matrix,nPattern,position2);
		}
		//adjust_scoreThresh_old(Matrix,nPattern,position1,position2,adjust);
		adjust_scoreThresh(Matrix,nPattern,position1,position2,length1,length2,p);
		get_frequency(Matrix,nPattern,position1,freqMatch1);
		get_frequency(Matrix,nPattern,position2,freqMatch2);
		if(p->progressFile){
			print_progress(iter,iter1,freqMatch1,freqMatch2,position1,position2,Matrix,nPattern,adjust,fp1);
		}
		if(iter == p->iterations) break;
		adjust_matrix(Matrix,nPattern,freqMatch1,freqMatch2,adjust,p->method);
	}
}
fprintf(fp,"Headers:\tName\tPattern\tThreshold");
if(p->headerIndex[3]>=0) fprintf(fp,"\tFreq");
if(p->headerIndex[4]>=0) fprintf(fp,"\tRatio");
fprintf(fp,"\tInfo\tScore");
if(p->headerIndex[7]>=0) fprintf(fp,"\tFDR");
if(p->headerIndex[8]>=0) fprintf(fp,"\tRepeat");
fprintf(fp,"\n");
for(im=0; im<nPattern; ++im){
	float sum, max, x1, entropy;
	int len, x, y, coef;
	MATRIX *M = &Matrix[im];
	len = M->len;
	M->info=0;
	M->pattern[0]='\0';
	M->patternRev[0]='\0';
	for(i=0; i<len; ++i){
		max=-1;
		sum=0;
		for(j=0; j<4; ++j){
			x1 = M->m[i][j];
			if(x1>max){ max=x1; }
			x1 = exp(x1*log10);
			sum += x1;
		}
		x=0;
		coef=1;
		entropy=0;
		for(j=0; j<4; ++j){
			y = 1;
			if(M->m[i][j] < max-0.301){ y=0; }
			x += y*coef;
			coef *= 2;
			x1 = exp(M->m[i][j]*log10)/sum;
			entropy -= x1*log(x1)/log2;
		}
		strcat(M->pattern,codes[x]);
		M->info += 2.0-entropy;
	}
	fprintf(fp,">%s\t%s\t%.4f",M->name,M->pattern,M->scoreThresh);
	if(p->headerIndex[3]>=0) fprintf(fp,"\t%d",M->freq);
	if(p->headerIndex[4]>=0) fprintf(fp,"\t%.4f",M->ratio);
	fprintf(fp,"\t%.4f\t%.4f",M->info,M->score);
	if(p->headerIndex[7]>=0) fprintf(fp,"\t%.4f",M->FDR);
	if(p->headerIndex[8]>=0) fprintf(fp,"\t%.2f",M->repeat);
	fprintf(fp,"\n");
	print_matrix(M, fp);
}
fclose(fp);
printf("patternTest done\n");
return(0);
}

