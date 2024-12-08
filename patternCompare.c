/*************************************
* patternCompare is a part of CisFinder software
* http://lgsun.grc.nia.nih.gov/CisFinder
*
* Function: compare motifs between 2 files based on PFM similarity.
* Similarity is measured by correlation of position-weight matrices (PWM) which
* are log-transformed PFMs.
* 
* Syntax:
* patternCompare -i1 motifFile1.txt -i2 motifFile2 -o outputFile [-match matchThreshold]
* 
* Comments:
* (a) matchThreshold = minimum similarity (correlation of PWMs) between motifs;
* default value = 0.75.
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
#define NHEADERS 9

int DEBUG=0;

typedef struct match_st{
	int ind;
	int dir;
	int off;
	float match;
}MATCH;
typedef struct matrix_st{
	int id;
	int len;
	int nMatch;
	MATCH *match;
	float m[MAXLENGTH][4];
	char pattern[MAXLENGTH];
	char patternRev[MAXLENGTH];
	char *name;
	char *info;
}MATRIX;
typedef struct param_st{
	char *inputFile1;
	char *inputFile2;
	char *outputFile;
	float matchThresh;
	int headerIndex[NHEADERS];
}PARAM;

static char *headerItems[NHEADERS]={"Name","Pattern","Threshold","Freq","Ratio","Info","Score","FDR","Repeat"};
static char *codes[16]={"O","A","C","M","G","R","S","V","T","W","Y","H","K","D","B","N"};
static char *revcodes[16]={"O","T","G","K","C","Y","S","B","A","W","R","D","M","H","V","N"};
int *stack, nStack=0;

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
		printf("ERROR: Out of memory\n");
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
for(i=len-1; i>=0; --i){
	for(j=3; j>=0; --j){
		M1->m[len-i-1][3-j] = M->m[i][j];
	}
}
return(M1);
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
void   print_matrix  (MATRIX *M, FILE *fp)
/*************************************/
{
int i,j;
float pp[4], sum;
for(i=0; i<M->len; ++i){
	fprintf(fp,"%d",i);
	sum = 0;
	for(j=0; j<4; ++j){
		pp[j] = exp(M->m[i][j]);
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
void  reverse_pattern (char *pattern, char *patternRev)
/***********************************************/
{
int i,j,len;

len = strlen(pattern);
for(i=0;i<len;++i){
	for(j=0;j<16;++j){
		if(pattern[i]==codes[j][0]) patternRev[len-1-i]=revcodes[j][0];
	}
}
patternRev[len]='\0';
return;
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
char *buffer, *ch, *shortname, *items[10];
FILE *fp;
MATRIX *M;
int N=0, i,pos,NNN;
float m[4];

check(buffer = (char*)malloc(3500*sizeof(char)));
fp = fopen(filename,"r");
shortname = truncate_filename(filename);
if(!fp){ printf("ERROR: Input file %s not found",shortname); exit(0); }
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
if(!N) error_message("No input motifs found!");
NNN = N;
check(*Mp = (MATRIX*)calloc(NNN,sizeof(MATRIX)));
M = *Mp;
N = 0;
while(fgets(buffer,3499,fp)){
	int done=0;
	float sd;
	while(buffer[0] != '>' && !done){
		if(!fgets(buffer,3499,fp)) done=1;
	}
	if(N>=NNN || done){ break; }
	ch = strchr(buffer,'\t');
	M[N].info = copy_string(ch);
	int n = split_string(&buffer[1],items,10);
	if(p->headerIndex[0]>=0) M[N].name = copy_string(items[p->headerIndex[0]]);
	M[N].name = copy_string(items[0]);
	pos=0;
	while(fgets(buffer,3499,fp)){
		float sum=0, max=0, sum1=0;
		if(strlen(buffer)<3){ break; }
		sscanf(buffer,"%d%f%f%f%f",&i,&m[0],&m[1],&m[2],&m[3]);
		if(i!=pos){ printf("ERROR: Wrong file format in line: %s",buffer); exit(0); }
		for(i=0; i<4; ++i){
			sum += m[i];
		}
		float pseudocount=sum*0.01;
		if(!pseudocount || pseudocount>1) pseudocount=1;
		sum = 0;
		for(i=0; i<4; ++i){
			if(m[i]<0) error_message("Negative value in PFM");
			m[i]+=pseudocount;
			if(max<m[i]) max=m[i];
			M[N].m[pos][i] = log(m[i]); 
			sum += M[N].m[pos][i];
		}
		int x=0;
		int coef=1;
		for(i=0; i<4; ++i){
			M[N].m[pos][i] -= sum/4;
			int y = 1;
			if(m[i] < max/2){ y=0; }
			x += y*coef;
			coef *= 2;
		}
		strcat(M[N].pattern,codes[x]);
		strcat(M[N].patternRev,revcodes[x]);
		++pos;
	}
	reverse_string(M[N].patternRev);
	M[N].len = pos;
	N++;
}
fclose(fp);
printf("File %s loaded. Input patterns: %d\n",shortname,N);
free(buffer);
return(NNN);
}

/******************************************/
void  compare_patterns (MATRIX *Matrix1, int nPattern1, MATRIX *Matrix2, int nPattern2, PARAM *p, FILE *fp)
/******************************************/
{
int i, ind1, ind2;
int offset, dir, nMatch;
float match, *score, *index;
MATCH *match_list;

check(match_list = (MATCH*)calloc(nPattern2,sizeof(MATCH)));
check(score = (float*)calloc(nPattern2,sizeof(float)));
check(index = (float*)calloc(nPattern2,sizeof(float)));
for(ind1=0; ind1<nPattern1; ++ind1){
	nMatch = 0;
	int offMin = 0;
	for(ind2=0; ind2<nPattern2; ++ind2){
		int nnn = align_matrixes_best(&Matrix1[ind1],&Matrix2[ind2],&match,&offset,&dir,p->matchThresh);
		if(match < p->matchThresh){ continue; }
		match_list[nMatch].match = match;
		match_list[nMatch].dir = dir;
		match_list[nMatch].off = offset;
		if(offMin > offset) offMin=offset;
		match_list[nMatch].ind = ind2;
		score[nMatch] = match+(1.0-match)*0.5*(1.0-6.0/(float)nnn);
		index[nMatch] = nMatch;
		++nMatch;
		//printf("%d %d %d %d %.4f %s %s\n",ind1,ind2,dir,offset,match,Matrix1[ind1].pattern,Matrix2[ind2].pattern);
	}
	fprintf(fp,">%s\t",Matrix1[ind1].name);
	for(i=0; i<-offMin; ++i) fprintf(fp," ");
	fprintf(fp,"%s\n",Matrix1[ind1].pattern);
	if(nMatch > 1)
		sortem(nMatch, score, 1, index, NULL, NULL, NULL, NULL, NULL, NULL);
	for(i=0; i<nMatch; ++i){
		int i1;
		int j = index[nMatch-1-i];
		ind2 = match_list[j].ind;
		dir = match_list[j].dir;
		fprintf(fp,"%s\t",Matrix2[ind2].name);
		for(i1=0; i1<match_list[j].off-offMin; ++i1) fprintf(fp," ");
		if(dir>0){
			fprintf(fp,"%s\t%.4f\t%d\t%d\n",Matrix2[ind2].pattern,match_list[j].match,dir,match_list[j].off);
		}else{
			fprintf(fp,"%s\t%.4f\t%d\t%d\n",Matrix2[ind2].patternRev,match_list[j].match,dir,match_list[j].off);
		}
	}
}
free(match_list);
free(score);
free(index);
return;
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
int iarg=1, score=-1;
static char *syntax = "patternCompare -i1 motifFile1.txt -i2 motifFile2 -o outputFile [-match matchThreshold]\n";

if(nargs < 3){ printf("%s\n",syntax); exit(0); }
check(p = (PARAM*)calloc(1,sizeof(PARAM)));
p->matchThresh=0.7;
while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i1") && iarg < nargs-1) p->inputFile1=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-i2") && iarg < nargs-1) p->inputFile2=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-match") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->matchThresh);
	else{
		printf("ERROR: Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}
if(!p->inputFile1[0] || !p->inputFile2[0] || !p->outputFile[0]){ printf("ERROR: %s\n",syntax); exit(0); }
if(p->matchThresh <0.5 || p->matchThresh>0.99) error_message("Match threshold should be from 0.5 to 0.99");
return(p);
}

/***********************************************/
int main (int argc, char **argv) 
/***********************************************/
{
MATRIX *Matrix, *Matrix1;
PARAM *p;
FILE *fp;
int nPattern, nPattern1, ind, ind1, ind2, nSimPairs=0, len, len1, i, j, nSeq=0;
float *sorted, *score;

p = read_parameters(argc, argv);
nPattern = read_file(p->inputFile1,&Matrix,p);
if(nPattern<=0) error_message("Input motif file #1 is empty");
nPattern1 = read_file(p->inputFile2,&Matrix1,p);
if(nPattern1<=0) error_message("Input motif file #2 is empty");
fp = fopen(p->outputFile,"w");
compare_patterns(Matrix,nPattern,Matrix1,nPattern1,p,fp);
fclose(fp);
return(0);
}
