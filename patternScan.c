/*************************************
* patternScan is a part of CisFinder software
* http://lgsun.grc.nia.nih.gov/CisFinder

* Function: finds sites in a sequence that match to motifs specified by PFM
* 
* Syntax:
* patternScan -i motifFile -f fastaFile -o outputFile [-cons conservationFile, 
* -maxlen maxSequenceLength, -strand strandOption, -n numberOfMotifs, 
* -fp falsePositives, -thresh matchThreshold, -userep, -one, -redund]
* 
* Comments:
* (a) conservation file = file with evolutionary conservation scores (see format below)
* (b) strandOption: 0 = search both strands, 1 = search foreward strand.
* (c) numberOfMotifs can be used to limit the number of input motifs (default=100)
* (d) falsePositives = number of expected false positives per 10000 bp in a random
* sequence (default = 5 if matrix-specific thresholds are not supplied).
* (e) matchThreshold = addition to matrix-specific thresholds. For example if
* it is equal to 0.7, then all match thresholds are incremented by 0.7 and
* search becomes more stringent.
* (f) userep: use repeats in sequence (lower-case in sequence)
* (g) one: consider not more than 1 motif occurrence per sequence (with highest 
* match score)
* (h) redund: include overlapping and redundant motifs. If this option is not selected
* then only the best matching motif out of over-lapping ones is recorded.
* (i) Matrix-specific match thresholds are automatically generated by the
* patternTest program, but they can be modified, added, or removed manually.
* The header for the threshold field should be "Threshold" (see file formats).
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
#define MAXSEQ 200000
#define MAXINT 65536
#define MINDIST 15
#define SCORE 2.00
#define NHEADERS 7
#define NSTACK 60

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
	float scoreThresh;
	float scoreThreshCore;
	float m[MAXLENGTH][4];
	char pattern[MAXLENGTH];
	char patternRev[MAXLENGTH];
	char name[MAXLENGTH];
	float multScore;
}MATRIX;
typedef struct position_st{
	int num;
	int alloc;
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
	char *conservFile;
	int use_repeats;
	int strand;
	int presence_option;
	int N_matrix;
	int interval;
	long max_length;
	int headerIndex[NHEADERS];
	float threshAdd;
	float falsePositives;
	int redundancy;
}PARAM;
static char *nucleotides[4]={"A","C","G","T"};
static char *headerItems[NHEADERS]={"Name","Pattern","Threshold","Freq","Ratio","Info","Score"};
static char *strands[2]={"-","+"};

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
long  get_random_sequence (int **sequence, long *seqlen, long *nseq, float falsePositives)
/***********************************************/
{
long i, len=1000, j;

*nseq = 500;
if(falsePositives < 1) *nseq /= falsePositives;
for(i=0; i<*nseq; ++i){
	seqlen[i] = len;
	check(sequence[i] = (int*)malloc(len*sizeof(int)));
	for(j=0; j<len; ++j){
		int x = rand()%4;
		sequence[i][j] = x;
	}
}
return(len*(*nseq));
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
int  split_string (char *string, char *items[], int num, char separator)
/***********************************************/
{
char *ch;
int i=0, len;

ch = strchr(string,'\n');
if(ch) *ch = '\0';
len = strlen(string);
while(len>0 && (string[len-1]==separator || string[len-1]==' ')) --len;
string[len] = '\0';
ch = string;
while(1){
	items[i] = ch;
	ch = strchr(ch,separator);
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

n = split_string(headers,items,10,'\t');
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
int N=0, i, pos, NNN, len;
float m[4];

check(buffer = (char*)malloc(3500*sizeof(char)));
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("ERROR: Input file %s not found",ch); exit(0); }
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
	float sumMax=0;
	while(buffer[0] != '>' && !done){
		if(!fgets(buffer,3499,fp)) done=1;
	}
	if(N>=NNN || done){ break; }
	int n = split_string(&buffer[1],items,10,'\t');
	if(p->headerIndex[0]>=0) strcpy(M[N].name,items[p->headerIndex[0]]);
	else strcpy(M[N].name,items[0]);
	if(p->headerIndex[1]>=0) strcpy(M[N].pattern,items[p->headerIndex[1]]);
	if(p->headerIndex[2]>=0) sscanf(items[p->headerIndex[2]],"%f",&M[N].scoreThresh);
	if(p->headerIndex[3]>=0) sscanf(items[p->headerIndex[3]],"%d",&M[N].freq);
	if(p->headerIndex[4]>=0) sscanf(items[p->headerIndex[4]],"%f",&M[N].ratio);
	if(p->headerIndex[5]>=0) sscanf(items[p->headerIndex[5]],"%f",&M[N].info);
	if(p->headerIndex[6]>=0) sscanf(items[p->headerIndex[6]],"%f",&M[N].score);
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
	if(M[N].scoreThresh && p->threshAdd){
		M[N].scoreThresh += p->threshAdd;
		if(M[N].scoreThresh > sumMax*0.9){ M[N].scoreThresh = sumMax*0.9; }
	}
	if(M[N].scoreThresh < 1.301) M[N].scoreThresh = 1.301;
	N++;
}
fclose(fp);
ch = truncate_filename(filename);
printf("File %s loaded. Input motifs: %d\n",ch,N);
free(buffer);
return(NNN);
}

/***********************************************/
int  *convert_to_numbers (char *x, int repeats, long *totalLen)
/***********************************************/
{
long len, i, pos;
int *y, xx;
char ch;

len = strlen(x);
check(y = (int*)malloc(len*sizeof(int)));
*totalLen += len;
if(repeats==0){
	for(i=0;i<len;++i){
		ch = x[i];
		if(ch=='A'){ xx=0; }
		else if(ch=='C'){ xx=1; }
		else if(ch=='G'){ xx=2; }
		else if(ch=='T'){ xx=3; }
		else{ xx=-1; --*totalLen; }
		y[i]=xx;
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
		y[i]=xx;
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
		y[i]=xx;
	}
}
return(y);
}

/***********************************************/
long read_file (char *filename, int rep, int **sequence, long *seqlength, long *nseq, char **seqNames)
/***********************************************/
{
char *seq, *buffer,*seqname,*seqname1;
long seqlen=0, len, n1=0, parts=0;
FILE *fp;
long totalLen=0;
char *ch;

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
			if(!parts) seqNames[n1] = copy_string(seqname);
			else{
				sprintf(seqname1,"%s-%d",seqname,++parts);
				seqNames[n1] = copy_string(seqname1);
			}
			sequence[n1] = convert_to_numbers(seq, rep, &totalLen);
			if(seqlen<10){ free(sequence[n1]); }
			else{ seqlength[n1++] = seqlen; }
			if(n1>=MAXSEQ-1) error_message("Too many sequences in input file");
		}
		strcpy(seqname,&buffer[1]);
		seq[0] = '\0';
		seqlen = 0;
		parts = 0;
	}else{
		seqlen += len;
		strcat(seq, buffer);
		if(seqlen+len >= 60000){
			parts++;
			sprintf(seqname1,"%s-%d",seqname,parts);
			seqNames[n1] = copy_string(seqname1);
			char a = seq[60000];
			seq[60000] = '\0';
			sequence[n1] = convert_to_numbers(seq, rep, &totalLen);
			seqlength[n1++] = 60000;
			if(n1>=MAXSEQ-1) error_message("Too many sequences in input file");
			seq[60000] = a;
			memmove(seq,&seq[60000],(MAXSEQ-60000)*sizeof(char));
			seq[0] = '\0';
			seqlen -= 60000;
		}
	}
}
if(!parts) seqNames[n1] = copy_string(seqname);
else{
	sprintf(seqname1,"%s-%d",seqname,++parts);
	seqNames[n1] = copy_string(seqname1);
}
sequence[n1] = convert_to_numbers(seq, rep, &totalLen);
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
long   get_sequence_index (char *seqName, long nseq, char **seqNames)
/***********************************************/
{
long iseq=0;
while(iseq<nseq && strcmp(seqName,seqNames[iseq])) ++iseq;
if(iseq==nseq) iseq=-1;
return(iseq);
}

/***********************************************/
int   read_conservation (char *filename, int **conservation, long nseq, char **seqNames, long *seqLength)
/***********************************************/
{
char *buffer,*seqname,*seqname1;
int *cons, interval=1;
FILE *fp;
char *ch, **items;
long seqlen=0, parts=0, i;

check(buffer = (char*)malloc(3100*sizeof(char)));
check(seqname = (char*)malloc(3000*sizeof(char)));
check(seqname1 = (char*)malloc(3000*sizeof(char)));
check(items = (char**)malloc(2000*sizeof(char*)));
check(cons = (int*)malloc(MAXINT*sizeof(int)));
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("ERROR: Conservation file %s not found",ch); exit(0); }
fgets(buffer,3099,fp);
ch = strstr(buffer,"interval");
if(!ch) ch = strstr(buffer,"Interval");
if(ch){
	sscanf(ch+9,"%d",&interval);
	if(!interval) interval=1;
}
if(buffer[0]=='>'){ rewind(fp); }
long maxLength = 60000/interval;
while(fgets(buffer,3099,fp)){
	int len = strlen(buffer);
	if(buffer[len-1]=='\n'){ buffer[len-1]='\0'; }
	if(buffer[0]=='>'){
		if(seqlen > 0){
			if(parts) sprintf(seqname1,"%s-%d",seqname,++parts);
			else strcpy(seqname1,seqname);
			int ind = get_sequence_index(seqname1,nseq,seqNames);
			if(ind>=0){
				long len2 = seqLength[ind]/interval+1;
				if(len2<seqlen){ len2=seqlen; }
				check(conservation[ind] = (int*)calloc(len2,sizeof(int)));
				memcpy(conservation[ind],cons,seqlen*sizeof(int));
			}
			seqlen = 0;
			parts = 0;
		}
		strcpy(seqname,&buffer[1]);
	}else{
		int len1 = split_string(buffer,items,200,',');
		for(i=0; i<len1; ++i){
			float x;
			sscanf(items[i],"%f",&x);
			cons[seqlen+i] = (int)x;
		}
		seqlen += len1;
		if(seqlen+len >= maxLength){
			parts++;
			sprintf(seqname1,"%s-%d",seqname,parts);
			int ind = get_sequence_index(seqname1,nseq,seqNames);
			if(ind>=0){
				check(conservation[ind] = (int*)malloc(maxLength*sizeof(int)));
				memcpy(conservation[ind],cons,maxLength*sizeof(int));
			}
			memmove(cons,&cons[maxLength],(MAXINT-maxLength)*sizeof(int));
			seqlen -= maxLength;
		}
	}
}
if(parts) sprintf(seqname1,"%s-%d",seqname,++parts);
else strcpy(seqname1,seqname);
int ind = get_sequence_index(seqname1,nseq,seqNames);
if(ind>=0){
	long len2 = seqLength[ind]/interval+1;
	if(len2<seqlen){ len2=seqlen; }
	check(conservation[ind] = (int*)calloc(len2,sizeof(int)));
	memcpy(conservation[ind],cons,seqlen*sizeof(int));
}
fclose(fp);
free(cons);
free(buffer);
free(seqname);
free(seqname1);
free(items);
return(interval);
}

/***********************************************/
void  matrix_offset (MATRIX *Matrix, int nPattern, PARAM *p)
/***********************************************/
{
float maxEntropy = 2, pp[4], entropy, *info, infoCum, sum, sumMean, sumVar, sumMax, sumMean1, sumVar1;
float *score;
int i, j, k, rep, NSIM = 10000;
MATRIX *M;

double log2 = log(2.0);
double log10 = log(10.0);
check(info = (float*)malloc(MAXLENGTH*sizeof(float)));
check(score = (float*)malloc(NSIM*sizeof(float)));
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

	//Monte-Carlo method to find threshold adjustment for the PFM core
	float adjust = 0;
	if(M->len > 8){
		for(rep=0; rep<NSIM; ++rep){
			sum = 0;
			for(j=0; j<M->len; ++j){
				if(j>=M->offset && j<M->offset+8) continue;
				int x = rand()%4;
				sum += M->m[j][x];
			}
			score[rep] = sum;
			//printf("%.4f\n",sum);
		}
		sortem(NSIM, score, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
		int nThr = NSIM*0.999;
		adjust = score[nThr]+0.1;
	}
	if(!M->scoreThresh && !p->falsePositives){ p->falsePositives=5; }
	if(p->falsePositives){
		M->scoreThresh = 1.301;
	}
	M->scoreThreshCore = M->scoreThresh - adjust;
	//printf("%.4f %.4f %.4f\n",M->scoreThresh,adjust,M->scoreThreshCore);
}
free(info);
free(score);
return;
}

/***********************************************/
void   clean_lookup_table(LOOKUP *lookupTable)
/***********************************************/
{
long ind;
for(ind=0; ind<MAXINT; ++ind){
	lookupTable[ind].num = 0;
	if(lookupTable[ind].alloc){
		free(lookupTable[ind].score);
		free(lookupTable[ind].imatr);
		free(lookupTable[ind].dir);
		lookupTable[ind].alloc = 0;
	}
}
return;
}

/***********************************************/
void  make_lookup_table (LOOKUP *lookupTable, MATRIX *Matrix, int nPattern, long *complement, int strand_option)
/***********************************************/
{
long ind, ind1;
int i, j, pos[8], k, im;
MATRIX *M;
LOOKUP *L;
float score, missing;

for(ind=0; ind<MAXINT; ++ind){
	lookupTable[ind].num=0;
}
for(ind=0; ind<MAXINT; ++ind){
	ind1 = ind;
	for(i=0; i<8; ++i){
		k = ind1%4;
		pos[i] = k;
		ind1 = (ind1-k)/4;
	}
	for(im=0; im<nPattern; ++im){
		M = &Matrix[im];
		score = 0;
		for(i=0; i<8; ++i){
			j = i + M->offset;
			if(j<M->len) score += M->m[j][pos[i]];
		}
		if(score >= M->scoreThreshCore){
			//printf("A %d %.4f\n",ind,score);
			L = &lookupTable[ind];
			if(!L->alloc){
				L->alloc = 10;
				check(L->score = (float*)malloc(L->alloc*sizeof(float)));
				check(L->imatr = (int*)malloc(L->alloc*sizeof(int)));
				check(L->dir = (int*)malloc(L->alloc*sizeof(int)));
			}else if(L->num==L->alloc){
				L->alloc += 20;
				check(L->score = (float*)realloc(L->score,L->alloc*sizeof(float)));
				check(L->imatr = (int*)realloc(L->imatr,L->alloc*sizeof(int)));
				check(L->dir = (int*)realloc(L->dir,L->alloc*sizeof(int)));
			}
			L->imatr[L->num] = im;
			L->score[L->num] = score;
			L->dir[L->num] = 1;
			L->num++;
			if(strand_option) continue;
			ind1 = complement[ind];
			if(ind1==ind) continue;
			L = &lookupTable[ind1];
			if(!L->alloc){
				L->alloc = 10;
				check(L->score = (float*)malloc(L->alloc*sizeof(float)));
				check(L->imatr = (int*)malloc(L->alloc*sizeof(int)));
				check(L->dir = (int*)malloc(L->alloc*sizeof(int)));
			}else if(L->num==L->alloc){
				L->alloc += 20;
				check(L->score = (float*)realloc(L->score,L->alloc*sizeof(float)));
				check(L->imatr = (int*)realloc(L->imatr,L->alloc*sizeof(int)));
				check(L->dir = (int*)realloc(L->dir,L->alloc*sizeof(int)));
			}
			L->imatr[L->num] = im;
			L->score[L->num] = score;
			L->dir[L->num] = -1;
			L->num++;
		}
	}
}
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
return;
}

/***********************************************/
int  remove_position (POSITION *P, int pos, int len)
/***********************************************/
{
int i, j, k;

i=P->num-1;
while(i>=0 && P->pos[i] > pos) --i; 
if(P->pos[i]==pos){
	for(j=i; j<P->num-1; ++j){
		P->dir[j] = P->dir[j+1];
		P->score[j] = P->score[j+1];
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
void  find_motifs (int **sequence, long *seqlen, long nseq, LOOKUP *lookupTable,
		MATRIX *Matrix, int nPattern, PARAM *p, FILE *fp, char **seqNames, int **conservation, float **scores, long *nhits)
/***********************************************/
{
long i, j, k, iseq, pos, ind, lastMissing=0;
int mult, *seq, dir, *pattern, N, x, im, *cons=NULL;
int nstack=0;
float score, score1, *stack;
MATRIX *M;
LOOKUP *L;
POSITION *position, *P;
char *buffer;

check(buffer = (char*)malloc(30300*sizeof(char)));
check(stack = (float*)malloc(NSTACK*sizeof(float)));
check(pattern = (int*)malloc(MAXLENGTH*sizeof(int)));
check(position = (POSITION*)calloc(nPattern,sizeof(POSITION)));
for(im=0; im<nPattern; ++im) nhits[im] = 0;
for(iseq=0; iseq<nseq; ++iseq){
	for(im=0; im<nPattern; ++im){
		position[im].num = 0;
	}
	seq = sequence[iseq];
	ind = 0;
	mult = 1;
	lastMissing = -1;
	nstack=0;
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
			for(j=nstack-1; j>=0 && p->redundancy==0; j-=3){
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
					if(!remove_position(&position[im1],pos1,M->len)) error_message("Deletion Failed");
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
				check(P->pos = (long*)malloc(P->alloc*sizeof(long)));
				check(P->dir = (int*)malloc(P->alloc*sizeof(int)));
				check(P->pattern = (int*)malloc(P->alloc*M->len*sizeof(int)));
			}else if(N == P->alloc){
				P->alloc += 100;
				check(P->score = (float*)realloc(P->score,P->alloc*sizeof(float)));
				check(P->pos = (long*)realloc(P->pos,P->alloc*sizeof(long)));
				check(P->dir = (int*)realloc(P->dir,P->alloc*sizeof(int)));
				check(P->pattern = (int*)realloc(P->pattern,P->alloc*M->len*sizeof(int)));
			}
			P->dir[N] = dir;
			P->score[N] = score;
			P->pos[N] = start;
			for(j=0; j<M->len; ++j){
				P->pattern[N*M->len+j] = pattern[j];
				//if(pattern[j]<0) P->pattern[N*M->len+j]=0;
				//if(pattern[j]>3) P->pattern[N*M->len+j]=0;
			}
			P->num++;
		}
	}
	if(conservation) cons = conservation[iseq];
	int found=0;
	for(im=0; im<nPattern; ++im){
		M = &Matrix[im];
		P = &position[im];
		int len = M->len;
		if(P->num) found=1;
		for(i=0; i<P->num; ++i){
			if(fp){  //print output
				for(j=0; j<len; ++j){
					buffer[j] = nucleotides[P->pattern[i*len+j]][0];
				}
				buffer[j]='\0';
				fprintf(fp,"%s\t%s\t%s\t%d\t%d\t%.4f\t%s",M->name,seqNames[iseq],strands[(P->dir[i]+1)/2],len,P->pos[i],P->score[i],buffer);
				if(cons){
					int meanCons=0, ncons=0;
					for(j=P->pos[i]/p->interval; j<=(P->pos[i]+len)/p->interval; ++j){
						meanCons += cons[j];
						++ncons;
					}
					if(ncons) meanCons /= ncons;
					fprintf(fp,"\t%d",meanCons);
				}
				fprintf(fp,"\n");
			}
			if(scores){
				scores[im][nhits[im]] = P->score[i];
			}
			nhits[im]++;
		}
	}
	if(!found && fp){
		fprintf(fp,"NONE\t%s\t0\t0\t0\t0\t0\n",seqNames[iseq]);
	}
}
for(im=0; im<nPattern; ++im){
	P = &position[im];
	if(P->alloc){
		free(P->dir);
		free(P->score);
		free(P->pos);
		free(P->pattern);
	}
}
free(position);
free(pattern);
free(stack);
free(buffer);
return;
}

/***********************************************/
PARAM *read_parameters (int nargs, char **argv)
/***********************************************/
{
PARAM *p;
int iarg=1, score=-1;
float max_length;

if(nargs < 4){ printf("patternScan -i motifFile -f fastaFile -o outputFile [-cons conservationFile, -maxlen maxSequenceLength, -strand strandOption, -n numberOfMotifs, -fp falsePositives, -thresh, matchThreshold, -userep, -one]"); exit(0); }
check(p = (PARAM*)calloc(1,sizeof(PARAM)));
p->max_length = 50000000;
p->N_matrix = 400;
while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1) p->inputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-f") && iarg < nargs-1) p->testFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-cons") && iarg < nargs-1) p->conservFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-maxlen") && iarg < nargs-1) sscanf(argv[++iarg],"%f",max_length);
	else if(!strcmp(argv[iarg],"-thresh") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->threshAdd);
	else if(!strcmp(argv[iarg],"-fp") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->falsePositives);
	else if(!strcmp(argv[iarg],"-userep")) p->use_repeats = 1;
	else if(!strcmp(argv[iarg],"-strand") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->strand);
	else if(!strcmp(argv[iarg],"-n") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->N_matrix);
	else if(!strcmp(argv[iarg],"-one")) p->presence_option = 1;
	else if(!strcmp(argv[iarg],"-redund")) p->redundancy = 1;
	else{
		printf("Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}
if(max_length){ p->max_length = max_length+1; }
return(p);
}

/***********************************************/
int main (int argc, char **argv) 
/***********************************************/
{
MATRIX *Matrix=NULL, *M;
PARAM *p;
FILE *fp, *fp1;
LOOKUP *lookupTable;
long *seqlen1, *complement, ind, ind1, ind2, n1;
int len, len1, i, j, im, nPattern, **sequence1, **conservation=NULL;
long length1=0, *nhits;
float x, y;
char **seqNames, *file_motif, *file_fasta;

p = read_parameters(argc, argv);
nPattern = read_input (p->inputFile,&Matrix,p);
if(nPattern<=0) error_message("Input motif file is empty");
check(seqlen1 = (long*)malloc(MAXSEQ*sizeof(long)));
check(sequence1 = (int**)malloc(MAXSEQ*sizeof(int*)));
check(lookupTable = (LOOKUP*)calloc(MAXINT,sizeof(LOOKUP)));
check(complement = (long*)malloc(MAXINT*sizeof(long)));
check(seqNames = (char**)malloc(MAXSEQ*sizeof(char*)));
check(nhits = (long*)malloc(nPattern*sizeof(long)));
length1 = read_file(p->testFile,p->use_repeats,sequence1,seqlen1,&n1,seqNames);
if(length1 < 10){ error_message("Test sequence is too short (<10)"); }
if(length1 < 5000){ printf("WARNING: Test sequence is too short (<5000)"); }
printf("Number of test sequences = %d. Total length = %d\n",n1,length1);
if(p->conservFile){
	check(conservation = (int**)calloc(n1,sizeof(int*)));
	p->interval = read_conservation(p->conservFile,conservation,n1,seqNames,seqlen1);
}
for(i=0; i<MAXINT; ++i){
	complement[i] = complementary(i);
}
matrix_offset(Matrix,nPattern,p);
if(p->falsePositives){
	int **seqRandom, im;
	long *seqlenRand, nseqRand=0;
	float **scores;
	check(seqRandom = (int**)malloc(MAXSEQ*sizeof(int*)));
	check(seqlenRand = (long*)malloc(MAXSEQ*sizeof(long)));
	check(scores = (float**)malloc(nPattern*sizeof(float*)));
	for(im=0; im<nPattern; ++im)
		check(scores[im] = (float*)malloc(MAXSEQ*sizeof(float)));
	long lengthRand = get_random_sequence(seqRandom,seqlenRand,&nseqRand,p->falsePositives);
	make_lookup_table(lookupTable,Matrix,nPattern,complement,p->strand);
	find_motifs(seqRandom,seqlenRand,nseqRand,lookupTable,Matrix,nPattern,p,NULL,seqNames,NULL,scores,nhits);
	clean_lookup_table(lookupTable);
	for(im=0; im<nPattern; ++im){
		M = &Matrix[im];
		long N = nhits[im];
		long expected = (double)lengthRand/10000*p->falsePositives;
		if(N >= expected || M->scoreThresh <= 1.301){
			float x = M->scoreThresh;
			if(N >= expected){
				sortem(N, scores[im], 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
				M->scoreThresh = scores[im][N-expected];
			}
			if(M->scoreThresh < 1.301) M->scoreThresh = 1.301;
			M->scoreThreshCore += M->scoreThresh - x;
			//printf("Done: %d %.4f\n",im,M->scoreThresh);
			continue;
		}
	}
	for(i=0; i<nseqRand; ++i) free(seqRandom[i]);
	for(im=0; im<nPattern; ++im) free(scores[im]);
	free(scores);
	free(seqRandom);
	free(seqlenRand);

}
make_lookup_table(lookupTable,Matrix,nPattern,complement,p->strand);
fp = fopen(p->outputFile,"w");
if(!fp) error_message("Output file not opened");
file_motif = truncate_filename(p->inputFile);
file_fasta = truncate_filename(p->testFile);
fprintf(fp,"Parameters:\tfile_motif=%s\tfile_fasta=%s\n",file_motif,file_fasta);
fprintf(fp,"Headers:\tMotifName\tSeqName\tStrand\tLen\tStart\tScore\tSequence");
if(p->conservFile) fprintf(fp,"\tConservation");
fprintf(fp,"\n");
find_motifs(sequence1,seqlen1,n1,lookupTable,Matrix,nPattern,p,fp,seqNames,conservation,NULL,nhits);
long N = 0;
for(im=0; im<nPattern; ++im){
	N += nhits[im];
}
if(!N) error_message("No motifs found! Try using less stringent conditions\n");
fclose(fp);
printf("patternScan done\n");
return(0);
}
