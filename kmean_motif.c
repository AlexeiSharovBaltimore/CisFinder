#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define MAXLEN 55
#define MAXINT 65536
#define MAX_ITERATIONS 30

//gcc kmean_motif.c -lm -o kmean_motif.exe
//kmean_motif -i motifs.txt -o kmean_out.txt -k 10 -n 1000 -meme clusters.meme
//This version reads MEME format and outputs MEME format, uses get_distance instead of get_distance2

int DEBUG=0;

typedef struct point_st{
	int ID;
	int clusterID;
	int len;
	int freq;
	float *m;
	char *name;
	char *pattern;
}POINT;
typedef struct pair_st{
	int ID1;
	int ID2;
	float dist;
}PAIR;
typedef struct param_st{
	char *inputFile;
	char *outputFile;
	char *clusterFile;
	int k;
	int n;
}PARAM;


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

/***********************************************/
char  *reverse_pattern(char *pattern, char *revpattern) 
/***********************************************/
{
static char *forward = "ACGTRYMKacgtrymk";
static char *reverse = "TGCAYRKMtgcayrkm";
int j=0;
for(int ii = strlen(pattern)-1; ii >= 0; ii--){
	char *ch1 = strchr(forward, pattern[ii]);
	if(ch1) revpattern[j] = reverse[ch1 - forward];
	else revpattern[j] = pattern[ii];
	j++;
}
revpattern[strlen(pattern)] = '\0';
return revpattern;
}

/***********************************************/
int  get_pattern(char *pattern, float *m, int len) 
/***********************************************/
{
static char *nuc_bg = "ACGT";
static char *nuc_sm = "acgt";
static char *nuc_pairsm = "mrwsyk";
float x[4];

pattern[0] = '\0';
for(int p1=0; p1<len; p1++){
	int imax,i1,i2,imin,i3,sum=0;
	float xmax=0,xmin=100000;
	for(int j1=0; j1<4; j1++){
		x[j1] = m[p1*4+j1];
		sum += x[j1];
		if(xmax < x[j1]){ imax=j1; xmax=x[j1]; }
		if(xmin > x[j1]){ imin=j1; xmin=x[j1]; }
	}
	if(sum==0) return(p1); //Correct the error
	i1=0; while(i1==imax || i1==imin) i1++;
	i2=i1+1;while(i2==imax || i2==imin) i2++;
	if(x[i1] < x[i2]){ i3=i2; i2=i1; i1=i3; }
	int code;
	if(imax==0){ code=i1; }				
	else if(i1==0){ code=imax; }				
	else if(imax==1 && i1==2 || imax==2 && i1==1){ code=3; }				
	else if(imax==1 && i1==3 || imax==3 && i1==1){ code=4; }				
	else{ code=5; }
	if(x[imax] >= 2*x[i1]){
		strncat(pattern, &nuc_bg[imax], 1);
	}else if(x[imax] >= 1.3*x[i1]){
		strncat(pattern, &nuc_sm[imax], 1);
	}else{
		if(x[i1] >= 2*(x[i2])){
			strncat(pattern, &nuc_pairsm[code], 1);
		}else{
			strncat(pattern, "_", 1);
		}
	}
}
return(len);
}

/*************************************/
float  get_distance  (POINT* P1, POINT* P2, int *nMatch, int *offset, int *direction)
/*************************************/
//P1 - cluster, P2 - point
{
int N1, N2, offset1, offset_start, offset_end, off, i, j, i1, i2, dir, Nmin;
float dist, dist1, dist2, thresh=75;
static float aver[4]={21.86,29.75,26.26,22.22};

if(!P1 || !P2){
	printf("get_distance null-pointers\n");
	return(100);
}
N1 = P1->len;
N2 = P2->len;
if(N1 > MAXLEN || N2 > MAXLEN){ printf("get_distance m length\n"); exit(0); }

dist = 1.0E20;
offset_start = 6-N2;
offset_end = N1-6;
Nmin = N2;
if(N2>N1){
	Nmin = N1;
}
for(off=offset_start; off<=offset_end; ++off){
	int n=0, n1=0, start=0, end;
	double sxy=0, sxy1=0, sxx=0, syy=0, syy1=0, sx=0, sy=0, sy1=0;
	if(off<0) start = off;
	end = N1;
	if(off+N2 > N1) end=off+N2;
	for(i1=start; i1<end; ++i1){
		i2 = i1-off;
		n++;
		for(j=0; j<4; ++j){
			float x, y, y1;
			if(i1>=0 && i1<N1) x = P1->m[i1*4+j];
			else x = aver[j];
			if(i2>=0 && i2<N2){
				y = P2->m[i2*4+j];
				y1 = P2->m[(N2-1-i2)*4+(3-j)];
			}else{
				y = aver[j];
				y1 = aver[j];
			}
			if(x > thresh){ x = thresh + 3*(x-thresh); }
			if(y > thresh){ y = thresh + 3*(y-thresh); }
			if(y1 > thresh){ y1 = thresh + 3*(y1-thresh); }
			sxy += x*y;
			sxy1 += x*y1;
			sxx += x*x;
			syy += y*y;
			syy1 += y1*y1;
			sx += x;
			sy += y;
			sy1 += y1;
			n1++; 
		}
	}
	dist1 = 1000;
	dist2 = 1000;
	dir = 0;
	if(n>5){
		float ssx = sxx-sx*sx/n1;
		dist1 = 100*(1.0-(sxy-sx*sy/n1)/(float)sqrt(ssx*(syy-sy*sy/n1)));
		dist2 = 100*(1.0-(sxy1-sx*sy1/n1)/(float)sqrt(ssx*(syy1-sy1*sy1/n1)));
		if(dist1 > dist2){
			dist1 = dist2;
			dir = 1;
		}
	}
	if(dist > dist1){
		*offset = off;
		*direction = dir;
		*nMatch = n;
		dist = dist1;
	}
	if(dist<0 && dist > -0.00001){ dist=0.001; }
	if(dist<0){ printf("D %.5f\n",dist); exit(0); }
}
return(dist);
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
POINT *read_file_meme (char *filename, int *nPoints, PARAM *p)
/***********************************************/
{
char *buffer, *ch, *items[600];
float aver[4], x[4], *m;
FILE *fp;
POINT *P;
int N=0, i,j,NNN,ii, len=0, nFreq=0;

check(buffer = (char*)malloc(3201*sizeof(char)));
fp = fopen(filename,"r");
ch = truncate_filename(filename);
if(!fp){ printf("Input file %s not found",ch); exit(0); }
while(fgets(buffer,3200,fp)){
	ch = strstr(buffer,"letter");
	if(ch != NULL){ N++; }
}
rewind(fp);
NNN = N;
for(i=0; i<4; i++) aver[i] = 0;
check(P = (POINT*)calloc(NNN,sizeof(POINT)));
ii=0;
while(fgets(buffer,3200,fp)){
	ch = strstr(buffer,"MOTIF");
	if(ch == NULL){ continue; }
	ch[strlen(ch)-1]='\0';
	P[ii].name = copy_string(ch + 6);
	while(fgets(buffer,3200,fp)){
		ch = strstr(buffer,"letter");
		if(ch){ break; }
	}
	P[ii].freq = 100;
	ch = strstr(buffer,"nsites");
	if(ch){ sscanf(ch+8,"%d",&P[ii].freq); }
	check(m = (float*)calloc(MAXLEN*4,sizeof(float)));
	i = 0;
	check(P[ii].m = (float*)calloc(MAXLEN*4,sizeof(float)));
	while(fgets(buffer,3200,fp)){
		if(i >= MAXLEN){printf("Max length err %d >= %d\n",len, MAXLEN); exit(0); }
		j = sscanf(buffer,"%f%f%f%f",&x[0],&x[1],&x[2],&x[3]);
		float sum = x[0]+x[1]+x[2]+x[3];
		if(j !=4 || sum<0.00001){ break; }
		for(j=0; j<4; j++){
			P[ii].m[i*4+j] = 100*x[j]/sum;
			aver[j] += P[ii].m[i*4+j];
		}
		i++;
	}
	nFreq += i;
	free(m);
	check(P[ii].pattern = (char*)calloc(MAXLEN,sizeof(char)));
	len = get_pattern(P[ii].pattern, P[ii].m, i);
	P[ii].len = len;
	if(len < 6){
		free(P[ii].name);
		free(P[ii].pattern);
		free(P[ii].m);
		continue;
	}
	ii++;
}
fclose(fp);
if(NNN > ii) NNN = ii;
if(NNN>p->n) NNN = p->n;
if(!NNN) error_message("No input motifs found!");
*nPoints = NNN;
for(i=0; i<4; i++) aver[i] /= nFreq;
printf("File %s loaded\nN input patterns = %d\nFrequency ACGT: %.2f %.2f %.2f %.2f\n",filename,NNN,aver[0],aver[1],aver[2],aver[3]);
free(buffer);
return(P);
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

if(nargs < 4){ printf("patternClusterNew -i inputMotifFile -o outputFile [-k Nclusters -meme ClusterOutput]\n"); exit(0); }
check(p = (PARAM*)calloc(1,sizeof(PARAM)));
p->k = 20;
p->n = 100000;
while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1) p->inputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-meme") && iarg < nargs-1) p->clusterFile=copy_string(argv[++iarg]);
	else if(!strcmp(argv[iarg],"-k") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->k);
	else if(!strcmp(argv[iarg],"-n") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->n);
	else{
		printf("Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}
if(p->k <3) error_message("Too small number of output clusters shown after -k option");
if(p->n <3) error_message("Too small number of input motifs shown after -n option");
if(!p->inputFile) error_message("Input file name needs to be specified after -i option");
if(!p->outputFile) error_message("Input file name needs to be specified after -o option");
return(p);
}

//***********************
float calc_total_distance(int n, int k, POINT *P, POINT *C, int *index)
//***********************
 {
int nMatch, offset, dir;
float tot_D = 0;
for (int ii = 0; ii < n; ii++){
	int icluster = index[ii];
	if (icluster >= 0){
		tot_D += get_distance(&C[icluster], &P[ii], &nMatch, &offset, &dir);
	}
}	  
return tot_D;
}

//***********************
void min_offset(int n, int k, int *index, int *offset, int *off_min)
//***********************
 {
int ii, jj, icluster;
for (int ii = 0; ii < k; ii++){ off_min[ii] = 30; }
for (int ii = 0; ii < n; ii++){
	int icluster = index[ii];
	if (icluster >= 0){
		int off = offset[ii*k+icluster];
		//if(abs(off) > 10) printf("min_offset offset = %d point %d\n",off,ii);
		if(off_min[icluster] > off){
			off_min[icluster] = off;
		}
	}
}	  
for (int ii = 0; ii < k; ii++){ 
	if(off_min[ii] == 30){
		//printf("min_offset cluster%d - off=%d\n",ii,off_min[ii]);
		off_min[ii] = 0;
	}
}
return;
}

//***********************
void calc_cluster_centroids(int n, int k, POINT *P, int *cluster_assignment, POINT *C, int *N_members, int *offset, int *direction, float *centroid, int *weight, float *dist)
//***********************
{
int ii, jj, i1;

for (ii = 0; ii < k; ii++){ N_members[ii] = 0; }
for (ii = 0; ii < k*MAXLEN*4; ii++){ centroid[ii] = 0; }
for (ii = 0; ii < k*MAXLEN; ii++){ weight[ii] = 0; }
for (ii = 0; ii < n; ii++){
	int icluster = cluster_assignment[ii];
	if(icluster < 0) continue;
	N_members[icluster]++;
}
for (ii = 0; ii < n; ii++){
	int icluster = cluster_assignment[ii];
	if(icluster < 0) continue;
	float d = dist[ii*k+icluster];
	if(N_members[icluster] > 10 && d > 25) continue; // do not use far members for calculating centroids
	int off = offset[ii*k+icluster];
	int dir = direction[ii*k+icluster];
	int len = P[ii].len;
	int wgt = P[ii].freq;
	if(wgt > 5000) wgt=5000;
	for (jj = 0; jj < len; jj++){
		int pos = icluster*MAXLEN + 6 + off + jj;
		weight[pos] += wgt;
		for (int j1 = 0; j1< 4; j1++){
			if(jj*4+j1 >= MAXLEN*4 || (len-jj-1)*4+(3-j1) >= MAXLEN*4 ||jj*4+j1 <0 || (len-jj-1)*4+(3-j1)<0) printf("ERR centroid\n");

			if(dir==0){
				centroid[pos*4+j1] += wgt*P[ii].m[jj*4+j1];
			}else{
				centroid[pos*4+j1] += wgt*P[ii].m[(len-jj-1)*4+(3-j1)];
			}
		}
	}
}
for (ii = 0; ii < k; ii++){
	//printf("Cl %d  %d\n",ii,N_members[ii]);
	if (N_members[ii] == 0){
		C[ii].freq = 0;
		for (jj = 0; jj < n; jj++){
			if(cluster_assignment[jj]==ii) cluster_assignment[jj]=-1;
		}
		continue;
	}
	int start=0, end=MAXLEN-1;
	int max = 0;
	for(jj = 0; jj < MAXLEN; jj++){
		int wgt = weight[ii*MAXLEN + jj];
		if(max < wgt){ max = wgt; }
		for (int j1 = 0; j1< 4 && wgt>0; j1++){
			centroid[(ii*MAXLEN + jj)*4+j1] /= wgt;
		}
	}
	if(max<=0){
		//printf("Centroid ERROR! %d freq=0\n",ii);
		C[ii].freq = 0;
		N_members[ii] == 0;
		for (i1 = 0; i1 < n; i1++){
			if(cluster_assignment[i1]==ii) cluster_assignment[i1]=-1;
		}
		continue;
	}
	while(start<MAXLEN){
		float *x = &centroid[(ii*MAXLEN + start)*4];
		float sum1 = x[0]+x[1]+x[2]+x[3];
		if(weight[ii*MAXLEN + start] < max*0.4 || sum1==0){
			start++;
			continue;
		}
		int ix[4]={0,1,2,3}; int sw;
		if(x[1] > x[0]){ ix[0]=1; ix[1]=0; }
		if(x[3] > x[2]){ ix[2]=3; ix[3]=2; }
		if(x[ix[2]] > x[ix[0]]){ sw=ix[2]; ix[2]=ix[0]; ix[0]=sw; }
		if(x[ix[3]] > x[ix[1]]){ sw=ix[3]; ix[3]=ix[1]; ix[1]=sw; }
		if(x[ix[2]] > x[ix[1]]){ sw=ix[2]; ix[2]=ix[1]; ix[1]=sw; }
		if(x[ix[0]] < 1.3*x[ix[1]] && x[ix[1]] < 2*x[ix[2]]){
			start++;
			continue;
		}
		break;
	}
	while(end > 0){
		float *x = &centroid[(ii*MAXLEN + end)*4];
		float sum1 = x[0]+x[1]+x[2]+x[3];
		if(weight[ii*MAXLEN + end] < max*0.4 || sum1==0){
			end--;
			continue;
		}
		int ix[4]={0,1,2,3}; int sw;
		if(x[1] > x[0]){ ix[0]=1; ix[1]=0; }
		if(x[3] > x[2]){ ix[2]=3; ix[3]=2; }
		if(x[ix[2]] > x[ix[0]]){ sw=ix[2]; ix[2]=ix[0]; ix[0]=sw; }
		if(x[ix[3]] > x[ix[1]]){ sw=ix[3]; ix[3]=ix[1]; ix[1]=sw; }
		if(x[ix[2]] > x[ix[1]]){ sw=ix[2]; ix[2]=ix[1]; ix[1]=sw; }
		if(x[ix[0]] < 1.3*x[ix[1]] && x[ix[1]] < 2*x[ix[2]]){
			end--;
			continue;
		}
		break;
	}
	C[ii].len = end-start+1;
	C[ii].freq = max;
	if(max==0 || C[ii].len<5){
		//printf("Centroid ERROR! %d freq=0\n",ii);
		N_members[ii] == 0;
		for (i1 = 0; i1 < n; i1++){
			if(cluster_assignment[i1]==ii) cluster_assignment[i1]=-1;
		}
		continue;
	}
	for(int jj = start; jj <= end; jj++){
		int wgt = weight[ii*MAXLEN + jj];
		int sum=0;
		for (int j1 = 0; j1< 4; j1++){
			sum += centroid[(ii*MAXLEN + jj)*4+j1];
		}
		for (int j1 = 0; j1< 4 && wgt>0; j1++){
			C[ii].m[(jj-start)*4 + j1] = centroid[(ii*MAXLEN + jj)*4+j1] *100/sum;
		}
	}
	int sum1=0;
	for (int j1 = 0; j1< 4; j1++){
		sum1 += C[ii].m[(C[ii].len-1)*4 + j1];
	}
	if(sum1==0){ printf("ERR: zeroes in cluster\n"); }
}
return;
}

//***********************
void copy_assignment_array(int n, int *src, int *tgt)
//***********************
{
for(int ii=0; ii<n; ii++)
	tgt[ii] = src[ii];
}

//***********************
int assignment_change_count(int n, int a[], int b[])
//***********************
{
int change_count = 0;
for (int ii = 0; ii < n; ii++)
	if (a[ii] != b[ii])
		change_count++;
	
return change_count;
}

//***********************
float  choose_all_clusters_from_distances(int n, int k, float *distance_array, int *cluster_assignment_index)
//***********************
{
float sum_dist = 0;
for (int ii = 0; ii < n; ii++){
	if(cluster_assignment_index[ii] < -1){ continue; }
	int best_index = -1;
	float closest_distance=1.0E20;
	for (int jj = 0; jj < k; jj++){
		float x = distance_array[ii*k + jj];
		if (jj==0 || x < closest_distance){
			best_index = jj;
			closest_distance = x;
		}
	}
	cluster_assignment_index[ii] = best_index;
	sum_dist += closest_distance;
}
return(sum_dist);
}

//***********************
void calc_all_distances(int n, int k, POINT *P, POINT *C, float *distance, int *offset, int *direction)
//***********************
{
for (int ii = 0; ii < n; ii++){ // for each point
	for (int jj = 0; jj < k; jj++){ // for each cluster
		 // calculate distance between point and cluster centroid
		int nMatch,off, dir;
		distance[ii*k + jj] = get_distance (&C[jj], &P[ii], &nMatch, &off, &dir);
		offset[ii*k + jj] = off;
		direction[ii*k + jj] = dir;
	}
}
return;
}

/******************************************/
void  get_linkage (PAIR *pairs, int nPairs, int **linked, int *nLinked, int k)
/******************************************/
{
int i, j, ind, ind1, ind2;
int *nalloc;

check(nalloc = (int*)calloc(k,sizeof(int)));
for(i=0; i < nPairs; ++i){
	ind  = pairs[i].ID1;
	ind1 = pairs[i].ID2;
	if(!nalloc[ind]){
		nalloc[ind] = 20;
		check(linked[ind] = (int*)malloc(20*sizeof(int)));
	}else if(nalloc[ind] <= nLinked[ind]){
		nalloc[ind] += 50;
		check(linked[ind] = (int*)realloc(linked[ind],nalloc[ind]*sizeof(int)));
	}
	if(!nalloc[ind1]){
		nalloc[ind1] = 20;
		check(linked[ind1] = (int*)malloc(20*sizeof(int)));
	}else if(nalloc[ind1] <= nLinked[ind1]){
		nalloc[ind1] += 50;
		check(linked[ind1] = (int*)realloc(linked[ind1],nalloc[ind1]*sizeof(int)));
	}
	linked[ind][nLinked[ind]++] = ind1;
	linked[ind1][nLinked[ind1]++] = ind;
}
free(nalloc);
return;
}

//***********************
void  rearrange_clusters (int n, int k, POINT *C, int *index, int *N_members, int *reserved, int iter)
//***********************
{
PAIR *pairs;
int nPairs=0, ii, jj, nMatch, off1, dir1, maxPairs;
int *nLinked;
int *stack, nStack=0, Nreserved=0, **linked;
float *N_members1, *icluster;

maxPairs = (int)(k/8)*k;
check(pairs = (PAIR*)calloc(maxPairs, sizeof(PAIR)));
for(ii=0; ii<k && nPairs<maxPairs; ii++){
	for(jj=ii+1; jj<k; jj++){
		float d=get_distance(&C[ii], &C[jj], &nMatch, &off1, &dir1);
		if(d > 3){ continue; }
		pairs[nPairs].ID1 = ii;
		pairs[nPairs].ID2 = jj;
		pairs[nPairs].dist = d;
		nPairs++;
		if(nPairs==maxPairs){
			printf("Max pairs %d is reached!\n",nPairs);
			break;		
		}
	}
}
check(linked = (int**)calloc(k, sizeof(int*)));
check(nLinked = (int*)calloc(k, sizeof(int)));
check(stack = (int*)calloc(MAXINT,sizeof(int)));
check(N_members1 = (float*)calloc(k, sizeof(float)));
check(icluster = (float*)calloc(k, sizeof(float)));
for(ii=0; ii< k; ++ii){
	N_members1[ii] = N_members[ii];
	icluster[ii] = ii;
}
//printf("N pairs %d\n",nPairs);

get_linkage(pairs, nPairs, linked, nLinked, k);
//Sort clusters by the number of members
sortem(k, N_members1, 1, icluster, NULL, NULL, NULL, NULL, NULL, NULL);
for(int m=0; m< k; ++m){
	ii = (int)icluster[k-m-1];
	if(nLinked[ii]==0){ //These are unique clusters
		if(iter>0 && N_members[ii] <= 0){
			//printf("A7 %d %d\n",ii,N_members[ii]);
			reserved[Nreserved++]=ii;
		}
		continue;
	} 
	int i1, found=0;
	for(i1=0;i1<Nreserved && !found;++i1){
		if(reserved[i1]==ii) found=1;
	}
	if(found) continue;

	stack[nStack++] = ii;
	while(nStack){
		jj = stack[--nStack];
		reserved[Nreserved++]=jj;
		//printf("A3 %d\n", jj);
		for(int i1=0; i1< n; ++i1){
			if(index[i1]==jj){ index[i1] = ii; }
		}
		for(int i=0; i<nLinked[jj]; ++i){
			int jj1 = linked[jj][i];
			for(i1=0;i1<Nreserved && !found;++i1){
				if(reserved[i1]==jj1) found=1;
			}
			for(i1=0;i1<nStack && !found;++i1){
				if(stack[i1]==jj1) found=1;
			}
			if(!found && N_members[jj1] > 1){
				stack[nStack++] = jj1;
			}
		}
	}
}
reserved[Nreserved++] = -1;
for(ii=0; ii< k; ++ii) if(nLinked[ii]) free(linked[ii]);
free(linked);
free(nLinked);
free(stack);
free(N_members1);
free(icluster);
return;
}

//***********************
void  add_remote_clusters(int n, int k, int *index, int *reserved, float *dist)
//***********************
{
float *dd, *ip;
int ii, icluster;

check(dd = (float*)calloc(n, sizeof(float*)));
check(ip = (float*)calloc(n, sizeof(float*)));
for(ii=0; ii<n; ii++){
	icluster = index[ii];
	if(icluster >= 0){
		dd[ii] = dist[ii*k+icluster];
	}
	ip[ii] = (float)ii;
}
sortem(n, dd, 1, ip, NULL, NULL, NULL, NULL, NULL, NULL);
int ires = 0;
ii = n - rand()%(n/10);
if(ii==n) ii--;
while(reserved[ires] >=0){
	//printf("Add cluster %d - %d %.2f\n",reserved[ires],(int)ip[ii],dd[ii]);
	index[(int)ip[ii--]] = reserved[ires++];
}
free(dd);
free(ip);
return;
}

//***********************
void  assign_cluster(int icl, int n, int k, POINT *P, POINT *C, int *cluster)
//***********************
{
int i1 = rand()%n;
while(cluster[i1] > -1 && i1<n-1) i1++;
if(cluster[i1]>-1){
	i1=0;
	while(cluster[i1] > -1 && i1<n-1) i1++;
}
cluster[i1] = icl;
C[icl].len = P[i1].len;
for(int jj=0; jj<P[i1].len*4; jj++){
	C[icl].m[jj] = P[i1].m[jj];
}
return;
}

//***********************
void kmeans(POINT *Point, int n, int k, POINT *C, int *cluster_assignment_final, int *offset, int *direction, float *dist, PARAM *p)
//***********************
{
int *cluster_assignment_cur  = (int *)malloc(sizeof(int) * n);
int *cluster_assignment_prev = (int *)malloc(sizeof(int) * n);
int *N_members = (int *)malloc(sizeof(int) * k);
int *reserved = (int*)malloc(k * sizeof(int));
float *centroid = (float*)malloc(k*MAXLEN*20*sizeof(float));
int *weight = (int*)malloc(k*MAXLEN*10* sizeof(int));
int *off_min = (int *)calloc(k, sizeof(int));
char pattern[MAXLEN], patternRev[MAXLEN];
FILE *fp;
int stop_iteration = MAX_ITERATIONS;

float prev_totD, totD;
int ii, jj;

if (!dist || !cluster_assignment_cur || !weight)
	error_message("Error allocating dist arrays");
	
// initial setup
for(ii=0; ii<n; ii++){ cluster_assignment_cur[ii] = -1; }
for(ii=0; ii<k; ii++){
	assign_cluster(ii, n, k, Point, C, cluster_assignment_cur);
	N_members[ii] = 1;	
}
rearrange_clusters(n, k, C, cluster_assignment_cur, N_members, reserved, 0);
int ires = 0;
while(reserved[ires]>=0){
	ii = reserved[ires++];
	assign_cluster(ii, n, k, Point, C, cluster_assignment_cur);
	N_members[ii] = 1;	
}
calc_all_distances(n, k, Point, C, dist, offset, direction);
prev_totD = choose_all_clusters_from_distances(n, k, dist, cluster_assignment_cur);

// BATCH UPDATE
int batch_iteration = 0;
printf("\nIteration N_changed totDist totDistChange\n");
while (batch_iteration < stop_iteration){
	int nMatch, off1, dir1, off2, dir2;
	copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);
	rearrange_clusters(n, k, C, cluster_assignment_cur, N_members, reserved, batch_iteration);
	add_remote_clusters(n, k, cluster_assignment_cur, reserved, dist);
	calc_cluster_centroids(n, k, Point, cluster_assignment_cur, C, N_members,offset,direction, centroid, weight, dist);
	calc_all_distances(n, k, Point, C, dist, offset, direction);
	totD = choose_all_clusters_from_distances(n, k, dist, cluster_assignment_cur);
	for (int ii = 0; ii < n; ii++){
		int icluster = cluster_assignment_cur[ii];
		if (icluster >= 0){
			if(icluster<0){ printf("kmeans1 cluster %d\n",icluster); continue; }
			float d=get_distance(&C[icluster],&Point[ii], &nMatch, &off1, &dir1);
			if(off1 != offset[ii*k+icluster]) printf("Err OFF %d %d %d\n",ii,off1,offset[ii*k+icluster]);
			if(dir1 != direction[ii*k+icluster]) printf("Err DIR %d %d %d\n",ii,dir1,direction[ii*k+icluster]);
			if(d != dist[ii*k+icluster]) printf("Err DIST %d %.2f %.2f\n",ii,d,dist[ii*k+icluster]);
		}
	}
	int change_count = assignment_change_count(n, cluster_assignment_cur, cluster_assignment_prev);
	printf("%d  %d %6.2f %6.2f\n", batch_iteration,change_count, totD, totD - prev_totD);
	fflush(stdout);
	if (totD > prev_totD && batch_iteration > 5){
		stop_iteration = batch_iteration+3;
	}else if(stop_iteration <= batch_iteration+3){
		printf("Last iteration completed\n");
		break;
	}
	if (change_count == 0){
		printf("No change made on this step - iteration completed \n");
		break;
	}
	prev_totD = totD;
	batch_iteration++;
}
if (batch_iteration == MAX_ITERATIONS){ printf("Max iterations completed \n"); }
copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_final);	
calc_all_distances(n, k, Point, C, dist, offset, direction);

int nMatch, off1, dir1;
for (int ii = 0; ii < n; ii++){
	int icluster = cluster_assignment_cur[ii];
	if (icluster >= 0){
		//printf("M1 %d %d\n", ii,icluster);
		float d=get_distance(&C[icluster],&Point[ii], &nMatch, &off1, &dir1);
		if(off1 != offset[ii*k+icluster]) printf("Err OFF %d %d %d\n",ii,off1,offset[ii*k+icluster]);
		if(dir1 != direction[ii*k+icluster]) printf("Err DIR %d %d %d\n",ii,dir1,direction[ii*k+icluster]);
		if(d != dist[ii*k+icluster]) printf("Err DIST %d %.2f %.2f\n",ii,d,dist[ii*k+icluster]);
	}
}
printf("Main output\n");
fp = fopen(p->outputFile,"w");
if(!fp) error_message("Output file not opened");
min_offset(n, k, cluster_assignment_cur, offset, off_min);
for(ii=0; ii<n; ii++){
	int dir=0, off=0;
	float d=0;
	jj = cluster_assignment_cur[ii];
	if(jj >= 0){
		dir = direction[ii*k+jj];
		d = dist[ii*k+jj];
		off = offset[ii*k+jj] - off_min[jj];
	}
	strcpy(pattern,	Point[ii].pattern);
	if(dir>0) reverse_pattern(Point[ii].pattern, pattern);
	reverse_pattern(pattern,patternRev);
	for(int i1=strlen(pattern); i1>=0; i1--){ pattern[i1+off] = pattern[i1]; }
	for(int i1=0; i1<off; i1++){ pattern[i1]=' '; }
	fprintf(fp,"%s\t%d\t%d\t%d\t%.2f\t%s\t%s\n",Point[ii].name,jj+1,Point[ii].freq,off,d,pattern,patternRev);
}
fclose(fp);
fp = fopen(p->clusterFile,"w");
fprintf(fp,"MEME version 4.5\nALPHABET= ACGT\n\n");
for(ii=0; ii<k; ii++){
	if(C[ii].freq==0){ continue; }
	char pattern[MAXLEN];
	int len1 = get_pattern(pattern, C[ii].m, C[ii].len);
	reverse_pattern(pattern,patternRev);
	if(len1 < C[ii].len) C[ii].len = len1;
	fprintf(fp,"MOTIF Cluster%d_%s_%s\n", ii+1, pattern, patternRev);
	fprintf(fp,"letter-probability matrix: alength= 4 w= 100 nsites= %d\n", C[ii].freq);
	for(int i1=0; i1<C[ii].len; i1++){
		fprintf(fp,"%.3f",C[ii].m[i1*4]);
		for(int i2=1; i2<4; i2++){
			fprintf(fp,"\t%.3f",C[ii].m[i1*4+i2]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}
fclose(fp);

free(cluster_assignment_cur);
free(cluster_assignment_prev);
free(N_members);
free(centroid);
free(weight);
free(off_min);
free(reserved);
return;
}		   

/***********************************************/
int main (int argc, char **argv) 
/***********************************************/
{
POINT *Point=NULL;
POINT *C;
PARAM *p;
int nPoints, i, j, k, i1, *cluster_assignment;
float m[4];
time_t t;
int *offset;
int *direction;
float *dist;

srand((unsigned) time(&t));
p = read_parameters(argc, argv);
k = p->k;
Point = read_file_meme(p->inputFile,&nPoints,p);
if(nPoints<2) error_message("Input file is empty");
check(cluster_assignment = (int*)calloc(nPoints, sizeof(int)));
check(offset = (int *)calloc(nPoints*k, sizeof(int)));
check(direction = (int *)calloc(nPoints*k, sizeof(int)));
check(dist = (float *)calloc(nPoints*k, sizeof(float)));

// Initialize centroids
check(C = (POINT*)calloc(k, sizeof(POINT)));
for(i=0; i<k; i++){
	check(C[i].m = (float*)calloc(MAXLEN*4, sizeof(float)));
}
kmeans(Point, nPoints, k, C, cluster_assignment, offset, direction, dist, p);
printf("Done\n");
return(0);
}


