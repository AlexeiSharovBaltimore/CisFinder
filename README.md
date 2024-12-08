# CisFinder
CisFinder software is installed at https://kolab.elixirgensci.com/cisfinder/.
CiFinder is a tool for for finding over-representing short DNA motifs (e.g., transcription factor binding
motifs) and estimates the Position Frequency Matrix (PFM) directly from word counts. It is designed for 
the analysis of ChIP-chip and ChIP-seq data and can process long input files (50 Mb). Identification 
and clustering of motifs takes 8 sec per 1 Mb. Other features: uses False discovery rate (FDR) and control 
sequences, finds motifs based on PFM.

Main components of the software are written in C, and the code is available from the installation website.
The UI is written in perl and html and the code is available in this projct.

DIRECTORY STRUSTURE

CisFinder(root) = cisfinder.ini configure file

CisFinder/cisfinder = html files

CisFinder/cisfinder/bin = cisfinder.cgi (CGI program)

CisFinder/cisfinder/images = images

CisFinder/cisfinder/output = writable directory for output files

CisFinder/cisfinder/download

CisFinder/cisfinderInfo/bin = other programs (C, perl), C-programs should be compiled here

CisFinder/cisfinder/info = writable directory for personal config files, login.txt

CisFinder/cisfinder/data = writable directory for data files

Compile C-programs in CisFinder/cisfinderInfo/bin:

gcc patternFind.c -lm -o patternFind

gcc patternCluster.c -lm -o patternCluster

gcc patternCompare.c -lm -o patternCompare

gcc patternDistrib.c -lm -o patternDistrib

gcc patternScan.c -lm -o patternScan

gcc patternTest.c -lm -o patternTest

gcc patternVar.c -lm -o patternVar

Configure two accounts to your name: administrator, public. Edit file login.txt and put your name and email address
Log into cisfinder and change your passwords for both account. Create your own working account in CisFinder

2. Program "patternFind"

Function: generates position frequency matrixes (PFM) for motifs over-represented in the
test DNA sequence compared to control sequence

Syntax:
patternFind -i inputFasta -o outputFile [-c controlFile, -pos positionFile,
-ratio minEnrichmentRatio, -FDR maxFDR, -maxlen maxSequenceLength, -len motifLength,
-strand strandOption, -score scoreOption, -userep, -getrep, -one, -cg, -brief,
-n numberOfMotifs]

Comments:
(a) If controlFile is not specified, then random sequence generated using 3rd order
Markov chain is used as a control.
(b) positionFile - stores position of each motif match in the test sequence file.
It can be later used for motif clustering on the basis of their co-occurrence.
(c) maxFDR = maximal False Discovery Rate (FDR) threshold. The program generates at least 
100 motifs even if they are not significant; additional motifs are included only 
if they are significant (i.e. FDR < FDR threshold)
(d) motifLength = 8 as a default. Possible values are 6, 8, and 10 only.
(e) minEnrichmentRatio = minimum enrichment ratio, default = 1.5
(f) strandOption: 0 = search both strands, 1 = search positive strand, 2 = search
positive strand and use negative strand as a control.

(g) scoreOption:

0 = use z-score (z) for motif over-representation for ordering motifs in the output file

1 = use z*(ratio-1) for sorting motifs, where ratio = over-representation ratio.

2 = use z*info for sorting motifs, where info = information content.

3 = use z*(ratio-1)*info

4 = use z*(1-selfsim) for sorting motifs, where selfsim = self-similarity of motif.

5 = use z*(ratio-1)*(1-selfsim)

6 = use z*info*(1-selfsim)

7 = use z*(ratio-1)*info*(1-selfsim)

(e) userep: use repeats in sequence (lower-case in sequence)

(f) getrep: generate repeat output file (motifs over-represented in repeats).

(g) one: consider not more than 1 motif occurrence per sequence

(h) cg: adjust motif abundance to C/G and CpG accurrence in test and control sequences

(i) brief: do not generate PFM.

(j) numberOfMotifs = maximum number of motifs to be generated (default = 500)

3. Program "patternCluster"

Function: clusters motifs with position frequency matrixes (PFM) 
based on PFM similarity and/or motif co-occurrence in the test sequence.
Similarity is measured by correlation of position-weight matrices (PWM) which
are log-transformed PFMs.

Syntax:
patternCluster -i inputMotifFile -o outputFile [-pos positionFile, 
-match matchThreshold, -n numberOfMotifs, -repeat maxRepeatEnrichment, -posonly]

Comments:
(a) positionFile - stores position of each motif match in the test sequence file.
It is used for motif clustering on the basis of their co-occurrence.
If positionFile is not specified, then clustering is done exclusively on the basis
of similarity of motifs, otherwise co-occurrence method is used at least for
clustering of motifs with high level of self-similarity (>0.5 on average). However
if "posonly" option is used, then clustering is based exclusively on co-occurrence
of motifs.
(b) posonly = motifs are clustered exclusively based on co-occurrence.
(c) matchThreshold = minimum similarity (correlation of PWMs) between motifs;
default value = 0.75.
(d) numberOfMotifs can be used to limit the number of input motifs
(e) maxRepeatEnrichment (ratio threshold) can be used to filter out motifs with 
high over-representation in repeats.

4. Program "patternCompare"

Function: compare motifs between 2 files based on PFM similarity.
Similarity is measured by correlation of position-weight matrices (PWM) which
are log-transformed PFMs.

Syntax:
patternCompare -i1 motifFile1.txt -i2 motifFile2 -o outputFile [-match matchThreshold]

Comments:
(a) matchThreshold = minimum similarity (correlation of PWMs) between motifs;
default value = 0.75.

5.  Program "patternTest"

Function: improve motifs using resampling method

Syntax:
patternTest -i motifFile -f fastaFile -o outputFile [-c controlFile, 
-prog progressFile, -maxlen maxSequenceLength, -strand strandOption, 
-n numberOfMotifs, -method resampleMethod, -iter numIteartions, 
-siter numSubiterations, -fp falsePositives, -userep, -one]

Comments:
(a) If controlFile is not specified, then random sequence generated using 3rd order
Markov chain is used as a control.
(b) progressFile keeps records for each iteration and subiteration.
(c) strandOption: 0 = search both strands, 1 = search positive strand, 2 = search
positive strand and use negative strand as a control.
(d) numberOfMotifs can be used to limit the number of input motifs (default=100)
(e) resampleMethod: 1=regression method (default), 2=simple resampling, 3=difference method.
(f) numIteartions = number of main iterations when a new set of motif matches 
is generated (default = 3)
(g) numSubiterations = number of sub-iterations which process the same set of motif
matches, and only the match score changes (default = 3).
(e) userep: use repeats in sequence (lower-case in sequence)
(g) one: consider not more than 1 motif occurrence per sequence (with highest 
match score)
(h) falsePositives = number of expected false positives per 10000 bp in the
control sequence (default = 5).

6. Program "patternScan"

Function: finds sites in a sequence that match to motifs specified by PFM

Syntax:
patternScan -i motifFile -f fastaFile -o outputFile [-cons conservationFile, 
-maxlen maxSequenceLength, -strand strandOption, -n numberOfMotifs, 
-fp falsePositives, -thresh matchThreshold, -userep, -one]

Comments:
(a) conservation file = file with evolutionary conservation scores (see format below)
(b) strandOption: 0 = search both strands, 1 = search foreward strand.
(c) numberOfMotifs can be used to limit the number of input motifs (default=100)
(d) falsePositives = number of expected false positives per 10000 bp in a random
sequence (default = 5 if matrix-specific thresholds are not supplied).
(e) matchThreshold = addition to matrix-specific thresholds. For example if
it is equal to 0.7, then all match thresholds are incremented by 0.7 and
search becomes more stringent.
(f) userep: use repeats in sequence (lower-case in sequence)
(g) one: consider not more than 1 motif occurrence per sequence (with highest 
match score)
(h) Matrix-specific match thresholds are automatically generated by the
patternTest program, but they can be modified, added, or removed manually.
The header for the threshold field should be "Threshold" (see file formats).

7. Program "patternDistrib"

Function: finds sites in a sequence that match to motifs specified by PFM

Syntax:
patternDistrib -i scanResults -f frequencyOutput -a abundanceOutput [-int intervalForFreq]

Comments:
(a) scanResults = file generated by "patternScan" program
(b) frequencyOutput = shows the frequency distribution of binding sites along 
sequences aligned by their starting position.
(c) abundanceOutput = shows the number of motif matches in each sequence. It can be used
for classifying sequences based on the composition of binding motifs.
(d) intervalForFreq = interval (bp) used for calculation of frequencyOutput 
file (default = 100bp).

8. File formats

All files used in CisFinder are in plain text format and use tabs for delimiting
fields within one line. Total there are 4 main types of files in CisFinder: sequences, motifs,
search results, and repeats. Motifs can be submited as PFMs and as degenerate
consensus sequences, which are converted to the PFM form. These 4 types of files
may have 2 optional lines of annotation that start with a keyword "Parameters:"
and "Headers:". Parameters may specify the origin of the file, e.g., algorithm
parameters for derivation, whereas headers specify the columns of tab-delimited lines.
Additional 3 kinds of files can be associated with a sequence file:
genomic coordinates, attributes (e.g., gene names), and conservation scores.<p>

Sequence file has a fasta format: sequence name is preceded by ">" sign
and the following lines contain the sequence. Example:
>PET070528
AACCCAAAGTATGATATGCTATGATAGATAACCAAAAGGTAATATTATGAAATTTTTATCAACTATAATTATATAACTTG
AAACTGTTTCCTAAATCCGCCCTAGAGCTTACACAAAGCTGAGGGAAGTTTGCTGGAAAGTTCAGGCTGAGTGGGATGTT
>PET070400
TACTATTGGCGCTTCAATCAGTATTCGTCTTTTATAATACAATAATGCTATTTTGGATAAGTAAGTTTCTATTCAAGGAC
ACGTGTGGGCAGCTGTAACACTAATAATGTCCCATAAATAAGCGAGCAGAGCACATACTGCTGAGACAGACATGTAAGAA

To extract DNA sequences use RSAT (http://rsat.scmbb.ulb.ac.be/rsat>RSAT)
or our PERL script "extract_genome_seq.pl" at http://lgsun.grc.nia.nih.gov/cisfinder/download.html.
We recommend using 200-bp sequences centered at the expected binding site of the TF.
ChIP-chip usually has a lower spatial resolution, thus the size of sequences can be
increased to 300 or even 400 bp.

Motif file is formatted as follows. The first line starts with a ">" sign
followed by motif name. The same line may contain additional tab-delimited fields:

Pattern 	(=consensus),
PatternRev 	(=reverse consensus),
Threshold 	(=threshold score),
Nmembers 	(=number of member motifs in motif cluster),
Freq 	(=number of motif matches in the test sequence),
Ratio 	(=enrichment ratio of motif matches in the test sequence),
Info 	(=information content of motif PFM),
Score 	(=motif score used for ordering),
FDR 	(=False Discpvery Rate),
Repeat 	(=enrichment ratio of this motif in repeat sequences),
Palindrome 	(=1 if palindrome, and =0 othewise),
Method 	(=method of motif clustering: 0 for similarity, and 1 for co-occurrence),
Species 	(=taxonomy of organism),

The PFM is formatted in 5 columns: position number (starting from 0), followed by frequency
of nucleotides A, C, G, and T respectively. There is an empty line at the end of each
matrix. Example of motif file:

Parameters:	matchThresh=0.8500	nucleotideOrder=A,C,G,T
Headers:	Name	Pattern	PatternRev	Freq	Ratio	Info	Score	p	FDR	Palindrome	Nmembers
>SOX9	CCWTTGTT	KAACAAWG	12813	4.5857	9.407	176.0950	0.0000	0.0000	0	151
0	2	72	13	13
1	3	72	2	23
2	46	1	2	51
3	1	2	1	96
4	1	1	1	96
5	2	4	93	1
6	3	2	2	93
7	6	26	13	55

>OCT	HATGCWAA	ATTWGCAT	404	3.1375	10.055	18.8350	0.0000	0.0000	0	4
0	93	3	3	2
1	1	1	2	96
2	1	1	97	1
3	3	76	10	11
4	61	3	1	36
5	80	2	17	2
6	96	1	2	1
7	3	15	14	69

If motifs are uploaded as a pattern (consensus), then the file is formatted as follows:

Parameters:	matchThresh=0.8500	nucleotideOrder=A,C,G,T
Headers:	Name	Pattern	TFgroup	Info
>MIT_001NRF1	RCGCANGCGY	NRF1	16
>MIT_002MYC	CACGTG	MYC	12
>MIT_003ELK	SCGGAAGY	ELK	14
>MIT_005NFY	GATTGGY	NFY	13
>MIT_006SP1	GGGCGGR	SP1	13
>MIT_007AP1	TGANTCA	AP1	12
>MIT_009ATF	TGAYRTCA	ATF	14
>MIT_010YY1	GCCATNTTG	YY1	16

Search results are formatted as a tab-delimited text file with the following
columns: MotifName, SeqName (sequence name), Strand, Len (motif length), Start (starting position counted from 0),
Score (matching score), Sequence (matching sequence), Conservation (evolutionary conservation score, from 0 to 100).
The "Parameters:" line should specify the name of the sequence (file_fasta) file that was searched
and the motif name that contain PFM (file_motif). Here is a example of search results:

Parameters:	file_motif=public-motif_pluripotent	file_fasta=public-P300_binding.fa
Headers:	MotifName	SeqName	Strand	Len	Start	Score	Sequence	Conservation
STAT3	Chen2008P300000002-900-1100	-	9	119	4.3932	TTCCCGGAA	40
TEF	Chen2008P300000005-900-1100	-	8	87	3.3398	AGGAATGC	0
NANOG	Chen2008P300000004-900-1100	+	9	116	3.5202	CCACTTCCT	1
KLF	Chen2008P300000004-900-1100	+	8	96	3.7154	CTCCACCC	1
KLF4	Chen2008P300000004-900-1100	+	9	109	4.0565	GCCACACCC	1
SOX9	Chen2008P300000003-900-1100	-	8	60	3.7513	CCATTGTT	28

Repeat file is formatted as a list of repeat motifs, each described with 2 lines.
The first line has the following tab-delimited fields:
Index for gap structure (shown at the left in Fig. 5)
Index for the 8-mer word (nucleotides A,C,G,T are encoded as 0,1,2,3; then the index is
calculated as c1+c2*4+c3*16+c4*64+..., where ci is the code of position i)
Enrichment ratio in comparison to non-repeats
Motif total length
Repeat name (optional)

The second line shows the full PFM (all lines concatenated), space- or tab-delimited.
Motif file can be generated by the "patternFind" program.

Sequence conservation file has a fasta format; sequence name is preceded by ">" 
sign and the following lines contain conservation scores for the sequence, comma separated.
Conservation score ranges from 0 to 100 and can be downloaded from UCSC web site (there it
ranges from 0 to 1, hence it needs to be multiplied by 100).
The first line in the file should specify parameters: interval, and genome. Normally we use
interval=5, which means that each conservation score correspond to 5 bp of the sequence. If
interval is not described it is assumed to =1. Here is an example of conservation file:
Parameters:	interval=5	genome=mm9
