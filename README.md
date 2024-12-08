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
