# CisFinder
CisFinder software is installed at https://kolab.elixirgensci.com/cisfinder/.
CiFinder is a tool for for finding over-representing short DNA motifs (e.g., transcription factor binding motifs) and estimates the Position Frequency Matrix (PFM) directly from word counts. It is designed for the analysis of ChIP-chip and ChIP-seq data and can process long input files (50 Mb). Identification and clustering of motifs takes 8 sec per 1 Mb. Other features: uses False discovery rate (FDR) and control sequences, finds motifs based on PFM.
Main components of the software are written in C, and the code is available from the installation website.
The UI is written in perl and html and the code is available in this projct.
