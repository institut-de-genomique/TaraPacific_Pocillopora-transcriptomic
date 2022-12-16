# RNAseq reads mapping
This directory contains the tools needed to calculate the number of RNAseq reads aligned on each gene of the *Pocillopora* coral and its *Cladocopium* symbiont. 
The command lines are below  

RNAseq reads from 102 *Pocillopora* coral colonies are available at the ebi under the project number : PRJEB52301  

The complete table of gene expression for *Pocillopora* and *Cladocopium* are available at https://doi.org/10.5281/zenodo.6341761  

## 1. *Cladocopium* and *Durusdinium*

The *Cladocopium* genome is available at http://symbs.reefgenomics.org/download/ (https://doi.org/10.1038/s42003-018-0098-3)  
The *Durusdinium* transcriptome is available at  https://doi.org/10.1186/1471-2148-12-217  

### A. Reference index
```bash
cat Cladocopium-CDS.fa Durusdinium-transcripts.fa > Symbiodiniaceae-transcriptome.fa
bwa index Symbiodiniaceae-transcriptome.fa
```
### B. Read alignement with BWA-mem v0.7.15 and filtering
For each metatranscriptomic readset of *Pocillopora* colony  
```bash
bwa mem -t 6 -M Symbiodiniaceae-transcriptome.fa READ1.fq READ2.fq | samtools view -b -@ 6 -F 4 /dev/stdin -o SAMPLE_Symbiodiniaceae.aln.bam
#BamFilters tool is available here : https://github.com/institut-de-genomique/bamFilters
bamFilters -i 98 -a 80 -r 75 -n 30 -b $a -o SAMPLE_Symbiodiniaceae.98i-80l.bam
samtools sort -@ 1 -o SAMPLE_Symbiodiniaceae.98i-80l.sort.bam SAMPLE_Symbiodiniaceae.98i-80l.bam
samtools index SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.sort.bam
```
### C. Calculation of read count per gene for *Cladocopium*
```bash
sam-carac -b SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.sort.bam -s allbest -a -p -o SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac
sort -k 2,2 SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac > SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac.sort
python Read-count-per-gene.py SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac.sort Transcripts-length.tab SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.readcount.tab
```

## 2. *Pocillopora* host

The *Pocillopora* genome is available at XXX

### A. Reference index
```bash
bwa index Pocillopora_meandrina_v3.fa
```
### B. Read alignement with BWA-mem v0.7.15 and filtering
For each metatranscriptomic readset of *Pocillopora* colony  
```bash
bwa mem -t 12 -M Pocillopora_meandrina_v3.fa READ1.fq READ2.fq | samtools view -b -@ 12 -F 4 /dev/stdin -o SAMPLE_Pocillopora_meandrina_v3.aln.bam
samtools sort -@ 1 -o SAMPLE_Pocillopora_meandrina_v3.aln.sort.bam SAMPLE_Pocillopora_meandrina_v3.aln.bam
#BamFilters tool is available here : https://github.com/institut-de-genomique/bamFilters
bamFilters -i 95 -a 50 -r 75 -n 30 -b $a -o SAMPLE_Pocillopora_meandrina_v3.aln.sort.95i-50l.bam
samtools index SAMPLE_Pocillopora_meandrina_v3.aln.sort.95i-50l.bam
```
