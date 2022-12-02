# RNAseq reads mapping
This file contains the tools and commands needed to calculate the number of RNAseq reads aligned on each gene of the *Pocillopora* coral and its *Cladocopium* symbiont.  

RNAseq reads from 102 *Pocillopora coral colonies* are available at the ebi under the project number : PRJEB52301

## 1. *Cladocopium* goreaui *Durusdinium*

The *Cladocopium* genome is available at http://symbs.reefgenomics.org/download/ (https://doi.org/10.1038/s42003-018-0098-3)  
The *Durusdinium* transcriptome is available at  https://doi.org/10.1186/1471-2148-12-217  

### A. Reference index
```bash
cat Cladocopium-CDS.fa Durusdinium-transcripts.fa > Symbiodiniaceae-transcriptome.fa
bwa index Symbiodiniaceae-transcriptome.fa
```
### A. Read alignement with BWA-mem v0.7.15 and filtering
For each *Pocillopora* readset
```bash
bwa mem -t 6 -M Symbiodiniaceae-transcriptome.fa READ1.fq READ2.fa | samtools view -b -@ 6 -F 4 /dev/stdin -o SAMPLE_Symbiodiniaceae.aln.bam;done
bamFilters -i 98 -a 80 -r 75 -n 30 -b $a -o SAMPLE_Symbiodiniaceae.98i-80l.bam
samtools sort -@ 1 -o SAMPLE_Symbiodiniaceae.98i-80l.sort.bam SAMPLE_Symbiodiniaceae.98i-80l.bam
samtools index SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.sort.bam
```
### C. Calculation of read count per gene for *Cladocopium*
For each *Pocillopora* readset
```bash
sam-carac -b SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.sort.bam -s allbest -a -p -o SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac
sort -k 2,2 SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac > SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac.sort
python Read-count-per-gene.py SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.carac.sort Transcripts-length.tab SAMPLE_Symbiodiniaceae.aln.sort.98i-80l.readcount.tab
```

## 2. *Pocillopora* host
