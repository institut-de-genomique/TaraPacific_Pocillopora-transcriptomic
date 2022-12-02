# RNAseq reads mapping
This file contains the tools and commands needed to calculate the number of RNAseq reads aligned on each gene of the *Pocillopora* coral and its *Cladocopium* symbiont.  

RNAseq reads from 102 *Pocillopora coral colonies* are available at the ebi under the project number : PRJEB52301

## 1. *Cladocopium* goreaui *Durusdinium*

The *Cladocopium* genome is available at http://symbs.reefgenomics.org/download/ (https://doi.org/10.1038/s42003-018-0098-3)
The *Durusdinium* transcriptome is available at  https://doi.org/10.1186/1471-2148-12-217

### A. Read alignement with BWA-mem v0.7.15
```bash
#index
cat Cladocopium-CDS.fa Durusdinium-transcripts.fa > Symbiodiniaceae-transcriptome.fa
bwa index Symbiodiniaceae-transcriptome.fa
#Alignement for each *Pocillopora* coral colony
bwa mem -t 6 -M Symbiodiniaceae-transcriptome.fa $b $c | samtools view -b -@ 6 -F 4 /dev/stdin -o SAMPLE_Symbiodiniaceae.aln.bam;done
samtools sort -@ 8 -o SAMPLE_Symbiodiniaceae.sort.bam SAMPLE_Symbiodiniaceae.aln.bam
```

### B. Read filtration and index
```bash
ls *_Symbiodiniaceae.aln.sort.bam | while read a ;do bamFilters -i 98 -a 80 -r 75 -n 30 -b $a -o ${a%.bam}.98i-80l.bam; done
ls *_Symbiodiniaceae.aln.sort.98i-80l.bam | while read a ; do "samtools index $a;done
```
### C. Calculation of read count per gene
```bash
ls *_Symbiodiniaceae.aln.sort.98i-80l.bam | while read a ;do sam-carac -b $a -s allbest -a -p -o ${a%.bam}.carac; done
ls *_Symbiodiniaceae.aln.sort.98i-80l.carac_perHit.stat | while read a; do sort -k 2,2 $a > $a.sort; done

ls *_Symbiodiniaceae.aln.sort.98i-80l.carac_perHit.stat.sort | while read a; do python program-jlehoang-Statistics-Alignment-Results-MetaT-Samples.py $a Size-Cladocopium-goreaui-C1-CDS-and-Durusdinium-D2-Ahyac-Transcriptome-Reference-Transcripts.txt ${a%_Cladocopium-C1-CDS-and-Durusdinium-D2-Transcriptome-DualT-Tara-DB.aln.sort.rg.98i-80l.carac_perHit.stat.sort}_Statistics-Alignment-Results-DualT-Samples-POC-C1-CDS-D2-transcriptome-98i-80l-Without-Rmdup.txt; done

## 2. *Pocillopora* host
