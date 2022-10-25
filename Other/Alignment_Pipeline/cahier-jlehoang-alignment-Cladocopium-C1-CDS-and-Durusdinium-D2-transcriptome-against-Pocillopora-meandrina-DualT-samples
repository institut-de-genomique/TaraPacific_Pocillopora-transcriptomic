# Dual-Transcriptomic alignment

We used the following pipeline to obtain read count of _Cladocopium_ and _Durusdinium_ in dual-transcriptomic samples. 

## 1. Indexation of the reference transcriptomes
bwa index Cladocopium-goreaui-C1-CDS-and-Durusdinium-D2-Ahyac-Transcriptome-with-taxonomic-affiliation-Tara-DB-Switch-Filter-90.fasta

## 2. Formated sample file as :
Sample_ID Path_Read1  Path_Read2

## 3. Alignment using BWA-MEM
cat DualT-Pocillopora-Samples-Readsets-List-Formated-Genoscope.txt | while read a b c; do bwa mem -t 6 -M Cladocopium-goreaui-C1-CDS-and-Durusdinium-D2-Ahyac-Transcriptome-with-taxonomic-affiliation-Tara-DB-Switch-Filter-90.fasta $b $c | samtools view -b -@ 6 -F 4 /dev/stdin -o ${a}_Cladocopium-C1-CDS-and-Durusdinium-D2-Transcriptome-DualT-Tara-DB.aln.bam;done

## 4. Samtools sort, to order reads
ls *.aln.bam | while read a ;do samtools sort -@ 8 -o ${a%.bam}.sort.bam $a --reference Cladocopium-goreaui-C1-CDS-and-Durusdinium-D2-Ahyac-Transcriptome-with-taxonomic-affiliation-Tara-DB-Switch-Filter-90.fasta; done

## 5. BamFilters, to filter reads
ls *.aln.sort.bam | while read a ;do bamFilters -i 98 -a 80 -r 75 -n 30 -b $a -o ${a%.bam}.98i-80l.bam; done

## 6. Bam index, to index file results
ls *.aln.sort.*i-*l.bam | while read a ; do "samtools index $a;done

## 7. Sam-carac, to have statistic view of reads
ls *.sort.*i-80l.bam | while read a ;do sam-carac -b $a -s allbest -a -p -o ${a%.bam}.carac; done

## 8. Statistic results: for each transcript the average identity (%); average coverage (%); count reads aligned
** STEP 1
ls *.aln.sort.rg.98i-80l.carac_perHit.stat | while read a; do sort -k 2,2 $a > $a.sort; done

** STEP 2
ls *_Cladocopium-C1-CDS-and-Durusdinium-D2-Transcriptome-DualT-Tara-DB.aln.sort.rg.98i-80l.carac_perHit.stat.sort | while read a; do python program-jlehoang-Statistics-Alignment-Results-MetaT-Samples.py $a Size-Cladocopium-goreaui-C1-CDS-and-Durusdinium-D2-Ahyac-Transcriptome-Reference-Transcripts.txt ${a%_Cladocopium-C1-CDS-and-Durusdinium-D2-Transcriptome-DualT-Tara-DB.aln.sort.rg.98i-80l.carac_perHit.stat.sort}_Statistics-Alignment-Results-DualT-Samples-POC-C1-CDS-D2-transcriptome-98i-80l-Without-Rmdup.txt; done
