# Assembly of Bacterial genome with PacBio and Illumina data

This is a pipeline to assemble a Bacterial Genome using PacBio and Illumina data. Here we describe the steps used to assemble the genome of [_Bacillus_ sp. RZ2MS9](https://www.ncbi.nlm.nih.gov/nuccore/CP049978), so you will see examples using this genome data thoughout this pipeline.

If you find this pipeline usefull, please cite [this article](XXX).


## Use [bamutil/1.0.14](https://github.com/statgen/bamUtil) to run _stats_ on _.bam_ file generated from PacBio reads.

`bam stats --in /path/to/PACBIO_FILE.bam --basic`


## Analyzing ZMW
_Questions: How many ZMW are in the .bam file?; How many subreads from the same ZMW?_

>ZMW Counts
First, use [Samtools/v1.9](https://github.com/samtools/samtools):

`samtools view /path/to/PACBIO_FILE.bam | awk -F '\t' '{print $1}' | awk -F '/' '{print $2}' | sort | uniq -c > ZMW_Counts.txt`

>This part was made with [R](https://www.google.com/search?client=safari&rls=en&q=r+software&ie=UTF-8&oe=UTF-8): 

`tb <- read.table("/paht/to/ZMW_Counts.txt",` `header=F,as.is=F)`

`sum(tb$V1)`

`dim(tb)`

`summary(tb$V1)`

`q ()` #to exit R in the terminal

**Example from _Bacillus_ sp. RZ2MS9 data:**

sum(tb$V1)

[1] 66292

dim(tb)

[1] 8274  2

summary(tb$V1)

Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  
  1.000   2.000   7.000   8.012  12.000  97.000
  
  
**Final Results**

|    Run.   |       Number of ZMW    | Average of subreads per ZMW |
|-----------|:-----------------------:|:------------------:|
|PACBIO_FILE |2 – 8274 Min: 1 Max: 97  |Average: 8.012086   |



## Load [smrtlink/8.0.0.80529](https://www.pacb.com/wp-content/uploads/SMRT-Link-User-Guide-v8.0.pdf) to run `ccs 4.0.0`

The [circular concensus sequencing](https://github.com/PacificBiosciences/ccs) or ccs considers the subreads of ZMW to correct PacBio reads. The PacBio reads must be [HIFI](https://www.pacb.com/smrt-science/smrt-sequencing/smrt-sequencing-modes/).

```
nohup ccs --reportFile=PACBIO_FILE.report \
/path/to/PACBIO_FILE.bam \
CCS.bam &
```


## Assembling the genome with [Canu v 1.8](https://canu.readthedocs.io/en/latest/index.html)

>"Canu specializes in assembling PacBio or Oxford Nanopore sequences. Canu operates in three phases: correction, trimming and assembly. The correction phase will improve the accuracy of bases in reads. The trimming phase will trim reads to the portion that appears to be high-quality sequence, removing suspicious regions such as remaining SMRTbell adapter. The assembly phase will order the reads into contigs, generate consensus sequences and create graphs of alternate paths.

First, convert the _.bam_ files into _.fq_ files. To that, I will load module `samtools/1.9`

`samtools fastq /path/to/CCS.bam > CCS.fq`

```
nohup canu -p canu_assemble -d ../Canu genomeSize=5.5m -pacbio-corrected path/to/CCS.fq &
```

**Note:** Please change the genome size to the expected genome size of your genome.

**You can use the [Bandage](https://rrwick.github.io/Bandage/) to visualize your Genome with the file canu_assemble.gfa.**


## Using `arrow/2.3.3` to correct the genome assemble [more information](https://github.com/PacificBiosciences/GenomicConsensus):

First, you will use [`pbalign/0.4.1`](https://github.com/PacificBiosciences/pbalign) to aligned PacBio raw data to my contig:

```
pbalign /path/to/PACBIO_FILE.fofn \
          /path/to/canu_assemble.fasta \
          canu_assemble-aligned-raw-data.bam 
```

**Note:** Information about a _.fofn_ file can be find [here](https://pb-falcon.readthedocs.io/en/latest/tutorial.html#create-fofn).

`samtools flagstat canu_assemble-aligned-raw-data.bam`

Then, we need to index the fasta file:

`samtools faidx /path/to/canu_assemble.fasta`

Running variantCaller:

```
nohup variantCaller /path/to/canu_assemble-aligned-raw-data.bam --algorithm=arrow \
                    -r /path/to/canu_assemble.fasta \
                    -o variants-raw.gff \
                    -o canu_assemble-arrowcorr.fasta &
```

## Correting Illumina reads with `htstream/1.1.0`

This pipeline is from [HTStream RNA preprocessing] (https://github.com/ucdavis-bioinformatics-training/HTStream_Training/blob/master/preproc.md) from Bioinformatics Core from UC Davis, with modifications.

```
hts_Stats -L GenomeHTSstats.log -1 /path/to/Illumina/Genome_S1_L001_R1_001.fastq.gz -2 /path/to/Illumina/Genome_S1_L001_R2_001.fastq.gz | \
hts_SeqScreener -A -L GenomeHTSstats.log | \
hts_SuperDeduper -A -L GenomeHTSstats.log | \
hts_Overlapper -n -A -L GenomeHTSstats.log | \
hts_QWindowTrim -n -A -L GenomeHTSstats.log | \
hts_NTrimmer -n -A -L GenomeHTSstats.log | \
hts_CutTrim -n -m 200 -A -L GenomeHTSstats.log | \
hts_Stats -A -L GenomeHTSstats.log -f Genome.htstream
```

## Running `pilon` on _arrow_ corrected data

First, you will need to align my Preprocessed Illumina reads against my corrected-contig. To that, you will need to create an index file to the _.fasta_ document.

`bwa index canu_assemble-arrowcorr.fasta`

To align:

`
bwa mem canu_assemble-arrowcorr.fasta /path/to/HTStream/Genome.htstream_SE.fastq  > canu_assemble-arrowcorr-illumina-alig-single.sam
`

`
bwa mem canu_assemble-arrowcorr.fasta /path/to/HTStream/Genome.htstream_R1.fastq /path/to/HTStream/Genome.htstream_R2.fastq > canu_assemble-arrowcorr-illumina-alig-paired.sam
`

_obs: do not run this command with nohup, it will generate a bad **.sam** file._


You will need to convert the _.sam_ file to its binary counterpart, _.bam_ format.

`samtools view -S -b canu_assemble-arrowcorr-illumina-alig-single.sam > canu_assemble-arrowcorr-illumina-alig-single.bam`


`samtools view -S -b canu_assemble-arrowcorr-illumina-alig-paired.sam > canu_assemble-arrowcorr-illumina-alig-paired.bam`


To report some stats of our alignment:


`samtools flagstat canu_assemble-arrowcorr-illumina-alig-single.bam`

You will need then to sort the files:

`
samtools sort canu_assemble-arrowcorr-illumina-alig-single.bam -o canu_assemble-arrowcorr-illumina-alig-single-sorted.bam
`

`
samtools sort canu_assemble-arrowcorr-illumina-alig-paired.bam -o canu_assemble-arrowcorr-illumina-alig-paired-sorted.bam
`

Then you will need to index these files.

`
samtools index canu_assemble-arrowcorr-illumina-alig-single-sorted.bam
`

`
samtools index canu_assemble-arrowcorr-illumina-alig-paired-sorted.bam
`

### Runing `pilon/1.23` on arrow previous corrected data.

Here you can use [pilon](https://github.com/broadinstitute/pilon/wiki) to polish the contig.

```
nohup pilon --genome /path/to/canu_assemble-arrowcorr.fasta \
--frags /path/to/canu_assemble-arrowcorr-illumina-alig-paired-sorted.bam \
--unpaired /path/to/canu_assemble-arrowcorr-illumina-alig-single-sorted.bam \
--output Final_Genome --changes --vcfqe --tracks &
```


## Final adjustments in your genome:

1. Since Bacterial Genome is circular, check if there is any overlap at the beginning and end of the Genome.

2. Also, using [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi), you can verify if your Genome is on the expected DNA frame (5'--- 3').

3. **Looking for the chromosomal replication initiator protein DnaA**. This will be the first base of your contig, if you are working with **Bacteria**. This is a gene of 1341 bp.
To find where DnaA is on your genome, you can use as reference the [chromosomal replication initiator protein DnaA](https://www.ncbi.nlm.nih.gov/nuccore/CP049978.1?from=1&to=1341) of _Bacillus_ sp. RZ2MS9.


### Using [`biopython`](https://biopython.org/wiki/Documentation) to do some modificiations based in what happened in the _Bacillus_ sp. RZ2MS9 genome assemble:

1. Cut the first 10.037 bp

2. Reverse complement it

3. Make base 4,522,510 the first one, it is the origin of replication

Using `python`: call `python3` on command line:

Once inside the python **>>>**:

- to know where the home directory is
>>> import os
>>> home = os.path.expanduser('~')
>>> home

- to see if I have `pip` to install `biophyton`
(command line) python3 -m ensurepip
(command line) python3 -m pip install biopython

- Setting work_directory
>>> import os
>>> os.chdir("/Users/mariabonatelli/Desktop")
>>> print("Current Working Directory " , os.getcwd())

_Running on BioPython_

```
from Bio import SeqIO

h = SeqIO.parse('Final_Genome.fasta', 'fasta')
s = next(h)

print(h)

_Trim off the left edge_:
sTrimmed = s[10037:]
print(sTrimmed)

_Reverse complement it_
sRC = sTrimmed.reverse_complement()
print(sRC)

print(sRC[0:5])
print(sTrimmed[-5:])

sRC.description = 'set_a_description'

SeqIO.write(sRC, handle='Final_Genome_trimmed_and_RC.fasta', format='fasta')
```


### Choosing the beginning of my sequence
I will use [bedtools2/2.27.0](https://github.com/arq5x/bedtools2) to determine the beginning of my sequence. First, I will need to construct a [_bedfile_](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).

The .bed format have three mandatory columns:

_chrom	chromStart	chromEnd_
set_a_description      0       843381
set_a_description      843381       5357194

*ps: chromStart considers the first base as 0.*

_Bacillus_ sp. RZ2MS9 contig has 5,357,194 bp.

One can use [`bedtools getfasta`](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html) to get the initial and final sections according to the bed files.

`
bedtools getfasta -fi Final_Genome_trimmed_and_RC.fasta -bed position1.bed -fo Final_Genome_trimmed_and_RC_beg.fasta
`

`
bedtools getfasta -fi Final_Genome_trimmed_and_RC.fasta -bed position2.bed -fo Final_Genome_trimmed_and_RC_end.fasta
`

Checking the sizes of the sequences -->

```
grep -v ">" Final_Genome_trimmed_and_RC_beg.fasta | wc | awk '{print $3-$1}'
843381

grep -v ">" Final_Genome_trimmed_and_RC_end.fasta | wc | awk '{print $3-$1}'
4513813

843381 + 4513813 = 5357194

grep -v ">" Final_Genome_trimmed_and_RC.fasta | wc | awk '{print $3-$1}'
5357194

You can manually join them.

grep -v ">" Final_Genome_trimmed_and_RC_total.fasta | wc | awk '{print $3-$1}'
5357194
```
