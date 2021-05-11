# RabbitIO-Ktrim
RabbitIO-Ktrim is an enhanced version of [Ktrim](https://github.com/hellosunking/Ktrim) based on [RabbitIO](https://github.com/RabbitBio/RabbitIO).

## Installation

you can re-compile the programs:

```bash

git https://github.com/RabbitBio/RabbitIO-Casestudy.git
cd RabbitIO-Casestudy/RabbitIO-Ktrim
mkdir build
cd build
cmake ..
make 
```

## Usage

```bash
Usage: RabbitIO-Ktrim [options] -1/-U Read1.fq [ -2 Read2.fq ] -o out.prefix

RabbitIO-Ktrim is an enhanced version of Ktrim based on RabbitIO.

Compulsory parameters:

  -1/-U Read1.fq  Specify the path to the files containing read 1
                  If your data is Paired-end, use '-1' and specify read 2 files using '-2' option
                  Note that if '-U' is used, specification of '-2' is invalid
                  If you have multiple files for your sample, use ',' to separate them
                  Gzipped files are supported from version 1.2.0

  -o out.prefix   Specify the prefix of the output files
                  Note that output files include trimmed reads in FASTQ format and statistics

Optional parameters:

  -2 Read2.fq     Specify the path to the file containing read 2
                  Use this parameter if your data is generated in paired-end mode
                  If you have multiple files for your sample, use ',' to separate them
                  and make sure that all the files are well paired in '-1' and '-2' options

  -t threads      Specify how many threads should be used (default: 1, single-thread)
                  You can set '-t' to 0 to use all threads (automatically detected)

  -p phred-base   Specify the baseline of the phred score (default: 33)
  -q score        The minimum quality score to keep the cycle (default: 20)
                  Note that 20 means 1% error rate, 30 means 0.1% error rate in Phred

  -w window       Set the window size for quality check (default: 5)
                  Ktrim will stop when all the bases in a consecutive window pass the quality threshold

                  Phred 33 ('!') and Phred 64 ('@') are the most widely used scoring system
                  Quality scores start from 35 ('#') in the FASTQ files is also common

  -s size         Minimum read size to be kept after trimming (default: 36; must be larger than 10)

  -k kit          Specify the sequencing kit to use built-in adapters
                  Currently supports 'Illumina' (default), 'Nextera', 'Transposase' and 'BGI'
  -a sequence     Specify the adapter sequence in read 1
  -b sequence     Specify the adapter sequence in read 2
                  If '-a' is set while '-b' is not, I will assume that read 1 and 2 use same adapter
                  Note that '-k' option has a higher priority (when set, '-a'/'-b' will be ignored)

  -m proportion   Set the proportion of mismatches allowed during index and sequence comparison
                  Default: 0.125 (i.e., 1/8 of compared base pairs)

  -h/--help       Show this help information and quit
  -v/--version    Show the software version and quit

Please refer to README.md file for more information (e.g., setting adapters).

RabbitIO-Ktrim: enhanced extra-fast and accurate adapter- and quality-trimmer based on Ktrim.
```

`Ktrim` contains built-in adapter sequences used by Illumina TruSeq kits, Nextera kits, Nextera transposase
adapters and BGI sequencing kits within the package. However, customized adapter sequences are also allowed
by setting '-a' (for read 1) and '-b' (for read 2; if it is the same as read 1, you can left it blank)
options. You may need to refer to the manual of your library preparation kit for the adapter sequences.
Note that in the current version of `Ktrim`, only 1 pair of adapters is allowed.

Here are the built-in adapter sequences (the copyright should belong to the corresponding companies):

```
Illumina TruSeq kits:
AGATCGGAAGAGC (for both read 1 and read 2)

Nextera kits (suitable for ATAC-seq data):
CTGTCTCTTATACACATCT (for both read 1 and read 2)

BGI adapters:
Read 1: AAGTCGGAGGCCAAGCGGTC
Read 2: AAGTCGGATCGTAGCCATGT
```

### Example 1

Your data is generated using Illumina TruSeq kit in Single-end mode, then you can run:

```bash
./RabbitIO-Ktrim -U /path/to/read1.fq -o /path/to/output/dir
```

### Example 2

Your data is generated using a kit with customized adapters; your data is composed of 3 lanes in Paired-end
mode and uses Phred scoring system starts from 35; you want to keep the high quality (Phred score >=30)
bases and reads longer than 50 bp after trimming; and you want to use 4 threads to speed-up the analysis,
then you can run:

```bash
./RabbitIO-Ktrim -1 /path/to/lane1.read1.fq.gz,/path/to/lane2.read1.fq.gz,/path/to/lane3.read1.fq \
                  -2 /path/to/lane1.read2.fq.gz,/path/to/lane2.read2.fq.gz,/path/to/lane3.read2.fq \
                  -t 4 -p 35 -q 30 -s 50 -o /path/to/output/dir \
                  -a READ1_ADAPTER_SEQUENCE -b READ2_ADAPTER_SEQUENCE
```

## Testing dataset

Under the `testing_dataset/` directory, a script named `simu.reads.pl` is provided to generate *in silico*
reads for testing purpose only. **Note that the results in the paper is based on the data generated by this
script.** Another script `check.accuracy.pl` is designed to evaluate the accuracies of the trimming tools.



## Outputs explanation

`Ktrim` outputs the trimmed reads in FASTQ format and key statistics (e.g., the numbers of reads that
contains adapters and the number of reads in the trimmed files).


## Notice

This is a case study of our RabbitIO project.
