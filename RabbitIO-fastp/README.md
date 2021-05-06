
# RabbitIO-fastp
RabbitIO-fastp is an enhanced version of [fastp](https://github.com/OpenGene/fastp) based on [RabbitIO](https://github.com/RabbitBio/RabbitIO).

# Build

we only provide the linux version of RabbitIO-fastp

```bash
cd RabbitIO-fastp && make
```
# Simple usage
* For single end data 
```
RabbitIO-fastp -w nthreads -i in.fq -o out.fq
```
* For paired end data
```
RabbitIO-fastp -w nthreads -i in.R1.fq -I in.R2.fq -o out.R1.fq -O out.R2.fq
```

# Options
For more help information, please refer to `RabbitIO-fastp -h`.

If `-w` opition is not specified, RabbitIO-fastp will set working thread number to total CPU cores - 2.

# Examples of report
`RabbitIO-fastp` creates reports in both HTML and JSON format.



<!--

## If you use RabbitIO-fastp for short read quality control please cite:

Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

## If you use RabbitIO-fastp for long read quality control please cite:

De Coster W, D’Hert S, Schultz D T, et al. NanoPack: visualizing and processing long-read sequencing data[J]. Bioinformatics, 2018, 34(15): 2666-2669.
-->
