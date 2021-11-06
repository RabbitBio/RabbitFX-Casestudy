
# RabbitFX-fastp
RabbitFX-fastp is an enhanced version of [fastp](https://github.com/OpenGene/fastp) based on [RabbitFX](https://github.com/RabbitBio/RabbitFX).

# Build

we only provide the linux version of RabbitFX-fastp

```bash
cd RabbitFX-fastp && make
```
# Simple usage
* For single end data 
```
RabbitFX-fastp -w nthreads -i in.fq -o out.fq
```
* For paired end data
```
RabbitFX-fastp -w nthreads -i in.R1.fq -I in.R2.fq -o out.R1.fq -O out.R2.fq
```

# Options
For more help information, please refer to `RabbitFX-fastp -h`.

If `-w` opition is not specified, RabbitFX-fastp will set working thread number to total CPU cores - 2.

# Examples of report
`RabbitFX-fastp` creates reports in both HTML and JSON format.
