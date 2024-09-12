# HLA-Inception

A repository for predicting MHC-I peptides using HLA-Inception

# Setting up HLA-Inception

## MacOS instructions

1. install [homebrew](https://brew.sh/)

2. install git using homebrew

```shell
brew install git
brew install git-lfs
```

3. install [Go](https://go.dev/doc/install)

## Linux instructions

1. install git

``` shell
sudo apt install git-all
```

2. install [Go](https://go.dev/doc/install)


## HLA-Inception Installation

1. Clone github repository

```shell
git clone git@github.com:eawilson-CompBio/HLA-Inception.git
```

2. go to peptide prediction directory
```shell
cd HLA-Inception
```

3. Run install script

``` shell
bash install.sh
```
*Look to make sure that all tests were passed*

# Using HLA-Inception

``` shell
./hla-inception -[flags]
```

option list
``` shell
Usage of ./HLA-Inception_pred:
  -P int
        input file type (1: peptides ; 0: fasta) default: 0
  -a string
        target MHC-I allele/s for prediction (default "A_02:01")
        * if predicting multiple alleles make sure that they are seperated by commas without spaces.
        * Alleles must be formatted to 4 digit resolution and without asterisks (EX. A_02:01) 
  -i string
        name of finput file  (default "example.fasta")
        * this can be fasta files or a list of peptides
  -l string
        peptide lengths for predictions (default 9mers)
        * peptides must be 8-15mers in length
  -o string
        output file for prediction (default "output.txt")
  -threshold
        threshold for binding predictions when doing fasta predictions (default "99.5")
  -w
        Length correction weights 
```

## example predictions  


### FASTA file predictions 

Predict A_02:01 9mer peptides with a binding percentile of 99.5 or higher from the test.fasta (SARS-CoV-2 proteome)
```shell
./hla-inception -i tests/test.fasta -l 9 -a A_02:01 -threshold 99.5 -o tests/single_allele_single_length.out 
```

Predict A_02:01 9 and 10 mer peptides with a binding percentile of 99.5 or higher from the test.fasta (SARS-CoV-2 proteome)
```shell
./hla-inception -i tests/test.fasta -l 9,10 -a A_02:01 -threshold 99.5 -o tests/single_allele_multiple_lengths.out 
```

Predict A_01:01 and A_02:01 9mer peptides with a binding percentile of 99.5 or higher from the test.fasta (SARS-CoV-2 proteome)
```shell
./hla-inception -i tests/test.fasta -l 9 -a A_01:01,A_02:01 -threshold 99.5 -o tests/multiple_alleles_single_length.out 
```

Predict A_01:01 and A_02:01 9 and 10 mer peptides with a binding percentile of 99.5 or higher from the test.fasta (SARS-CoV-2 proteome)
```shell
./hla-inception -i tests/test.fasta -l 9,10 -a A_01:01,A_02:01 -threshold 99.5 -o tests/multiple_alleles_multiple_lengths.out
```

### Peptide predictions 

Predict binding affinity of peptides in peps.test to A_02:01
```shell
./hla-inception -i tests/peps.test -P 1 -a A_02:01 -o tests/peptide_prediction.out
```

*note: that peptide predictions will return the affinity of all peptides. Therefore, no threshold need be specified.*
