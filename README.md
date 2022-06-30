# HLA-Inception
A repository for predicting MHC-I peptides using HLA-Inception


# Setting up HLA-Inception

## MacOS instructions

1. install [homebrew](https://brew.sh/)

2. install git using homebrew

```shell
brew install git
		
```
   
3. install [Go](https://go.dev/doc/install)

## Linux instructions

1. install git

``` shell
sudo apt install git-all
```

2. install [Go](https://go.dev/doc/install)


## Install HLA-Inception 

1. Clone github repository

```shell
git clone git@github.com:eawilson-CompBio/HLA-Inception.git
```

2. Build HLA inception

``` shell
cd HLA-Inception
go build scripts/
```

## Prediction from protein fasta

``` shell
./HLA-Inception_pred -t -p 
```

option list
``` shell
Usage of ./predict_peptidesHLA-Inception_pred:
  -P int
        input file type (1: peptides ; 0: fasta) default: 0
  -a string
        target MHC-I allele for prediction (default "A_02:01")
  -i string
        name of finput file  (default "example.fasta")
  -max int
        maximum peptide length (default 9)
  -min int
        minimum peptide length (default 9)
  -neg int
        keep peptides with negative score? (1: yes ; 0: no) default: 0
  -o string
        output file for prediction (default "output.txt")
  -w int
        use length weights. (1: yes ; 0: no). requires -w flag to be set. default: 0
  -weights_file string
        length correction weights (default "length_weights.txt")
```

