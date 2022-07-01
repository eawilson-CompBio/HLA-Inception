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

2. go to peptide prediction directory
```shell
cd HLA-Inception
cd HLA-Inception_pred/
```

3. Build HLA-Inception peptide prediction algorithm

``` shell
go build scripts/HLA-Inception_pred.go
```

4. uncompress data files 

``` shell
cd data
tar -xzvf LO_matrices.tgz 
tar -xvzf dist_cutoffs.tgz
cd ..
```

5. set data path to files (change the path to the location of the data directory)

``` shell
export HI_PRED_DATA=path/to/data/directory
```


## Using HLA-Inception to predict peptides

``` shell
./HLA-Inception_pred -[flags] 
```

option list
``` shell
Usage of ./HLA-Inception_pred:
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

* some of the functionality is still be updated. Only use this for prediction of 9mers for the time being.

## Testing Installation

``` shell
cd test
../HLA-Inception_pred -a A_02:01 -i SARS_CoV_2.fasta -o test.out

```
