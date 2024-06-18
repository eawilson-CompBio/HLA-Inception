# run tests of the HLA-inception algorithm 


# cd to directory holding this script 
cd "$(dirname "$0")"

if ! command -v go &> /dev/null
then
    echo "go could not be found"
    echo "please install go and try again"
    exit
fi


#check if env variable has been set 
if [[ -z "${HI_PRED_PATH}" ]]; then
    echo "HI_PRED_PATH is not set"
    echo "Do you want to add HLA-inception directory path to .bashrc? [y/n]"
    read add_path
    if [[ $add_path == "y" || $add_path == "Y" ]]; then
        #check shell type 
        if [[ $SHELL == *"bash"* ]]; then
            echo -e "\nexport HI_PRED_PATH=`pwd`\n" >> ~/.bashrc
            echo "added to .bashrc"
        elif [[ $SHELL == *"zsh"* ]]; then
            echo -e "\nexport HI_PRED_PATH=`pwd`\n" >> ~/.zshrc
            echo "added to .zshrc"
        else 
            echo "shell type not supported"
            echo "In the future, you will need to set HI_PRED_PATH manually everytime you open a new terminal"
        fi
    else 
        echo "not added to .bashrc"
        echo "In the future, you will need to set HI_PRED_PATH manually everytime you open a new terminal"
    fi
    export HI_PRED_PATH=`pwd`
    
fi
# ask user if they want to add directory to  .bashrc 

#extract data files in data/ directory and delete files after
echo "extracting data files"
tar -xzf data/Dist.json.tgz -C data/
tar -xzf data/LO.json.tgz -C data/

# compile golang binaray 
echo "compile golang binary"
go mod init HLA-Inception
go build -o hla-inception ./src/


# run tests

# function to check if two files are the same allowing for differences in row orders
function CheckDiff () {
    RED='\033[0;31m'
    NC='\033[0m'
    GREEN='\033[0;32m'
    if [[ $(sort $1) == $(sort $2) ]]; then
        echo "${GREEN}test passed${NC}\n"
        return 1
    else
        echo "${RED}test failed${NC}\n"
        return 0
    fi
}
export -f CheckDiff
echo "running tests"

echo -e "running single allele single length test"
./hla-inception -i tests/test.fasta -l 9 -a A_02:01 -o tests/single_allele_single_length.out > /dev/null
CheckDiff tests/single_allele_single_length.out tests/single_allele_single_length.out.expected


echo -e "running single allele multiple lengths test"
./hla-inception -i tests/test.fasta -l 9,10 -a A_02:01 -o tests/single_allele_multiple_lengths.out > /dev/null
CheckDiff tests/single_allele_multiple_lengths.out tests/single_allele_multiple_lengths.out.expected


echo -e "running multiple alleles single length test"
./hla-inception -i tests/test.fasta -l 9 -a A_01:01,A_02:01 -o tests/multiple_alleles_single_length.out > /dev/null
CheckDiff tests/multiple_alleles_single_length.out tests/multiple_alleles_single_length.out.expected

echo -e "running multiple alleles multiple length test"
./hla-inception -i tests/test.fasta -l 9,10 -a A_01:01,A_02:01 -o tests/multiple_alleles_multiple_lengths.out > /dev/null
CheckDiff tests/multiple_alleles_multiple_lengths.out  tests/multiple_alleles_multiple_lengths.out.expected

echo -e "Testing allele alignment search"
./hla-inception -i tests/test.fasta -l 9 -a A_01:02 -o tests/allele_alignment_search.out > /dev/null
CheckDiff tests/allele_alignment_search.out tests/allele_alignment_search.out.expected

echo -e "checking peptide prediction test"
./hla-inception -i tests/peps.test -P 1 -a A_02:01 -o tests/peptide_prediction.out > /dev/null
CheckDiff tests/peptide_prediction.out tests/peptide_prediction.out.expected

