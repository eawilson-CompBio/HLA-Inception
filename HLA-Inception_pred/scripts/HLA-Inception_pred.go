package main

import (
        "bufio"
        "fmt"
        "log"
        "os"
	"strings"
	"sync"
	"runtime"
	"flag"
	"encoding/csv"
	"strconv"
)


func readfasta(input string) ([]string , []string) {
	f, err := os.Open(input)

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()

	scanner := bufio.NewScanner(f)
	name := make([]string, 0)
	seqs := make([]string , 0)
	out := ""
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] == '>' {
			
			name = append(name,line[1:])
			if out != "" {
				seqs = append(seqs,out)
				out = ""
				
			}
			continue
		}
		out += scanner.Text()
	}

	seqs = append(seqs,out)
	
	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
	return name,seqs

}

func ReadPeptides(input string) ([]string) {
	f, err := os.Open(input)

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()

	scanner := bufio.NewScanner(f)
	seqs := make([]string , 0)
	for scanner.Scan() {
		seqs = append(seqs,scanner.Text())

	}

	
	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
	return seqs

}


func readLength_weights(input string) []float64 {
	f, err := os.Open(input)

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()

	scanner := bufio.NewScanner(f)
	weight := make([]float64 , 0)
	for scanner.Scan() {
		num,_ := strconv.ParseFloat(scanner.Text(),64)
		weight = append(weight,num)
	}

	return weight

}



func gen_pep(seq string , name string) {
	prot_len := len(seq)
	s := strings.TrimSpace(seq)
	for i:= 8 ; i < 15 ; i++ {
		num_peps := prot_len - i
		for j := 1 ; j <= num_peps; j++ {
			end := j + i
			fmt.Println(s[j:end],name)
		}
	}
}

func createLO_matrix(data [][]string) ([]map[string]interface{}, int) {

	AA := [20]string{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}

	LO := make([]map[string]interface{},1)
	ent :=make([]float64,9)
	for i, line := range data {
		AA_dist := make(map[string]interface{})
		for j, field := range line {
			if j <= 19 {
				AA_dist[AA[j]] = field
			} else {
				ent[i],_ = strconv.ParseFloat(field,64)
			}
			
		}
	    	LO = append(LO, AA_dist)			

	}
	return LO[1:10], WhichMax(ent[0:8])
}


func WhichMax (vec []float64) int {
	max := 0.00
	index := 0
	for i,j := range vec {
		if j >= max {
			index = i
			max = j
		}
	}

	return index
}

func CheckNonStandard(AA string) bool {
	switch AA {
	case "X":
		return true
	case "J":
		return true
	case "O":
		return true
	case "B":
		return true
	case "Z":
		return true
	case "U":
		return true
	default:
		return false

	}
}


func score_peptide(peptide string, motif_mat []map[string]interface{}, allele string ,protein string, weight []float64, thres float64) string {

	pep_len := len(peptide)
	len_weight := pep_len - 8
	if (pep_len < 9) {
		out := 0.00
		for i, j := range motif_mat {
			if (i <= 6) {
				next := i+1
				if CheckNonStandard(peptide[i:next]) {
					out += 0.00
				} else {
					
					num,_ := strconv.ParseFloat(j[peptide[i:next]].(string),64)
					out += num
				}
			} else if (i == 8) {
				if CheckNonStandard(peptide[7:8]) {
					out += 0.00
				} else {
					
					num,_ := strconv.ParseFloat(j[peptide[7:8]].(string),64)
					out += num
				}

			}
		}
		//		fmt.Println(peptide,allele,protein,out/float64(pep_len))
		if out >= thres {
			//fmt.Fprintln(o,peptide,allele,protein,out*weight[len_weight])
			out_string := peptide+" "+allele+" "+protein+" "+fmt.Sprintf("%f",out*weight[len_weight])
			return out_string
		}
		return ""

	} else if pep_len == 9 {
		out := 0.00
		for i, j := range motif_mat {
			
			next := i+1
			if CheckNonStandard(peptide[i:next]) {
				out += 0.00
			} else {
				num,_ := strconv.ParseFloat(j[peptide[i:next]].(string),64)
				out += num
			}
		}
		//fmt.Println(peptide,allele,protein,out/float64(pep_len))
		if out >= thres {
			//fmt.Fprintln(o,peptide,allele,protein,out*weight[len_weight])
			out_string := peptide+" "+allele+" "+protein+" "+fmt.Sprintf("%f",out*weight[len_weight])
			return out_string			
		}
		return ""
	} else {
		out := 0.00
		for i, j := range motif_mat {
			if i <= 7 {
				
				next := i+1
				if CheckNonStandard(peptide[i:next]) {
					out += 0.00
				} else {

					num,_ := strconv.ParseFloat(j[peptide[i:next]].(string),64)
					out += num
				}
			} else  {
				i1 := pep_len -1
				i2 := pep_len
				if CheckNonStandard(peptide[i1:i2]) {
					out += 0.00
				} else {

					num,_ := strconv.ParseFloat(j[peptide[i1:i2]].(string),64)
					out += num
				}
			}
		}
		//		fmt.Println(peptide,allele,protein,out/float64(pep_len))
		if out >= thres {
			//fmt.Fprintln(o,peptide,allele,protein,out*weight[len_weight])
			out_string := peptide+" "+allele+" "+protein+" "+fmt.Sprintf("%f",out*weight[len_weight])
			return out_string
		}
		return ""
		
	}
}


func ReadThreshold (Threshold float64, threshold_path string ) (float64) {
	f, err := os.Open(threshold_path)

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()
	thres := 0.00
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.Split(scanner.Text(), " ")
		s,_ := strconv.ParseFloat(line[0],64)
		if s == Threshold {
			thres,_ = strconv.ParseFloat(line[1],64)
			break
		}
	}

	return thres

}

func main() {
	//	var name, seqs = readfasta("test.fasta")
	fastaPtr := flag.String("i", "example.fasta", "name of finput file ")
	weightPtr := flag.String("weights_file", "length_weights.txt", "length correction weights")
	minPtr := flag.Int("min", 9, "minimum peptide length")
	maxPtr := flag.Int("max", 9, "maximum peptide length")
	thresPtr := flag.Float64("threshold", 99.5, "threshold percentile for predictions")
	allelePtr := flag.String("a", "A_02:01", "target MHC-I allele for prediction")
	FilePtr := flag.String("o", "output.txt", "output file for prediction")
	PeptidePtr := flag.Int("P", 0, "input file type (1: peptides ; 0: fasta) default: 0")
	UseWeights := flag.Int("w", 0, "use length weights. (1: yes ; 0: no). requires -w flag to be set. default: 0")
	flag.Parse()	
	LO_matrix_path := os.Getenv("HI_PRED_DATA")+"LO_matrices/"+*allelePtr+"_LO_matrix.txt"
	threshold_path := os.Getenv("HI_PRED_DATA")+"/dist_cutoffs/"+*allelePtr+".dist"

	if _, err := os.Stat(LO_matrix_path); os.IsNotExist(err) {
		fmt.Printf("ERROR: unable to find log odd matrix for"+*allelePtr+ "  \n");
		os.Exit(1)
	}


	f, err := os.Open(LO_matrix_path)
	if err != nil {
		log.Fatal(err)
	}

	// remember to close the file at the end of the program
	defer f.Close()

	// read csv values using csv.Reader
	csvReader := csv.NewReader(f)
	data, err := csvReader.ReadAll()
	if err != nil {
		log.Fatal(err)
	}

	out, err := os.Create(*FilePtr)
	if err != nil {
		fmt.Println(err)
		out.Close()
		return
	}

	weights := []float64 {1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00}
	if (*UseWeights == 1) {
		weights = readLength_weights(*weightPtr)
	} 
	// convert records to array of structs
	LO_matrix,_ := createLO_matrix(data)
	thres := ReadThreshold(*thresPtr,threshold_path)
	fmt.Println(thres)
	if *PeptidePtr == 0 {
		if _, err := os.Stat(*fastaPtr); os.IsNotExist(err) {
			fmt.Printf("ERROR: fasta file not found!\n");
			os.Exit(1)
		} 
		var names , seqs = readfasta(*fastaPtr)

		runtime.GOMAXPROCS(runtime.NumCPU())
		var wg sync.WaitGroup
		numFasta :=len(names)
		wg.Add(numFasta)
		for k := 0 ; k < numFasta; k ++ {
			go func(k int) {
				defer wg.Done()
				prot_len := len(seqs[k])
				s := strings.TrimSpace(seqs[k])
				for i:= *minPtr ; i <= *maxPtr ; i++ {
					num_peps := prot_len - i
					for j := 0 ; j <= num_peps; j++ {
						end := j + i
						pep_string := score_peptide(s[j:end],LO_matrix,*allelePtr,names[k], weights, thres)
						if pep_string != "" {
							fmt.Fprintln(out,pep_string)
							//fmt.Println(pep_string)
						}
					}
				}


			}(k)
		}

		wg.Wait()
	} else {
		if _, err := os.Stat(*fastaPtr); os.IsNotExist(err) {
			fmt.Printf("ERROR: peptide file not found!\n");
			os.Exit(1)
		} 
		var seqs = ReadPeptides(*fastaPtr)
		runtime.GOMAXPROCS(runtime.NumCPU())
		var wg sync.WaitGroup
		numFasta :=len(seqs)
		wg.Add(numFasta)
		for k := 0 ; k < numFasta; k ++ {
			go func(k int) {
				defer wg.Done()
				s := strings.TrimSpace(seqs[k])
				pep_string := score_peptide(s,LO_matrix,*allelePtr,"peptide", weights, thres)
				if pep_string != "" {
					fmt.Fprintln(out,pep_string)
					//fmt.Println(pep_string)
				}

			}(k)
		}

		wg.Wait()


	}
	err = out.Close()
	if err != nil {
		fmt.Println(err)
		return
	}
}
