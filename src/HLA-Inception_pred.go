package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"sync"
)

type LO_Matrix struct {
	HLA  string      `json:"HLA"`
	POS1 [20]float64 `json:"Pos1"`
	POS2 [20]float64 `json:"Pos2"`
	POS3 [20]float64 `json:"Pos3"`
	POS4 [20]float64 `json:"Pos4"`
	POS5 [20]float64 `json:"Pos5"`
	POS6 [20]float64 `json:"Pos6"`
	POS7 [20]float64 `json:"Pos7"`
	POS8 [20]float64 `json:"Pos8"`
	POS9 [20]float64 `json:"Pos9"`
}

type LO_Score struct {
	Scores [9]map[string]interface{}
	Thres  float64
	Dist   LO_Dist
	Allele string
}

type LO_Dist struct {
	HLA   string        `json:"HLA"`
	Score [1000]float64 `json:"Score"`
	Perc  [1000]float64 `json:"Perc"`
}

func readfasta(input string) ([]string, []string) {
	f, err := os.Open(input)

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()

	scanner := bufio.NewScanner(f)
	name := make([]string, 0)
	seqs := make([]string, 0)
	out := ""
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] == '>' {

			name = append(name, line[1:])
			if out != "" {
				seqs = append(seqs, out)
				out = ""

			}
			continue
		}
		out += scanner.Text()
	}

	seqs = append(seqs, out)

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
	return name, seqs

}

func ReadPeptides(input string) []string {
	f, err := os.Open(input)

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()

	scanner := bufio.NewScanner(f)
	seqs := make([]string, 0)
	for scanner.Scan() {
		seqs = append(seqs, scanner.Text())

	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
	return seqs

}

func readLengthWeights(input string) map[string]float64 {
	f, err := os.Open(input)

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()
	weight := make(map[string]float64)
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.Split(scanner.Text(), "\t")
		num, _ := strconv.ParseFloat(line[1], 64)
		weight[line[0]] = num
	}

	return weight

}

func gen_pep(seq string, name string) {
	prot_len := len(seq)
	s := strings.TrimSpace(seq)
	for i := 8; i < 15; i++ {
		num_peps := prot_len - i
		for j := 1; j <= num_peps; j++ {
			end := j + i
			fmt.Println(s[j:end], name)
		}
	}
}

func FindBestMatch(allele string) string {
	f, err := os.Open(os.Getenv("HI_PRED_PATH") + "/data/AlleleAlignments.txt")

	if err != nil {
		log.Fatal(err)
	}

	defer f.Close()
	var Match string
	var selAllele string
	var score float64
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		line := strings.Split(scanner.Text(), "\t")
		selAllele = line[0]
		Match = line[1]
		score, _ = strconv.ParseFloat(line[2], 64)
		if selAllele == allele {
			break
		}

	}
	fmt.Println("LOG: best match allele for ", allele, ": ", Match, " score: ", score)
	if score < 0.5 {
		fmt.Println("Warning: Best match allele score is < 0.5 which is some cases indicate a poor alignment.")
		fmt.Println("Proceed with caution.")
	}
	return Match
}

func createLO_matrix(allele string, Thres float64, LOMats []LO_Matrix, Dists []LO_Dist) []LO_Score {
	// set amino acids names
	AA := [20]string{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}
	Alleles := strings.Split(allele, ",")
	LO_outs := make([]LO_Score, len(Alleles))

	for i, a := range Alleles {
		var SEL_ALLELE LO_Matrix
		var SEL_DIST LO_Dist
		var ThresScore float64
		var LO LO_Score
		LO.Allele = a

		// this is where the best match allele function will be called
		BestMatchAllele := FindBestMatch(a)

		NumMats := len(LOMats)
		for i := 0; i < NumMats; i++ {
			if LOMats[i].HLA == BestMatchAllele {
				SEL_ALLELE = LOMats[i]
				break
			} else if i == (NumMats - 1) {
				fmt.Println("Error: Allele not supported")
				os.Exit(1)
			}
		}

		NumDists := len(Dists)
		for i := 0; i < NumDists; i++ {
			if Dists[i].HLA == BestMatchAllele {
				SEL_DIST = Dists[i]
				for index, p := range Dists[i].Perc {
					if p >= Thres {
						ThresScore = Dists[i].Score[index]
						break
					}
				}
				break
			} else if i == (NumMats - 1) {
				fmt.Println("Error: Allele not supported")
				os.Exit(1)
			}
		}

		for pos := 1; pos <= 9; pos++ {
			AA_dist := make(map[string]interface{})
			LO_VEC := reflect.ValueOf(SEL_ALLELE).Field(pos)

			for i := 0; i < LO_VEC.Len(); i++ {

				AA_dist[AA[i]] = LO_VEC.Index(i).Interface()
			}

			LO.Scores[pos-1] = AA_dist
		}
		LO.Thres = ThresScore
		LO.Dist = SEL_DIST
		LO_outs[i] = LO

	}

	return LO_outs
}

func WhichMax(vec []float64) int {
	max := 0.00
	index := 0
	for i, j := range vec {
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

func GetPercentile(Score float64, Dist LO_Dist) float64 {
	Perc := 0.00
	for i, S := range Dist.Score {
		if S > Score || i == len(Dist.Score)-1 {
			Perc = Dist.Perc[i-1]
			break
		} else if S == Score {
			Perc = Dist.Perc[i]
			break
		}
	}
	return Perc
}

// format score string
func FormatScore(out float64, peptide string, allele string, protein string, thres float64, Len int, Start int, End int, SelThres float64, Dist LO_Dist, upstream string, downstream string, LenCorOut float64) string {
	b := [...]string{
		peptide,
		strconv.Itoa(Len),
		allele,
		protein,
		strconv.Itoa(Start),
		strconv.Itoa(End),
		fmt.Sprintf("%f", SelThres),
		upstream,
		downstream,
		fmt.Sprintf("%f", out),
		fmt.Sprintf("%f", GetPercentile(out, Dist)),
		fmt.Sprintf("%f", LenCorOut),
		fmt.Sprintf("%f", GetPercentile(LenCorOut, Dist))}
	return strings.Join(b[:], ";")

}

func FormatScorePeptideMode(out float64, peptide string, allele string, Len int, Dist LO_Dist, LenCorOut float64) string {
	b := [...]string{
		peptide,
		strconv.Itoa(Len),
		allele,
		"PEPTIDE",
		fmt.Sprintf("%f", out),
		fmt.Sprintf("%f", GetPercentile(out, Dist)),
		fmt.Sprintf("%f", LenCorOut),
		fmt.Sprintf("%f", GetPercentile(LenCorOut, Dist))}
	return strings.Join(b[:], ";")

}

// initalize output file
func CreateOutput(file string) *os.File {
	out, err := os.Create(file)
	if err != nil {
		fmt.Println(err)
		out.Close()

	} else {
		b := [...]string{
			"Peptide",
			"Length",
			"Allele",
			"Protein",
			"Start",
			"End",
			"SelectedThres",
			"upstream",
			"downstream",
			"RawScore",
			"RawScorePercentile",
			"LengthCorrectedScore",
			"LengthCorrectedScorePercentile"}

		fmt.Fprintln(out, strings.Join(b[:], ";"))

	}
	return out

}

func CreateOutputPeptideMode(file string) *os.File {
	out, err := os.Create(file)
	if err != nil {
		fmt.Println(err)
		out.Close()

	} else {
		b := [...]string{
			"Peptide",
			"Length",
			"Allele",
			"Protein",
			"RawScore",
			"RawScorePercentile",
			"LengthCorrectedScore",
			"LengthCorrectedScorePercentile"}

		fmt.Fprintln(out, strings.Join(b[:], ";"))

	}
	return out

}

// ! this uses the fast method for different lengths which is different then what was used in the manuscript
// TODO: I will need to come back to this to make it match to the paper
// note: this is where is screens all possible combinations and picks the highest scoring peptide
func score_peptide(peptide string, motif_mat [9]map[string]interface{}, allele string, protein string, weight map[string]float64, thres float64, Len int, Start int, End int, SelThres float64, Dist LO_Dist, upstream string, downstream string) string {

	pep_len := len(peptide)
	len_weight := strconv.Itoa(pep_len)
	if pep_len < 9 {
		out := 0.00
		for i, j := range motif_mat {
			if i <= 6 {
				next := i + 1
				if CheckNonStandard(peptide[i:next]) {
					out += 0.00
				} else {

					num, _ := j[peptide[i:next]].(float64)
					out += num
				}
			} else if i == 8 {
				if CheckNonStandard(peptide[7:8]) {
					out += 0.00
				} else {

					num, _ := j[peptide[7:8]].(float64)
					out += num
				}

			}
		}
		LenCorOut := out + weight[len_weight]
		if LenCorOut >= thres {
			out_string := FormatScore(out, peptide, allele, protein, thres, Len, Start, End, SelThres, Dist, upstream, downstream, LenCorOut)
			return out_string
		}
		return ""

	} else if pep_len == 9 {
		out := 0.00
		for i, j := range motif_mat {

			next := i + 1
			if CheckNonStandard(peptide[i:next]) {
				out += 0.00
			} else {
				num, _ := j[peptide[i:next]].(float64)
				out += num
			}
		}

		LenCorOut := out + weight[len_weight]
		if LenCorOut >= thres {
			out_string := FormatScore(out, peptide, allele, protein, thres, Len, Start, End, SelThres, Dist, upstream, downstream, LenCorOut)
			return out_string
		}
		return ""
	} else {
		out := 0.00
		for i, j := range motif_mat {
			if i <= 7 {

				next := i + 1
				if CheckNonStandard(peptide[i:next]) {
					out += 0.00
				} else {

					num, _ := j[peptide[i:next]].(float64)
					out += num
				}
			} else {
				i1 := pep_len - 1
				i2 := pep_len
				if CheckNonStandard(peptide[i1:i2]) {
					out += 0.00
				} else {

					num, _ := j[peptide[i1:i2]].(float64)
					out += num
				}
			}
		}

		LenCorOut := out + weight[len_weight]
		if LenCorOut >= thres {
			out_string := FormatScore(out, peptide, allele, protein, thres, Len, Start, End, SelThres, Dist, upstream, downstream, LenCorOut)
			return out_string
		}
		return ""

	}
}

func score_peptide_peptide_mode(peptide string, motif_mat [9]map[string]interface{}, allele string, weight map[string]float64, Dist LO_Dist) string {

	pep_len := len(peptide)
	len_weight := strconv.Itoa(pep_len)
	if pep_len < 9 {
		out := 0.00
		for i, j := range motif_mat {
			if i <= 6 {
				next := i + 1
				if CheckNonStandard(peptide[i:next]) {
					out += 0.00
				} else {

					num, _ := j[peptide[i:next]].(float64)
					out += num
				}
			} else if i == 8 {
				if CheckNonStandard(peptide[7:8]) {
					out += 0.00
				} else {

					num, _ := j[peptide[7:8]].(float64)
					out += num
				}

			}
		}
		LenCorOut := out + weight[len_weight]

		out_string := FormatScorePeptideMode(out, peptide, allele, pep_len, Dist, LenCorOut)
		return out_string

	} else if pep_len == 9 {
		out := 0.00
		for i, j := range motif_mat {

			next := i + 1
			if CheckNonStandard(peptide[i:next]) {
				out += 0.00
			} else {
				num, _ := j[peptide[i:next]].(float64)
				out += num
			}
		}

		LenCorOut := out + weight[len_weight]
		out_string := FormatScorePeptideMode(out, peptide, allele, pep_len, Dist, LenCorOut)
		return out_string

	} else {
		out := 0.00
		for i, j := range motif_mat {
			if i <= 7 {

				next := i + 1
				if CheckNonStandard(peptide[i:next]) {
					out += 0.00
				} else {

					num, _ := j[peptide[i:next]].(float64)
					out += num
				}
			} else {
				i1 := pep_len - 1
				i2 := pep_len
				if CheckNonStandard(peptide[i1:i2]) {
					out += 0.00
				} else {

					num, _ := j[peptide[i1:i2]].(float64)
					out += num
				}
			}
		}

		LenCorOut := out + weight[len_weight]
		out_string := FormatScorePeptideMode(out, peptide, allele, pep_len, Dist, LenCorOut)
		return out_string

	}
}

func GetFlank(s string, start int, end int) (upstream string, peptide string, downstream string) {
	FlankBegin := (start - 3)
	FlankEnd := (end + 3)
	// upstream:=""
	// peptide:=""
	// downstream:=""
	if FlankBegin < 0 {
		upstream = strings.Repeat("_", FlankBegin*-1) + s[0:start+1]
	} else {
		upstream = s[FlankBegin : start+1]
	}
	if FlankEnd >= len(s) {
		downstream = s[end-1:len(s)] + strings.Repeat("_", FlankEnd-len(s))
	} else {
		downstream = s[end-1 : FlankEnd]
	}
	peptide = s[start:end]
	return upstream, peptide, downstream
}

func main() {
	//	var name, seqs = readfasta("test.fasta")
	fastaPtr := flag.String("i", "example.fasta", "name of finput file ")
	weightPtr := flag.String("w", "", "length correction weights")
	LenPtr := flag.String("l", "9", "minimum peptide length")
	thresPtr := flag.Float64("threshold", 99.5, "threshold percentile for predictions")
	allelePtr := flag.String("a", "A_02:01", "target MHC-I allele for prediction")
	FilePtr := flag.String("o", "output.txt", "output file for prediction")
	PeptidePtr := flag.Int("P", 0, "input file type (1: peptides ; 0: fasta) default: 0")
	flag.Parse()
	LO_matrix_data := os.Getenv("HI_PRED_PATH") + "/data/LO.json"
	LO_dist_data := os.Getenv("HI_PRED_PATH") + "/data/Dist.json"

	// check the LenPtr flag to make sure it is a valid length and
	// split the string into a slice if multiple lengths are given
	var out *os.File
	var Lens []int
	if *PeptidePtr == 0 {
		out = CreateOutput(*FilePtr)
		for _, l := range strings.Split(*LenPtr, ",") {
			num, _ := strconv.Atoi(l)
			if num >= 8 || num <= 15 {
				Lens = append(Lens, num)
			}
		}
	} else {
		out = CreateOutputPeptideMode(*FilePtr)
	}

	file, err := os.Open(LO_matrix_data)
	if err != nil {
		fmt.Println("Error opening file:", err)
		return
	}
	defer file.Close()

	var LOMatObjects []LO_Matrix
	decoder := json.NewDecoder(file)
	err = decoder.Decode(&LOMatObjects)
	if err != nil {
		fmt.Println("Error decoding JSON:", err)
		return
	}

	file2, err := os.Open(LO_dist_data)
	if err != nil {
		fmt.Println("Error opening file:", err)
		return
	}
	defer file2.Close()

	var LO_Dists []LO_Dist
	decoder2 := json.NewDecoder(file2)
	err = decoder2.Decode(&LO_Dists)
	if err != nil {
		fmt.Println("Error decoding JSON:", err)
		return
	}

	weights := map[string]float64{}
	if *weightPtr != "" {
		weights = readLengthWeights(*weightPtr)
	} else {
		weights = readLengthWeights(os.Getenv("HI_PRED_PATH") + "/data/DefaultLengthWeights.txt")
	}

	// convert records to array of structs
	LO_matrices := createLO_matrix(*allelePtr, *thresPtr, LOMatObjects, LO_Dists)

	if *PeptidePtr == 0 {

		if _, err := os.Stat(*fastaPtr); os.IsNotExist(err) {
			fmt.Printf("ERROR: fasta file not found!\n")
			os.Exit(1)
		}
		var names, seqs = readfasta(*fastaPtr)

		runtime.GOMAXPROCS(runtime.NumCPU())
		for _, LO := range LO_matrices {
			var wg sync.WaitGroup
			numFasta := len(names)
			wg.Add(numFasta)
			LO_matrix := LO.Scores
			thres := LO.Thres
			SEL_DIST := LO.Dist
			allele := LO.Allele
			fmt.Println("LOG: running: ", allele, " in fastamode with a threshold of:", *thresPtr, "(", thres, ")", " and a length/s of:", Lens)
			for k := 0; k < numFasta; k++ {
				go func(k int) {
					defer wg.Done()
					prot_len := len(seqs[k])
					s := strings.TrimSpace(seqs[k])
					// ? This is a bit of a legacy line. I orginally wanted the ability to do a range of lengths, but no I am going to let individual lengths get batchced togethee
					for _, PepLen := range Lens {
						num_peps := prot_len - PepLen
						for j := 0; j <= num_peps; j++ {
							end := j + PepLen
							upstream, peptide, downstream := GetFlank(s, j, end)
							pep_string := score_peptide(peptide, LO_matrix, allele, names[k], weights, thres, PepLen, j+1, end, *thresPtr, SEL_DIST, upstream, downstream)
							if pep_string != "" {
								fmt.Fprintln(out, pep_string)
								//fmt.Println(pep_string)
							}
						}
					}

				}(k)
			}
			wg.Wait()
		}
	} else {
		if _, err := os.Stat(*fastaPtr); os.IsNotExist(err) {
			fmt.Printf("ERROR: peptide file not found!\n")
			os.Exit(1)
		}
		seqs := ReadPeptides(*fastaPtr)

		runtime.GOMAXPROCS(runtime.NumCPU())
		for _, LO := range LO_matrices {
			var wg sync.WaitGroup
			numSeqs := len(seqs)
			wg.Add(numSeqs)
			LO_matrix := LO.Scores
			SEL_DIST := LO.Dist
			allele := LO.Allele
			for k := 0; k < numSeqs; k++ {
				go func(k int) {
					defer wg.Done()
					pep_string := score_peptide_peptide_mode(seqs[k], LO_matrix, allele, weights, SEL_DIST)
					fmt.Fprintln(out, pep_string)
				}(k)
			}

			wg.Wait()
		}
	}
	err = out.Close()
	if err != nil {
		fmt.Println(err)
		return
	}
}
