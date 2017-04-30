package main

import (
	"compress/gzip"   // pour compresser les fichiers
	"fmt"             // pour print
	"io/ioutil"       // pour lire un fichier
	"os"              // pour créer, ouvrir, fermer les fichiers
	"strings"	  // pour les chaines de caractères
	"math"		  // pour calculer la racine carrée
	"strconv"	  // pour convertir une chaine de caractères en nombre (ASCII)
)

func openFile(fileToOpen string) (*os.File, error) {
    return os.OpenFile(fileToOpen, openFileOptions, openFilePermissions)
}

func closeFile(handle *os.File) {
    if handle == nil {
        return
    }

    err := handle.Close()
    if err != nil {
        fmt.Println("[ERROR] Closing file:", err)
    }
}

const openFileOptions int = os.O_CREATE | os.O_RDWR
const openFilePermissions os.FileMode = 0660
const nbFiles int = 15

var fnaNames = []string{"Bacillus_subtilis", "Bacillus_amyloliquefaciens_FZB42", "Bacillus_pumilus_SAFR_032",
			"Bacillus_thuringiensis_BMB171", "Bacillus_cereus_03BB102", "Bacillus_anthracis_Ames",
			"Bacillus_coagulans_2_6", "Bacillus_atrophaeus_1942", "Bacillus_licheniformis_ATCC_14580",
			"Escherichia_coli_K_12_substr__MG1655", "Pseudomonas_aeruginosa_LESB58", "Rhodobacter_sphaeroides_ATCC_17025",
			"Streptomyces_flavogriseus_ATCC_33331", "Micrococcus_luteus_NCTC_2665_uid59033", "Lactococcus_lactis_Il1403"}

var fnaReferences = []string{"NC_000913.fna", "NC_000964.fna", "NC_002662.fna", "NC_003997.fna", "NC_006322.fna",
			     "NC_009428.fna", "NC_009725.fna", "NC_009848.fna", "NC_011770.fna", "NC_012472.fna",
			     "NC_012803.fna", "NC_014171.fna", "NC_014639.fna", "NC_015634.fna", "NC_016114.fna"}

var R_clustering = []string{"NC_000964", "NC_009725", "NC_009848", "NC_014171", "NC_012472",
			    "NC_003997", "NC_015634", "NC_014639", "NC_006322", "NC_000913", 
			    "NC_011770", "NC_009428", "NC_016114", "NC_012803", "NC_002662"}

func append(a string, b string){
	file1,_:=ioutil.ReadFile("fnaFiles/"+a)
	file2,_:=ioutil.ReadFile("fnaFiles/"+b)
	s1 := string(file1)
	s2 := string(file2)
	s3 := s1 + s2
	
	filename := strings.Trim(a, ".fna") + "_" + strings.Trim(b, ".fna") + ".fna"
	os.Create("append/"+filename)
	ioutil.WriteFile("append/"+filename, []byte(s3), 0644)
}

func compress(fnaFilePath string, archivesPath string, fna string) float64 {
	comp := fna
	var zipFile = archivesPath+strings.Trim(comp, ".fna")+".gz"
		
	handle, err := openFile(zipFile)
    
	if err != nil {
        	fmt.Println("[ERROR] Opening file:", err)
    	}	

	zipWriter, err := gzip.NewWriterLevel(handle, 9)
	if err != nil {
        	fmt.Println("[ERROR] New gzip writer:", err)
    	}	
    	file2, _ := ioutil.ReadFile(fnaFilePath+fna)
    	numberOfBytesWritten, err := zipWriter.Write([]byte(file2))
    	if err != nil {
        	fmt.Println("[ERROR] Writing:", err)
    	}
    	err = zipWriter.Close()
    	if err != nil {
        	fmt.Println("[ERROR] Closing zip writer:", err)
    	}
    	fmt.Println("[INFO] Number of bytes written:", numberOfBytesWritten)

    	closeFile(handle)
	
	stat, _ := os.Stat(zipFile)
	size1 := stat.Size()
	
	return float64(size1)
}

func size(zipTarget string) float64{
	stat, _ := os.Stat(zipTarget)
	size := stat.Size()
	
	return float64(size)
}

func distance(z_A float64,z_B float64,z_AB float64,z_AA float64,z_BB float64) float64{
	var d_AB float64
	
	op1 := z_AB/(z_A+z_B)
	op2 := z_AA/(4*z_A)
	op3 := z_BB/(4*z_B)
	d_AB = op1 - op2 - op3
	
	return math.Abs(d_AB)
}

func min_matrix(matrix [][]float64) (int, int, float64) {
	min := matrix[0][1]
	index_i := 0
	index_j := 1
	for i := 0; i < nbFiles; i++ {
		for j := i; j < nbFiles; j++ {
			if matrix[i][j] < min && i != j {
				min = matrix[i][j]
				index_i = i
				index_j = j
			}
		}
	}

	return index_i, index_j, min
}

func minimum_maximum(x, y int) (int, int) {
	if x < y {
		return x, y
	} else {
		return y, x
	}
}

func fmini(x, y float64) float64 {
	if x < y {
		return x
	} else {
		return y
	}
}

func singleLinkage(matrix [][]float64) {
	chaine := ""
	var i_index int
	var j_index int
	var min float64

	nb_colonne := nbFiles

	for nb_colonne != 1 {
		i_index, j_index, min = min_matrix(matrix)
		ind_ensemble, ind_supprime := minimum_maximum(i_index, j_index)
		fmt.Printf("%g, %d, %d\n", min, i_index, j_index)
		
		chaine = "(" + fnaNames[i_index] + "," + fnaNames[j_index] + ")"

		fnaNames[ind_ensemble] = chaine
		fnaNames[ind_supprime] = ""

		for i := 0; i < nbFiles; i++ {
			if i != i_index || i != j_index {
				matrix[i][ind_ensemble] = fmini(matrix[i][i_index], matrix[i][j_index])
				matrix[ind_ensemble][i] = matrix[i][ind_ensemble]
			}
		}

		for i := 0; i < nbFiles; i++ {
			matrix[i][ind_supprime] = 90
			matrix[ind_supprime][i] = 90
		}

		fmt.Println("")
		nb_colonne = nb_colonne - 1
	}

	chaine = chaine + ";"

	content := []byte(chaine)
	ioutil.WriteFile("Tree.txt", content, 0640)
}

func distanceFileFillIn(matrix [][]float64) {
	var res string
	f, _ := os.Create("Distances_R.txt")
	fmt.Fprint(f, "References"+";")
	
	for i := 0; i < nbFiles; i++ {
		fmt.Fprint(f, R_clustering[i], ";")
	}

	fmt.Fprint(f, "\n")
	
	for i := 0; i < nbFiles; i++ {
		fmt.Fprint(f, R_clustering[i], ";")
		for j := 0; j < nbFiles; j++ {
			res = res + strconv.FormatFloat(matrix[i][j], 'f', 6, 64) + " "
			fmt.Fprint(f, matrix[i][j], ";")
		}
		res = res + "\n"
		f.Sync()
		fmt.Fprint(f, "\n")
	}
	
	content := []byte(res)
	ioutil.WriteFile("Distances.txt", content, 0640)
}

func main() {
	
	fnaFiles := make([]string, nbFiles)
	appendedFnaFiles := make([]string, 120)
	gzFiles := make([]string, 135)
	tab := make([]float64, 135)
	
	fmt.Println("\n Creating New Directory for FNA Files....")
	os.Mkdir("fnaFiles",0777)
	
	fmt.Println("\n\n Creating New Directory for Appended FNA Files....")
	os.Mkdir("append",0777)
	
	fmt.Println("\n\n Creating New Directory for Archived FNA Files....")
	os.Mkdir("archives",0777)
	
	for i := 0; i < nbFiles; i++ {
		os.Rename(fnaReferences[i], "fnaFiles/"+fnaReferences[i])
	}
	
	fmt.Println("\n\n Reading FNA Files....")
	dirScan, _ := ioutil.ReadDir("fnaFiles/")

 	for i, file := range dirScan {
		fnaFiles[i]=file.Name()
	}
	
	fmt.Printf("\n\n Appending....")
	for k := 0; k < nbFiles; k++ {
		append(fnaFiles[k],fnaFiles[k])
	}
	fmt.Printf("Done!\n\n")
	
	for i := 0; i < nbFiles; i++ {
		for j := i+1; j < nbFiles; j++ {
			append(fnaFiles[i],fnaFiles[j])
		}
	}
	
	dirScan1, _ := ioutil.ReadDir("append/")
		
	for i, file := range dirScan1 {
		appendedFnaFiles[i]=file.Name()
	}
	
	for i := 0; i < 120; i++ {
		fmt.Println(appendedFnaFiles[i])
	}
	
	fmt.Println("\n\n Compressing....")
	for i := 0; i < nbFiles; i++ {
		tab[i]=compress("fnaFiles/","archives/",fnaFiles[i])
	}
	
	for i := 0; i < 120; i++ {
		tab[i+nbFiles]=compress("append/","archives/",appendedFnaFiles[i])
	}
	fmt.Printf("\nDone!\n\n")
	
	dirScan2, _ := ioutil.ReadDir("archives/")
		
	for i, file := range dirScan2 {
		gzFiles[i]=file.Name()
		fmt.Println(file.Name())
	}
	
	for i := 0; i < 135; i++ {
		tab[i]=size("archives/"+gzFiles[i])
	}
	
	var distMatrix [][]float64
	distMatrix = make([][]float64, nbFiles, nbFiles)

	for i := 0; i < nbFiles; i++ {
		distMatrix[i] = make([]float64, nbFiles)
	}
	
	var k,l,increment,increment_i,increment_j,ij  int
	l = 0
	k = 0
	increment_i = 0 
	increment_j = 0
	increment = nbFiles+1
	
	for i := 0; i < nbFiles; i++ {
		
		ij = increment_i+1
		l=0
		increment_j = increment_i
		
		for j := i; j < nbFiles; j++ {
			
			distMatrix[i][j]=distance(tab[increment_i], tab[increment_j], tab[ij], tab[increment_i+1], tab[increment_j+1])
			distMatrix[j][i]=distMatrix[i][j]
			ij++
			
			increment_j += increment-l
			l++
		}
		
		increment_i += (increment)
		k++
		increment--
	}
	
	fmt.Println("\n\n Saving distances....")
	distanceFileFillIn(distMatrix)
	fmt.Printf("\nDone!\n\n")
	
	fmt.Println("\n\n Single Linkage Algorithm\n\n")
	singleLinkage(distMatrix)
	fmt.Printf("\nDone!\n\n")
}
