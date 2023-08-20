# DNAMarkMaker

## Table of contents
 - [Introduction of DNAMarkMaker](#Introduction-of-DNAMarkMaker)
 - [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [Installation using bioconda](#Installation-using-bioconda)
 - [Usage](#Usage)
  + [command : target_SNP_selection](#command-:-target_SNP_selection)
  + [command : ARMS_preparation](#command-:-ARMS_preparation)
  + [command : tri_ARMS](#command-:-tri_ARMS)
  + [command : tetra_ARMS](#command-:-tetra_ARMS)
  + [command : CAPS](#command-:-CAPS)
 - [The input file format](#The-input-file-format)
 - [The example of execution](#The-example-of-format-execution)
 - [The output file format](#The-output-file-format)

## Introduction of DNAMarkMaker
  <img src="https://github.com/YuSugihara/QTL-seq/blob/master/images/1_logo.png" width=200>
  
  DNAMarkMaker is a tool to develop ARMS and/or CAPS markers that target SNPs between cultivars/lines, utilizing data from next-generation sequencing.
  
## Installation
### Dependencies
   - [samtools](https://github.com/samtools/samtools)
   - [primer3-py](https://github.com/libnano/primer3-py)

### Installation using bioconda
  You can install DNAMarkMaker using bioconda.
  ```
  conda install -c bioconda DNAMarkMaker
  ```
  Alternatively, if you want to create DNAMarkMaker specific environment.
  ```
  conda create -n DNAMarkMaker -c bioconda DNAMarkMaker
  conda activate DNAMarkMaker 
  ```
  
## Usage
  'DNAMarkMaker' offers five commands. Users can specify which command to run using the '-w' option. Here is a summary of each command.

### command : target_SNP_selection
  Purpose: Identify SNPs from the BAM file, which contains alignment data of two breeds.
  Description: This command should be the first one executed when using DNAMarkMaker, as it lays the foundation for subsequent marker design by selecting target SNPs.

  Required options
  ```
  -Abam ABAM                                  Full path of A bam
  -Bbam BBAM                                  Full path of B bam
  -reference REFERENCE                        Full path of reference fasta
  -position POSITION                          Target chromosome position [chr:start:end]
  -o OUTPUT_DIR                               Output directory
  ```
  
  Additional options
  ```
  -Aname ANAME                                A name (A)
  -Bname BNAME                                B name (B)
  -min_depth MIN_DEPTH                        Minimum depth of target SNP (10)
  -max_depth MAX_DEPTH                        Maximum depth of target SNP (99)
  -minMQ MINMQ                                Minimum mapping quality detected from bam (0)
  -minBQ MINBQ                                Minimum base quality detected from bam (13)
  -Bhetero BHETERO                            Whether to target heterozygous SNP in B (no)
  -Bsim BSIM                                  B simulation file
  -Cbam CBAM                                  Full path of C bam
  -Csim CSIM                                  C simulation file
  ```


### command : ARMS_preparation
   Purpose: Design breed-specific primers for the creation of ARMS markers.
   Description: This command utilizes the inter-cultivar SNP information previously identified to design primers specific to each breed or line, preparing for the next stage of ARMS marker development.

  Required options
  ```
  -o OUTPUT_DIR                              Output directory
  ```
  
  Additional options
  ```
  -recipe RECIPE                             Full path of primer recipe file
  ```

### command : tri_ARMS
   Purpose: Develop tri-ARMS markers.
   Description: Leveraging the primers designed in the `ARMS_preparation` phase, this command facilitates the development of tri-ARMS markers, which involve three primers for amplification.

    Required options
    ```
    -o OUTPUT_DIR                             Output directory
    ```

    Additional options
    ```
    -recipe RECIPE                           Full path of primer recipe file
    -PCR_max_size PCR_MAX_SIZE               Maximum size of PCR product (700)
    -PCR_min_size PCR_MIN_SIZE               Minimum size of PCR product (100)
    -SNP_dist SNP_DIST                       Target SNP distance (100:300)
    -make_html MAKE_HTML                     Whether to html file (yes)
    ```

### command : tetra_ARMS
   Purpose: Develop tetra-ARMS markers.
   Description: Similarly to the `tri_ARMS` command, this command uses the primers from `ARMS_preparation` but develops tetra-ARMS markers, which utilize four primers in the amplification process.

    Required options
    ```
    -o OUTPUT_DIR                 Output directory
    ```

    Additional options
    ```
    -first_size FIRST_SIZE                   Size of first band (100:500)
    -second_size SECOND_SIZE                 Size of second band (600:1000) 
    -make_html MAKE_HTML                     Whether to html file (yes)
    ```

### command : CAPS
   Purpose: Create CAPS markers.
   Description: This command supports the development of CAPS markers by introducing specific restriction enzymes to the identified inter-cultivar SNP information, resulting in markers that are identified based on the cleavage patterns of these enzymes on the DNA fragments.

    Required options
    ```
    -o OUTPUT_DIR                            Output directory
    -restriction_enzyme RESTRICTION_ENZYME   Full path of restriction enzyme file
    ```

    Additional options
    ```
    -recipe RECIPE                           Full path of primer recipe file
    -PCR_max_size PCR_MAX_SIZE               Maximum size of PCR product (1000)
    -PCR_min_size PCR_MIN_SIZE               Minimum size of PCR product (500)
    -fragment_min_size FRAGMENT_MIN_SIZE     Minimum fragment size of restricted PCR product (200)
    -make_html MAKE_HTML                     Whether to html file (yes)
    ```


Band size of PCR product amplified by developped marker
XXXXXXXXXX

## The input file format
BAM File
 A sorted binary format file that stores NGS reads alignment data bai file.

Reference 
A fasta format file for alignmnet.

Simulation Files
A txt file contains space-separated arbitrary simulation confidence intervals for each depth. Users can obtain this from the provided URL (https://github.com/SegawaTenta/DNAMarkMaker_manual/tree/main/simulation) and modify the values if required for specific analyses.

Recipe File
A txt file contains the options used by the primer design tool, primer3. The default version of this file can be downloaded from the provided URL (https://github.com/SegawaTenta/DNAMarkMaker_manual/blob/main/primer_recipe/primer_recipe.txt), but users can also tailor the values inside to suit specific needs or experimental conditions.

Restriction enzyme file
A txt file contains space-separated arbitrary the name of restriction enzyme and recognizing sequense (https://github.com/SegawaTenta/DNAMarkMaker_manual/blob/main/restriction_enzyme/restriction_enzyme.txt) and modify the contents if required for specific analyses.

## The example of execution

target_SNP_selection
ex1) Homozygous plant
```
 DNAMarkMaker -w target_SNP_selection \
              -reference Full/path/to/referance.fasta \
              -Abam Full/path/to/A.bam \
              -Bbam Full/path/to/B.bam \
              -position chr1:10000:50000 \
              -o example1
```                           

ex2) Autotetraploid plants (terget of simplex in B)
```
  DNAMarkMaker -w target_SNP_selection \
               -reference Full/path/to/referance.fasta \
               -Abam Full/path/to/A.bam \
               -Bbam Full/path/to/B.bam \
               -position chr1:10000:50000 \
               -o example3 \
               -min_depth 59 \
               -max_depth 300 \
               -Bhetero yes\
               -Bsim downloaded/sim_simplex_AAAa_95.txt
```

ex3) Heterozygous plants
```
 DNAMarkMaker -w target_SNP_selection \
              -reference Full/path/to/referance.fasta \
              -Abam Full/path/to/A.bam \
              -Bbam Full/path/to/B.bam \
			        -Cbam Full/path/to/F1.bam \
              -position chr1:10000:50000 \
              -o example3 \
```

ARMS_preparation
ex)
```
 DNAMarkMaker -w ARMS_preparation \
              -o example1
```

tri_ARMS
ex)
```
 DNAMarkMaker -w tri_ARMS \
              -o example1
```

tetra_ARMS
ex)
```
 DNAMarkMaker -w tetra_ARMS \
              -o example1
```

CAPS
ex)
```
 DNAMarkMaker -w CAPS \
		          -restriction_enzyme downloaded/restriction_enzyme.txt\
              -o example1
```

## The output file format
The formats and contents of the output files align with those delineated in the manual for the GUI version, available at URL(https://github.com/SegawaTenta/DNAMarkMaker_manual/blob/main/Manual.pdf).
