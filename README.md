# Learned-Rank-Select-ALENEX21
This repository contains the code to reproduce the experiments in the ALENEX21 paper:

> Antonio Boffa, Paolo Ferragina, and Giorgio Vinciguerra. A "learned" approach to quicken and compress rank/select dictionaries. In Proceedings of the Symposium on Algorithm Engineering and Experiments (ALENEX). SIAM, 2020.

In brief, the paper follows a recent line of research on the so-called learned data structures. It provides the first “learned”
scheme for implementing a rank/select dictionary over compressed space. In particular, it introduces a novel lossless
compressed storage scheme for the input dictionary S which turns this problem into the one of approximating
a set of points in the Cartesian plane via segments, so that the storage of S can be defined by means of a compressed 
encoding of these segments and the “errors” they do in approximating the input integers. Proper algorithms and data 
structures are then added to this compressed storage scheme to support fast rank and select operations.

## Build and run

Please clone this repository using  the flag `--recurse-submodules`. To run the experiments you need CMake 3.8+, and a compiler with support for C++17, Boost, and OpenMP.
To compile the executables, issue the following commands:

    mkdir build && cmake . -B build -DCMAKE_BUILD_TYPE=Release
    cd build && make
    cp build/Comparison ./Comparison
    cp build/Generate_datasets ./Generate_datasets

The latter commands compile both the program for the benchmark (`comparison.cpp`) and the program for the datasets
generation (`generate_datasets.cpp`). To get the already manipulated datasets you can download them [here](https://drive.google.com/drive/folders/1K78tr9maRMPBhjx0Uo_SklogPY9wC7S5?usp=sharing).
Put them in the directory `datasets`. 

The experiments can be run with the following script, which will populate a `result` directory with a csv file:

    bash run_all.sh
    
The experiments require at least 32 GB of RAM and may take quite some time to finish (approximately 9 hours on our machine, whose specs are detailed below). 

## Tests

Since this benchmark deals with very different data structures implementations that have slightly different operations
there is a bunch of tests that check the actual correctness of every implementation and every wrapper. 
They are in the directory `tests`. To run them: 
    
    tests/my_tests 

## Datasets 

They are generated in the following way: 

### GOV2

Download the full GOV2 Inverted Index [here](https://drive.google.com/drive/folders/1K78tr9maRMPBhjx0Uo_SklogPY9wC7S5?usp=sharing). Then run:

    Generate_datasets write_integers_from_GOV2 <gov2_file>

The latter command creates the directories: 

* `datasets/GOV2/100K-1M` of all the inverted lists with a number of elements between 100 thousand an 1 million.
* `datasets/GOV2/1M-10M` of all the inverted lists with a number of elements between 1 million an 10 million.
* `datasets/GOV2/10M-` of all the inverted lists with a number of elements grater than 10 million.

All these files containing inverted lists are ready for the benchmark.

### URL 

Download the dataset of URLs used from Classification Task available [here](https://www.kaggle.com/shawon10/url-classification-dataset-dmoz) and extract the second column: 
 
    awk -F "\"*,\"*" '{print $2}' URL_Classification.csv > URL_Classification_justURL.csv 
    
Then download the dataset of URLs of News available [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ILAT5B) and extract the third column of both files: 

    awk -F "\"*,\"*" '{print $3}' news-week-18aug24.csv > news-week_justURL.csv
    awk -F "\"*,\"*" '{print $3}' news-week-17aug24.csv >> news-week_justURL.csv
    
    
Then download the dataset of URLs of DOI available [here](https://archive.org/details/doi-urls) and extract the second column:
    
    awk -F "\"*,\"*" '{print $2}' 2013.csv > 2013_just_URL.csv 
    
Merge and sort all the URLs:
    
    cat URL_Classification_justURL.csv 2013_just_URL.csv news-week_just_URL.csv >>  my_url_all.txt
    sort -o my_url_all.txt my_url_all.txt 
        
Compile the Burrows–Wheeler transform tool available [here](https://people.unipmn.it/manzini/bwtdisk/index.html) and run: 

    bwte my_url_all.txt -m 8162
    
The latter command creates the file `my_url_all.txt.bwt`
In the end extract the positions of the most common character using: 

    Generate_datasets write_integers_from_BWT my_url_all.txt.bwt
    
### 5GRAM

We use more than 100 GB of 5GRAM from [Google corpus](https://storage.googleapis.com/books/ngrams/books/datasetsv3.html) Version 20090715. 
The exact list of the used files is in [this](datasets/list_5gram_used) file.
To aggregate all the 5GRAMs using run: 

    Generate_datasets write_google_ngram_v1 googlebooks-eng-all-5gram-20090715-xxx.csv ... googlebooks-eng-all-5gram-20090715-yyy.csv
   
This creates a file called `ALL_google-eng-all-5gram.txt` that contains all the different 5grams in the files.
Compile the Burrows–Wheeler transform tool available [here](https://people.unipmn.it/manzini/bwtdisk/index.html) and run:

    bwte ALL_google-eng-all-5gram.txt -m 8162
    
The latter command creates the file `ALL_google-eng-all-5gram.txt.bwt`. In the end extract the positions of the most common character using: 

    Generate_datasets write_integers_from_BWT ALL_google-eng-all-5gram.txt.bwt
    
### DNA

Download the full human genome file from [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39), 
file type Genomic FASTA (.fna) and extract it. Then get the just "ACGT" version of the genome:  
      
    grep -v '^>' GCF_000001405.39_GRCh38.p13_genomic.fna > genomemic_tmp.unmasked
    awk '{if(/^[^>]/)$0=toupper($0);print $0}' genomemic_tmp.unmasked  > tmp.txt && mv tmp.txt genomemic_tmp.unmasked
    tr -d 'N' < genomemic_tmp.unmasked > tmp.txt && mv tmp.txt genomemic_tmp.unmasked
    tr -d 'B' < genomemic_tmp.unmasked > tmp.txt && mv tmp.txt genomemic_tmp.unmasked
    tr -d 'M' < genomemic_tmp.unmasked > tmp.txt && mv tmp.txt genomemic_tmp.unmasked
    tr -d '\n' < genomemic_tmp.unmasked > tmp.txt && mv tmp.txt genomemic.unmasked
    rm genomemic_tmp.unmasked
    
The latter commands return the file `genomemic.unmasked` of 3.2GB, truncate the file to the first 1GB:

    truncate -s 1000000000 genomemic.unmasked
    
In the end to extract the positions of all the nucleotides run: 

    Generate_datasets write_integers_from_genome genomemic.unmasked

## Analyse the results

The output files can be analysed using the following Jupyter notebooks:

- `plots/plot_bar.ipynb` generates the Figure 3 of the paper. The output is the file `plot_bars.tex`.
- `plots/plot_space_time.ipynb` generates the Figure 4 and the Figure 5 of the paper. The outputs are the files `plot_rank.tex` and `plot_select.tex`.


## Test environment

The code was tested on the following machine:

| Component | Specs                                     |
|-----------|-------------------------------------------|
| CPU       | Intel(R) Xeon(R) CPU E5-2407 v2 @ 2.40GHz |
| RAM       | 41 GB                                     |
| L1 cache  | 32 KB (data) 32 KB (instructions)         |
| L2 cache  | 256 KB                                    |
| L3 cache  | 1 MB                                      |
| OS        | Ubuntu 16.04.6 LTS                        |
| Compiler  | gcc 9.2.1                                 |
| CMake     | version 3.13.2                            |

## License

This project is licensed under the terms of the GNU General Public License v3.0.

If you use this code for your research, please cite:

```
@inproceedings{Boffa:2021,
	Author = {Boffa, Antonio and Ferragina, Paolo and Vinciguerra, Giorgio},
	Booktitle = {Proceedings of the Symposium on Algorithm Engineering and Experiments (ALENEX)},
	Publisher = {{SIAM}},
	Title = {A ``learned'' approach to quicken and compress rank/select dictionaries},
	Year = {2021}}
```