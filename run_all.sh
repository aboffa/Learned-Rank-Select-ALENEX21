#!/bin/bash

mkdir -p results

time {

  ./Comparison datasets/DNA_1 datasets/DNA_2 datasets/DNA_3 \
  datasets/5GRAM_1 datasets/5GRAM_2 datasets/5GRAM_3 \
  datasets/URL_1 datasets/URL_2 datasets/URL_3  \
  >> results/DNA_5GRAM_URL_comparison.csv \
  2>> results/comparison.err

./Comparison \
  datasets/GOV2/100K-1M/* \
  >> results/GOV2_100K-1M_comparison.csv \
  2>> results/comparison.err

./Comparison \
  datasets/GOV2/1M-10M/* \
  >> results/GOV2_1M-10M_comparison.csv \
  2>> results/comparison.err

./Comparison \
  datasets/GOV2/10M-/* \
  >> results/GOV2_10M-_comparison.csv \
  2>> results/comparison.err

  printf "AVERAGE_GOV2_10M-" >> results/GOV2_averages.csv
  for i in {2..114}
  do
    printf "," >> results/GOV2_averages.csv
    awk -v var="$i" -F ',' '{ total += $var } END { printf total/NR }' results/GOV2_10M-_comparison.csv >> results/GOV2_averages.csv
  done
  printf "\n" >> results/GOV2_averages.csv

  printf "AVERAGE_GOV2_1M-10M" >> results/GOV2_averages.csv
  for i in {2..114}
  do
    printf "," >> results/GOV2_averages.csv
    awk -v var="$i" -F ',' '{ total += $var } END { printf total/NR }' results/GOV2_1M-10M_comparison.csv >> results/GOV2_averages.csv
  done
  printf "\n" >> results/GOV2_averages.csv

  printf "AVERAGE_GOV2_100K-1M" >> results/GOV2_averages.csv
  for i in {2..114}
  do
    printf "," >> results/GOV2_averages.csv
    awk -v var="$i" -F ',' '{ total += $var } END { printf total/NR }' results/GOV2_100K-1M_comparison.csv >> results/GOV2_averages.csv
  done
  printf "\n" >> results/GOV2_averages.csv

  cat results/DNA_5GRAM_URL_comparison.csv > results/comparison.csv
  cat results/GOV2_averages.csv >> results/comparison.csv
  sed -n '6p' <  results/GOV2_10M-_comparison.csv >> results/comparison.csv

}