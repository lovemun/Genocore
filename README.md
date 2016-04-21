# GenoCore : a simple and fast algorithm for core subset selection from large genotype datasets
Usage : Rscript run_genocore.R input_file -cv coverage -d difference -o result_file_name

Genocore is available at https://github.com/lovemun/Genocore. Source code was written in R language and supported on windows and linux platform. 

# Example

$ git clone https://github.com/lovemun/Genocore

$ cd Genocore

$ Rscript run_genocore.R wheat_subset.csv -cv 100 -d 0.001 -o example &

# Contact

lovemun@hanyang.ac.kr
