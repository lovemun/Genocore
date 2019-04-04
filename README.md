# GenoCore : a simple and fast algorithm for core subset selection from large genotype datasets
Usage : Rscript run_genocore.R input_file -cv coverage -d difference -o result_file_name

Genocore is available at https://github.com/lovemun/Genocore. Source code was written in R language and supported on windows and linux platform. 

## Requirement
- python : rst2pdf modules
  - install command in command line : pip install rst2pdf
- R : argparse package
  - install command in R : install.packages("argparse")

## Example

$ git clone https://github.com/lovemun/Genocore

$ cd Genocore

$ Rscript run_genocore.R wheat_subset.csv -cv 99 -d 0.001 -o example &

## Contact

lovemun@kribb.re.kr

## Citation

Jeong S, Kim JY, Jeong SC, Kang ST, Moon JK, et al. (2017) GenoCore: A simple and fast algorithm for core subset selection from large genotype datasets. PLOS ONE 12(7): e0181420. https://doi.org/10.1371/journal.pone.0181420
