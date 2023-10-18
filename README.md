Agreement between different scores for the estimation of cardiovascular risk in Brazil, 2013

- `rcv_2013.R` - run this to reproduce results
- `functions.R` - ancillary functions
- `tests.R` - used in development to ensure functions work as intended
- `rcv_2013.Rproj` - open this project file with RStudio for your convenience; otherwise make sure to set its directory (folder) as the working directory in R
- `LICENSE` - MIT license for both the original analysis code and that under `contrib/`
- `README.md` - this file ;)
- `contrib/` - grabbed from the [GitHub package](https://github.com/boyercb/globorisk) so you don't have to; see `LICENSE`
  - `globorisk.R` - function `globorisk()`
  - `sysdata.rda` - coefficients
- `data/` - microdata from the [National Health Survey website](https://www.pns.icict.fiocruz.br/bases-de-dados/)
  - `dicionario_de_variaveis_exames_pns_2013_05052023.xlsx` - data dictionary, downloaded as a zip file also containing the next file
  - `EXAMES-PNS-2013-FINAL_05052023.xlsx` - original data, downloaded as a zip file also containing the previous file
  - `pns2013_exames_2023-05-05.csv` - semicolon-separated values file, obtained by opening the original data with LibreOffice Calc 7.6 and saving as a CSV file
  - `data_dictionary.csv` - comma-separated values file describing the variables used in the analysis
  - `Questionario-PNS-2013.pdf` - the questionnaire (in Portuguese) used in the National Health Survey, downloaded from [its website](https://www.pns.icict.fiocruz.br/questionarios/). Please note that, compared to the questionnaire, variable names in the data have an additional zero between the letter and the digits