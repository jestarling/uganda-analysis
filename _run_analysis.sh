## This script runs the analysis in Starling 2020, 
## "Monotone function estimation in the presence of extreme 
## data coarsening: Analysis of preeclampsia and birth weight in urban Uganda""
## and creates all figures and tables in the ./output
## directory.  This file should be run from the project
## root, with the three R scripts in a "code" folder and
## the data file in a "data" folder.

## Create output directory.
rm -r output
mkdir output

### Run Early Medical Abortion analysis.
Rscript ./code/01-sim-study-mse-inflation.R
Rscript ./code/02-sim-study-monotonicity.R
Rscript ./code/03-mulago-analysis.R


