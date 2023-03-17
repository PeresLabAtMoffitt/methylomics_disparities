# methylomics_disparities

# Data
The pipeline needs the data "cleaned_07082022.rda" as a start.

# RE data preparation
I created 3 Rmd files to prep the 3 different RE.  
You should run first "01.LTR_retrotransposons.Rmd" as the grooming step is done in this file only.  
Then you can run "02.L1_retrotransposons.Rmd" and "03.Alu_retrotransposons.Rmd".
Within each file, I ran the prediction, trimming and aggregation steps. After each of the step, I saved a rds file for a faster run of the file. 
You will have to un-comment the code to create the rds the first time.

# Clustering

# Statistical analyses
