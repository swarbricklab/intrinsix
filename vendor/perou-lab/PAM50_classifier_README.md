# PAM50 classifier
Original PAM50 code from Parker et al., JCO 2009, PMID: 19204204

## Repository contents
1. `PAM50_R/`: Contains the functions and fixed parameters for running the PAM50 classifier. These contents include annotation files and the PAM50 centroids file (`pam50_centroids.txt`). The contents of this folder should not need to be changed.
   
2. `PAM50_example/`: Example input/output for running the PAM50 classifier.
   
   a. `PAM50_subtypeAssignment.R`: An example parameter file for running the PAM50 classifier. This file should be modified to point to the correct files and to run with the desired parameters. Note that setting the `calibrationParameters` parameter to NA will force median gene centering within the input test set, and setting this parameter to -1 will not perform any adjustment. This parameter should be set to -1 if the user has already adjusted the input test set (e.g. if [ER/HER2-based centering](https://github.com/perou-lab/PAM50_ER_HER2_normalization) was performed).
   
   b. `sampleInputFile_200subjects.txt`: An example input gene expression file. Rows represent genes and columns represent samples.
   
   c. `sampleOutput/`: Example output files. To test the PAM50 classifier, the parameter file can be run as-is with the example input gene expression file. The output of this run should match this example set of output files. The correlations to PAM50 centroids and PAM50 calls for each sample can be found in the output file ending with `_pam50scores.txt`.

## Dependencies
Before running the PAM50 classifier, the `impute`, `ctc`, and `heatmap.plus` R packages should first be installed.
