# Matlab code
This repository contains the codes used in the paper 'Revealing static and dynamic biomarkers from postprandial metabolomics data through coupled matrix and tensor factorizations'.

All the implementations are tested on MacOS version 10.15.3.

## Packages needed for implementation
*  Brett W. Bader, Tamara G. Kolda and others. MATLAB Tensor Toolbox, Version 3.1. Available online at https://www.tensortoolbox.org, 2020
*  Evrim Acar,  Evangelos E Papalexakis et al, “Structure-revealing data fusion”, 2014. CMTF Toolbox Available online at https://github.com/eacarat/CMTF_Toolbox
*  Daniel M. Dunlavy, Tamara G. Kolda, and Evrim Acar, “Poblano v1.0: A Matlab Toolbox for Gradient-Based Optimization”, 2010. Available online at https://github.com/sandialabs/poblano_toolbox
* Eigenvector Research, DataSet Object, available online at https://eigenvector.com/software/dataset-object/

All the packages should be installed as a subfolder under Matlab path


## Descriptions of the codes 

*  The file named 'static_dynamic_acmtf.m' is the code for fitting the fasting and T0-corrected dynamic metabolomics data with an ACMTF model using the CMTF toolbox with multiple initialisations, and also plotting the best factors by the ACMTF model 
*  The file named 'static_dynamic_acmtf_replicability.m' is an example code for checking the replicability of an ACMTF model
     
   ### Descriptions of the codes under the folder 'functions'
   The file 'show_spread.m' and 'bar_wrange.m' are used to plot the spread of the weights by an ACMTF model in matrix and tensor, where the weights are taken from specific runs with the smallest objective function values
   The file 'check_spread_only.m' exports the factors of an ACMTF model from specific runs which give the smallest objective function values
   The file 'plot_metab_set.m' is used to generate indexes and plotting markers for metabolite classess according to their density and size
   The file 'preprocess.m' is used to preprocess a tensor/matrix by centering across the first mode, and scaling within the third mode (second mode if it is for matrix)
