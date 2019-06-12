# PBLR
An effective tool for imputing scRNA-seq data by considering cell heterogeneity and prior expression level of dropouts
===============

Overview
--------

PBLR (cell sub-Population based Bounded Low-Rank method) is an effective tool for scRNA-seq data imputation, which can
--------
(1) recover transcriptomic level and dynamics masked by dropouts, 

(2) improve low-dimensional representation, 

(3) restore the gene-gene co-expression relationship, 

(4) and also detect accurate and robust cell subpopulations automatically, shedding light its flexibility and generality for scRNA-seq data analysis. 


Usage
-----
Download the source codes and unzip the MATLAB package. Change the current directory in MATLAB to the folder containing the scripts.

We provide an example with MATLAB live scripts to show how to use PBLR. Please see the /example_data folder and the following details.

Example
-----
In this Matlab live script (generated by Matlab R2018b), we provide an example workflow ([Walkthrough](https://github.com/amsszlh/PBLR/blob/master/example_data/Walkthrough_of_PBLR_on_example_data.pdf)) that outlines the key steps of PBLR. Input data are gene expression data matrix (rows are genes and columns are cells). 

Help
-----
If you have any problem or question using the package please contact zhanglihua@amss.ac.cn

