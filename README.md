# PBLR
An accurate tool for imputing scRNA-seq data by considering cell heterogeneity and prior expression level of dropouts
===============

Overview
--------

PBLR (cell sub-Population based Bounded Low-Rank method) is an effective tool for scRNA-seq data imputation, which can
--------
(1) recover transcriptomic level and dynamics masked by dropouts, 

(2) improve low-dimensional representation, 

(3) restore the gene-gene co-expression relationship, 

(4) and also detect accurate and robust cell subpopulations automatically, shedding light its flexibility and generality for scRNA-seq data analysis. 

Systems Requirements
--------------------
This Package has been tested using MATLAB 2018a on Mac OS/64-bit Windows. 


Usage
-----

Unzip the package. Change the current directory in Matlab to the folder containing the scripts.

PBLR_demo.m -- an example run of PBLR on a synthetic dataset

Please refer to PBLR_demo.m for instructions on how to use this code.
Input Data are gene expression data matrix (rows are genes and columns are cells). 

If you have any problem or question using the package please contact zhanglihua@amss.ac.cn

