% This is a demo showing how to running PBLR using the synthetic dataset 2
clear;clc
addpath(genpath('./'))
addpath PROPACK;
addpath Rank
%% load data
iniData = readtable('raw.txt','Delimiter','\t','ReadRowNames',true,'ReadVariableNames',true);% M is raw data in table form, rows are genes and columns are cells.
%% processing or not 
minGenes = 0; minCells = 0; libararyflag = 0; logNormalize = 1;  % if you don't want do libarary normalization, please set it to 0.
proData = preprocessing(iniData, minCells, minGenes, libararyflag,logNormalize);
M = proData.data; % data matrix
clear iniData
%% select informative genes for clustering
id = gene_selection(M); % id is the index of the selected informative genes; Generally, the number of informative genes can be 1000-2000.
%% clustering
M0 = M(id,:); % submatrix with selected genes
K = 3; % the number of desired clusters given by user. One can also infer the number of clusters based on clustering stability index coph, see function clusteing_NMFs.m and our paper 
numCores = 3; % defined by user
[group,coph] = clusteing(M0,K,numCores);
dlmwrite('cell_group.txt',group)
%% boundary selection by visually checking the estimated boundary. by_default = sophisticated
boundary_selection(M);
%% run PBLR with selected boundary function
% 1: exponetial function; 2: simple piecewise function; 3: sophisticated piecewise function
boundary_function = 3; %  3 is default 
imputation_all = true; % determine impute all genes or not(false by default)
accelate = true;
tic;
X = PBLR_main(M,id,group,boundary_function,imputation_all,numCores,accelate); % return the imputed data matrix
toc;
if isequal(imputation_all,true)
    id = 1:size(M,1);
end
T = array2table(X,'VariableNames',proData.cells,'RowNames',proData.genes(id));
writetable(T,'PBLR_impute.txt','Delimiter','\t','WriteRowNames',1);
