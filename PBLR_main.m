function Xi = PBLR_main(M,id,group,boundary_function,imputing_all,numCores,accelate)
% M: the data matrix
% group: the label information is know in advance or estimated by clustering
% boundary_function: the estimated boundary type. 1: exponetial function; 2: simple piecewise function; 3: sophisticated piecewise function
% imputing_all: true or false, determine impute all genes or not(false by default)
% numCores: the number of cores used for parallel
% accelate: true or false. true represents adopting accelated algorithm:
% truncked SVD, false represents adopting full svd
if ~exist('imputing_all','var') || isempty(imputing_all)
    imputing_all = false;
end
poolobj = gcp('nocreate');
if exist('numCores','var') && isempty(poolobj)
    parpool('local',numCores);    
end
%% run PBLR on each group
M1 = M(id,:); [m,n] = size(M); nid = setdiff(1:m,id); M2 = M(nid,:); group2 = ones(1,n);
disp('imputation on submatrix with selected genes: ')
Xi1 = PBLR(M1,group,boundary_function,accelate); 
if imputing_all
    disp('imputation on the remaining submatrix: ')
    Xi2 = PBLR(M2,group2,boundary_function,accelate);
    Xi = zeros(m,n);
    Xi(id,:) = Xi1; Xi(nid,:) = Xi2;
else
    Xi = Xi1;
end
