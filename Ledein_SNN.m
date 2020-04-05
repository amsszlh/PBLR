function [group,coph] = Ledein_SNN(Data,system_used)
if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac';
end

% use graph based leiden algorithm
% Replace the following line by the appropriate path for Rscript
if strcmp(system_used,'Windows')
    Rscript = '"C:\Program Files\R\R-3.5.1\bin\Rscript"'; % for 64-bit windows
elseif strcmp(system_used,'Mac')
    Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
end

filefolder = 'intermediateFiles';
if ~isfolder(filefolder)
    mkdir intermediateFiles
end
writetable(Data,fullfile(filefolder,'data_temporal.txt'),'Delimiter','\t','WriteRowNames',1);
% Calling R
RscriptFileName = ' ./identify_clusters_fast.R ';
eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);
res = 0.1:0.1:0.5; n = length(res);
identity_res_record = cell(1,n);
for j = 1:n
    identity = readtable(fullfile(filefolder,['identity_clustering_res',num2str(res(j)),'.txt']),'Delimiter','\t','ReadRowNames',1);
    identity_res_record{1,j} = identity.Var1;
end
[group,coph] = build_consensus(identity_res_record);
