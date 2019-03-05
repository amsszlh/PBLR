function Xi = PBLR(M,group,boundary_function,accelate)
% <INPUT>
% M: the data matrix
% group: the label information is know in advance or estimated by clustering
% boundary_function: the estimated boundary type. 1: exponetial function; 2: simple piecewise function; 3: sophisticated piecewise function
% accelate: true or false. true represents adopting accelated algorithm:
% truncked SVD, false represents adopting full svd
% <OUTPUT>
% Xi: the imputed data matrix

if ~exist('accelate','var') || isempty(accelate)
    accelate = true;
end

K = max(group);  X_groups = cell(1,K); U_record = cell(1,K);
for i = 1:K
    idxCluster = group == i;
    x = M(:,idxCluster);
    X_groups{1,i} = x;
 
% compute the dropout on each group
    num_nonzeros = sum(x~=0,2);
    avgs = sum(x,2)./num_nonzeros; avgs(isnan(avgs)) = 0; 
    zero_ratio = 1-num_nonzeros/size(x,2);
    ix = intersect(find(zero_ratio > 0),find(zero_ratio < 1));
    %ix = zero_ratio~=0 & zero_ratio~=1;% do not consider all zeros or all non-zeros
    avgsx = avgs(ix); zrs = zero_ratio(ix);
    % three situations
    if boundary_function == 1 %exponetial function
        % order avgl
        [avgso,I] = sort(avgsx); zrso = zrs(I);
        mavg = mean(avgso); mzr = mean(zrso);
        a0 = -log(mzr+1)/(mavg^2+eps);
        ft = fittype('exp(-a*x^2)');% exponetial function
        fitobject = fit(avgso,zrso,ft,'StartPoint',a0);
        p = predint(fitobject,avgso,0.95);% 95% confidece intervel
        bound_values = exp(-p(2)*zrs.^2);
    elseif boundary_function == 2 % simple piecewise function
        bound_values = zeros(length(ix),1);
        for k = 1:length(ix)
            % between 0.05
            id = intersect(find(zrs> zrs(k)-0.05), find(zrs < zrs(k)+0.05));
            % corresponding avg values
            va_id = avgsx(id);
            % if the rate high than 0.8 take min as upper bound
            if zrs(k) >= 0.8
                bound_values(k) = min(va_id);
            else
                bound_values(k) = max(va_id);
            end
        end
        
    elseif boundary_function == 3 % sophisticated piecewise function
        bound_values = zeros(length(ix),1);
        for k = 1:length(ix)
            % between 0.05
            id = abs(zrs-zrs(k)) < 0.05;
            % corresponding avg values
            va_id = avgsx(id);
            if zrs(k) >= 0.8
                bound_values(k) = min(va_id);
            elseif zrs(k) >= 0.6
                bound_values(k) = median(va_id(va_id<median(va_id)));%1/4
            elseif zrs(k)>= 0.4
                bound_values(k) = median(va_id(va_id>median(va_id)));% 3/4
            else
                bound_values(k) = max(va_id);
            end
        end
    end
    U = avgs;  U(ix) = bound_values; U_record{1,i} = repmat(U,1,size(x,2));     
end
%% impute by PBLR
Xi_record = cell(1,K); delta = 0;
parfor i = 1:K
    x = X_groups{1,i}; U = U_record{1,i};
    xi = MC_ADMM(x,delta,U,accelate);
    Xi_record{1,i} = xi;
end
% integrate to one matrix
Xi = M;
for i = 1:K
    idxCluster = group == i;
    Xi(:,idxCluster) = Xi_record{1,i};
end

