function Z = SVD_thresholding(L,tau,accelate)
% accelate: 1 or 0. 1 represents adopting accelated algorithm:
% truncked SVD, false represents adopting full svd
if isequal(accelate, 0)
    [U,S,V] = svd(L);
    s = diag(S); r = sum(s>tau);
    St = diag(s(1:r)-tau); U = U(:,1:r); V = V(:,1:r);
else
    %     %% select k on columns
    %     if size(L,2) > 2000
    %         s = randperm(size(L,2));
    %         L0 = L(:,s(1:1000));
    %         [~,S0,~] = svd(L0); S0 = diag(S0);
    %         % compute the effective number
    %         k0 = min(1000/size(S0,1),size(S0,1)/1000);
    %         k1 = optimal_SVHT_coef(k0, 0)*median(S0);
    %         k1 = min(round(k1),min(size(L)));
    %     else
    %         k1 = min(size(L));
    %     end
    %     %% select k on rows
    %     if size(L,1) > 5000
    %         s = randperm(size(L,1));
    %         L0 = L(s(1:1000),:);
    %         [~,S0,~] = svd(L0); S0 = diag(S0);
    %         % compute the effective number
    %         k0 = min(1000/size(S0,2),size(S0,2)/1000);
    %         k2 = optimal_SVHT_coef(k0, 0)*median(S0);
    %         k2 = min(round(k2),min(size(L)));
    %     else
    %         k2 = min(size(L));
    %     end
    %     k = min(k1,k2);
    % SVD thresholding on L
    %     if k > 0.05*min(size(L))
    %         [U,S,V] = svd(L);  s = diag(S); r = sum(s>tau);
    %         St = diag(s(1:r)-tau);
    %         U = U(:,1:r); V = V(:,1:r);
    %     else
    %         [U,S,V] = lansvd(L,k,'L');
    %         St = diag(diag(S)-tau);
    %     end
    K = min(100,min(size(L)));
    [U,S,V] = rsvd(L,K);
    St = diag(diag(S)-tau);
end
Z = U*St*V';


