function X = MC_ADMM(M,delta,U,accelate,varargin)
%imputation matrix with ADMM method
% agrmin ||X||_*, 
% s.t. ||X_omega-M||_F^2 <= delta
%      0 <= X_omega' <= U 
% <Inputs>:
% M:     The observed data matrix
% delta: Parameter for relaxation model
% U: the upper boundary
% accelate: 1 or 0. 1 represents adopting accelated algorithm:
% truncked SVD, false represents adopting full svd

%        (Below are optional arguments: can be set by providing name-value pairs)
%        X_INIT:     The initial value of decision matrix X       
%        Z_INIT:     The initial value of Lagrange multiplier
%        ITER_MAX:   The maximum iterations
%        TOL:        The predefined presision


% <Outputs>:
% X :     Record the decision matrix of each iteration
%% initialization
[m,n] = size(M); Omega = find(M == 0); X_init = zeros(m,n); Z_init = zeros(m,n); iter_max = 200; tol = 10^(-6); 
if (rem(length(varargin),2) == 1)
    error('Error:Optional parameters should always go by pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'X_INIT',              X_init = varargin{i+1};
            case 'Z_INIT',              Z_init = varargin{i+1};
            case 'ITER_MAX',            iter_max = varargin{i+1};
            case 'TOL',                 tol = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end

X = X_init;  Z = Z_init; 
%beta = 0.3/(lansvd(M,1,'L'));
%r_s = length(Omega) / (m * n); r = 1.1 + 2.5 * r_s;
r = 1.6; beta = 2.5/sqrt(m*n);
stop_value = zeros(1,iter_max); 
for iter = 1:iter_max
    % update Y    
    B = X-1/beta*Z; 
    % project to the range
    L = zeros(size(B));
    d1 = B-L; index1 = d1 < 0; B(index1) = L(index1);
    d2 = U-B; index2 = d2 < 0; B(index2) = U(index2);   
    b_m = B-M; b_m(Omega) = 0; para_y = min(delta/norm(b_m,'fro'),1);
    Y = (para_y - 1)*b_m + B; 
    A = Y + 1/beta*Z;  
    %% update X
    X = SVD_thresholding(A,1/beta,accelate); 
    % update Z
    Z = Z-r*beta*(X-Y); 
    % stop control
    if iter > 1
        stop_value(iter) = norm(X-X1,'fro')/norm(X1,'fro');
        if stop_value(iter) < tol
            break,
        end
    end
    X1 = X;
    %% update beta
    %beta = beta*r;
end
X(X < 10^(-4)) = 0;
