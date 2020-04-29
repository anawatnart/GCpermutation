function M = gen_sparseVAR(PARAMETER,opts,varargin)
% gen_VAR_param generates multiple sparse VAR parameters [A_1, ..., A_p]
% y_k(t) = A_k_1*y_k(t-1) + ... + A_k_p*y_k(t-p) + u(t)
% Where 'k' is k th model
% Input parameters are
% 'PARAMETER' = [n p k]
% 'opts' : a struct that contains field
%          '.diffdensity : differential level on similar structure, all
%          similar model have same differential level
%          '.density' : nonzero density of VAR coefficient
%          '.type' : 'common', 'similar', 'differential' or manually
%          specify Kx2 vector where type(i,:) = [a,b] refers to
%           index (k=i) has relation type 'b' with index (k=j) (j<i)
%                    'b' = 2, 1 ,0 if no relation,  similar edges, common sparsity
%                    respectively
% The output are
% 'M' : M contains k struct with field
%        '.ind_nz' : nonzero indices
%        '.A' : VAR parameter with size [n,n,p,k]
%        '.info' : contains cells ind_nz
%% Function init
n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);
M = struct();
M.info.dim = [n,p,K];
if numel(varargin) > 0
    DMAT = varargin{1};
end
if strcmp(opts.type,'common')
    tmp = [1 2;repmat([1 0],K-1,1)];
elseif strcmp(opts.type,'similar')
    K=K+1;
    tmp = [1 2;repmat([1 1],K-1,1)];
    
elseif strcmp(opts.type,'differential')
    tmp = [1 2;repmat([1 1],K-1,1)];
else
    tmp = opts.type;
end
M.A = zeros(n,n,p,K);
M.info.density = opts.density;
M.info.type = opts.type;

if isfield(opts,'diffdensity')
    M.info.diffdensity = opts.diffdensity;
else
    opts.diffdensity = 0;
    M.info.diffdensity = 0;
end


for kk=1:K
    if (kk==1) && ~isequal(tmp(1,:),[1,2])
        error('First model must be different, set opts.type(1) = [1,2]')
    end
    if (tmp(kk,1)>=kk) && (kk~=1)
        error('parameter "a" in opts.type is in incorrect order, must be the index that model that is already generated')
    end
    if exist('DMAT','var')
        [M.info.ind_nz{kk},M.A(:,:,:,kk)] = subfunc(n,p,opts.density,opts.diffdensity,1,DMAT(:,:,:,kk));
    else
        [M.info.ind_nz{kk},M.A(:,:,:,kk)] = subfunc(n,p,opts.density,opts.diffdensity,tmp(kk,2),M.A(:,:,:,tmp(kk,1)));

    end
end
        if strcmp(opts.type,'similar')
            M.A = M.A(:,:,:,2:end);
        end
end


function [ind_nz,A] = subfunc(n,p,density,diffdensity,opts,S) %No. Groups, Similarity
% This function generates a sparse vector autoregressive model
% y(t) = A1*y(t-1) + A2*y(t-2) + ... + Ap*y(t-p) + u(t)
%
% 'A' represents AR coefficients A1,A2,...,Ap and is stored as a p-dimensional array
%
% The input arguments are
% 'n': dimension of AR coefficient matrices
% 'p': AR model order
% 'noise_var': variance of u(t) (noise)
% 'density': the fraction of nonzero entries in AR coefficients
% 'Num': number of data points in time series
% 'S' : Original AR coefficients
% 'opts' : 1 will return common zero pattern with S,0 will return
%  similar coefficients with S
%
% The AR coefficients are sparse with a common sparsity pattern. The
% indices of nonzero entries are saved in 'ind_nz'.
% if p = 0, 'y' is simply a random variable. In this case, A is the
% covariance matrix of u with sparse inverse.
%% Common
if (opts == 2)
    S = sprand(n,n,density)+eye(n);
    opts = 0;
end

%% Static case
if (p==0)
    S = sparse(2*eye(n)+sign(sprandsym(n,density)));
    [i,j]=find(S);
    S = S+sparse(ceil(max(0,-min(eig(S))))*eye(n));
    A = S\eye(n); % covariance matrix with sparse inverse
    %             R = chol(phi);
    %             y = R'*randn(n,Num); % y reduces to a random variable with covariance 'phi'
    ind_nz = sub2ind([n n],i,j);
    %             figure;plot_spy(ind_nz,n,'image');
    %             title('correct sparsity');
    return;
end

%% Randomize AR coefficients
if (opts==0) % Generate common structure with given matrix 'S'
    MAX_EIG = 1;
    diag_ind = find(eye(n));
    k = length(diag_ind);
    diag_ind3D = kron(n^2*(0:p-1)',ones(k,1))+kron(ones(p,1),diag_ind);
    off_ind3D = (1-repmat(eye(n),1,1,p));
    A = zeros(n,n,p);
    ii = 0;
    while MAX_EIG
        ii = ii+1;
        %                 disp(ii)
        for k=1:p
            A(:,:,k) = 0.1*sprandn(S(:,:,1));
        end
        
        poles = -0.7+2*0.7*rand(n,p); % make the poles inside the unit circle
        characeq = zeros(n,p+1);
        for jj=1:n
            characeq(jj,:) = poly(poles(jj,:)); % each row is [1 -a1 -a2 ... -ap]
        end
        aux = -characeq(:,2:end);
        A(diag_ind3D) = aux(:); % replace the diagonal entries with stable polynomial
        %                 size(A)
        %                 A(find(off_ind3D));
        A(((abs(A)<0.2)&(A~=0))&(off_ind3D))=0.4.*sign(A(((abs(A)<0.2)&(A~=0))&(off_ind3D)));
        AA = zeros(n,n*p);
        for k=1:p
            AA(1:n,k*n-(n-1):k*n) = A(:,:,k);
        end
        AA = sparse([AA ; [eye(n*(p-1)) zeros(n*(p-1),n)]]);
        %                 disp(max(abs(eigs(AA))))
        if (max(abs(eigs(AA))) < 1)
            MAX_EIG = 0;
        end
    end
elseif (opts==1) % Generate Similar model
    MAX_EIG = 1;
    off_ind3D = (1-repmat(eye(n),1,1,p));
    diag_ind = (1:n+1:n^2)';
    A = zeros(n,n,p);
    R = sprand(n,n,diffdensity);
    R(diag_ind) = 0;
    R(S(:,:,1)~=0) = 0;
    ii = 0;
    while MAX_EIG
        ii = ii+1;
        for k=1:p
            A(:,:,k) = S(:,:,k)+ 0.2*sprandn(R); % add some off-diagonal entries
        end
        A(((abs(A)<0.2)&(A~=0))&(off_ind3D))=0.2.*sign(A(((abs(A)<0.2)&(A~=0))&(off_ind3D)));
        AA = zeros(n,n*p);
        for k=1:p
            AA(1:n,k*n-(n-1):k*n) = A(:,:,k);
        end
        AA = sparse([AA ; [speye(n*(p-1)) zeros(n*(p-1),n)]]);
        if max(abs(eigs(AA))) < 1
            MAX_EIG = 0;
        end
    end
end
ind_nz = find(A(:,:,1));
end