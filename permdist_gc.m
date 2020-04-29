%==========================================================================
%
%	This function calculated a permutation distribution of GC matrix from
%	the input time series data by partitioning the data into windows with
%	length wsize, then the permutation of these windows are done for each
%	data channel seperately. For channel j, we obtain the distribution for
%	Fij for i = 1,...,n
%
%	INPUTS
%           data        =   input time series data
%           wsize       =   size of windows for permutation
%           nperm       =   number of permutation performed
%           nstate      =   number of states for subspace identification
%	OUTPUTS
%           permdist    =   permutation samples of GC matrix
%
%   Developper: Anawat Nartkulpat
%
%========================================================================== 
function permdist = permdist_gc(data,wsize,nperm,nstate)
[nc,nt] = size(data); nw = floor(nt/wsize);
data = data(:,1:nw*wsize);
permdist = zeros(nc,nc,nperm);
for j = 1:nc
    parfor i = 1:nperm
        datagc = data;
        dataj = reshape(data(j,:),1,wsize,nw);
        dataj = dataj(:,:,randperm(nw));
        dataj = reshape(dataj,1,nw*wsize);
        datagc(j,:) = dataj;

        [A,~,C,~,K,R] = subid(datagc,[],ceil(2*nstate/nc),nstate,[],'CVA',1);
        F = calgcss(A,C,K*R*K',R,K*R);
        
        permdist(:,j,i) = F(:,j);
    end
end
end

