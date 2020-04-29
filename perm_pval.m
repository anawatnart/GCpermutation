%==========================================================================
%
%	This function calculates p-value from the permutation distribution of
%	the GC matrix under the null hypothesis and the estimated GC matrix.
%
%	INPUTS
%           F           =   estimated GC matrix
%           Fdist       =   permutation samples of GC matrix
%	OUTPUTS
%           Fpval       =   p-value matrix corresponding to GC matrix
%
%========================================================================== 
function Fpval = perm_pval(F,Fdist)
[n,~,np] = size(Fdist);
Fpval = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            xdist = reshape(Fdist(i,j,1:np),1,np,1);
            
            I = sum(xdist >= F(i,j));
            Fpval(i,j) = (I)/(length(xdist));
        end
    end
end
end