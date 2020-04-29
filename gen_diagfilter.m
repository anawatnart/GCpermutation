%========================================================================== 
%
%   This function generate a stable, invertible diagonal filer with
%   its ii element equal to
%                   Gii(z)  = pii(z)/qii(z)
%
%                              (z-z_1)(z-z_2)(...)(z-z_2p)
%                           =  ---------------------------
%                              (z-p_1)(z-p_2)(...)(z-p_nq)
%
%   where the order of the polynomial pii(z) and qii(z) are np and nq
%   respectively.
%
%   INPUTS
%           np      =	order of the numerator polynomial
%           nq      =   order of the denominator polynomial
%           N       =   filter matrix size
%           option  =   'same' for using the same poles on all entries
%                       'diff' for random a new set of poles for each entry
%                       no input will uses 'same' as the default
%   OUTPUTS
%           Gz      =   filter's transfer function matrix
%
%========================================================================== 
function Gz = gen_diagfilter(np,nq,N,option)
if nargin < 4
    option = 'same';
end
%   Choose 2*floor(np/4) poles and 2*floor(nq/4) zeros be complex
npc = floor(np/4); nqc = floor(nq/4);
%   Generate real and complex zeros
Zr = (2*randi([0,1],np-2*npc,N)-1).*(0.05 + 0.75*rand(np-2*npc,N));
Zc = (0.05 + 0.75*rand(npc,N)).*exp(1j*(pi/10+(8*pi/10)*rand(npc,N)));
%   Generate real and complex poles
if strcmp(option,'diff')
    Pr = (2*randi([0,1],nq-2*nqc,N)-1).*(0.05 + 0.75*rand(nq-2*nqc,N));
    Pc = (0.05 + 0.75*rand(nqc,N)).*exp(1j*(pi/10+(8*pi/10)*rand(nqc,N)));
elseif strcmp(option,'same')
    Pr = (2*randi([0,1],nq-2*nqc,1)-1).*(0.05 + 0.75*rand(nq-2*nqc,1));
    Pc = (0.05 + 0.75*rand(nqc,1)).*exp(1j*(pi/10+(8*pi/10)*rand(nqc,1)));
end
num = repmat({[0]},N,N); den = repmat({[1]},N,N);
Gz = tf(num,den,1);
for i = 1:N
    Z = [Zr(:,i); Zc(:,i); conj(Zc(:,i))];
    if strcmp(option,'diff')
        P = [Pr(:,i); Pc(:,i); conj(Pc(:,i))];
    elseif strcmp(option,'same')
        P = [Pr(:,1); Pc(:,1); conj(Pc(:,1))];
    end
    ZP = zpk(Z,P,1,1);
    Gz(i,i) = tf(ZP);
end
end

