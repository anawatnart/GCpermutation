function [est_model,sys_n4sid,sys_pem,fit_B,fit_noB,fit_pem] = ssiden_noinput(y,n)
%
%  SSIDEN_NOINPUT estimates state-space model equations without deterministic input term using n4sid.
% 
%  MATLAB n4sid uses the following state-space equation
%
%       x(t+1)  = Ax(t) + Bu(t)+ Ke(t)        
%       y(t)    = Cx(t) + Du(t) + e(t)          
%   
% however, when we estimate time series data, u(t) is zero, hence B is forced to be zero
%
%    We perform an initial guess by "n4sid" estimation with enforce stability option 
%   then set structure of model and refine estimated parameters by using "pem"
%
% INPUT PARAMETERS:
% 	y           =   output measurement of size "DIMENSION  x TIME POINTS"
%  	n           =   dimension of state variable in the estimated model
%
% OUTPUT PARAMETERS:
%  	est_model   =   estimated model in MATLAB format based on the equation
%
%	        z(t+1)      =   Az(t) + w(t)
%	        y(t)        =   Hz(t) + e(t)
%
% Developper: Nattaporn Plubin, Anawat, Jitkomut 
% 
% N. Plub-in and J. Songsiri, "State-space Model of EEG Time Series
% for Learning Granger Causality of Brain Connectivity"
%
%%========================================================================== 
    [m,Num] = size(y);

    %------------------------------ Initial Guess ------------------------- 
    id_y = iddata(y',zeros(Num,1));     % set iddata
    n4opt =  n4sidOptions('EnforceStability',true);
    tic;
    sys_n4sid = n4sid(id_y,n,n4opt);    % subspace identification w/ n4sid algorithm
    toc;
%     sys_n4sid = n4sid(id_y,n);
    
    %---------------- Set Structure for estimated model ------------------- 
    [n,~] = size(sys_n4sid.A);

    % A C K x0 set to be estimated parameters
    As(1:n,1:n) = NaN;
    Cs(1:m,1:n) = NaN;
    Ks(1:n,1:m) = NaN;
    x0s(1:n,1) = NaN;

    % B D set to be zeros
    Bs = zeros(n,1);
    Ds = zeros(m,1);

    sys_init = sys_n4sid;                	% initial guess for pem
    setstruc(sys_init,As,Bs,Cs,Ds,Ks,x0s); 	% set initial guess structure
    set(sys_init,'Ts',1);                  	% set initial guess to be discrete

    % pem estimation
    tic;
    sys_pem = pem(id_y,sys_init);
    toc;
        
    % collect estimated model
    est_model.A(:,:) = sys_pem.A;
    est_model.H(:,:) = sys_pem.C;
    est_model.S(:,:) = sys_pem.K*sys_pem.NoiseVariance;
    est_model.W(:,:) = sys_pem.K*sys_pem.NoiseVariance*sys_pem.K';
    est_model.E(:,:) = sys_pem.NoiseVariance;
    est_model.K(:,:) = sys_pem.K;
    
    figure;
    [~,fit_B,~] = compare(id_y,sys_n4sid,1);
    compare(id_y,sys_n4sid,1);
    figure;
    sys_n4sid.B = zeros(size(sys_n4sid.B));
    [~,fit_noB,~] = compare(id_y,sys_n4sid,1);
    compare(id_y,sys_n4sid,1);
    figure;
    [~,fit_pem,~] = compare(id_y,sys_pem,1);
    compare(id_y,sys_pem,1);
%     figure;
%     predict(sys_pem,id_y);    figure;
%     predict(sys_pem,id_y,inf);
% figure;
% plot(y'); hold on;
% YP=predict(sys_n4sid,id_y,1); 
% figure;
% compare(id_y,YP);
end

