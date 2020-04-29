%==========================================================================
%
%   INPUTS
%       Fest    =   The estimated F
%    	Fhat    =   must have size(Fhat) = [m,m,k,n] where m is the
%                   dimension of the system, n is the number of F used in
%                   one averaging and k is the number of averaged samples. 
%                   The diagonal entries of Fhat(:,:,i,j) are always
%                   ignored and never used in calculation.
%       auto    =   1 for automatically choosing the order of GMM
%                   0 for choosing the order of GMM manually
%       maxmode =   number of the maximum mode of the estimated GMM. All
%                   GMM with all mode up to maxmode is estimated for
%                   selection later.
%
%   OUTPUTS
%       F_gmm   =   Granger causality matrix via GMM methods
%
%
% Developper: Nattaporn Plubin, Anawat Nartkulpat
% 
% N. Plub-in and J. Songsiri, "State-space Model of EEG Time Series
% for Learning Granger Causality of Brain Connectivity"
%
%==========================================================================
function [Fgmm,GMM] = gmm_gc(Fest,Fhat,auto,maxmode)
    [m,~,n,~] = size(Fhat);
    Fbar = mean(Fhat,4);
    vecFbar_train = reshape(Fbar,m*m,n);
    vecFbar_train(1:m+1:m*m,:) = [];     % Remove diagonal entries
    vecFbar_train = reshape(vecFbar_train,size(vecFbar_train,1)*size(vecFbar_train,2),1);
    Fmin = min(vecFbar_train,[],'all');
    Fmax = max(vecFbar_train,[],'all');

    % GMM estimation
    GMModels = cell(1,maxmode);
    options = statset('MaxIter',1000);  % set the number of iteration

    for k = 1:maxmode
        % initial point by k-means++
        S = 'plus';

        % Replicate with 10 initial values (algorithm from 'Start') then return the model with the largest likelihood.
        num_replicate = 10;
        GMModels{k} = fitgmdist(vecFbar_train,k,'Options',options,'Start',S,'Replicates',num_replicate);

        % collect AIC BIC logL
        AIC(k)= GMModels{k}.AIC;
        BIC(k)= GMModels{k}.BIC;
    end
    if auto == 0
        figure; subplot(211);
        histogram(vecFbar_train,100,'Normalization','pdf'); hold on;
        axis([Fmin Fmax 0 150]);
        legend('Histogram of Fij'); 
        subplot(212); hold on;
        for k = 1:maxmode
            y=pdf(GMModels{k},[Fmin:0.005:Fmax]');
            plot([Fmin:0.005:Fmax]',y,'LineWidth',1.25,'DisplayName',['order ',num2str(k)]);
            axis([Fmin Fmax 0 150]);
        end
        
        %======================================================
        figure;
        histogram(vecFbar_train,100,'Normalization','pdf'); hold on;
        x = [0:0.001:0.5]';
        y = pdf(GMModels{3},x);
        plot(x,y,'LineWidth',1.25);
        xlabel('Magnitude of F_{ij}'); 
        legend('Histogram','pdf of fitted GMM');
        title('An example of fitting a GMM to data with 3 modes');
        %======================================================
        
        legend show;
        figure; subplot(211); plot(AIC); ylabel('AIC');
        subplot(212); plot(BIC); ylabel('BIC');
        n_chosen = input('Enter the number of GMM mode: ');
    elseif auto == 1
        % auto select model order
        [~,AIC_chosen] = min(AIC); % select model from AIC score
        [~,BIC_chosen] = min(BIC); % select model from BIC score
        n_chosen = max(round(mean([AIC_chosen,BIC_chosen])));
    elseif auto == 2
        % auto select model order
        n_chosen = 2;
    end 

    %% GMM selectionF_result
    GMM = GMModels{n_chosen};
    [~,sort_ind_rBIC] = sort(GMM.mu);
    
    % clustering by GMM
    [vecCluster_result] = cluster(GMM,reshape(Fest,m*m,1));
    F_sig = zeros(size(Fest,1),size(Fest,2),size(Fest,3));
    F_sig(logical(ones(size(Fest,1),size(Fest,2),size(Fest,3)))) = (vecCluster_result~=sort_ind_rBIC(1));

    % collect F significance result
    Fgmm = F_sig;
end

