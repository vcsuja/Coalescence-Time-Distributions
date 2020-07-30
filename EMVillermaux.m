function [mu,ratios] = EMVillermaux(data,K,itr)
% This function estimates the parameters of K Villermaux distributions that describe the data using the
% Estimation-Maximization algorithm. Here mu is the scale parameter 
% (Equivalent to lambda^(4/3) of the Weibull distribution)

%itr - Number of iterations. (Default:1000)
%K - Number of Villermaux distributions 

% Handle optional arguments
if nargin<3
    itr=1000;
end

%Prealocate the variables in the calculation:
mixRatio = zeros(K,itr+1);
muEstimate  = zeros(K,itr+1);
log_likelihood = zeros(1,itr+1);
ClusterDistance = zeros(K, length(data));
ClusterLikelihood = zeros(K, length(data));

% Initialize the variables in the calculation:
mixRatio(:,1)=1/K;
muEstimate(:,1) = logspace(log10(max(min(data),0.1)), log10(max(data)),K);
  
%EM algorithm:
for i =1:itr
    sumDistance = 0;
    % Compute expectation:
    for j =1:K
        ClusterDistance(j,:)=weibull_dist(data,muEstimate(j,i),4/3)+eps^6;
        sumDistance = sumDistance + mixRatio(j,i)*ClusterDistance(j,:);
    end
    log_likelihood(i)=sum(log(sumDistance));
    for j=1:K
        ClusterLikelihood(j,:) = (mixRatio(j,i)*ClusterDistance(j,:))./sumDistance;
        % Maximization Step
        muEstimate(j,i+1) = sum(ClusterLikelihood(j,:).*data.^(4/3))/sum(ClusterLikelihood(j,:));
        mixRatio(j,i+1)=sum(ClusterLikelihood(j,:))/length(data);
    end
end

%Return variables:
mu = muEstimate(:,end);
ratios = mixRatio(:,end);

% Function to calculate the Rayleigh distance:
function [y] = weibull_dist(x,mu,k)
%Rayleigh_DIST function for Rayleigh distribution
    y=(k/mu)*(x).^(k-1).*exp(-(x.^k)/mu);
end

end