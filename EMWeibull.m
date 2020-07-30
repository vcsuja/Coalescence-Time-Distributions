function [mu,k,ratios] = EMWeibull(data,K,itr)
% This function estimates the parameters of K Weibull distributions that describe the data using the
% Estimation-Maximization algorithm. Here mu is the scale parameter 
% (wikipedia denotes it by lambda^k) of the distribution and k is the shape parameter.

%itr - Number of iterations. (Default:1000)
%K - Number of Weibull distributions 

% Handle optional arguments
if nargin<3
    itr=1000;
end

%Prealocate the variables in the calculation:
mixRatio = zeros(K,itr+1);
muEstimate  = zeros(K,itr+1);
kEstimate  = zeros(K,itr+1);
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
        ClusterDistance(j,:)=weibull_dist(data,muEstimate(j,i),kEstimate(j,i))+eps^6;
        sumDistance = sumDistance + mixRatio(j,i)*ClusterDistance(j,:);
    end
    log_likelihood(i)=sum(log(sumDistance));
    for j=1:K
        ClusterLikelihood(j,:) = (mixRatio(j,i)*ClusterDistance(j,:))./sumDistance;
        % Maximization Step
        % Numerically obtain kEstimate(j,i+1)               
        if i == 1
            kEstimate(j,1)=muEstimate(j,1)/sqrt(var(data.*ClusterLikelihood(j,:)));  % The guess for the shape parameter
        end
        kEstimate(j,i+1) =  abs(fzero(@(k) calculatek(k,data,ClusterLikelihood(j,:)), kEstimate(j,i)));
        muEstimate(j,i+1) = sum(ClusterLikelihood(j,:).*data.^( kEstimate(j,i+1)))/sum(ClusterLikelihood(j,:));
        mixRatio(j,i+1)=sum(ClusterLikelihood(j,:))/length(data);
    end
end

%Return variables:
mu = muEstimate(:,end);
k =  kEstimate(:,end);
ratios = mixRatio(:,end);

% Function to calculate the Rayleigh distance:
function [y] = weibull_dist(x,mu,k)
%Rayleigh_DIST function for Rayleigh distribution
    y=(k/mu)*(x).^(k-1).*exp(-(x.^k)/mu);
end

function y = calculatek(k,x,p)
y = k*sum(p.*x.^k.*log(x))./sum(p.*x.^k)-1 -k*sum(p.*log(x))/sum(p);
end

end