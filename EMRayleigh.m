function [mu,ratios] = EMRayleigh(data,K,itr)
% This function estimates the parameters of K Rayleigh distributions that describe the data using the
% Estimation-Maximization algorithm. Here mu is the scale parameter
% (wikipedia denotes it by sigma^2) of the distribution.

%itr - Number of iterations. (Default:1000)
%K - Number of Rayleigh distributions 

% Handle optional arguments
if nargin<3
    itr=1000;
end

% Prealocate the variables in the calculation:
mixRatio = zeros(K,itr+1);
muEstimate  = zeros(K,itr+1);
log_likelihood = zeros(1,itr+1);
ClusterDistance = zeros(K, length(data));
ClusterLikelihood = zeros(K, length(data));

% Initialize the variables in the calculation:
mixRatio(:,1)=1/K;
muEstimate(1:2,1)=[(min(data)+ max(data))/2,max(data)];
for i=3:K  %Might need to initialize the scale parameter better:
    muEstimate(i,1)=(i-2)* (min(data)+ max(data))/(K-1);
end

%EM algorithm:
for i =1:itr
    sumDistance = 0;
    % Compute expectation:
    for j =1:K
        ClusterDistance(j,:)=raleigh_dist(data,muEstimate(j,i))+eps^6;
        sumDistance = sumDistance + mixRatio(j,i)*ClusterDistance(j,:);
    end
    log_likelihood(i)=sum(log(sumDistance));
    for j=1:K
        ClusterLikelihood(j,:) = (mixRatio(j,i)*ClusterDistance(j,:))./sumDistance;
        % Maximization Step
        muEstimate(j,i+1) = 0.5*sum(ClusterLikelihood(j,:).*data.^2)/sum(ClusterLikelihood(j,:));
        mixRatio(j,i+1)=sum(ClusterLikelihood(j,:))/length(data);
    end
end

%Return variables:
mu = muEstimate(:,end);
ratios = mixRatio(:,end);

% Function to calculate the Rayleigh distance:
function [ y ] = raleigh_dist(x,mu)
%Rayleigh_DIST function for Rayleigh distribution
    y=(x/mu).*exp((-x.^2)/(2*mu));
end



end