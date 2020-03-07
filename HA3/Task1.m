clc
close all
clear
addpath('Data')
load coal_mine_disasters.mat

% Number of breakpoints (not including start and end)
d = 5;
% Hyperparameter psi, chosen by us
psi = 25;
% Initial breakpoints (together with start- and endpoint)
t = linspace(1658, 1980, d+2)';
% Initialize lambda
cond_lambda = 5;
    
steps = 100000;
t_tracker = zeros(d+2, steps);
accidents = zeros(d+1, 1);
for step = 1:steps
    % Get number of accidents in each interval
    startpoints = t(1:end-1);
    endpoints = t(2:end);
    for i = 1:length(startpoints)
        accidents(i) = sum(length(T(T > startpoints(i) & T < endpoints(i))));
    end
    
    % Draw (theta|lambda, t, T) 
    cond_theta = gamrnd((2*d+2)*ones(d+1,1), (1/(psi+sum(cond_lambda)))*ones(d+1,1));
    % Draw (lambda|theta, t, T)  
    cond_lambda = gamrnd((accidents+2), 1./((endpoints-startpoints)+cond_theta));
    % Calculate (t|theta, lambda, T) with a Metropolis-Hastings sampler  
    t = MCMC_MH(cond_lambda, t, T);
    % Save the current state of the random walk
    t_tracker(:,step) = t;
    
    if(st
    
end

% Plot the random walks
figure
hold on
for i = 1:d+2
    plot(t_tracker(i,:))
end
ylim([1600 2000])
