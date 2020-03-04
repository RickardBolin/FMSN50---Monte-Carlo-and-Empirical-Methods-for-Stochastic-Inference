clc
clear
load coal_mine_disasters.mat

% Number of breakpoints (not including start and end)
d = 5;

% Hyperparameter chosen by us
psi = 25;

% Initial breakpoints (together with start- and endpoint)
t = linspace(1658, 1980, d+2)';

% Draw initial values from the priors. With d breakpoints, we get d+1 intervals
theta = gamrnd(2*ones(d+1,1), ones(d+1,1)*1./psi);
lambda = gamrnd(2*ones(d+1,1), 1./theta);
cond_lambda = lambda; % For later

t_tracker = [];
accidents = zeros(d+1, 1);

for steps = 1:100000
    % Get number of accidents in each interval
    startpoints = t(1:end-1);
    endpoints = t(2:end);
    for i = 1:length(startpoints)
        accidents(i) = sum(length(T(find(T > startpoints(i) & T < endpoints(i)))));
    end

    cond_theta = gamrnd((2*d+2)*ones(d+1,1), (1/(psi+sum(cond_lambda)))*ones(d+1,1));
    cond_lambda = gamrnd((accidents+2), 1./((endpoints-startpoints)+cond_theta));
    t = MCMC_MH(cond_lambda, t, T);
    t_tracker = [t_tracker, t];
end

figure
hold on
for i = 2:d+1
    plot(t_tracker(i,:))
end
ylim([1658 2000])
