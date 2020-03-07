clc
close all
clear
addpath('Data')
load coal_mine_disasters.mat

% Number of breakpoints (not including start and end)
d = 5;
% Hyperparameter psi, chosen by us
psi = 15;
% Initial breakpoints (together with start- and endpoint)
t = linspace(1658, 1980, d+2)';
% Initialize lambda
cond_lambda = 5;
    
steps = 1e5;
burn_in = 3000;
jump = 1;
t_tracker = zeros(d+2, (steps-burn_in)/jump);

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
    % Save the current state of the random walk every 50 steps to get 
    % relatively independent samples
    if(mod(step, jump) == 0 && step > burn_in)
        t_tracker(:,(step-burn_in)/jump) = t;
    end
    
    if(mod(step,floor(steps/100)) == 0)
       clc
       disp([num2str(100*step/steps) '% |' char(ones(1,floor(50*step/steps))*'=') char(ones(1, ceil(50 - 50*step/steps))*' ') '|'])       
    end
    
end
%%
% Plot the random walks
figure
hold on
title(['Chains of ' num2str(d+1) ' breakpoints obtained from hybrid Metropolis-Hastings sampler'])
xlabel('Step') 
ylabel('Year') 
for i = 1:d+2
    plot(t_tracker(i,:))
end
ylim([1600 2000])
%%
deriv = zeros(d+1,1);
figure
plot(T,cumsum(T > 0), 'r')
hold on
title('Spline interpolation of lines with corresponding lambda as gradients')
xlabel('Year') 
ylabel('Accumulated number of accidents')
c_map = parula(12);
start = 0;
for i = 1:d+1
    plot([t(i) (t(i))],[0 800], 'Color', c_map((mod(i,2)+1)*3,:))
    y = find(T > t(i) & T < t(i+1));
    y1 = y(1);
    y2 = y(end);
    deriv(i) = (y2-y1)/(t(i+1)-t(i));
    line = cond_lambda(i)*linspace(0,t(i+1)-t(i)) + start;
    start = line(end);
    plot(linspace(t(i), t(i+1)), line, 'Color', c_map((mod(i,2)+1)*3,:))
end
    plot([1980 1980],[0 800])

%% Calculate autocorrelation
figure
acf(t_tracker(6,:)', 550);
