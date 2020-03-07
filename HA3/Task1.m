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
    
steps = 10000;
t_tracker = zeros(d+2, steps);
accidents = zeros(d+1, 1);


rhos = linspace(0,0.1,50);
nrhos = length(rhos);

psis = 25;%linspace(0,50,50);
npsis = length(psis);

acceptance_rate = zeros(nrhos,npsis);
theta_on_psi = zeros(npsis, d + 1);
lambda_on_psi = zeros(npsis, d + 1);
meant_on_psi = zeros(npsis, d + 2);


for psi_index = 1:npsis
    psi = psis(psi_index);
    for rho_index = 1:nrhos
        rho = rhos(rho_index);
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
            [t, accepted] = MCMC_MH(rho ,cond_lambda, t, T);

            acceptance_rate(rho_index,psi_index) = acceptance_rate(rho_index,psi_index) + accepted;
            % Save the current state of the random walk
            t_tracker(:,step) = t;

            if(mod(step,floor(steps/100)) == 0)
               q = step/steps;
               clc
               disp([num2str(rho_index + (psi_index - 1)*nrhos) '/' num2str(nrhos*npsis) ' ' num2str(100*q) '% |' char(ones(1,floor(50*q))*'=') char(ones(1, ceil(50 - 50*q))*' ') '|'])       
            end
        end
    end
    theta_on_psi(psi_index,:) = cond_theta;
    lambda_on_psi(psi_index,:) = cond_lambda;
    meant_on_psi(psi_index,:) = mean(t_tracker,2);
end
acceptance_rate = acceptance_rate./steps;
%% plots


figure
plot(psis, mean(theta_on_psi,2))
title('Mean theta dependent on Psi')
ylabel('Mean theta')
xlabel('Psi')

figure
plot(psis, lambda_on_psi)
title('Lambda parameters dependent on Psi')
ylabel('Lambda')
xlabel('Psi')


figure
plot(psis, meant_on_psi)
title('Mean t parameters dependent on Psi')
ylabel('t')
xlabel('Psi')


%plot acceptance rate
figure
scatter(rhos, acceptance_rate(:,1))
title('Acceptance rate dependent on rho')
ylabel('Acceptance rate')
xlabel('Rho')

figure
scatter(psis, acceptance_rate(1,:))
title('Acceptance rate dependent on psi')
ylabel('Acceptance rate')
xlabel('Psi')


% Plot the random walks
figure
hold on
for i = 1:d+2
    plot(t_tracker(i,:))
end
hold off
ylim([1600 2000])

deriv = zeros(d+1,1);
figure
plot(T,cumsum(T > 0))
hold on

start = 0;
for i = 1:d+1
    y = find(T > t(i) & T < t(i+1));
    y1 = y(1);
    y2 = y(end);
    deriv(i) = (y2-y1)/(t(i+1)-t(i));
    line = cond_lambda(i)*linspace(0,t(i+1)-t(i)) + start;
    start = line(end);
    plot(linspace(t(i), t(i+1)), line)
end


