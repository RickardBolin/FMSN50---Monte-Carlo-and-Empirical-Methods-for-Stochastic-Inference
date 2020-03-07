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
    
steps = 5000;
accidents = zeros(d+1, 1);
burn_in = 3000;
jump = 1;
t_tracker = zeros(d+2, (steps-burn_in)/jump);
theta_tracker = zeros(d+1, (steps-burn_in)/jump);
lambda_tracker = zeros(d+1, (steps-burn_in)/jump);

rhos = 0.055;
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
            % Save the current state of the random walk every 50 steps to get 
            % relatively independent samples
            if(mod(step, jump) == 0 && step > burn_in)
                t_tracker(:,(step-burn_in)/jump) = t;
                theta_tracker(:,(step-burn_in)/jump) = cond_theta;
                lambda_tracker(:,(step-burn_in)/jump) = cond_lambda;
            end
            
            if(mod(step,floor(steps/100)) == 0 )
               q = step/steps;
               clc
               disp([num2str(rho_index + (psi_index - 1)*nrhos) '/' num2str(nrhos*npsis) ' ' num2str(100*q) '% |' char(ones(1,floor(50*q))*'=') char(ones(1, ceil(50 - 50*q))*' ') '|'])       
            end
        end
    end
    theta_on_psi(psi_index,:) = mean(theta_tracker,2);
    lambda_on_psi(psi_index,:) = mean(lambda_tracker,2);
    meant_on_psi(psi_index,:) = mean(t_tracker,2);
end
acceptance_rate = acceptance_rate./steps;
%% plots
close all

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
plot(t_tracker')
title(['Chains of ' num2str(d+1) ' breakpoints obtained from hybrid Metropolis-Hastings sampler'])
xlabel('Step') 
ylabel('Year') 
ylim([1600 2000])

deriv = zeros(d+1,1);
figure
plot(T,cumsum(T > 0), 'r')
hold on
title('Spline interpolation of lines with corresponding lambda as gradients')
xlabel('Year') 
ylabel('Accumulated number of accidents')
c_map = parula(12);
start = 0;
cond_lambda = mean(lambda_tracker,2);
t = mean(t_tracker,2);
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

% Calculate autocorrelation
figure
acf(t_tracker(6,:)', 550);
