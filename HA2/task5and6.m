clear
close all

N = 10000; % Number of particles
steps = 10 + 1; 
dims = 2; % Dimensions
X = zeros(steps,dims,N);
w = zeros(steps,N);

% Set all starting weights to 1
w(1,:) = 1;

for k = 2:steps
   weight_distro = cumsum(w(k - 1,:)./sum(w(k - 1,:)));
   for i = 1:N
       
        % Resample
        parent_index = find((weight_distro > rand) == 1);
        parent_index = parent_index(1);
        X(1:k - 1,:,i) = X(1:k - 1,:,parent_index);        
       
        % Sample possible directions
        [oneHotDir, nfree] = g(X(1:k - 1,:,i));
        X(k,:,i) = X(k - 1,:,i) + oneHotDir;
        
        
        %1Sn, is it self avoiding?
        z = 1;
        if(length(unique(X(1:k,:,i),'row')) < k)
            z = 0;
        end
        
        % Set weight
        w(k,i) = (z/(1/(nfree))); 
   end 
   % Print progress
   if mod(k,5) == 0
      disp("Progress: Step " + k + "/" + (steps - 1)) 
   end
end

% Calculate approximation of number of self avoiding walks
cn = cumprod(mean(w(2:end,:),2));

%%
figure
plot(cn)
title('Approximate number of self avoiding random walks')
xlabel('Walk length')
ylabel('Approximate number of SAW') 

%%
% Find the number of surviving particles at each step
survivors = zeros(steps, 1);
for i = 1:steps
    survivors(i) = N - sum(length(find(w(i,:)==0)));
end

%% Regression 
Y = log(cn) + log(1:steps-1)';
X_reg = [ones(steps-1,1) (1:steps-1)' log(1:steps-1)'];

beta = X_reg\Y;
expbeta = exp(beta)';

A2_reg = expbeta(1);
mu2_reg = expbeta(2);
gamma2_reg = beta(3);

cn_reg = A2_reg*(mu2_reg.^(1:steps-1)).*((1:steps-1).^(gamma2_reg - 1));
