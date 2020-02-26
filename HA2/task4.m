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
   X(k,:,:) = X(k-1,:,:);
   for i = 1:N
        [oneHotDir, nfree] = g(X(1:k,:,i));
        % Move in a valid direction
        X(k,:,i) = X(k,:,i) + oneHotDir;

        z = 1;
        % Is it self avoiding?
        if(length(unique(X(1:k,:,i),'row')) < k)
            z = 0;
        end
        % Set weight
        w(k,i) = w(k-1,i) * (z/(1/(nfree))); 
   end
   % Print progress
   if mod(k,5) == 0
      disp("Progress: Step " + k + "/" + (steps - 1) ) 
   end
end

% Calculate approximation of number of self avoiding walks
cn = mean(w(2:end,:),2);

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
