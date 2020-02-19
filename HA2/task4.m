
clear
close all

N = 1000; %nbr of particles
steps = 20; 
dims = 2; %dimensions
X = zeros(steps,dims,N);
w = zeros(steps,N);



for i = 1:N
    for d = 1:dims
         X(1,d,i) = 0;
    end
    w(1,i) = 1;
end


for k = 2:steps
   X(k,:,:) = X(k-1,:,:);
   for i = 1:N
        [oneHotDir, nfree] = g(X(1:k,:,i));
        X(k,:,i) = X(k,:,i) + oneHotDir;

        z = 1;
        if(length(unique(X(1:k,:,i),'row')) < k)
            z = 0;
        end
        w(k,i) = w(k-1,i) * (z/(1/(nfree))); 
   end
end

probs = mean(w,2);

plot(probs)

