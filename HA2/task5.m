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
   weight_distro = cumsum(w(k - 1,:)./sum(w(k - 1,:)));
   for i = 1:N
       
        %Resample
        parent_index = find((weight_distro > rand) == 1);
        parent_index = parent_index(1);
        X(1:k - 1,:,i) = X(1:k - 1,:,parent_index);        
       
        %Sample possible directions
        [oneHotDir, nfree] = g(X(1:k - 1,:,i));
        X(k,:,i) = X(k - 1,:,i) + oneHotDir;
        
        
        %1Sn, is it self avoiding?
        z = 1;
        if(length(unique(X(1:k,:,i),'row')) < k)
            z = 0;
        end
        
        %set weight
        w(k,i) = (z/(1/(nfree))); 
   end 
end

probs = cumprod(mean(w,2));

plot(probs)

