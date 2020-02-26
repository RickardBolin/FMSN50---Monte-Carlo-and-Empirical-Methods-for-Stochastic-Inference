clear
close all

N = 1000; %nbr of particles
steps = 10; 
dims = 5; %dimensions
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
       % if(length(unique(X(1:k,:,i),'row')) < k)
       %     z = 0;
       % end
        
        %set weight
        w(k,i) = (z/(1/(nfree))); 
   end 
   if mod(k,5) == 0
      k 
   end
end

cn = cumprod(mean(w,2));
%%
mu2lim = nthroot(cn', 1:steps);
mu2 = mu2lim(end);
gamma2 = 43/32;
mu2nngamma2minus1 = (mu2.^(1:steps)).*((1:steps).^(gamma2-1));

A2 = cn./mu2nngamma2minus1';
%%
Y = log(cn) + log(1:steps)';
X_reg = [ones(steps,1) (1:steps)' log(1:steps)'];

beta = X_reg\Y;

beta = (X_reg'*X_reg)\X_reg'*Y;


expbeta = exp(beta)';

A_reg = expbeta(1);
mu_reg = expbeta(2);
gamma_reg = 1;%beta(3);

cn_reg = A_reg*(mu_reg.^(1:steps)).*((1:steps).^(gamma_reg - 1));

e = cn_reg' - cn;

mu_theo = 2*dims - 1 - 1/(2*dims) - 3/((2*dims)^2) - 16/((2*dims)^2);

