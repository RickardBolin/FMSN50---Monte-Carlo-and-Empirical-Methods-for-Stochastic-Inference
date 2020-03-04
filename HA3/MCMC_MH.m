% Pseudo-code can be found in slide 5 in lecture 10:
% http://www.maths.lth.se/matstat/kurser/fmsn50masm11/2020/material/L10.pdf
function t = MCMC_MH(lambda, t, T)
rho = 0.01;

for i = 2:length(t)-1
    R = rho*(t(i+1)-t(i-1));
  
    % Keep drawing until we get an X that is inside the interval
    X = t(i) + R*(2*rand(1)-1);
    
    while(X < t(i-1) || X > t(i+1))
        X = t(i) + R*(2*rand(1)-1);
    end
    
    new_t = t;
    new_t(i) = X;
    alpha = min(1, f(lambda, new_t, T)/f(lambda, t, T));
        
    if rand(1) <= alpha
        t(i) = X;
    end
    
end


