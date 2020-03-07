% Pseudo-code can be found in slide 5 in lecture 10:
% http://www.maths.lth.se/matstat/kurser/fmsn50masm11/2020/material/L10.pdf
function [t, accepted] = MCMC_MH(rho,lambda, t, T)
accepted = 0;
% Selecting one t_i at random
i = randi(length(t) - 2) + 1;
%for i = 2:length(t)-1
    % Calculate the maximum step size based on the interval width
    
    R = rho*(t(i+1)-t(i-1));
  
    X = -1;
    % Keep drawing until we get an X that is inside the interval
    while(X < t(i-1) || X > t(i+1))
        X = t(i) + R*(2*rand(1)-1);
    end
    
    % Create a temporary new break point list with the new value
    % from the random walk 
    new_t = t;
    new_t(i) = X;
    
    % Check if we should accept or reject the new value
    % f := f(t|theta, lambda, T) 
    if rand(1) <= f(lambda, new_t, T)/f(lambda, t, T)
        t = new_t;
        accepted = 1;
    end
%end
end
