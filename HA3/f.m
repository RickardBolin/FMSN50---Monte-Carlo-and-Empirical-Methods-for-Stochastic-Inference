% f(t|theta, lambda, T) 
function fX = f(lambda,t,T)
% Get number of accidents in each interval
startpoints = t(1:end-1);
endpoints = t(2:end);
accidents = zeros(length(t)-1, 1);
for j = 1:length(startpoints)
    accidents(j) = sum(length(T(T > startpoints(j) & T < endpoints(j))));
end

% Calculate X_{k+1} 
fX = exp(sum(log(lambda).*accidents + log(t(2:end)-t(1:end-1)) - lambda.*(t(2:end)-t(1:end-1))));
end

