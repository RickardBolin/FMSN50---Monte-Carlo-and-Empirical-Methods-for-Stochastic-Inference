clear
clc
close all

load powercurve_V112.mat

lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k =      [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];

avgLambda = mean(lambda);
avgK = mean(k);

avgWind = 8.1;

x = 0:0.1:40;

%%




%  Tror vi missat grejer h?r! Vi m?ste v?l ha med 0-3 och >25 ocks? p?
%  n?got s?tt?



weibulDist = wblpdf(x,lambda(1),k(1));

figure
plot(x, weibulDist)
figure
plot(x, P(x))
N = 100000;
powers = zeros(N,12);
means = zeros(1,12);
stds = zeros(1,12);
for i = 1:12
    randomWeibuls = sort(wblrnd(lambda(i), k(i), N, 1)); % Sorting to get a nice plot
    powers(:,i) = P(randomWeibuls);
    means(i) = mean(powers(:,i));
    stds(i) = std(powers(:,i));
end

figure
scatter(randomWeibuls, powers(:,1))

figure
scatter(1:12, means)
figure
scatter(1:12, stds)

%%
confIntervals = [means - 1.96*(stds/sqrt(N)); means + 1.96*(stds/sqrt(N))];
width = abs(confIntervals(1,:)-confIntervals(2,:));

figure
scatter(1:12, width)

% Konfidensintervallet minskar n?r samples ?kar. Kanske rimligt, men bara
% till en viss gr?ns antar jag?

%%

a = 3;
b = 25;
N = 1000000;
U = rand(1,N);

powers = zeros(N,12);
means = zeros(1,12);
stds = zeros(1,12);

for i = 1:12
    F = @(x) wblcdf(x,lambda(i),k(i));
    c = F(b) - F(a); %integral(f, a, b);
    Finv = @(x) wblinv(x, lambda(i), k(i));
    FCondInv = @(x) Finv(x*c + F(a));  
    powers(:,i) = P(FCondInv(U));
    means(i) = mean(powers(:,i));
    stds(i) = std(powers(:,i));
end

figure
scatter(1:12, means)
figure
scatter(1:12, stds)

%% 

clc
close all
lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];
a = 3;
b = 25;
% Importance sampling: 
% http://www.inferencelab.com/importance-sampling-matlab-demo/

%x = 0:0.1:40;
%figure
%plot(x, wblpdf(x,lambda(1),k(1))/(b-a))
%hold on
%plot(x, normpdf(x, 15, 3));

N = 1000000;

normNbrs = normrnd((b-a)/2, 3, N, 1);

powers = zeros(N,12);
means = zeros(1,12);
stds = zeros(1,12);

for i = 1:12
    powers(:,i) = P(normNbrs).*wblpdf(normNbrs, lambda(i), k(i))./normpdf(normNbrs,(b-a)/2, 3); % Divide by value at x in normal dist.
end

means = mean(powers);
stds = std(powers);
figure
scatter(1:12, means)
figure
scatter(1:12, stds)



%%
close all
a = 3;
b = 25;
N = 1000000;
U = rand(1,N/2);
Uinv = 1 - U;

powers = zeros(N/2,12);
means = zeros(1,12);
stds = zeros(1,12);

for i = 1:12
    Finv = @(x) wblinv(x, lambda(i), k(i));
    powers(:,i) = (P(Finv(U)) +  P(Finv(Uinv)))/2;
    means(i) = mean(powers(:,i));
    stds(i) = std(powers(:,i));
end

figure
scatter(1:12, means)
figure
scatter(1:12, stds)

%%

%Capacity?
mean(means)/(3.075*10^6)

%Availability?
avg_integral


