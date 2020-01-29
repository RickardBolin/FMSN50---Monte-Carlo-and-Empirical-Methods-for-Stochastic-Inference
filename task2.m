clear
clc

load powercurve_V112.mat

lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k =      [2.0  2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0  2.0];

avgLambda = mean(lambda);
avgK = mean(k);

avgWind = 8.1;

t = 0:0.1:40;
weibulDist = wblpdf(t,lambda(1),k(1));

figure
plot(t, weibulDist)
figure
plot(t, P(t))

randomWeibuls = wblrnd(lambda, k, 1000);
powers = P(randomWeibuls);
means = mean(p);
variances = var(powers);

figure
scatter(1:12, means)
figure
scatter(1:12, variances)