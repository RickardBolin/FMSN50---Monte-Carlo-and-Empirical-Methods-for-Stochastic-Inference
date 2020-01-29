clear
clc

load powercurve_V112.mat

lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k =      [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];

avgLambda = mean(lambda);
avgK = mean(k);

avgWind = 8.1;

t = 0:0.1:40;
weibulDist = wblpdf(t,lambda(1),k(1));

%figure
%plot(t, weibulDist)
%figure
%plot(t, P(t))
samples = 100000;
powers = zeros(samples,12);
means = zeros(1,12);
variances = zeros(1,12);
for i = 1:12
    randomWeibuls = wblrnd(lambda(i), k(i), samples, 1);
    powers(:,i) = P(randomWeibuls);
    means(i) = mean(powers(:,i));
    variances(i) = std(powers(:,i));
end

figure
scatter(1:12, means)
figure
scatter(1:12, variances)