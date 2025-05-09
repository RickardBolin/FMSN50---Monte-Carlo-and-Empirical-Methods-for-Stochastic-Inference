clear
clc
close all

load powercurve_V112.mat

lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k =      [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};


N = 1000000;

%Standard Monte Carlo
powers = zeros(N,12);
for i = 1:12
    randomWeibuls = wblrnd(lambda(i), k(i), N, 1);
    powers(:,i) = P(randomWeibuls);
end
standard_means = mean(powers);
standard_stds = std(powers);

standard_conf = [standard_means - 1.96*(standard_stds/sqrt(N)); standard_means + 1.96*(standard_stds/sqrt(N))];
standard_width = abs(standard_conf(1,:)-standard_conf(2,:));

figure
scatter(1:12, standard_means, 'rx')
hold on 
m1 = plot(1:12, standard_means, 'r-');
scatter(1:12,standard_conf(2,:), 'b^')
scatter(1:12, standard_conf(1,:), 'bV')
s1 = plot(1:12, standard_conf(2,:), 'b-');
plot(1:12, standard_conf(1,:), 'b-')
hold off
xticks(1:12)
xticklabels(months)
xlim([1 12])
legend([m1 s1],{'Mean','Standard'})


%Truncated inverse sampling
a = 3;
b = 25;
U = rand(1,N);

powers = zeros(N,12);
for i = 1:12
    F = @(x) wblcdf(x,lambda(i),k(i));
    c = F(b) - F(a); %integral(f, a, b);
    Finv = @(x) wblinv(x, lambda(i), k(i));
    FCondInv = @(x) Finv(x*c + F(a));
    powers(:,i) = c*P(FCondInv(U));
end
trunc_means = mean(powers);
trunc_stds = std(powers);

trunc_conf = [standard_means - 1.96*(trunc_stds/sqrt(N)); standard_means + 1.96*(trunc_stds/sqrt(N))];
trunc_width = abs(trunc_conf(1,:)-trunc_conf(2,:));


figure
scatter(1:12, standard_means, 'rx');
hold on 
m1 = plot(1:12, standard_means, 'r-');

scatter(1:12,standard_conf(2,:), 'b^')
scatter(1:12, standard_conf(1,:), 'bV')
s1 = plot(1:12, standard_conf(2,:), 'b-');
plot(1:12, standard_conf(1,:), 'b-')

scatter(1:12, trunc_conf(2,:), 'g^')
scatter(1:12, trunc_conf(1,:), 'gV')
s2 = plot(1:12, trunc_conf(2,:), 'g-');
plot(1:12, trunc_conf(1,:), 'g-')
hold off
xticks(1:12)
xticklabels(months)
xlim([1 12])
legend([m1 s1 s2],{'Mean','Standard','Truncated'})

%% 
%Importance Sampling

mu = 11.5;
sigma = 5;
x = 0:0.1:40;
g = @(x) normpdf(x, mu, sigma);
gx = g(x);
phi = max(gx)*0.5e-5*P(x).*wblpdf(x,mean(lambda(1)),mean(k(1)))';
figure
plot(x, phi)
hold on
plot(x, gx);
hold off
legend({'phi(x)','g(x)'})

normNbrs = normrnd(mu, sigma, N, 1);

powers = zeros(N,12);
for i = 1:12
    powers(:,i) = P(normNbrs).*wblpdf(normNbrs, lambda(i), k(i))./g(normNbrs); % Divide by value at x in normal dist.
end

imp_means = mean(powers);
imp_stds = std(powers);

imp_conf = [standard_means - 1.96*(imp_stds/sqrt(N)); standard_means + 1.96*(imp_stds/sqrt(N))];
imp_width = abs(imp_conf(1,:)-imp_conf(2,:));

figure
scatter(1:12, standard_means, 'rx');
hold on 
m1 = plot(1:12, standard_means, 'r-');


scatter(1:12,standard_conf(2,:), 'b^')
scatter(1:12, standard_conf(1,:), 'bV')
s1 = plot(1:12, standard_conf(2,:), 'b-');
plot(1:12, standard_conf(1,:), 'b-')

scatter(1:12, trunc_conf(2,:), 'g^')
scatter(1:12, trunc_conf(1,:), 'gV')
s2 = plot(1:12, trunc_conf(2,:), 'g-');
plot(1:12, trunc_conf(1,:), 'g-')

scatter(1:12, imp_conf(2,:), 'c^')
scatter(1:12, imp_conf(1,:), 'cV')
s3 = plot(1:12, imp_conf(2,:), 'c-');
plot(1:12, imp_conf(1,:), 'c-')
hold off
xticks(1:12)
xticklabels(months)
xlim([1 12])
legend([m1 s1 s2 s3],{'Mean','Standard','Truncated','Importance Sampling'})

%%
%Antithetic sampling

U = rand(1,N/2);
Uinv = 1 - U;

powers = zeros(N/2,12);

for i = 1:12
    Finv = @(x) wblinv(x, lambda(i), k(i));
    powers(:,i) = (P(Finv(U)) +  P(Finv(Uinv)))/2;

end

anti_means = mean(powers);
anti_stds = std(powers);

anti_conf = [standard_means - 1.96*(anti_stds/sqrt(N)); standard_means + 1.96*(anti_stds/sqrt(N))];
anti_width = abs(anti_conf(1,:)-anti_conf(2,:));

figure
scatter(1:12, standard_means, 'rx');
hold on 
m1 = plot(1:12, standard_means, 'r-');

scatter(1:12,standard_conf(2,:), 'b^')
scatter(1:12, standard_conf(1,:), 'bV')
s1 = plot(1:12, standard_conf(2,:), 'b-');
plot(1:12, standard_conf(1,:), 'b-')

scatter(1:12, trunc_conf(2,:), 'g^')
scatter(1:12, trunc_conf(1,:), 'gV')
s2 = plot(1:12, trunc_conf(2,:), 'g-');
plot(1:12, trunc_conf(1,:), 'g-')

scatter(1:12, imp_conf(2,:), 'c^')
scatter(1:12, imp_conf(1,:), 'cV')
s3 = plot(1:12, imp_conf(2,:), 'c-');
plot(1:12, imp_conf(1,:), 'c-')

scatter(1:12, anti_conf(2,:), 'm^')
scatter(1:12, anti_conf(1,:), 'mV')
s4 = plot(1:12, anti_conf(2,:), 'm-');
plot(1:12, anti_conf(1,:), 'm-')

hold off
xticks(1:12)
xticklabels(months)
xlim([1 12])
legend([m1 s1 s2 s3 s4],{'Mean','Standard','Truncated','Importance Sampling','Antithetic'})



%%
%Power coefficient
rho = 1.225;
d = 112;
avgPtot = 0;
powers = zeros(N,12);
means = zeros(1,12);

for i = 1:12
    Ptot = @(v) ((rho*pi*(d.^2).*(v.^3))/8).*wblpdf(v,lambda(i), k(i));
    pow = integral(Ptot, 0, 100);
    avgPtot = avgPtot + pow/12;
    randomWeibuls = wblrnd(lambda(i), k(i), N, 1);
    powers(:,i) = P(randomWeibuls);
    means(i) = mean(powers(:,i));
end

[sortedPowers,I] = sort(powers(:,12));
avgPowerCoeff = mean(means/avgPtot);

stds = std(means/avgPtot);
conf = [avgPowerCoeff - 1.96*stds/sqrt(N); avgPowerCoeff + 1.96*stds/sqrt(N)];


%%
%Capacity
capacities = zeros(1,12);
for i = 1:12
    capacities(i) = means(i)/3.075e6;
end
figure
scatter(1:12, capacities(1:12))
mean(capacities)

%%
figure
scatter(1:12, standard_means, 'rx');


hold on 
scatter(1:12, trunc_means, 'gx');
scatter(1:12, imp_means, 'cx');
scatter(1:12, anti_means, 'mx');
m1 = plot(1:12, standard_means, 'r-');
m2 = plot(1:12, trunc_means, 'g-');
m3 = plot(1:12, imp_means, 'c-');
m4 = plot(1:12, anti_means, 'm-');


hold off
xticks(1:12)
xticklabels(months)
xlim([1 12])
legend([m1 s1 s2 s3 s4],{'Standard','Truncated','Importance Sampling','Antithetic'})