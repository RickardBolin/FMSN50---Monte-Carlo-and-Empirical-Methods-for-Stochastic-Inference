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
    f = @(x) wblpdf(x,lambda(i),k(i));
    c = wblcdf(b, lambda(i), k(i)) - wblcdf(a, lambda(i), k(i)); %integral(f, a, b);
    Finv = @(x) wblinv(x, lambda(i), k(i));
    FCondInv = @(x) Finv(x*c + f(a));  
    %plot(x, FCondInv(x))    
    powers(:,i) = P(sort(FCondInv(U)));
    means(i) = mean(powers(:,i));
    stds(i) = std(powers(:,i));
end

figure
scatter(sort(FCondInv(U)), powers(:,1))  %Kolla denna, skumt att vindarna bara g?r till typ 20? Borde v?l g? till 25?
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
    powers(:,i) = P(normNbrs).*wblpdf(normNbrs, lambda(i), k(i))./normpdf(normNbrs,(b-a)/2, 3); x% Divide by value at x in normal dist.
end

means = mean(powers);
stds = std(powers);
figure
scatter(1:12, means)
figure
scatter(1:12, stds)

%%


%L?gg in antithetic sampling h?r!




%%
N = 1000000;
a = 3;
b = 25;

box_height = 0.1;
box_width = b-a;
x = a + (b-a)*sort(rand(1,N),1,'descend');
y = box_height*sort(rand(1,N),1,'descend');
integral_sum = 0;
for i = 1:12
    weibulDist = @(x) wblpdf(x,lambda(i),k(i));
    bool = y < weibulDist(x);
    area = sum(bool)/N;
    integral_approx = area*box_width*box_height;
    integral_sum = integral_sum + integral_approx;
end
avg_integral = integral_sum/12;

%%

%Capacity?
mean(means)/(3.075*10^6)

%Availability?
avg_integral

%%
%Combined Weibul cumulative distribution function 
k = 1.96;
lambda = 9.13;
F = @(x) wblcdf(x,lambda,k);

%Combined
alpha = 0.638;
p = 3; 
q = 1.5;
Fcomb = @(v1,v2) F(v1)'*F(v2).*(1+alpha*((1-F(v1).^p)).^q.*(1-F(v2).^p).^q);

v1 = -1:0.3:25;
v2 = -1:0.3:25;
z = Fcomb(v1, v2);

surf(v1,v2,z)

%%
%Combined Weibul density function 
f = @(x) wblpdf(x,lambda,k);
fcomb = @(v1, v2) f(v1)'*f(v2).*(1 + alpha*((1 - F(v1).^p).^(q - 1)).*((1 - F(v2).^p).^(q - 1)).*((F(v1).^p)*(1 + p*q) - 1).*((F(v2).^p)*(1 + p*q) - 1));

v1 = -1:0.3:25;
v2 = -1:0.3:25;
z = fcomb(v1, v2);

surf(v1,v2,z)
%%
mu = [1 1]*(b-a)/2;
sigma = eye(2)*30;

x1 = -10:1:30;
x2 = -10:1:30;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

y = mvnpdf(X,mu,sigma);
y = reshape(y,length(x2),length(x1));

surf(x1,x2,y)
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density')
%%

%normNbrs = mvnrnd(mu, sigma, N, 2);


newP = @(x1,x2) P(x1).*fcomb(x1,x2);
surf(x1,x2,newP(x1,x2))




%%

a = 3;
b = 25;
N = 1000000;
U = rand(1,N/2);
Uinv = 1 - U;

powers = zeros(N,12);
means = zeros(1,12);
stds = zeros(1,12);

for i = 1:12
    f = @(x) wblpdf(x,lambda(i),k(i));
    c = integral(f, a, b);
    Finv = @(x) wblinv(x, lambda(i), k(i));
    FCondInv = @(x) Finv(x*c + f(a));  
    %plot(x, FCondInv(x))    
    powers(1:N/2,i) = P(sort(Finv(U)));
    powers(N/2 + 1:end, i) = P(sort(Finv(Uinv)));
    means(i) = mean(powers(:,i));
    stds(i) = std(powers(:,i));
end



