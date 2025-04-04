mu = [1 1]*11.5;
sigma = eye(2)*30;
load powercurve_V112.mat
N = 10000000;

%Combined Weibul cumulative distribution function 
k = 1.96;
lambda = 9.13;
F = @(x) wblcdf(x,lambda,k);

%Combined
alpha = 0.638;
p = 3; 
q = 1.5;
Fcomb = @(v1,v2) F(v1)'*F(v2).*(1+alpha*((1-F(v1).^p)).^q.*(1-F(v2).^p).^q);

v1 = -5:1:30;
v2 = -5:1:30;
z = Fcomb(v1, v2);

surf(v1,v2,z)

%% Combined Weibul density function 
f = @(x) wblpdf(x,lambda,k);
fcomb = @(v1, v2) f(v1).*f(v2).*(1 + alpha*((1 - F(v1).^p).^(q - 1)).*((1 - F(v2).^p).^(q - 1)).*((F(v1).^p)*(1 + p*q) - 1).*((F(v2).^p)*(1 + p*q) - 1));

[V1,V2] = meshgrid(v1,v2);

z = fcomb(V1, V2);

surf(v1,v2,z)
%% Plot of the probability density

x1 = -5:1:30;
x2 = -5:1:30;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

y = mvnpdf(X,mu,sigma);
y = reshape(y,length(x1),length(x2));

surf(x1,x2,y)
alpha 0.5
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
xlabel('v1')
ylabel('v2')
zlabel('Probability Density')
%% Plot of P(v1,v2)*f(v1,v2)
P2d = @(x,y) P(x)+P(y);
z = P2d(X1,X2);
z = reshape(z,length(x1),length(x2));

surf(x1, x2, z)
fcmb = fcomb(X1,X2);
figure
Z = fcmb.*z;
surf(x1,x2,Z)

%% Density function and P(v1,v2)*f(v1,v2) in the same figure for comparison
 
surf(x1,x2,y)
hold on
alpha 0.1
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])

Z = fcmb.*z;
surf(x1,x2,max(max(y))*Z/max(max(Z)))


%% Importance sampling on the 2D power function newP
N = 1000000;

normNbrs = mvnrnd(mu, sigma, N);
g = mvnpdf(normNbrs,mu,sigma);
pofv1 = P(normNbrs(:,1)).*fcomb(normNbrs(:,1), normNbrs(:,2))./g;
pofv2 = P(normNbrs(:,2)).*fcomb(normNbrs(:,1), normNbrs(:,2))./g;
power = pofv1 + pofv2;
Pstd = std(power);
Pmean = mean(power);
cov(pofv1, pofv2)

%% Calculate probabilities that the combined power is above or below
%% a certain threshold as well as  their confidence intervals 
pi1 = sum(pofv1 + pofv2 > 3.075e6)/N;
pi2 = sum(pofv1 + pofv2 < 3.075e6)/N;
interval = sqrt(pi1*(1 - pi1)/N);
[pi1-(1.96*interval), pi1+(1.96*interval)];

