
%Combined Weibul cumulative distribution function 
k = 1.96;
lambda = 9.13;
F = @(x) wblcdf(x,lambda,k);

%Combined
alpha = 0.638;
p = 3; 
q = 1.5;
Fcomb = @(v1,v2) F(v1)'*F(v2).*(1+alpha*((1-F(v1).^p)).^q.*(1-F(v2).^p).^q);

v1 = -1:0.5:40;
v2 = -1:0.5:40;
z = Fcomb(v1, v2);

surf(v1,v2,z)

%%
%Combined Weibul density function 
f = @(x) wblpdf(x,lambda,k);
fcomb = @(v1, v2) f(v1).*f(v2).*(1 + alpha*((1 - F(v1).^p).^(q - 1)).*((1 - F(v2).^p).^(q - 1)).*((F(v1).^p)*(1 + p*q) - 1).*((F(v2).^p)*(1 + p*q) - 1));

v1 = -1:0.3:25;
v2 = -1:0.3:25;
[V1,V2] = meshgrid(v1,v2);

z = fcomb(V1, V2);

surf(v1,v2,z)
%%
b = 25;
a = 3;
mu = [1 1]*(b-a)/2;
sigma = eye(2)*30;

x1 = -10:1:30;
x2 = -10:1:30;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

y = mvnpdf(X,mu,sigma);
y = reshape(y,length(x1),length(x2));

surf(x1,x2,y)
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density')
%%

%normNbrs = mvnrnd(mu, sigma, N, 2);
P2d = @(x,y) P(x)+P(y);
z = P2d(X1,X2);
z = reshape(z,length(x1),length(x2));

surf(x1, x2, z)
newP = @(x1,x2) P2d(x1,x2).*fcomb(x1,x2);
fcmb = fcomb(X1,X2);
figure
Z = fcmb.*z;
surf(x1,x2,Z)

%%
% Importance sampling on the 2D power function newP
mu = [1 1]*(b-a)/2;
sigma = eye(2)*30;

N = 1000000;

normNbrs = mvnrnd(mu, sigma, N);
g = mvnpdf(normNbrs,mu,sigma);

power = newP(normNbrs(:,1), normNbrs(:,2))./g; % Divide by value at x in normal dist.

Pmean = mean(power);
Pstd = std(power);

%%

N = 1000000;

normNbrs = mvnrnd(mu, sigma, N);
g = mvnpdf(normNbrs,mu,sigma);

pofv1 = P(normNbrs(:,1)).*fcomb(normNbrs(:,1), normNbrs(:,2))./g;
pofv2 = P(normNbrs(:,2)).*fcomb(normNbrs(:,1), normNbrs(:,2))./g;
cov(pofv1, pofv2)

%%
var(pofv1+pofv2)
std(pofv1+pofv2)

%%
pi1 = sum(pofv1 + pofv2 > 3.075e6)/N;
pi2 = sum(pofv1 + pofv2 < 3.075e6)/N;
interval = sqrt(pi1*(1 - pi1)/N);
[pi1-1.96*interval, pi1+1.96*interval]