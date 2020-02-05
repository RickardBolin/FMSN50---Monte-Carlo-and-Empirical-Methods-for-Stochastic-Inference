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
