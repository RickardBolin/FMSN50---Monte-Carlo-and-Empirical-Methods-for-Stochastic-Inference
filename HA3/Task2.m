clc
clear
close all
addpath('Data')
f = fopen('atlantic.txt');
data = textscan(f,'%s');
fclose(f);
atlantic = str2double(data{1});
n = length(atlantic);

[beta,mu] = est_gumbel(atlantic);

Finv = @(x, mu1, beta1) mu1 - beta1*log(-log(x));

nsamples = 1000;
bs_mus = zeros(nsamples,1);
bs_betas = zeros(nsamples,1);
bs_retwave = zeros(nsamples,1);

T = 3*14*100;

for i = 1:nsamples
    u = rand(n,1);
    fake_atlantic = Finv(u,mu,beta);
    [bs_beta, bs_mu] = est_gumbel(fake_atlantic);
    
    bs_retwave(i) = Finv(1 - 1/T, bs_mu, bs_beta);
    bs_mus(i) = bs_mu;
    bs_betas(i) = bs_beta;
    
    if(mod(i,floor(nsamples/100)) == 0)
       q = i/nsamples;
       clc
       disp([num2str(100*q) '% |' char(ones(1,floor(50*q))*'=') char(ones(1, ceil(50 - 50*q))*' ') '|'])       
    end
end
conf_mu = mu + 1.96*std(bs_mus)*[-1,1]
conf_beta = beta + 1.96*std(bs_betas)*[-1,1]
conf_retwave = mean(bs_retwave) + 1.64*std(bs_retwave)*[-1,1]





