clc
clear
close all
addpath('Data')
f = fopen('atlantic.txt');
data = textscan(f,'%s');
fclose(f);
atlantic = str2double(data{1});
n = length(atlantic);
%Slajd 33
%http://www.maths.lth.se/matstat/kurser/fmsn50masm11/2020/material/L13.pdf
[beta,mu] = est_gumbel(atlantic);

Finv = @(x, mu1, beta1) mu1 - beta1*log(-log(x));

nsamples = 1000;
bs_mus = zeros(nsamples,1);
bs_betas = zeros(nsamples,1);
bs_retwave = zeros(nsamples,1);

T = 3*14*100;
big_wave = Finv(1 - 1/T, mu, beta);


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

delta_wave =  big_wave - bs_retwave;
delta_mu = mu - bs_mus;
delta_beta = beta - bs_betas;
delta_mu = sort(delta_mu);
delta_beta = sort(delta_beta);
delta_wave = sort(delta_wave);
alpha = 0.05;

conf_mu =  mu - [delta_mu(ceil((1 - alpha/2)*nsamples)), delta_mu(ceil(alpha*nsamples/2))]
conf_beta =  beta - [delta_beta(ceil((1 - alpha/2)*nsamples)), delta_beta(ceil(alpha*nsamples/2))]
conf_wave = big_wave + delta_wave(ceil((1 - alpha)*nsamples))

