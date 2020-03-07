clc
clear
close all
addpath('Data')
f = fopen('atlantic.txt');
data = textscan(f,'%s');
fclose(f);
atlantic = str2double(data{1});

[beta,mu] = est_gumbel(atlantic);