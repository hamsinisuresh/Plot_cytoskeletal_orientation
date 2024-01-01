%% function applying eLoG kernel convolution on images; based on outline in Zemel et al (2010)
%
% Stefan Marien, 11/02/2016 and Hamsini Suresh, 01/02/2018
%
function eLoG()
clc
clear 
close all

index = 1;
strpaths = {'SLB1','SLB2','SLB3','SLC1','SLC2',...
    'SLC3','SLC4','SLE1','SLE2','Homog'};

str2 = strcat('Plot_cytoskeletal_orientation\raw_images\',strpaths(index));
str3 = char(str2);
addpath(str3);
numcells = [45; 60; 48; 48; 57; 48; 50; 61; 57; 52];
L = [0,-1,0; -1,4,-1; 0,-1,0];


for ty=1:1:numcells(index)    
strfile = sprintf('idcs_adapthisteq_%d.tif',ty);
I = imread(strfile);


% ------ simple eLoG ------ %

% The amplitude of the noise:
% noise=0.1;

% Parameters of the Gaussian filter:
phi_tr = 0:pi/25:(pi-pi/25);
sigma1=1; sigma2=1;
I3 = im2double(I);
f2old = zeros(size(I3,1)+2,size(I3,2)+2);

for tr=1:25
theta=phi_tr(tr);
n1=2*ceil(2*sigma1)+1; n2=2*ceil(2*sigma2)+1;

filter1=d2gauss(n1,sigma1,n2,sigma2,theta);
f1=conv2(I3,filter1,'same');
f2 = conv2(f1,L);
f3 = max(f2old, f2);

f2old = f3;
end

strfile2 = sprintf('image_eLoG_%d.tif',ty);
imwrite(f3, strfile2)
end

end

% ---------------------- Function "d2gauss.m" ---------------------- %
% This function returns a 2D Gaussian filter with size n1*n2; theta is 
% the angle that the filter rotated counter clockwise; and sigma1 and sigma2
% are the standard deviation of the gaussian functions.
function h = d2gauss(n1,std1,n2,std2,theta)
r=[cos(theta) -sin(theta);
   sin(theta)  cos(theta)];
h = zeros(n2,n1);
for i = 1 : n2 
    for j = 1 : n1
        u = r * [j-(n1+1)/2 i-(n2+1)/2]';
        h(i,j) = gauss(u(1),std1)*gauss(u(2),std2);
    end
end
h = h / sqrt(sum(sum(h.*h)));
end
% Function "gauss.m":
function y = gauss(x,std)
y = exp(-x^2/(2*std^2)) / (std*sqrt(2*pi));
end