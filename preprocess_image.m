%% script to pre-process microscopy images before applying LoG kernel
%
% Hamsini Suresh, 23/01/2018
%

clc
clear
close all

% Import image
for image_id=1:1:10
strfile = sprintf('cells_%d.tif',image_id);
I1 = imread(strfile);
img_dims = size(I1);
if length(img_dims)>2
    I = rgb2gray(I1);
else
    I = I1;
end
figure(1); clf
imshow(I);

% adjust contrast
I2 = I;
Ish = imsharpen(I2, 'Radius', 1, 'Amount', 1.8, 'Threshold', 0);
% figure;
% imshowpair(I2,Ish,'montage')
% title('Sharpened image'); 
strfile1 = sprintf('sharpened_%d.tif',ty);
imwrite(Ish,strfile1)

% enhance contrast
I4 = adapthisteq(Ish);
figure(2); clf
imshow(I4)
rows2 = size(I4,1);
cols2 = size(I4,2);
I5 = I4;
strfile2 = sprintf('contrast_%d.tif',ty);
imwrite(I5, strfile2)
end