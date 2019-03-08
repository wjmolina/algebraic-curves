clear;
clc;

n = 4;
load('Image01.mat');

K = PowerBaseAlg(ConvolveAndDownsample(I, n + ceil(n / 2), 1, 1), n);
J = PowerBaseImg(K, n, size(I, 1), size(I, 2));

figure;
imshow(I);
figure;
imshow(J);
figure;
imshow(abs(I - J));