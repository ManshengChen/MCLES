close all;
load('MSRC-v1.mat');
%X = {gist',hog',lbp'};
gt = Y;
for i=1:4
    temp = X{i};  %Xv(dv*n)
    X{i} = temp';
end

maxIters = 30;
alpha = 0.8;
beta = 0.4;
d = 70; %lower dimension
gamma = 0.004;

result = MCLES(X, alpha, beta, d, gamma, maxIters, gt)

