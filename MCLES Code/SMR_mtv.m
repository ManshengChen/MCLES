function [H] =  SMR_mtv(M,P,S,alpha)
N = size(S,2); % number of data points
I = eye(N);

A_syl = P'*P;
B_syl = alpha*(I-S)*(I-S)';
C_syl = -P'*M;

%X = lyap(A,B,C) solves the Sylvester equation,AX+XB+C=0 so here C is negative
H = lyap(A_syl,B_syl,C_syl);

end