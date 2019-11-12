function [matout]=normcol_lessequal(matin)
% solve the proximal problem
% matout = argmin||matout-matin||_F^2, s.t. matout(:,i)<=1
     matout = matin./repmat(max(1,sqrt(sum(matin.*matin,1)+eps)),size(matin,1),1);

