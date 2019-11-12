function [Z]=UpdateS(K,F,alpha,beta)
%"Twin Learning for Similarity  and Clustering: A Unified Kernel Approach,
%Zhao Kang; Chong Peng; Qiang Cheng, AAAI-17." 
addpath('./qpc');
[~,n]=size(K);
%options = optimset( 'Algorithm','interior-point-convex','Display','off','TolFun',10^-3);

for ij=1:n
    for ji=1:n
        all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
    end
    H=2*alpha*eye(n)+2*K;
    H=(H+H')/2;
    ff=beta/2*all'-2*K(:,ij);
    
% we use the free package to solve quadratic equation: http://sigpromu.org/quadprog/index.html
    
    [Z(:,ij),err,lm] = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));
    %Z(:,ij)=quadprog(H,(beta/2*all'-2*K(:,ij))',[],[],ones(1,n),1,zeros(n,1),ones(n,1),Z(:,ij),options);
    
end
    

