function [result] = MCLES(X, alpha, beta, d, gamma, maxIters, gt)
%C: number of classes
C = size(unique(gt),1);
V = size(X,2);
N = size(X{1},2); % number of data points
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans

for i=1:V
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);  %normalized
end

for i=1:V
    D{i} = size(X{i},1); % dimension of each view
end
SD = 0;
M = [];
for i=1:V
    SD = SD + D{i};
    M = [M;X{i}]; %X
end
P = zeros(SD,d);
H = rand(d,N);
S = zeros(N);
F = rand(N,C);

for it=1:maxIters
    
    %------update P--------
    P = UpdateP(H,M,P);
    
    %------update H--------
    H = SMR_mtv(M,P,S,alpha);
    
    %------update S--------
    S = UpdateS(H'*H,F,beta/alpha,gamma/alpha);
    
    %------update F--------
    Z = S;
    Z= (Z+Z')/2;
    D = diag(sum(Z));
    L = D-Z;
    [F, temp, ev]=eig1(L, C, 0);
    
    %------print OBJ-------
    Obj(it) = norm((M-P*H),'fro')^2+alpha*norm((H-H*S),'fro')^2+beta*norm(S,'fro')^2+gamma*trace(F'*L*F);
    if (it>1 && (abs(Obj(it)-Obj(it-1))/Obj(it-1)) < 10^-2)
        break;
    end
end

%actual_ids= kmeans(F, C, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
l = kmeans(F,C,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[ACC, NMI, PUR] = ClusteringMeasure(gt,l);
result = [ACC NMI PUR];

