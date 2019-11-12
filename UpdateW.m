function [D_Mat] = UpdateP(Coef, Data, D_Mat)
%Coef means H
%Data means M
%Input D_Mat means the previous P

Imat= eye(size(Coef,1)); %number of rows of Coef

TempCoef = Coef;
TempData = Data;
rho = 1;
rate_rho = 1.2;
TempS = D_Mat;
TempT = zeros(size(TempS));
previousD = D_Mat;
Iter = 1;ERROR=1;
while(ERROR>1e-6&&Iter<100)
    
    TempD   = (rho*(TempS-TempT)+TempData*TempCoef')/(rho*Imat+TempCoef*TempCoef');
    TempS   = normcol_lessequal(TempD+TempT);
    TempT   = TempT+TempD-TempS;
    rho     = rate_rho*rho;
    ERROR = mean(mean((previousD- TempD).^2));
    previousD = TempD;
    Iter=Iter+1;
end
D_Mat = TempD;
